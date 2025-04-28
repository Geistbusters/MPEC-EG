c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine updateg(c0,c1,c2,c3,ih0mo,ioooo,ionex,jmats,kmats,i3x, 1d9s18
     $    ih02,ioooo2,idoub,iacto,noc,nvirtc,nbasdwsc,ignew,iden,
     $j2den,itmat,iumat,itotm,rmssz,iter1,idwsdeb,ios,ior,numpro,bc,ibc)11d9s22
      implicit real*8 (a-h,o-z)
      integer*8 ipack,i8a,i8b,i8c,i8d,i8e,i8f,i8g,i8h,i8i,i8j,i8k
      dimension ioooo2(*),idoub(*),iacto(*),noc(*),ignew(*),nvirtc(*),
     $     idoit(5),nbasdwsc(*),iden(*),ioooo(*),ionex(*),jmats(*),     1d9s18
     $     kmats(*),i3x(*),j2den(*),itmat(*),iumat(*),itotm(*),rmssz(8) 4d2s18
      include "common.store"
      include "common.hf"
      include "common.print"                                            2d13s20
      logical lprint                                                    4d5s18
      dimension ios(numpro,idbk),ior(numpro,idbk),ioonn(idbk),          2d26s19
     $     iovnn(idbk),                                                 2d26s19
     $     iokx(idbk),i1xa(idbk),i1xb(idbk),i4ok(idbk),ijnn(idbk),      3d16s18
     $     ionex2(idbk),ionex3(idbk),i3xa(4,idbk),i3xb(2,idbk),
     $     iok3x(3,idbk)                                                4d13s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      equivalence (pack,ipack)
      data idoit/5*1/
      common/unitcm/iunit                                               11d9s17
      data icall/0/
      save
      lprint=.false.                                                     4d10s18
      if(idwsdeb.gt.1)lprint=.true.                                     5d2s18
      if(iprtr(14).ne.0)lprint=.true.                                   2d13s20
c
c     transform integrals to new basis so we can compute the new
c     gradient and energy, then compute the new gradient.
c
c     general scheme: transform indicies local to each processor, then
c     reshuffle across processors so that the remaining two indices are
c     now local so that we can complete the transformation.
c     For K or 3x, we need to preload so that we have nbasdwsc rather
c     than just nvirtc indices.
c
c     form product of act-doub and virt-occ rotations.
c
      icall=icall+1
      xnan=3.14d0
      if(lprint)then                                                    2d14s20
       write(6,*)('in updateg ...'),icall,xnan                          2d14s20
       write(6,*)('starting h0: ')                                      2d14s20
       jh0mo=ih0mo                                                      2d14s20
       do isb=1,nsymb                                                   2d14s20
        if(noc(isb).gt.0)then                                           2d14s20
         write(6,*)('for symmetry block '),isb                          2d14s20
         call prntm2(bc(jh0mo),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         jh0mo=jh0mo+nbasdwsc(isb)**2                                   2d14s20
        end if                                                          2d14s20
       end do                                                           2d14s20
       write(6,*)('starting oooo: ')
       do is=1,nsdlk                                                    2d14s20
        if(isblk(1,is).eq.isblk(2,is))then                              2d14s20
         nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 2d14s20
        else                                                            2d14s20
         nrow=noc(isblk(1,is))*noc(isblk(2,is))                         2d14s20
        end if                                                          2d14s20
        if(isblk(3,is).eq.isblk(4,is))then                              2d14s20
         ncol=(noc(isblk(3,is))*(noc(isblk(3,is))+1))/2                 2d14s20
        else                                                            2d14s20
         ncol=noc(isblk(3,is))*noc(isblk(4,is))                         2d14s20
        end if                                                          2d14s20
        if(min(nrow,ncol).gt.0)then                                     2d14s20
         write(6,*)('integral type '),(isblk(j,is),j=1,4)
         call prntm2(bc(ioooo(is)),nrow,ncol,nrow)                      2d14s20
        end if                                                          2d14s20
       end do                                                           2d14s20
      end if                                                            2d14s20
      do isb=1,nsymb                                                    1d10s18
       if(nbasdwsc(isb).gt.0.and.noc(isb).gt.0)then                     2d25s19
       if(lprint)then                                                   2d14s20
        write(6,*)('for symmetry block '),isb
        write(6,*)('starting umat: ')
        call prntm2(bc(iumat(isb)),nbasdwsc(isb),noc(isb),nbasdwsc(isb))2d14s20
        write(6,*)('starting tmat: ')
        call prntm2(bc(itmat(isb)),noc(isb),noc(isb),noc(isb))
       end if                                                           2d14s20
       call dgemm('n','n',nbasdwsc(isb),noc(isb),noc(isb),1d0,
     $      bc(iumat(isb)),nbasdwsc(isb),bc(itmat(isb)),noc(isb),0d0,   1d10s18
     $      bc(itotm(isb)),nbasdwsc(isb),                               1d10s18
     d' updateg.  1')
       end if                                                           2d25s19
       iad1=iumat(isb)+nbasdwsc(isb)*noc(isb)-1                         1d10s18
       iad2=itotm(isb)+nbasdwsc(isb)*noc(isb)-1                         1d10s18
       do i=1,nvirtc(isb)                                               1d10s18
        do j=1,nbasdwsc(isb)                                            1d10s18
         bc(iad2+j)=bc(iad1+j)                                          1d10s18
        end do                                                          1d10s18
        iad2=iad2+nbasdwsc(isb)                                         1d10s18
        iad1=iad1+nbasdwsc(isb)                                         1d10s18
       end do                                                           1d10s18
      end do                                                            1d10s18
c
c     make room for temporary storage of oooo2: integrals for           1d10s18
c     new energy calculation                                            1d10s18
c
      ntot2=0d0                                                         1d10s18
      do is=1,nsdlk                                                     1d10s18
       ioooo2(is)=ibcoff                                                1d10s18
       if(isblk(1,is).eq.isblk(2,is))then
        nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                  1d10s18
       else
        nrow=noc(isblk(1,is))*noc(isblk(2,is))                          1d10s18
       end if
       need=nrow*noc(isblk(3,is))*noc(isblk(4,is))                      1d10s18
       ibcoff=ibcoff+need                                               1d10s18
       ntot2=ntot2+need                                                 1d10s18
      end do                                                            1d10s18
      call enough('updateg.  1',bc,ibc)
      do i=0,ntot2-1                                                    1d10s18
       bc(ioooo2(1)+i)=0d0                                              3d26s18
      end do                                                            1d10s18
c
c     ih02 will have h0mo transformed to new basis.                     1d16s17
c     however we will also need h0mo half transformed by tmat and       1d16s17
c     fully transformed by umat. this will be stored in ih02h.          1d16s17
c
      ih02=ibcoff                                                       1d16s17
      do isb=1,nsymb                                                    1d16s17
       ibcoff=ibcoff+nbasdwsc(isb)**2                                   1d16s17
      end do                                                            1d16s17
      ih02h=ibcoff                                                      1d16s17
      do isb=1,nsymb                                                    1d16s17
       ibcoff=ibcoff+nbasdwsc(isb)*noc(isb)                             1d16s17
      end do                                                            1d16s17
      call enough('updateg.  2',bc,ibc)
      jh0mo=ih0mo                                                       1d16s17
      jh02=ih02                                                         1d16s17
      jh02h=ih02h                                                       1d16s17
      do isb=1,nsymb                                                    1d16s17
       itmp=ibcoff                                                      1d16s17
       ibcoff=itmp+nbasdwsc(isb)**2                                     1d16s17
       call enough('updateg.  3',bc,ibc)
       if(nbasdwsc(isb).gt.0)then                                       2d25s19
       call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),1d0,1d16s17
     $      bc(jh0mo),nbasdwsc(isb),bc(itotm(isb)),nbasdwsc(isb),0d0,   1d16s17
     $      bc(itmp),nbasdwsc(isb),                                     1d16s17
     d' updateg.  2')
       end if                                                           2d25s19
       itmpt=ibcoff                                                     1d16s17
       ibcoff=itmpt+nbasdwsc(isb)**2                                    1d16s17
       call enough('updateg.  4',bc,ibc)
       do i=0,nbasdwsc(isb)-1                                           1d16s17
        do j=0,nbasdwsc(isb)-1                                          1d16s17
         ji=itmp+j+nbasdwsc(isb)*i                                      1d16s17
         ij=itmpt+i+nbasdwsc(isb)*j                                     1d16s17
         bc(ij)=bc(ji)                                                  1d16s17
        end do                                                          1d16s17
       end do                                                           1d16s17
       if(nbasdwsc(isb).gt.0.and.noc(isb).gt.0)then                     2d25s19
       call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),1d0,1d16s17
     $       bc(itmpt),nbasdwsc(isb),bc(itotm(isb)),nbasdwsc(isb),0d0,  1d16s17
     $      bc(jh02),nbasdwsc(isb),                                     1d16s17
     d' updateg.  3')
       call dgemm('n','n',noc(isb),nbasdwsc(isb),nbasdwsc(isb),1d0,     4d4s18
     $      bc(itmpt),nbasdwsc(isb),bc(iumat(isb)),nbasdwsc(isb),0d0,   4d4s18
     $      bc(itmp),noc(isb),                                          4d4s18
     d' updateg.  4')
       do i=0,nbasdwsc(isb)-1                                           4d4s18
        do j=0,noc(isb)-1                                               4d4s18
         ji=itmp+j+noc(isb)*i                                           4d4s18
         ij=jh02h+i+nbasdwsc(isb)*j                                     4d4s18
         bc(ij)=bc(ji)                                                  4d4s18
        end do                                                          4d4s18
       end do                                                           4d4s18
       end if                                                           2d25s19
       ibcoff=itmp                                                      1d16s17
       jh0mo=jh0mo+nbasdwsc(isb)**2                                     1d16s17
       jh02=jh02+nbasdwsc(isb)**2                                       1d16s17
       jh02h=jh02h+nbasdwsc(isb)*noc(isb)                               1d16s17
      end do                                                            1d16s17
c
c     ioonn points to space for the half transformed integrals.         1d11s18
c
      do is=1,nsdlk                                                     1d10s18
       if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.          1d10s18
     $      noc(isblk(3,is)).ne.0.and.noc(isblk(4,is)).ne.0)then        1d10s18
        call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,         1d10s18
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            1d10s18
        ncol=ih+1-il                                                    1d10s18
        nrow=noc(isblk(3,is))*noc(isblk(4,is))                          1d10s18
        ioonn(is)=ibcoff                                                1d10s18
        ibcoff=ioonn(is)+nrow*ncol                                      1d11s18
       end if                                                           1d11s18
      end do                                                            1d10s18
      if(c1.ne.0d0)then                                                 1d26s18
c
c     iovnn points to space for the half transformed integrals.         1d11s18
c
       do is=1,nsdlk1                                                   1d26s18
        if(noc(isblk1(1,is)).ne.0.and.noc(isblk1(2,is)).ne.0.and.       1d26s18
     $      noc(isblk1(3,is)).ne.0.and.nbasdwsc(isblk1(4,is)).ne.0)then 1d27s18
         call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,      1d27s18
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            1d10s18
         ncol=ih+1-il                                                   1d26s18
         nrow=noc(isblk1(3,is))*nvirtc(isblk1(4,is))                    1d27s18
         iovnn(is)=ibcoff                                               1d26s18
         ibcoff=iovnn(is)+nrow*ncol                                     1d26s18
        end if                                                          1d26s18
       end do                                                           1d26s18
      end if                                                            1d26s18
      call enough('updateg.  5',bc,ibc)
      if(c2.ne.0d0)then                                                 2d22s18
c                                                                       3d16s18
c     iknn, ijnn points to space for half transformed kmats jmats       3d16s18
c
       do is=1,nsdlk1                                                   1d26s18
        if(noc(isblk1(1,is)).ne.0.and.noc(isblk1(2,is)).ne.0.and.       1d26s18
     $      noc(isblk1(3,is)).ne.0.and.nbasdwsc(isblk1(4,is)).ne.0)then 1d27s18
         call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,      1d27s18
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            1d10s18
         ncol=ih+1-il                                                   1d26s18
         nrow=noc(isblk1(3,is))*nvirtc(isblk1(4,is))                    1d27s18
         ionex2(is)=ibcoff                                              3d16s18
         ibcoff=ionex2(is)+nrow*ncol                                    3d16s18
         ionex3(is)=ibcoff                                              3d29s18
         n12=noc(isblk1(1,is))*noc(isblk1(2,is))                        3d29s18
         ibcoff=ionex3(is)+nrow*n12                                     3d29s18
         call enough('updateg.  6',bc,ibc)
         do i=0,nrow*(ncol+n12)-1                                       3d29s18
          bc(ionex2(is)+i)=0d0                                          3d16s18
         end do                                                         3d16s18
        end if                                                          1d26s18
       end do                                                           1d26s18
       do is=1,nsdlk                                                    3d16s18
        if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.         3d16s18
     $       nbasdwsc(isblk(3,is)).ne.0.and.nbasdwsc(isblk(4,is)).ne.0) 3d16s18
     $       then                                                       3d16s18
         call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,        3d16s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           3d16s18
         ncol=ih+1-il                                                   3d16s18
         nrow=nbasdwsc(isblk(3,is))*nbasdwsc(isblk(4,is))               3d16s18
         ijnn(is)=ibcoff                                                3d16s18
         ibcoff=ijnn(is)+nrow*ncol                                      3d16s18
        end if                                                          3d16s18
       end do                                                           3d16s18
c
c     for transformation of K and 3x we need additional information
c     so get that now.
c
c     space for results
c
       do is=1,nsdlkk                                                   3d16s18
        iokx(is)=ibcoff                                                 2d22s18
c
c     recall K_nm^ab is (nb|ma) with ab distributed across procs.
c     ab are virtual indices, but we need occ as well for the
c     transformation, so we need (no|ma), (nb|mo), (no|mp) with o p occ
c     distributed over ao, ob, and mp, respectively.
c
        nrow=nbasdwsc(isblkk(3,is))*nbasdwsc(isblkk(4,is))              3d16s18
        call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,       3d16s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           3d16s18
        nhere=ih+1-il                                                   3d16s18
        ibcoff=ibcoff+nrow*nhere                                        3d16s18
       end do                                                           2d22s18
       if(c3.ne.0d0)then                                                4d13s18
        do is3=1,nsdlk1                                                 4d13s18
         nn=noc(isblk1(3,is3))*nbasdwsc(isblk1(4,is3))                  4d13s18
         if(nn.ne.0)then                                                4d13s18
          call ilimts(nvirtc(isblk1(1,is3)),noc(isblk1(2,is3)),         4d13s18
     $        mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                 4d13s18
          nhere=ih+1-il                                                 4d13s18
          iok3x(1,is3)=ibcoff                                           4d13s18
          ibcoff=ibcoff+nn*nhere                                        4d13s18
          call enough('updateg.  7',bc,ibc)
          do i=0,nn*nhere-1
           bc(iok3x(1,is3)+i)=xnan
          end do
          call ilimts(noc(isblk1(1,is3)),noc(isblk1(2,is3)),            4d13s18
     $        mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                 4d13s18
          nhere=ih+1-il                                                 4d13s18
          iok3x(2,is3)=ibcoff                                           4d13s18
          ibcoff=ibcoff+nn*nhere                                        4d13s18
          call enough('updateg.  8',bc,ibc)
          do i=0,nn*nhere-1
           bc(iok3x(2,is3)+i)=xnan
          end do
          if(isblk1(1,is3).ne.isblk1(2,is3))then                        4d13s18
           call ilimts(noc(isblk1(1,is3)),nvirtc(isblk1(2,is3)),         4d13s18
     $        mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                 4d13s18
           nhere=ih+1-il                                                 4d13s18
           iok3x(3,is3)=ibcoff                                           4d13s18
           ibcoff=ibcoff+nn*nhere                                        4d13s18
           call enough('updateg.  9',bc,ibc)
           do i=0,nn*nhere-1
            bc(iok3x(3,is3)+i)=xnan
           end do
          end if                                                        4d13s18
         end if                                                         4d13s18
        end do                                                          4d13s18
        call enough('updateg. 10',bc,ibc)
       end if                                                           4d13s18
c
c     for debugging...
c
c
c     ipass = 1 means count, ipass = 2 means do
c
       isend=ibcoff                                                     2d22s18
       irecv=isend+mynprocg                                             2d22s18
       nsend=irecv+mynprocg                                             2d22s18
       nrecv=nsend+mynprocg                                             2d22s18
       ibcoff=nrecv+mynprocg                                            2d22s18
       call enough('updateg. 11',bc,ibc)
       do i=0,mynprocg-1                                                 1d10s18
        bc(nsend+i)=0                                                    1d10s18
        bc(nrecv+i)=0                                                    1d10s18
        do j=1,idbk
         ios(i+1,j)=0
         ior(i+1,j)=0
        end do                                                           1d10s18
       end do                                                            1d10s18
       do ipass=1,2                                                     2d22s18
c
c     for sending
c
        me2mea=0
        me2meb=0
        me2mec=0
        me2med=0
        if(c3.ne.0d0)then                                               4d10s18
         do is3=1,nsdlk1                                                4d10s18
          call ilimts(noc(isblk1(3,is3)),noc(isblk1(4,is3)),            4d10s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               4d10s18
          ncolb=noc(isblk1(3,is3))*noc(isblk1(4,is3))                    4d10s18
          nhere0b=ncolb/mynprocg                                          4d10s18
          nleftb=ncolb-nhere0b*mynprocg                                    4d10s18
          nhereb=nhere0b+1                                                4d10s18
          nhere0b=max(1,nhere0b)                                          4d10s18
          ncutb=nleftb*nhereb                                           4d10s18
          call ilimts(noc(isblk1(3,is3)),nvirtc(isblk1(4,is3)),         4d10s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               4d10s18
          ncol=noc(isblk1(3,is3))*nvirtc(isblk1(4,is3))                 4d10s18
          nhere0=ncol/mynprocg                                          4d10s18
          nleft=ncol-nhere0*mynprocg                                    4d10s18
          nhere=nhere0+1                                                4d10s18
          nhere0=max(1,nhere0)                                          4d10s18
          ncut=nleft*nhere                                              4d10s18
c
c     to transform 1st 2 indices, need (nm|ov), (v'm|ov) and (nv'|ov)   4d10s18
c     and to transform the last 2 indices, need (vm|op), (nm|op),       4d10s18
c     (nv|op), and (vv'|op). (nm|op) is usually already in place, but
c     we need to account for it being distributed (nm|po).              4d13s18
c     (nm|ov) is already in place.                                      4d13s18
c
          do is1=1,nsdlk1                                               4d10s18
           if(noc(isblk1(4,is3)).ne.0.and.noc(isblk1(2,is3)).ne.0)then  4d10s18
c     (vm|op)
           if(isblk1(4,is1).eq.isblk1(1,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(2,is3).and.                     4d10s18
     $          isblk1(1,is1).eq.isblk1(3,is3))then                     4d10s18
            call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            if(isblk1(1,is1).eq.isblk1(2,is1))then                      4d10s18
             nrow1=(noc(isblk1(1,is1))*(noc(isblk1(1,is1))+1))/2        4d10s18
             iswitch=0                                                  4d10s18
            else                                                        4d10s18
             nrow1=noc(isblk1(1,is1))*noc(isblk1(2,is1))                4d10s18
             iswitch=1                                                  4d10s18
            end if                                                      4d10s18
            do i4=0,noc(isblk1(2,is1))-1                                4d10s18
             do i3=0,noc(isblk1(1,is1))-1                               4d10s18
              icol=i3+noc(isblk1(1,is1))*i4                             4d10s18
              if(icol.lt.ncutb)then                                      3d12s18
               ip=icol/nhereb                                            3d12s18
              else                                                      3d12s18
               ip=nleftb+(icol-ncutb)/nhere0b                              3d12s18
              end if                                                    3d12s18
              if(ip.lt.0.or.ip.ge.mynprocg)then
               write(6,*)('ip error no. 1: '),ip
               write(6,*)('icol = '),icol,noc(isblk1(1,is1)),
     $              noc(isblk1(2,is1))
               write(6,*)('ncolb: '),ncolb,ncutb,nhereb,nleftb,
     $              nhere0b
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ipass.eq.1)then                                        4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                        4d10s18
              else                                                      4d10s18
               ix=max(i3,i4)                                            4d10s18
               in=min(i3,i4)                                            4d10s18
               ieq=((ix*(ix+1))/2)+in                                   4d10s18
               inot=i3+noc(isblk1(1,is1))*i4                            4d10s18
               i1x=ionex(is1)+(inot-ieq)*iswitch+ieq                    4d10s18
               i10=i1s                                                     4d10s18
               i1n=noc(isblk1(3,is1))                                   4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(i1x)                              4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                 i1x=i1x+nrow1                                           4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           else if(isblk1(4,is1).eq.isblk1(1,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(2,is3).and.                     4d10s18
     $          isblk1(2,is1).eq.isblk1(3,is3))then                     4d10s18
            call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            nrow1=noc(isblk1(1,is1))*noc(isblk1(2,is1))                 4d10s18
            do i4=0,noc(isblk1(2,is1))-1                                4d10s18
             do i3=0,noc(isblk1(1,is1))-1                               4d10s18
              icol=i4+noc(isblk1(2,is1))*i3                             4d10s18
              if(icol.lt.ncutb)then                                      3d12s18
               ip=icol/nhereb                                            3d12s18
              else                                                      3d12s18
               ip=nleftb+(icol-ncutb)/nhere0b                           4d10s18
              end if                                                    3d12s18
              if(ip.lt.0.or.ip.ge.mynprocg)then
               write(6,*)('ip error no. 2: '),ip
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ipass.eq.1)then                                        4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                        4d10s18
              else                                                      4d10s18
               i1x=ionex(is1)+i3+noc(isblk1(1,is1))*i4                  4d13s18
               i10=i1s                                                     4d10s18
               i1n=noc(isblk1(3,is1))                                   4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(i1x)                              4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                 i1x=i1x+nrow1                                           4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           end if                                                       4d10s18
           end if                                                       4d10s18
c     (nv|op)
           if(noc(isblk1(1,is3)).ne.0.and.noc(isblk1(4,is3)).ne.0)then  4d10s18
           if(isblk1(4,is1).eq.isblk1(2,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(1,is3).and.                     4d10s18
     $          isblk1(1,is1).eq.isblk1(3,is3))then                     4d10s18
            call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            if(isblk1(1,is1).eq.isblk1(2,is1))then                      4d10s18
             nrow1=(noc(isblk1(1,is1))*(noc(isblk1(1,is1))+1))/2        4d10s18
             iswitch=0                                                  4d10s18
            else                                                        4d10s18
             nrow1=noc(isblk1(1,is1))*noc(isblk1(2,is1))                4d10s18
             iswitch=1                                                  4d10s18
            end if                                                      4d10s18
            do i4=0,noc(isblk1(2,is1))-1                                4d10s18
             do i3=0,noc(isblk1(1,is1))-1                               4d10s18
              icol=i3+noc(isblk1(1,is1))*i4                             4d10s18
              if(icol.lt.ncutb)then                                      3d12s18
               ip=icol/nhereb                                            3d12s18
              else                                                      3d12s18
               ip=nleftb+(icol-ncutb)/nhere0b                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then                                        4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                        4d10s18
              else                                                      4d10s18
               ix=max(i3,i4)                                            4d10s18
               in=min(i3,i4)                                            4d10s18
               ieq=((ix*(ix+1))/2)+in                                   4d10s18
               inot=i3+noc(isblk1(1,is1))*i4                            4d10s18
               i1x=ionex(is1)+(inot-ieq)*iswitch+ieq                    4d10s18
               i10=i1s                                                     4d10s18
               i1n=noc(isblk1(3,is1))                                   4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(i1x)                              4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                 i1x=i1x+nrow1                                           4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           else if(isblk1(4,is1).eq.isblk1(2,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(1,is3).and.                     4d10s18
     $          isblk1(2,is1).eq.isblk1(3,is3))then                     4d10s18
            call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            nrow1=noc(isblk1(1,is1))*noc(isblk1(2,is1))                 4d10s18
            do i4=0,noc(isblk1(2,is1))-1                                4d10s18
             do i3=0,noc(isblk1(1,is1))-1                               4d10s18
              icol=i4+noc(isblk1(2,is1))*i3                             4d10s18
              if(icol.lt.ncutb)then                                      3d12s18
               ip=icol/nhereb                                            3d12s18
              else                                                      3d12s18
               ip=nleftb+(icol-ncutb)/nhere0b                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then                                        4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                        4d10s18
              else                                                      4d10s18
               i1x=ionex(is1)+i3+noc(isblk1(1,is1))*i4                     4d10s18
               i10=i1s                                                     4d10s18
               i1n=noc(isblk1(3,is1))                                   4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(i1x)                              4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                 i1x=i1x+nrow1                                           4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           end if                                                       4d10s18
           end if                                                       4d10s18
          end do                                                        4d10s18
          if(noc(isblk1(4,is3)).ne.0)then                               4d10s18
          do isj=1,nsdlk                                                4d10s18
c     (vv'|op)
           if(isblk(3,isj).eq.isblk1(1,is3).and.                        4d10s18
     $        isblk(4,isj).eq.isblk1(2,is3).and.                        4d10s18
     $          isblk(1,isj).eq.isblk1(3,is3))then                        4d10s18
            call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            if(isblk(1,isj).eq.isblk(2,isj))then                        4d10s18
             nrowj=(noc(isblk(1,isj))*(noc(isblk(1,isj))+1))/2          4d10s18
             iswitch=0                                                  4d10s18
            else                                                        4d10s18
             nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                  4d10s18
             iswitch=1                                                  4d10s18
            end if                                                      4d10s18
            do i4=0,noc(isblk(2,isj))-1                                 4d10s18
             do i3=0,noc(isblk(1,isj))-1                                4d10s18
              icol=i3+noc(isblk(1,isj))*i4                              4d10s18
              if(icol.lt.ncutb)then                                      3d12s18
               ip=icol/nhereb                                            3d12s18
              else                                                      3d12s18
               ip=nleftb+(icol-ncutb)/nhere0b                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then                                        4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                        4d10s18
              else                                                      4d10s18
               ix=max(i3,i4)                                            4d10s18
               in=min(i3,i4)                                            4d10s18
               ieq=((ix*(ix+1))/2)+in                                   4d10s18
               inot=i3+noc(isblk(1,isj))*i4                             4d10s18
               jm=jmats(isj)+(inot-ieq)*iswitch+ieq                     4d10s18
               i10=i1s                                                     4d10s18
               i1n=nvirtc(isblk(3,isj))                                 4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(jm)                               4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                 jm=jm+nrowj                                             4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           else if(isblk(4,isj).eq.isblk1(1,is3).and.                   4d10s18
     $        isblk(3,isj).eq.isblk1(2,is3).and.                        4d10s18
     $          isblk(1,isj).eq.isblk1(3,is3))then                      4d10s18
            call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                   4d10s18
            do i4=0,noc(isblk(2,isj))-1                                 4d10s18
             do i3=0,noc(isblk(1,isj))-1                                4d10s18
              icol=i3+noc(isblk(1,isj))*i4                              4d10s18
              if(icol.lt.ncutb)then                                      3d12s18
               ip=icol/nhereb                                            3d12s18
              else                                                      3d12s18
               ip=nleftb+(icol-ncutb)/nhere0b                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then                                        4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                        4d10s18
              else                                                      4d10s18
               jm=jmats(isj)+i3+noc(isblk(1,isj))*i4                    4d10s18
               i10=i1s                                                     4d10s18
               i1n=nvirtc(isblk(3,isj))                                 4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(jm)                               4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                 jm=jm+nrowj                                             4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           else if(isblk(4,isj).eq.isblk1(1,is3).and.                   4d10s18
     $           isblk(3,isj).eq.isblk1(2,is3).and.                        4d10s18
     $           isblk(2,isj).eq.isblk1(3,is3))then                      4d10s18
            call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $           mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                   4d10s18
            do i4=0,noc(isblk(2,isj))-1                                 4d10s18
             do i3=0,noc(isblk(1,isj))-1                                4d10s18
              icol=i4+noc(isblk(2,isj))*i3                              4d10s18
              if(icol.lt.ncutb)then                                      4d10s18
               ip=icol/nhereb                                             4d10s18
              else                                                       4d10s18
               ip=nleftb+(icol-ncutb)/nhere0b                               4d10s18
              end if                                                     4d10s18
              if(ipass.eq.1)then                                         4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                          4d10s18
              else                                                        4d10s18
               jm=jmats(isj)+i3+noc(isblk(1,isj))*i4                      4d10s18
               i10=i1s                                                    4d10s18
               i1n=nvirtc(isblk(3,isj))                                   4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(jm)                                 4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                            4d10s18
                 jm=jm+nrowj                                               4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                      4d10s18
             end do                                                      4d10s18
            end do                                                       4d10s18
           else if(isblk(3,isj).eq.isblk1(1,is3).and.                   4d10s18
     $           isblk(4,isj).eq.isblk1(2,is3).and.                     4d13s18
     $           isblk(2,isj).eq.isblk1(3,is3))then                      4d10s18
            call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $           mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            mhere=jh+1-jl                                               4d10s18
            nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                   4d10s18
            do i4=0,noc(isblk(2,isj))-1                                 4d10s18
             do i3=0,noc(isblk(1,isj))-1                                4d10s18
              icol=i4+noc(isblk(2,isj))*i3                              4d10s18
              if(icol.lt.ncutb)then                                      4d10s18
               ip=icol/nhereb                                             4d10s18
              else                                                       4d10s18
               ip=nleftb+(icol-ncutb)/nhere0b                               4d10s18
              end if                                                     4d10s18
              if(ipass.eq.1)then                                         4d10s18
               ibc(nsend+ip)=ibc(nsend+ip)+mhere                          4d10s18
              else                                                        4d10s18
               jm=jmats(isj)+i3+noc(isblk(1,isj))*i4                      4d10s18
               i10=i1s                                                    4d10s18
               i1n=nvirtc(isblk(3,isj))                                   4d10s18
               do i2=i2s,i2e                                               4d10s18
                if(i2.eq.i2e)i1n=i1e                                       4d10s18
                do i1=i10,i1n                                              4d10s18
                 bc(ibc(isend+ip))=bc(jm)                                 4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                            4d10s18
                 jm=jm+nrowj                                               4d10s18
                end do                                                     4d10s18
                i10=1                                                      4d10s18
               end do                                                      4d10s18
              end if                                                      4d10s18
             end do                                                      4d10s18
            end do                                                       4d10s18
           end if                                                       4d10s18
c
c     (nm|op) if distributed (nm|po)
c
           if(isblk(3,isj).eq.isblk1(3,is3).and.                        4d13s18
     $          isblk(4,isj).eq.isblk1(4,is3))then
            iok=1
           else if(isblk(3,isj).eq.isblk1(4,is3).and.
     $           isblk(4,isj).eq.isblk1(3,is3))then                     4d13s18
            if(isblk(1,isj).eq.isblk1(1,is3))then                       4d13s18
             call ilimts(noc(isblk(3,isj)),noc(isblk(4,isj)),           4d13s18
     $           mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
             mhere=jh+1-jl                                               4d10s18
             nrowo=noc(isblk(1,isj))*noc(isblk(2,isj))                  4d13s18
             i10=i1s                                                    4d13s18
             i1n=noc(isblk(3,isj))                                      4d13s18
             i4o=ioooo(isj)                                             4d13s18
             do i2=i2s,i2e                                              4d13s18
              i2m=i2-1                                                  4d13s18
              if(i2.eq.i2e)i1n=i1e                                      4d13s18
              do i1=i10,i1n                                             4d13s18
               i1m=i1-1
               icol=i2m+noc(isblk(4,isj))*i1m                            4d13s18
               if(icol.lt.ncutb)then                                      4d10s18
                ip=icol/nhereb                                             4d10s18
               else                                                       4d10s18
                ip=nleftb+(icol-ncutb)/nhere0b                               4d10s18
               end if                                                     4d10s18
               if(ipass.eq.1)then                                         4d10s18
                ibc(nsend+ip)=ibc(nsend+ip)+nrowo                       4d13s18
               else                                                        4d10s18
                do i12=0,nrowo-1                                        4d13s18
                 bc(ibc(isend+ip))=bc(i4o+i12)                          4d13s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d13s18
                end do                                                  4d13s18
               end if
               i4o=i4o+nrowo
              end do                                                    4d13s18
              i10=1                                                     4d13s18
             end do                                                     4d13s18
            else if(isblk(2,isj).eq.isblk1(1,is3))then                  4d13s18
             call ilimts(noc(isblk(3,isj)),noc(isblk(4,isj)),           4d13s18
     $           mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
             mhere=jh+1-jl                                               4d10s18
             nrowo=noc(isblk(1,isj))*noc(isblk(2,isj))                  4d13s18
             i10=i1s                                                    4d13s18
             i1n=noc(isblk(3,isj))                                      4d13s18
             i4o=ioooo(isj)                                             4d13s18
             do i2=i2s,i2e                                              4d13s18
              i2m=i2-1                                                  4d13s18
              if(i2.eq.i2e)i1n=i1e                                      4d13s18
              do i1=i10,i1n                                             4d13s18
               i1m=i1-1
               icol=i2m+noc(isblk(4,isj))*i1m                            4d13s18
               if(icol.lt.ncutb)then                                      4d10s18
                ip=icol/nhereb                                             4d10s18
               else                                                       4d10s18
                ip=nleftb+(icol-ncutb)/nhere0b                               4d10s18
               end if                                                     4d10s18
               if(ipass.eq.1)then                                         4d10s18
                ibc(nsend+ip)=ibc(nsend+ip)+nrowo                       4d13s18
               else                                                        4d10s18
                do i4=0,noc(isblk(2,isj))-1                                4d13s18
                 do i3=0,noc(isblk(1,isj))-1                            4d13s18
                  i21=i4o+i4+noc(isblk(2,isj))*i3                       4d13s18
                  bc(ibc(isend+ip))=bc(i21)                             4d13s18
                  ibc(isend+ip)=ibc(isend+ip)+1                          4d13s18
                 end do                                                 4d13s18
                end do                                                  4d13s18
               end if                                                   4d13s18
               i4o=i4o+nrowo                                            4d13s18
              end do                                                    4d13s18
              i10=1                                                     4d13s18
             end do                                                     4d13s18
            end if                                                      4d13s18
           end if                                                       4d13s18
          end do                                                        4d10s18
          end if                                                        4d10s18
          do isk=1,nsdlkk                                               4d10s18
           if(noc(isblk1(2,is3)).ne.0)then                              4d10s18
c     (v'm|ov)=K_{mo}^{vv'}=K_{om}^{v'v}
           if(isblkk(1,isk).eq.isblk1(2,is3).and.                       4d10s18
     $        isblkk(2,isk).eq.isblk1(3,is3).and.                       4d10s18
     $          isblkk(3,isk).eq.isblk1(4,is3))then                       4d10s18
            call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            nrowk=noc(isblkk(1,isk))*noc(isblkk(2,isk))                 4d10s18
            i10=i1s                                                     4d10s18
            i1n=nvirtc(isblkk(3,isk))                                   4d10s18
            km=kmats(isk)                                               4d10s18
            do i2=i2s,i2e                                               4d10s18
             if(i2.eq.i2e)i1n=i1e                                       4d10s18
             do i1=i10,i1n                                              4d10s18
              i1m=i1-1                                                  4d10s18
              do i4=0,noc(isblkk(2,isk))-1                              4d10s18
               icol=i4+noc(isblkk(2,isk))*i1m                           4d10s18
               if(icol.lt.ncut)then                                      3d12s18
                ip=icol/nhere                                            3d12s18
               else                                                      3d12s18
                ip=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(ipass.eq.1)then                                       4d10s18
                ibc(nsend+ip)=ibc(nsend+ip)+noc(isblkk(1,isk))          4d10s18
               else                                                     4d10s18
                do i3=0,noc(isblkk(1,isk))-1                            4d10s18
                 iad1=km+i3+noc(isblkk(1,isk))*i4                       4d10s18
                 bc(ibc(isend+ip))=bc(iad1)                             4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                end do                                                  4d10s18
               end if                                                   4d10s18
              end do                                                    4d10s18
              km=km+nrowk                                                4d10s18
             end do                                                     4d10s18
             i10=1                                                      4d10s18
            end do                                                      4d10s18
           else if(isblkk(2,isk).eq.isblk1(2,is3).and.                       4d10s18
     $        isblkk(1,isk).eq.isblk1(3,is3).and.                       4d10s18
     $          isblkk(4,isk).eq.isblk1(4,is3))then                       4d10s18
            call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            nrowk=noc(isblkk(1,isk))*noc(isblkk(2,isk))                 4d10s18
            i10=i1s                                                     4d10s18
            i1n=nvirtc(isblkk(3,isk))                                   4d10s18
            km=kmats(isk)                                               4d10s18
            do i2=i2s,i2e                                               4d10s18
             i2m=i2-1                                                   4d10s18
             if(i2.eq.i2e)i1n=i1e                                       4d10s18
             do i1=i10,i1n                                              4d10s18
              do i3=0,noc(isblkk(1,isk))-1                              4d10s18
               icol=i3+noc(isblkk(1,isk))*i2m                           4d10s18
               if(icol.lt.ncut)then                                      3d12s18
                ip=icol/nhere                                            3d12s18
               else                                                      3d12s18
                ip=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(ipass.eq.1)then                                       4d10s18
                ibc(nsend+ip)=ibc(nsend+ip)+noc(isblkk(2,isk))          4d10s18
               else                                                     4d10s18
                do i4=0,noc(isblkk(2,isk))-1                            4d10s18
                 iad1=km+i3+noc(isblkk(1,isk))*i4                       4d10s18
                 bc(ibc(isend+ip))=bc(iad1)                             4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                end do                                                  4d10s18
               end if                                                   4d10s18
              end do                                                    4d10s18
              km=km+nrowk                                                4d10s18
             end do                                                     4d10s18
             i10=1                                                      4d10s18
            end do                                                      4d10s18
           end if
           end if                                                       4d10s18
c     (nv'|ov)=K_{no}^{vv'}=K_{on}^{v'v}
           if(noc(isblk1(1,is3)).ne.0)then                              4d10s18
           if(isblkk(1,isk).eq.isblk1(1,is3).and.                       4d10s18
     $        isblkk(2,isk).eq.isblk1(3,is3).and.                       4d10s18
     $          isblkk(3,isk).eq.isblk1(4,is3))then                       4d10s18
            call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            nrowk=noc(isblkk(1,isk))*noc(isblkk(2,isk))                 4d10s18
            i10=i1s                                                     4d10s18
            i1n=nvirtc(isblkk(3,isk))                                   4d10s18
            km=kmats(isk)                                               4d10s18
            do i2=i2s,i2e                                               4d10s18
             if(i2.eq.i2e)i1n=i1e                                       4d10s18
             do i1=i10,i1n                                              4d10s18
              i1m=i1-1                                                  4d10s18
              do i4=0,noc(isblkk(2,isk))-1                              4d10s18
               icol=i4+noc(isblkk(2,isk))*i1m                           4d10s18
               if(icol.lt.ncut)then                                      3d12s18
                ip=icol/nhere                                            3d12s18
               else                                                      3d12s18
                ip=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(ipass.eq.1)then                                       4d10s18
                ibc(nsend+ip)=ibc(nsend+ip)+noc(isblkk(1,isk))          4d10s18
               else                                                     4d10s18
                do i3=0,noc(isblkk(1,isk))-1                            4d10s18
                 iad1=km+i3+noc(isblkk(1,isk))*i4                       4d10s18
                 bc(ibc(isend+ip))=bc(iad1)                             4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                end do                                                  4d10s18
               end if                                                   4d10s18
              end do                                                    4d10s18
              km=km+nrowk                                                4d10s18
             end do                                                     4d10s18
             i10=1                                                      4d10s18
            end do                                                      4d10s18
           else if(isblkk(2,isk).eq.isblk1(1,is3).and.                  4d10s18
     $        isblkk(1,isk).eq.isblk1(3,is3).and.                       4d10s18
     $          isblkk(4,isk).eq.isblk1(4,is3))then                       4d10s18
            call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $          mynprocg,mynowprog,jl,jh,i1s,i1e,i2s,i2e)               4d10s18
            nrowk=noc(isblkk(1,isk))*noc(isblkk(2,isk))                 4d10s18
            i10=i1s                                                     4d10s18
            i1n=nvirtc(isblkk(3,isk))                                   4d10s18
            km=kmats(isk)                                               4d10s18
            do i2=i2s,i2e                                               4d10s18
             i2m=i2-1                                                   4d10s18
             if(i2.eq.i2e)i1n=i1e                                       4d10s18
             do i1=i10,i1n                                              4d10s18
              do i3=0,noc(isblkk(1,isk))-1                              4d10s18
               icol=i3+noc(isblkk(1,isk))*i2m                           4d10s18
               if(icol.lt.ncut)then                                      3d12s18
                ip=icol/nhere                                            3d12s18
               else                                                      3d12s18
                ip=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(ipass.eq.1)then                                       4d10s18
                ibc(nsend+ip)=ibc(nsend+ip)+noc(isblkk(2,isk))          4d10s18
               else                                                     4d10s18
                do i4=0,noc(isblkk(2,isk))-1                            4d10s18
                 iad1=km+i3+noc(isblkk(1,isk))*i4                       4d10s18
                 bc(ibc(isend+ip))=bc(iad1)                             4d10s18
                 ibc(isend+ip)=ibc(isend+ip)+1                          4d10s18
                end do                                                  4d10s18
               end if                                                   4d10s18
              end do                                                    4d10s18
              km=km+nrowk                                               4d13s18
             end do                                                     4d10s18
             i10=1                                                      4d10s18
            end do                                                      4d10s18
           end if
           end if                                                       4d10s18
          end do                                                        4d10s18
         end do                                                         4d10s18
        end if                                                          4d10s18
        do is=1,nsdlkk                                                  2d22s18
         call ilimts(nvirtc(isblkk(3,is)),noc(isblkk(4,is)),mynprocg,   2d22s18
     $       mynowprog,il1,ih1,i1s,i1e,i2s,i2e)                           2d22s18
         nhere1=ih1+1-il1                                                 2d22s18
         call ilimts(noc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,   2d22s18
     $       mynowprog,il2,ih2,i1s,i1e,i2s,i2e)                           2d22s18
         nhere2=ih2+1-il2                                                 2d22s18
         call ilimts(noc(isblkk(3,is)),noc(isblkk(4,is)),mynprocg,      2d22s18
     $       mynowprog,il3,ih3,i1s,i1e,i2s,i2e)                           2d22s18
         nhere3=ih3+1-il3                                                 2d22s18
         do is1x=1,nsdlk1                                               2d22s18
          if(isblk1(1,is1x).eq.isblkk(1,is)
     $         .and.isblk1(3,is1x).eq.isblkk(2,is)                      3d12s18
     $         .and.isblk1(2,is1x).eq.isblkk(4,is).and.                 3d13s18
     $         noc(isblk1(2,is1x)).gt.0)then                            3d13s18
           call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               2d22s18
           ncol=noc(isblk1(2,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(1,nhere0)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           if(isblk1(1,is1x).eq.isblk1(2,is1x))then                     2d22s18
            nrowx=(noc(isblk1(1,is1x))*(noc(isblk1(1,is1x))+1))/2       3d14s18
            iswitch=0                                                   2d22s18
           else                                                         2d22s18
            nrowx=noc(isblk1(1,is1x))*noc(isblk1(2,is1x))               2d22s18
            iswitch=1                                                   2d22s18
           end if                                                       2d22s18
           i10=i1s                                                      2d22s18
           i1n=noc(isblk1(3,is1x))
           i1x=ionex(is1x)                                              2d22s18
           do i2=i2s,i2e                                                2d22s18
            if(i2.eq.i2e)i1n=i1e                                        2d22s18
            i2m=i2-1
            do i1=i10,i1n                                               2d22s18
             i1m=i1-1
             do io=0,noc(isblk1(2,is1x))-1                              3d13s18
              icol=io+i2m*noc(isblk1(2,is1x))                           3d12s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ip.lt.0.or.ip.ge.mynprocg)then
               write(6,*)('error computing ip: '),ip
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk1(1,is1x))
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk1(1,is1x))
              else                                                      2d23s18
               do in=0,noc(isblk1(1,is1x))-1
                inot=in+noc(isblk1(1,is1x))*io                          3d12s18
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                           3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                iad=i1x+(inot-ieq)*iswitch+ieq                          3d12s18
                bc(ibc(isend+ip))=bc(iad)                               3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i1x=i1x+nrowx
            end do
            i10=1
           end do                                                       2d22s18
          else if(isblk1(2,is1x).eq.isblkk(1,is)
     $         .and.isblk1(1,is1x).eq.isblkk(4,is)
     $         .and.isblk1(3,is1x).eq.isblkk(2,is).and.                 3d13s18
     $          noc(isblk1(1,is1x)).gt.0)then                           3d13s18
           call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               2d22s18
           ncol=noc(isblk1(1,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(1,nhere0)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           if(isblk1(1,is1x).eq.isblk1(2,is1x))then                     3d12s18
            nrowx=(noc(isblk1(1,is1x))*(noc(isblk1(2,is1x))+1))/2       3d12s18
            iswitch=0                                                   3d12s18
           else                                                         3d12s18
            nrowx=noc(isblk1(1,is1x))*noc(isblk1(2,is1x))                3d12s18
            iswitch=1
           end if                                                       3d12s18
           i10=i1s                                                      2d22s18
           i1n=noc(isblk1(3,is1x))
           i1x=ionex(is1x)                                              2d22s18
           do i2=i2s,i2e                                                2d22s18
            if(i2.eq.i2e)i1n=i1e                                        2d22s18
            i2m=i2-1
            do i1=i10,i1n                                               3d14s18
             i1m=i1-1
             do io=0,noc(isblk1(1,is1x))-1                              3d13s18
              icol=io+i2m*noc(isblk1(1,is1x))                           3d12s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ip.lt.0.or.ip.ge.mynprocg)then
               write(6,*)('error computing ip: '),ip
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk1(2,is1x))          3d12s18
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk1(2,is1x))            3d12s18
              else                                                      2d23s18
               do in=0,noc(isblk1(2,is1x))-1
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                          3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                inot=io+noc(isblk1(1,is1x))*in                          3d12s18
                iad1=i1x+(inot-ieq)*iswitch+ieq                         3d12s18
                bc(ibc(isend+ip))=bc(iad1)                              3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i1x=i1x+nrowx
            end do
            i10=1
           end do                                                       2d22s18
          end if
          if(isblk1(1,is1x).eq.isblkk(2,is)                             3d13s18
     $         .and.isblk1(3,is1x).eq.isblkk(1,is)                      3d13s18
     $         .and.isblk1(2,is1x).eq.isblkk(3,is).and.                 3d13s18
     $         noc(isblk1(2,is1x)).gt.0)then                            3d13s18
           call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               2d22s18
           ncol=noc(isblk1(2,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(nhere0,1)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           if(isblk1(1,is1x).eq.isblk1(2,is1x))then                     2d22s18
            nrowx=(noc(isblk1(1,is1x))*(noc(isblk1(1,is1x))+1))/2       3d14s18
            iswitch=0                                                   2d22s18
           else                                                         2d22s18
            nrowx=noc(isblk1(1,is1x))*noc(isblk1(2,is1x))               2d22s18
            iswitch=1                                                   2d22s18
           end if                                                       2d22s18
           i10=i1s                                                      2d22s18
           i1n=noc(isblk1(3,is1x))
           i1x=ionex(is1x)                                              2d22s18
           do i2=i2s,i2e                                                2d22s18
            if(i2.eq.i2e)i1n=i1e                                        2d22s18
            i2m=i2-1
            do i1=i10,i1n                                               3d14s18
             i1m=i1-1
             do io=0,noc(isblk1(2,is1x))-1                              3d13s18
              icol=io+i2m*noc(isblk1(2,is1x))                           3d12s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ip.lt.0.or.ip.ge.mynprocg)then
               write(6,*)('error computing ip: '),ip
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk1(1,is1x))
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk1(1,is1x))
              else                                                      2d23s18
               do in=0,noc(isblk1(1,is1x))-1
                inot=in+noc(isblk1(1,is1x))*io                          3d12s18
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                           3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                iad=i1x+(inot-ieq)*iswitch+ieq                          3d12s18
                bc(ibc(isend+ip))=bc(iad)                               3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i1x=i1x+nrowx
            end do
            i10=1                                                       3d14s18
           end do                                                       2d22s18
          else if(isblk1(2,is1x).eq.isblkk(2,is)                        3d13s18
     $         .and.isblk1(1,is1x).eq.isblkk(3,is)                      3d13s18
     $         .and.isblk1(3,is1x).eq.isblkk(1,is).and.                 3d13s18
     $          noc(isblk1(1,is1x)).gt.0)then                           3d13s18
           call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               2d22s18
           ncol=noc(isblk1(1,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(1,nhere0)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           if(isblk1(1,is1x).eq.isblk1(2,is1x))then                     3d12s18
            nrowx=(noc(isblk1(1,is1x))*(noc(isblk1(2,is1x))+1))/2       3d12s18
            iswitch=0                                                   3d12s18
           else                                                         3d12s18
            nrowx=noc(isblk1(1,is1x))*noc(isblk1(2,is1x))                3d12s18
            iswitch=1
           end if                                                       3d12s18
           i10=i1s                                                      2d22s18
           i1n=noc(isblk1(3,is1x))
           i1x=ionex(is1x)                                              2d22s18
           do i2=i2s,i2e                                                2d22s18
            if(i2.eq.i2e)i1n=i1e                                        2d22s18
            i2m=i2-1
            do i1=i10,i1n                                               3d14s18
             i1m=i1-1
             do io=0,noc(isblk1(1,is1x))-1                              3d13s18
              icol=io+i2m*noc(isblk1(1,is1x))                           3d12s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ip.lt.0.or.ip.ge.mynprocg)then
               write(6,*)('error computing ip: '),ip
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk1(2,is1x))          3d12s18
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk1(2,is1x))            3d12s18
              else                                                      2d23s18
               do in=0,noc(isblk1(2,is1x))-1
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                          3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                inot=io+noc(isblk1(1,is1x))*in                          3d12s18
                iad1=i1x+(inot-ieq)*iswitch+ieq                         3d12s18
                bc(ibc(isend+ip))=bc(iad1)                              3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i1x=i1x+nrowx
            end do
            i10=1                                                       3d14s18
           end do                                                       2d22s18
          end if
         end do                                                         2d22s18
         do iso=1,nsdlk                                                 3d14s18
          if(isblk(1,iso).eq.isblkk(1,is).and.                          3d15s18
     $       isblk(2,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(3,iso).eq.isblkk(2,is).and.noc(isblk(2,iso)).gt.0    3d15s18
     $       .and.noc(isblk(4,iso)).gt.0)then                           3d14s18
           call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
           ncol=noc(isblk(4,iso))*noc(isblk(2,iso))                     3d15s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(1,nhere0)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           i10=i1s
           i1n=noc(isblk(3,iso))
           if(isblk(1,iso).eq.isblk(2,iso))then
            nrow=(noc(isblk(1,iso))*(noc(isblk(1,iso))+1))/2
            iswitch=0
           else
            nrow=noc(isblk(1,iso))*noc(isblk(2,iso))
            iswitch=1
           end if
           i4o=ioooo(iso)
           do i2=i2s,i2e
            if(i2.eq.i2e)i1n=i1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             do io=0,noc(isblk(2,iso))-1
              icol=i2m+noc(isblk(4,iso))*io                             3d15s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk(1,iso))            3d14s18
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk(1,iso))              3d14s18
              else                                                      2d23s18
               do in=0,noc(isblk(1,iso))-1                              3d14s18
                inot=in+noc(isblk(1,iso))*io                            3d14s18
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                           3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                iad=i4o+(inot-ieq)*iswitch+ieq                          3d14s18
                bc(ibc(isend+ip))=bc(iad)                               3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i4o=i4o+nrow
            end do
            i10=1
           end do
          else if(isblk(2,iso).eq.isblkk(1,is).and.                     3d15s18
     $       isblk(1,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(3,iso).eq.isblkk(2,is).and.noc(isblk(1,iso)).gt.0    3d15s18
     $       .and.noc(isblk(4,iso)).gt.0)then                           3d14s18
           call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
           ncol=noc(isblk(1,iso))*noc(isblk(4,iso))                     3d14s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(1,nhere0)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           i10=i1s
           i1n=noc(isblk(3,iso))
           if(isblk(1,iso).eq.isblk(2,iso))then
            nrow=(noc(isblk(1,iso))*(noc(isblk(1,iso))+1))/2
            iswitch=0
           else
            nrow=noc(isblk(1,iso))*noc(isblk(2,iso))
            iswitch=1
           end if
           i4o=ioooo(iso)
           do i2=i2s,i2e
            if(i2.eq.i2e)i1n=i1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             do io=0,noc(isblk(1,iso))-1
              icol=i2m+noc(isblk(4,iso))*io                             3d15s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk(2,iso))
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk(2,iso))
              else                                                      2d23s18
               do in=0,noc(isblk(2,iso))-1                              3d14s18
                inot=io+noc(isblk(1,iso))*in                            3d14s18
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                           3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                iad=i4o+(inot-ieq)*iswitch+ieq                          3d14s18
                bc(ibc(isend+ip))=bc(iad)                               3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i4o=i4o+nrow
            end do
            i10=1
           end do
          else if(isblk(2,iso).eq.isblkk(1,is).and.                     3d15s18
     $       isblk(1,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(4,iso).eq.isblkk(2,is).and.noc(isblk(1,iso)).gt.0    3d15s18
     $       .and.noc(isblk(3,iso)).gt.0)then                           3d14s18
           call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
           ncol=noc(isblk(1,iso))*noc(isblk(3,iso))                     3d15s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(1,nhere0)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           i10=i1s
           i1n=noc(isblk(3,iso))
           if(isblk(1,iso).eq.isblk(2,iso))then
            nrow=(noc(isblk(1,iso))*(noc(isblk(1,iso))+1))/2
            iswitch=0
           else
            nrow=noc(isblk(1,iso))*noc(isblk(2,iso))
            iswitch=1
           end if
           i4o=ioooo(iso)
           do i2=i2s,i2e
            if(i2.eq.i2e)i1n=i1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             do io=0,noc(isblk(1,iso))-1
              icol=i1m+noc(isblk(3,iso))*io                             3d15s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk(2,iso))
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk(2,iso))
              else                                                      2d23s18
               do in=0,noc(isblk(2,iso))-1                              3d14s18
                inot=io+noc(isblk(1,iso))*in                            3d14s18
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                           3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                iad=i4o+(inot-ieq)*iswitch+ieq                          3d14s18
                bc(ibc(isend+ip))=bc(iad)                               3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i4o=i4o+nrow
            end do
            i10=1
           end do
          else if(isblk(1,iso).eq.isblkk(1,is).and.                     3d15s18
     $       isblk(2,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(4,iso).eq.isblkk(2,is).and.noc(isblk(2,iso)).gt.0    3d15s18
     $       .and.noc(isblk(3,iso)).gt.0)then                           3d14s18
           call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
           ncol=noc(isblk(2,iso))*noc(isblk(3,iso))                     3d15s18
           nhere0=ncol/mynprocg                                         3d12s18
           nleft=ncol-nhere0*mynprocg                                   3d12s18
           nhere=nhere0+1                                               3d12s18
           nhere0=max(1,nhere0)                                         3d14s18
           ncut=nleft*nhere                                             3d12s18
           i10=i1s
           i1n=noc(isblk(3,iso))
           if(isblk(1,iso).eq.isblk(2,iso))then
            nrow=(noc(isblk(1,iso))*(noc(isblk(1,iso))+1))/2
            iswitch=0
           else
            nrow=noc(isblk(1,iso))*noc(isblk(2,iso))
            iswitch=1
           end if
           i4o=ioooo(iso)
           do i2=i2s,i2e
            if(i2.eq.i2e)i1n=i1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             do io=0,noc(isblk(2,iso))-1
              icol=i1m+noc(isblk(3,iso))*io                             3d15s18
              if(icol.lt.ncut)then                                      3d12s18
               ip=icol/nhere                                            3d12s18
              else                                                      3d12s18
               ip=nleft+(icol-ncut)/nhere0                              3d12s18
              end if                                                    3d12s18
              if(ipass.eq.1)then
               ibc(nsend+ip)=ibc(nsend+ip)+noc(isblk(1,iso))
               ios(ip+1,is)=ios(ip+1,is)+noc(isblk(1,iso))
              else                                                      2d23s18
               do in=0,noc(isblk(1,iso))-1                              3d14s18
                inot=in+noc(isblk(1,iso))*io                            3d14s18
                ix=max(in,io)                                           3d12s18
                inx=min(in,io)                                           3d12s18
                ieq=((ix*(ix+1))/2)+inx                                 3d12s18
                iad=i4o+(inot-ieq)*iswitch+ieq                          3d14s18
                bc(ibc(isend+ip))=bc(iad)                               3d12s18
                ibc(isend+ip)=ibc(isend+ip)+1                           2d23s18
               end do
              end if
             end do
             i4o=i4o+nrow
            end do
            i10=1
           end do
          end if
         end do                                                         3d14s18
        end do                                                          2d22s18
        if(ipass.eq.2)then                                               1d10s18
         jbufs=0                                                         1d10s18
         jbufr=0
         do ip=0,mynprocg-1                                              1d10s18
          iad1=jbufs+ibufs
          ibc(isend+ip)=jbufs                                            1d10s18
          jbufs=jbufs+ibc(nsend+ip)                                      1d10s18
          ibc(irecv+ip)=jbufr                                           3d13s18
          jbufr=jbufr+ibc(nrecv+ip)                                     3d13s18
          n4=ibc(nsend+ip)
         end do                                                          1d10s18
         call dws_all2allvb8(bc(ibufs),ibc(nsend),ibc(isend),bc(ibufr),   1d10s18
     $       ibc(nrecv),ibc(irecv))                                     1d10s18
         jbufr=ibufr                                                     1d10s18
         do ip=0,mynprocg-1                                              1d10s18
          n4=ibc(nrecv+ip)
          ibc(irecv+ip)=jbufr                                            1d10s18
          jbufr=jbufr+ibc(nrecv+ip)                                      1d10s18
         end do                                                          1d10s18
        end if                                                           1d10s18
c
c     for receiving
c
        me2mea=0
        me2meb=0
        me2mec=0
        me2med=0
        do ip=0,mynprocg-1                                              2d23s18
         if(c3.ne.0d0)then                                              4d10s18
          do is3=1,nsdlk1                                                4d10s18
           call ilimts(noc(isblk1(3,is3)),nvirtc(isblk1(4,is3)),         4d10s18
     $         mynprocg,mynowprog,il3,ih3,i1s,i1e,i2s,i2e)              4d10s18
           ncol=noc(isblk1(3,is3))*nvirtc(isblk1(4,is3))                 4d10s18
           nhere0=ncol/mynprocg                                          4d10s18
           nleft=ncol-nhere0*mynprocg                                    4d10s18
           nhere=nhere0+1                                                4d10s18
           nhere0=max(1,nhere0)                                          4d10s18
           ncut=nleft*nhere                                              4d10s18
           il3=il3-1                                                    4d10s18
           nhere3m=ih3-il3                                               4d10s18
           call ilimts(noc(isblk1(3,is3)),noc(isblk1(4,is3)),            4d10s18
     $         mynprocg,mynowprog,il0,ih0,i1s,i1e,i2s,i2e)              4d10s18
           il0=il0-1                                                    4d10s18
           nhere0m=ih0-il0                                               4d10s18
           ncolb=noc(isblk1(3,is3))*noc(isblk1(4,is3))                    4d10s18
           nhere0b=ncolb/mynprocg                                          4d10s18
           nleftb=ncolb-nhere0b*mynprocg                                    4d10s18
           nhereb=nhere0b+1                                                4d10s18
           nhere0b=max(1,nhere0b)                                          4d10s18
           ncutb=nleftb*nhereb                                              4d10s18
           do is1=1,nsdlk1                                               4d10s18
            if(noc(isblk1(4,is3)).ne.0.and.noc(isblk1(2,is3)).ne.0)then  4d10s18
c     (vm|op)
             if(isblk1(4,is1).eq.isblk1(1,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(2,is3).and.                     4d10s18
     $          isblk1(1,is1).eq.isblk1(3,is3))then                     4d10s18
              call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk1(2,is1))-1                                4d10s18
               do i3=0,noc(isblk1(1,is1))-1                               4d10s18
                icol=i3+noc(isblk1(1,is1))*i4                             4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                            3d12s18
                else                                                      3d12s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                              3d12s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                                4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                        4d10s18
                 else                                                      4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=noc(isblk1(3,is1))                                   4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                             4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                            4d10s18
                    iad=i3xa(1,is3)+icol-il0                            4d10s18
     $                   +nhere0m*(i2m+nvirtc(isblk1(1,is3))*i1m)       4d12s18
                    bc(iad)=bc(ibc(irecv+ip))                            4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                        4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if                                                  4d10s18
               end do                                                     4d10s18
              end do                                                      4d10s18
             else if(isblk1(4,is1).eq.isblk1(1,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(2,is3).and.                     4d10s18
     $          isblk1(2,is1).eq.isblk1(3,is3))then                     4d10s18
              call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk1(2,is1))-1                                4d10s18
               do i3=0,noc(isblk1(1,is1))-1                               4d10s18
                icol=i4+noc(isblk1(2,is1))*i3                             4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                            3d12s18
                else                                                      3d12s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                              3d12s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                                4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                        4d10s18
                 else                                                      4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=noc(isblk1(3,is1))                                   4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                             4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                            4d10s18
                    iad=i3xa(1,is3)+icol-il0                            4d10s18
     $                   +nhere0m*(i2m+nvirtc(isblk1(1,is3))*i1m)       4d10s18
                    bc(iad)=bc(ibc(irecv+ip))                            4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                        4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if                                                  4d10s18
               end do                                                     4d10s18
              end do                                                      4d10s18
             end if
            end if                                                       4d10s18
c     (nv|op)
            if(noc(isblk1(1,is3)).ne.0.and.noc(isblk1(4,is3)).ne.0)then  4d10s18
             if(isblk1(4,is1).eq.isblk1(2,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(1,is3).and.                     4d10s18
     $          isblk1(1,is1).eq.isblk1(3,is3))then                     4d10s18
              call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk1(2,is1))-1                                4d10s18
               do i3=0,noc(isblk1(1,is1))-1                               4d10s18
                icol=i3+noc(isblk1(1,is1))*i4                             4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                            3d12s18
                else                                                      3d12s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                      4d10s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                               4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                        4d10s18
                 else                                                   4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=noc(isblk1(3,is1))                                   4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                              4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                             4d10s18
                    iad=i3xa(2,is3)+icol-il0                             4d10s18
     $                   +nhere0m*(i1m+noc(isblk1(1,is3))*i2m)           4d10s18
                    bc(iad)=bc(ibc(irecv+ip))                            4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                          4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if                                                  4d10s18
               end do                                                     4d10s18
              end do                                                      4d10s18
             else if(isblk1(4,is1).eq.isblk1(2,is3).and.                       4d10s18
     $          isblk1(3,is1).eq.isblk1(1,is3).and.                     4d10s18
     $          isblk1(2,is1).eq.isblk1(3,is3))then                     4d10s18
              call ilimts(noc(isblk1(3,is1)),nvirtc(isblk1(4,is1)),       4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk1(2,is1))-1                                4d10s18
               do i3=0,noc(isblk1(1,is1))-1                               4d10s18
                icol=i4+noc(isblk1(2,is1))*i3                             4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                      4d10s18
                else                                                    4d10s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                      4d10s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                              4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                        4d10s18
                 else                                                      4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=noc(isblk1(3,is1))                                   4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                             4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                            4d10s18
                    iad=i3xa(2,is3)+icol-il0                            4d10s18
     $                   +nhere0m*(i1m+noc(isblk1(1,is3))*i2m)          4d10s18
                    bc(iad)=bc(ibc(irecv+ip))                           4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                       4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if
               end do                                                     4d10s18
              end do                                                      4d10s18
             end if                                                       4d10s18
            end if
           end do                                                       4d10s18
           if(noc(isblk1(4,is3)).ne.0)then                              4d10s18
            do isj=1,nsdlk                                               4d10s18
c     (vv'|op)
             if(isblk(3,isj).eq.isblk1(1,is3).and.                        4d10s18
     $        isblk(4,isj).eq.isblk1(2,is3).and.                        4d10s18
     $          isblk(1,isj).eq.isblk1(3,is3))then                        4d10s18
              call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk(2,isj))-1                                 4d10s18
               do i3=0,noc(isblk(1,isj))-1                                4d10s18
                icol=i3+noc(isblk(1,isj))*i4                              4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                      4d10s18
                else                                                      3d12s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                      4d10s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                              4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                     4d10s18
                 else                                                      4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=nvirtc(isblk(3,isj))                                 4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                             4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                            4d10s18
                    iad=i3xa(3,is3)+icol-il0                            4d10s18
     $                   +nhere0m*(i1m+nvirtc(isblk1(1,is3))*i2m)       4d10s18
                    bc(iad)=bc(ibc(irecv+ip))                           4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                          4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if
               end do                                                     4d10s18
              end do                                                      4d10s18
             else if(isblk(4,isj).eq.isblk1(1,is3).and.                   4d10s18
     $        isblk(3,isj).eq.isblk1(2,is3).and.                        4d10s18
     $          isblk(1,isj).eq.isblk1(3,is3))then                      4d10s18
              call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk(2,isj))-1                                 4d10s18
               do i3=0,noc(isblk(1,isj))-1                                4d10s18
                icol=i3+noc(isblk(1,isj))*i4                              4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                      4d10s18
                else                                                      3d12s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                      4d10s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                              4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                     4d10s18
                 else                                                      4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=nvirtc(isblk(3,isj))                                 4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                             4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                            4d10s18
                    iad=i3xa(3,is3)+icol-il0                            4d10s18
     $                   +nhere0m*(i2m+nvirtc(isblk1(1,is3))*i1m)       4d10s18
                    bc(iad)=bc(ibc(irecv+ip))                           4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                          4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if
               end do                                                     4d10s18
              end do                                                      4d10s18
             else if(isblk(4,isj).eq.isblk1(1,is3).and.                   4d10s18
     $           isblk(3,isj).eq.isblk1(2,is3).and.                        4d10s18
     $           isblk(2,isj).eq.isblk1(3,is3))then                      4d10s18
              call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk(2,isj))-1                                 4d10s18
               do i3=0,noc(isblk(1,isj))-1                                4d10s18
                icol=i4+noc(isblk(2,isj))*i3                              4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                      4d10s18
                else                                                      3d12s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                      4d10s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                              4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                     4d10s18
                 else                                                      4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=nvirtc(isblk(3,isj))                                 4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                             4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                            4d10s18
                    iad=i3xa(3,is3)+icol-il0                            4d10s18
     $                   +nhere0m*(i2m+nvirtc(isblk1(1,is3))*i1m)       4d10s18
                    bc(iad)=bc(ibc(irecv+ip))                           4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                          4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if
               end do                                                     4d10s18
              end do                                                      4d10s18
             else if(isblk(3,isj).eq.isblk1(1,is3).and.                   4d10s18
     $           isblk(4,isj).eq.isblk1(2,is3).and.                     4d13s18
     $           isblk(2,isj).eq.isblk1(3,is3))then                      4d10s18
              call ilimts(nvirtc(isblk(3,isj)),nvirtc(isblk(4,isj)),      4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              mhere=jh+1-jl                                               4d10s18
              do i4=0,noc(isblk(2,isj))-1                                 4d10s18
               do i3=0,noc(isblk(1,isj))-1                                4d10s18
                icol=i4+noc(isblk(2,isj))*i3                              4d10s18
                if(icol.lt.ncutb)then                                      3d12s18
                 iptry=icol/nhereb                                      4d10s18
                else                                                      3d12s18
                 iptry=nleftb+(icol-ncutb)/nhere0b                      4d10s18
                end if                                                    3d12s18
                if(iptry.eq.mynowprog)then                              4d10s18
                 if(ipass.eq.1)then                                        4d10s18
                  ibc(nrecv+ip)=ibc(nrecv+ip)+mhere                     4d10s18
                 else                                                      4d10s18
                  i10=i1s                                                     4d10s18
                  i1n=nvirtc(isblk(3,isj))                                 4d10s18
                  do i2=i2s,i2e                                               4d10s18
                   i2m=i2-1                                             4d10s18
                   if(i2.eq.i2e)i1n=i1e                                       4d10s18
                   do i1=i10,i1n                                              4d10s18
                    i1m=i1-1                                            4d10s18
                    iad=i3xa(3,is3)+icol-il0                            4d10s18
     $                   +nhere0m*(i1m+nvirtc(isblk1(1,is3))*i2m)       4d10s18
                    bc(iad)=bc(ibc(irecv+ip))                           4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                          4d10s18
                   end do                                                     4d10s18
                   i10=1                                                      4d10s18
                  end do                                                      4d10s18
                 end if                                                    4d10s18
                end if
               end do                                                     4d10s18
              end do                                                      4d10s18
             end if
c
c     (nm|op) if distributed (nm|po)
c
             if(isblk(3,isj).eq.isblk1(3,is3).and.                      4d13s18
     $          isblk(4,isj).eq.isblk1(4,is3))then
              iok=1
             else if(isblk(3,isj).eq.isblk1(4,is3).and.
     $           isblk(4,isj).eq.isblk1(3,is3))then                     4d13s18
              if(isblk(1,isj).eq.isblk1(1,is3).or.                      4d13s18
     $           isblk(2,isj).eq.isblk1(1,is3))then                     4d13s18
               call ilimts(noc(isblk(3,isj)),noc(isblk(4,isj)),           4d13s18
     $           mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                     4d13s18
               mhere=jh+1-jl                                               4d10s18
               nrowo=noc(isblk(1,isj))*noc(isblk(2,isj))                  4d13s18
               i10=i1s                                                    4d13s18
               i1n=noc(isblk(3,isj))                                      4d13s18
               do i2=i2s,i2e                                              4d13s18
                i2m=i2-1                                                  4d13s18
                if(i2.eq.i2e)i1n=i1e                                      4d13s18
                do i1=i10,i1n                                             4d13s18
                 i1m=i1-1
                 icol=i2m+noc(isblk(4,isj))*i1m                            4d13s18
                 if(icol.lt.ncutb)then                                      4d10s18
                  iptry=icol/nhereb                                             4d10s18
                 else                                                       4d10s18
                  iptry=nleftb+(icol-ncutb)/nhere0b                               4d10s18
                 end if                                                     4d10s18
                 if(iptry.eq.mynowprog)then                             4d13s18
                  if(ipass.eq.1)then                                         4d10s18
                   ibc(nrecv+ip)=ibc(nrecv+ip)+nrowo                    4d13s18
                  else                                                        4d10s18
                   do i12=0,nrowo-1                                        4d13s18
                    iad=i3xa(4,is3)+icol-il0+nhere0m*i12                4d13s18
                    bc(iad)=bc(ibc(irecv+ip))                           4d13s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                          4d13s18
                   end do                                                  4d13s18
                  end if
                 end if                                                 4d13s18
                end do                                                    4d13s18
                i10=1                                                     4d13s18
               end do                                                     4d13s18
              end if                                                      4d13s18
             end if                                                       4d13s18
            end do                                                       4d10s18
           end if                                                       4d10s18
           do isk=1,nsdlkk                                              4d10s18
            if(noc(isblk1(2,is3)).ne.0)then                              4d10s18
c     (v'm|ov)=K_{mo}^{vv'}=K_{om}^{v'v}
             if(isblkk(1,isk).eq.isblk1(2,is3).and.                       4d10s18
     $        isblkk(2,isk).eq.isblk1(3,is3).and.                       4d10s18
     $          isblkk(3,isk).eq.isblk1(4,is3))then                       4d10s18
              call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              i10=i1s                                                     4d10s18
              i1n=nvirtc(isblkk(3,isk))                                   4d10s18
              do i2=i2s,i2e                                               4d10s18
               i2m=i2-1                                                 4d10s18
               if(i2.eq.i2e)i1n=i1e                                       4d10s18
               do i1=i10,i1n                                              4d10s18
                i1m=i1-1                                                  4d10s18
                do i4=0,noc(isblkk(2,isk))-1                              4d10s18
                 icol=i4+noc(isblkk(2,isk))*i1m                           4d10s18
                 if(icol.lt.ncut)then                                      3d12s18
                  iptry=icol/nhere                                            3d12s18
                 else                                                      3d12s18
                  iptry=nleft+(icol-ncut)/nhere0                              3d12s18
                 end if                                                    3d12s18
                 if(iptry.eq.mynowprog)then                             4d10s18
                  if(ipass.eq.1)then                                       4d10s18
                   ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblkk(1,isk))          4d10s18
                  else                                                     4d10s18
                   do i3=0,noc(isblkk(1,isk))-1                            4d10s18
                    iad1=i3xb(1,is3)+icol-il3                           4d10s18
     $                   +nhere3m*(i2m+nvirtc(isblk1(1,is3))*i3)        4d10s18
                    bc(iad1)=bc(ibc(irecv+ip))                          4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                       4d10s18
                   end do                                                  4d10s18
                  end if                                                   4d10s18
                 end if                                                 4d10s18
                end do                                                    4d10s18
               end do                                                     4d10s18
               i10=1                                                      4d10s18
              end do                                                      4d10s18
             else if(isblkk(2,isk).eq.isblk1(2,is3).and.                       4d10s18
     $        isblkk(1,isk).eq.isblk1(3,is3).and.                       4d10s18
     $          isblkk(4,isk).eq.isblk1(4,is3))then                       4d10s18
              call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $          mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                      4d10s18
              i10=i1s                                                     4d10s18
              i1n=nvirtc(isblkk(3,isk))                                   4d10s18
              do i2=i2s,i2e                                               4d10s18
               if(i2.eq.i2e)i1n=i1e                                       4d10s18
               i2m=i2-1                                                 4d10s18
               do i1=i10,i1n                                              4d10s18
                i1m=i1-1                                                  4d10s18
                do i3=0,noc(isblkk(1,isk))-1                              4d10s18
                 icol=i3+noc(isblkk(1,isk))*i2m                           4d10s18
                 if(icol.lt.ncut)then                                      3d12s18
                  iptry=icol/nhere                                            3d12s18
                 else                                                      3d12s18
                  iptry=nleft+(icol-ncut)/nhere0                              3d12s18
                 end if                                                    3d12s18
                 if(iptry.eq.mynowprog)then                             4d10s18
                  if(ipass.eq.1)then                                       4d10s18
                   ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblkk(2,isk))          4d10s18
                  else                                                     4d10s18
                   do i4=0,noc(isblkk(2,isk))-1                            4d10s18
                    iad1=i3xb(1,is3)+icol-il3                           4d10s18
     $                   +nhere3m*(i1m+nvirtc(isblk1(1,is3))*i4)        4d10s18
                    bc(iad1)=bc(ibc(irecv+ip))                          4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                       4d10s18
                   end do                                                  4d10s18
                  end if                                                   4d10s18
                 end if                                                 4d10s18
                end do                                                    4d10s18
               end do                                                     4d10s18
               i10=1                                                      4d10s18
              end do                                                      4d10s18
             end if
            end if
c     (nv'|ov)=K_{no}^{vv'}=K_{on}^{v'v}
            if(noc(isblk1(1,is3)).ne.0)then                              4d10s18
             if(isblkk(1,isk).eq.isblk1(1,is3).and.                       4d10s18
     $           isblkk(2,isk).eq.isblk1(3,is3).and.                       4d10s18
     $           isblkk(3,isk).eq.isblk1(4,is3))then                       4d10s18
              call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $           mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                     4d10s18
              i10=i1s                                                     4d10s18
              i1n=nvirtc(isblkk(3,isk))                                   4d10s18
              do i2=i2s,i2e                                               4d10s18
               if(i2.eq.i2e)i1n=i1e                                       4d10s18
               i2m=i2-1                                                 4d10s18
               do i1=i10,i1n                                              4d10s18
                i1m=i1-1                                                  4d10s18
                do i4=0,noc(isblkk(2,isk))-1                              4d10s18
                 icol=i4+noc(isblkk(2,isk))*i1m                           4d10s18
                 if(icol.lt.ncut)then                                      3d12s18
                  iptry=icol/nhere                                            3d12s18
                 else                                                      3d12s18
                  iptry=nleft+(icol-ncut)/nhere0                              3d12s18
                 end if                                                    3d12s18
                 if(iptry.eq.mynowprog)then                             4d10s18
                  if(ipass.eq.1)then                                       4d10s18
                   ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblkk(1,isk))          4d10s18
                  else                                                     4d10s18
                   do i3=0,noc(isblkk(1,isk))-1                            4d10s18
                    iad1=i3xb(2,is3)+icol-il3                           4d10s18
     $                   +nhere3m*(i3+noc(isblk1(1,is3))*i2m)           4d10s18
                    bc(iad1)=bc(ibc(irecv+ip))                          4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                       4d10s18
                   end do                                                  4d10s18
                  end if                                                   4d10s18
                 end if                                                 4d10s18
                end do                                                     4d10s18
               end do                                                      4d10s18
               i10=1                                                      4d10s18
              end do                                                    4d10s18
             else if(isblkk(2,isk).eq.isblk1(1,is3).and.                  4d10s18
     $         isblkk(1,isk).eq.isblk1(3,is3).and.                       4d10s18
     $         isblkk(4,isk).eq.isblk1(4,is3))then                       4d10s18
              call ilimts(nvirtc(isblkk(3,isk)),nvirtc(isblkk(4,isk)),    4d10s18
     $           mynprocg,ip,jl,jh,i1s,i1e,i2s,i2e)                     4d10s18
              i10=i1s                                                     4d10s18
              i1n=nvirtc(isblkk(3,isk))                                   4d10s18
              do i2=i2s,i2e                                               4d10s18
               if(i2.eq.i2e)i1n=i1e                                       4d10s18
               i2m=i2-1                                                 4d10s18
               do i1=i10,i1n                                              4d10s18
                i1m=i1-1                                                  4d10s18
                do i3=0,noc(isblkk(1,isk))-1                              4d10s18
                 icol=i3+noc(isblkk(1,isk))*i2m                           4d10s18
                 if(icol.lt.ncut)then                                      3d12s18
                  iptry=icol/nhere                                            3d12s18
                 else                                                      3d12s18
                  iptry=nleft+(icol-ncut)/nhere0                              3d12s18
                 end if                                                    3d12s18
                 if(iptry.eq.mynowprog)then                             4d10s18
                  if(ipass.eq.1)then                                       4d10s18
                   ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblkk(2,isk))          4d10s18
                  else                                                     4d10s18
                   do i4=0,noc(isblkk(2,isk))-1                            4d10s18
                    iad1=i3xb(2,is3)+icol-il3                           4d10s18
     $                   +nhere3m*(i4+noc(isblk1(1,is3))*i1m)           4d10s18
                    bc(iad1)=bc(ibc(irecv+ip))                          4d10s18
                    ibc(irecv+ip)=ibc(irecv+ip)+1                       4d10s18
                   end do                                                  4d10s18
                  end if                                                   4d10s18
                 end if                                                 4d10s18
                end do                                                     4d10s18
               end do                                                      4d10s18
               i10=1                                                      4d10s18
              end do                                                    4d10s18
             end if
            end if
           end do                                                       4d10s18
          end do                                                        4d10s18
         end if                                                         4d10s18
         do is=1,nsdlkk                                                 2d23s18
          do is1x=1,nsdlk1                                              3d13s18
           if(isblk1(1,is1x).eq.isblkk(1,is).and.                       3d13s18
     $          isblk1(3,is1x).eq.isblkk(2,is).and.                     3d13s18
     $          isblk1(2,is1x).eq.isblkk(4,is).and.                     3d13s18
     $          noc(isblk1(2,is1x)).gt.0)then                           3d13s18
            call ilimts(noc(isblk1(2,is1x)),nvirtc(isblk1(4,is1x)),     3d13s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                      3d13s18
            il=il-1                                                     3d14s18
            call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,ip,ilq,ihq,i1s,i1e,i2s,i2e)                      3d13s18
            ncol=noc(isblk1(2,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(nhere0,1)                                        3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s                                                      2d22s18
            i1n=noc(isblk1(3,is1x))
            do i2=i2s,i2e                                                2d22s18
             if(i2.eq.i2e)i1n=i1e                                        2d22s18
             i2m=i2-1
             do i1=i10,i1n                                              3d14s18
              i1m=i1-1
              do io=0,noc(isblk1(2,is1x))-1                             3d14s18
               icol=io+i2m*noc(isblk1(2,is1x))                          3d14s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                        3d13s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                          3d14s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d13s18
                if(ipass.eq.1)then                                      3d13s18
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk1(1,is1x))        3d13s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk1(1,is1x))          3d13s18
                else                                                    3d13s18
                 do in=0,noc(isblk1(1,is1x))-1                          3d13s18
                  iad=i1xa(is)+in+noc(isblk1(1,is1x))*(i1m              3d13s18
     $                 +noc(isblk1(3,is1x))*(icol-il))                  3d13s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d13s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d13s18
                 end do                                                 3d13s18
                end if                                                  3d13s18
               end if                                                   3d13s18
              end do                                                    3d13s18
             end do                                                     3d13s18
             i10=1                                                      3d13s18
            end do                                                      3d13s18
           else if(isblk1(2,is1x).eq.isblkk(1,is).and.                  3d13s18
     $          isblk1(3,is1x).eq.isblkk(2,is).and.                     3d13s18
     $          isblk1(1,is1x).eq.isblkk(4,is).and.                     3d13s18
     $           noc(isblk1(1,is1x)).gt.0)then                          3d13s18
            call ilimts(noc(isblk1(1,is1x)),nvirtc(isblk1(4,is1x)),     3d13s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                      3d13s18
            il=il-1                                                     3d14s18
            call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,ip,ilq,ihq,i1s,i1e,i2s,i2e)                      3d13s18
            ncol=noc(isblk1(1,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(nhere0,1)                                        3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s                                                      2d22s18
            i1n=noc(isblk1(3,is1x))
            do i2=i2s,i2e                                                2d22s18
             if(i2.eq.i2e)i1n=i1e                                        2d22s18
             i2m=i2-1
             do i1=i10,i1n                                              3d14s18
              i1m=i1-1
              do io=0,noc(isblk1(1,is1x))-1                             3d14s18
               icol=io+i2m*noc(isblk1(1,is1x))                          3d14s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                        3d13s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                          3d14s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d13s18
                if(ipass.eq.1)then                                      3d13s18
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk1(2,is1x))        3d13s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk1(2,is1x))          3d13s18
                else                                                    3d13s18
                 do in=0,noc(isblk1(2,is1x))-1                          3d13s18
                  iad=i1xa(is)+in+noc(isblk1(2,is1x))*(i1m              3d13s18
     $                 +noc(isblk1(3,is1x))*(icol-il))                  3d13s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d13s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d13s18
                 end do                                                 3d13s18
                end if                                                  3d13s18
               end if                                                   3d13s18
              end do                                                    3d13s18
             end do                                                     3d13s18
             i10=1                                                      3d13s18
            end do                                                      3d13s18
           end if                                                       3d13s18
           if(isblk1(1,is1x).eq.isblkk(2,is).and.                       3d13s18
     $          isblk1(3,is1x).eq.isblkk(1,is).and.                     3d13s18
     $          isblk1(2,is1x).eq.isblkk(3,is).and.                     3d13s18
     $          noc(isblk1(2,is1x)).gt.0)then                           3d13s18
            call ilimts(noc(isblk1(2,is1x)),nvirtc(isblk1(4,is1x)),     3d13s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                      3d13s18
            il=il-1                                                     3d14s18
            call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,ip,ilq,ihq,i1s,i1e,i2s,i2e)                      3d13s18
            ncol=noc(isblk1(2,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(nhere0,1)                                        3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s                                                      2d22s18
            i1n=noc(isblk1(3,is1x))
            do i2=i2s,i2e                                                2d22s18
             if(i2.eq.i2e)i1n=i1e                                        2d22s18
             i2m=i2-1
             do i1=i10,i1n                                              3d14s18
              i1m=i1-1
              do io=0,noc(isblk1(2,is1x))-1                             3d14s18
               icol=io+i2m*noc(isblk1(2,is1x))                          3d14s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                        3d13s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                          3d13s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d13s18
                if(ipass.eq.1)then                                      3d13s18
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk1(1,is1x))        3d13s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk1(1,is1x))          3d13s18
                else                                                    3d13s18
                 do in=0,noc(isblk1(1,is1x))-1                          3d13s18
                  iad=i1xb(is)+i1m+noc(isblk1(3,is1x))*(in              3d26s18
     $                 +noc(isblk1(1,is1x))*(icol-il))                  3d26s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d13s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d13s18
                 end do                                                 3d13s18
                end if                                                  3d13s18
               end if                                                   3d13s18
              end do                                                    3d13s18
             end do                                                     3d13s18
             i10=1                                                      3d14s18
            end do                                                      3d13s18
           else if(isblk1(2,is1x).eq.isblkk(2,is).and.                  3d13s18
     $          isblk1(3,is1x).eq.isblkk(1,is).and.                     3d13s18
     $          isblk1(1,is1x).eq.isblkk(3,is).and.                     3d13s18
     $          noc(isblk1(1,is1x)).gt.0)then                           3d13s18
            call ilimts(noc(isblk1(1,is1x)),nvirtc(isblk1(4,is1x)),     3d13s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                      3d13s18
            il=il-1                                                     3d14s18
            call ilimts(noc(isblk1(3,is1x)),nvirtc(isblk1(4,is1x)),      2d22s18
     $          mynprocg,ip,ilq,ihq,i1s,i1e,i2s,i2e)                      3d13s18
            ncol=noc(isblk1(1,is1x))*nvirtc(isblk1(4,is1x))              3d12s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(nhere0,1)                                        3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s                                                      2d22s18
            i1n=noc(isblk1(3,is1x))
            do i2=i2s,i2e                                                2d22s18
             if(i2.eq.i2e)i1n=i1e                                        2d22s18
             i2m=i2-1
             do i1=i10,i1n                                              3d14s18
              i1m=i1-1
              do io=0,noc(isblk1(1,is1x))-1                             3d14s18
               icol=io+i2m*noc(isblk1(1,is1x))                          3d14s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                        3d13s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                          3d13s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d13s18
                if(ipass.eq.1)then                                      3d13s18
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk1(2,is1x))        3d13s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk1(2,is1x))          3d13s18
                else                                                    3d13s18
                 do in=0,noc(isblk1(2,is1x))-1                          3d13s18
                  iad=i1xb(is)+i1m+noc(isblk1(3,is1x))*(in              3d26s18
     $                 +noc(isblk1(2,is1x))*(icol-il))                  3d26s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d13s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d13s18
                 end do                                                 3d13s18
                end if                                                  3d13s18
               end if                                                   3d13s18
              end do                                                    3d13s18
             end do                                                     3d13s18
             i10=1                                                      3d14s18
            end do                                                      3d13s18
           end if                                                       3d13s18
          end do                                                        3d13s18
          do iso=1,nsdlk                                                3d14s18
           if(isblk(1,iso).eq.isblkk(1,is).and.                         3d15s18
     $       isblk(2,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(3,iso).eq.isblkk(2,is).and.noc(isblk(2,iso)).gt.0    3d15s18
     $       .and.noc(isblk(4,iso)).gt.0)then                           3d14s18
            call ilimts(noc(isblk(4,iso)),noc(isblk(2,iso)),mynprocg,   3d15s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
            il=il-1
            call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         ip,ilx,ihx,i1s,i1e,i2s,i2e)                              3d14s18
            ncol=noc(isblk(2,iso))*noc(isblk(4,iso))                    3d15s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(1,nhere0)                                         3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s
            i1n=noc(isblk(3,iso))
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             i2m=i2-1
             do i1=i10,i1n
              i1m=i1-1
              do io=0,noc(isblk(2,iso))-1
               icol=i2m+noc(isblk(4,iso))*io                            3d15s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                            3d12s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d14s18
                if(ipass.eq.1)then
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk(1,iso))            3d14s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk(1,iso))              3d14s18
                else                                                      2d23s18
                 do in=0,noc(isblk(1,iso))-1                              3d14s18
                  iad=i4ok(is)+in+noc(isblk(1,iso))*(i1m                3d15s18
     $                 +noc(isblk(3,iso))*(icol-il))                    3d15s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d14s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d14s18
                 end do
                end if
               end if
              end do
             end do
             i10=1
            end do
           else if(isblk(2,iso).eq.isblkk(1,is).and.                    3d15s18
     $       isblk(1,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(3,iso).eq.isblkk(2,is).and.noc(isblk(1,iso)).gt.0    3d15s18
     $       .and.noc(isblk(4,iso)).gt.0)then                           3d14s18
            call ilimts(noc(isblk(4,iso)),noc(isblk(1,iso)),mynprocg,    3d14s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
            il=il-1
            call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         ip,ilx,ihx,i1s,i1e,i2s,i2e)                              3d14s18
            ncol=noc(isblk(1,iso))*noc(isblk(4,iso))                    3d15s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(1,nhere0)                                         3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s
            i1n=noc(isblk(3,iso))
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             i2m=i2-1
             do i1=i10,i1n
              i1m=i1-1
              do io=0,noc(isblk(1,iso))-1
               icol=i2m+noc(isblk(4,iso))*io                            3d15s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                            3d12s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d14s18
                if(ipass.eq.1)then
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk(2,iso))            3d14s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk(2,iso))              3d14s18
                else                                                      2d23s18
                 do in=0,noc(isblk(2,iso))-1                              3d14s18
                  iad=i4ok(is)+in+noc(isblk(2,iso))*(i1m                3d15s18
     $                 +noc(isblk(3,iso))*(icol-il))                    3d15s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d14s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d14s18
                 end do
                end if
               end if
              end do
             end do
             i10=1
            end do
           else if(isblk(2,iso).eq.isblkk(1,is).and.                    3d15s18
     $       isblk(1,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(4,iso).eq.isblkk(2,is).and.noc(isblk(1,iso)).gt.0    3d15s18
     $       .and.noc(isblk(3,iso)).gt.0)then                           3d14s18
            call ilimts(noc(isblk(3,iso)),noc(isblk(1,iso)),mynprocg,   3d15s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
            il=il-1
            call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         ip,ilx,ihx,i1s,i1e,i2s,i2e)                              3d14s18
            ncol=noc(isblk(1,iso))*noc(isblk(3,iso))                    3d15s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(1,nhere0)                                         3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s
            i1n=noc(isblk(3,iso))
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             i2m=i2-1
             do i1=i10,i1n
              i1m=i1-1
              do io=0,noc(isblk(1,iso))-1
               icol=i1m+noc(isblk(3,iso))*io                            3d15s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                            3d12s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d14s18
                if(ipass.eq.1)then
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk(2,iso))            3d14s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk(2,iso))              3d14s18
                else                                                      2d23s18
                 do in=0,noc(isblk(2,iso))-1                              3d14s18
                  iad=i4ok(is)+in+noc(isblk(2,iso))*(i2m                3d15s18
     $                 +noc(isblk(4,iso))*(icol-il))                    3d15s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d14s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d14s18
                 end do
                end if
               end if
              end do
             end do
             i10=1
            end do
           else if(isblk(1,iso).eq.isblkk(1,is).and.                    3d15s18
     $       isblk(2,iso).eq.isblkk(4,is).and.                          3d15s18
     $       isblk(4,iso).eq.isblkk(2,is).and.noc(isblk(2,iso)).gt.0    3d15s18
     $       .and.noc(isblk(3,iso)).gt.0)then                           3d14s18
            call ilimts(noc(isblk(3,iso)),noc(isblk(2,iso)),mynprocg,   3d15s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d14s18
            il=il-1
            call ilimts(noc(isblk(3,iso)),noc(isblk(4,iso)),mynprocg,    3d14s18
     $         ip,ilx,ihx,i1s,i1e,i2s,i2e)                              3d14s18
            ncol=noc(isblk(2,iso))*noc(isblk(3,iso))                    3d15s18
            nhere0=ncol/mynprocg                                         3d12s18
            nleft=ncol-nhere0*mynprocg                                   3d12s18
            nhere=nhere0+1                                               3d12s18
            nhere0=max(1,nhere0)                                         3d14s18
            ncut=nleft*nhere                                             3d12s18
            i10=i1s
            i1n=noc(isblk(3,iso))
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             i2m=i2-1
             do i1=i10,i1n
              i1m=i1-1
              do io=0,noc(isblk(2,iso))-1
               icol=i1m+noc(isblk(3,iso))*io                            3d15s18
               if(icol.lt.ncut)then                                      3d12s18
                iptry=icol/nhere                                            3d12s18
               else                                                      3d12s18
                iptry=nleft+(icol-ncut)/nhere0                              3d12s18
               end if                                                    3d12s18
               if(iptry.eq.mynowprog)then                               3d14s18
                if(ipass.eq.1)then
                 ibc(nrecv+ip)=ibc(nrecv+ip)+noc(isblk(1,iso))            3d14s18
                 ior(ip+1,is)=ior(ip+1,is)+noc(isblk(1,iso))              3d14s18
                else                                                      2d23s18
                 do in=0,noc(isblk(1,iso))-1                              3d14s18
                  iad=i4ok(is)+in+noc(isblk(1,iso))*(i2m                3d15s18
     $                 +noc(isblk(4,iso))*(icol-il))                    3d15s18
                  bc(iad)=bc(ibc(irecv+ip))                             3d14s18
                  ibc(irecv+ip)=ibc(irecv+ip)+1                         3d14s18
                 end do
                end if
               end if
              end do
             end do
             i10=1
            end do
           end if
          end do                                                        3d14s18
         end do                                                         2d23s18
        end do                                                          2d23s18
        if(ipass.eq.1)then
         nst=0                                                          3d13s18
         nrt=0                                                          3d13s18
         do ip=0,mynprocg-1                                             3d13s18
          nst=nst+ibc(nsend+ip)                                         3d13s18
          nrt=nrt+ibc(nrecv+ip)                                         3d13s18
         end do                                                         3d13s18
         ibufs=ibcoff                                                   3d13s18
         ibufr=ibufs+max(nst,nrt)                                       3d13s18
         ibcoff=ibufr+nrt                                               3d13s18
         do i=0,max(nst,nrt)-1                                          4d12s18
          bc(ibufs+i)=0d0
         end do                                                         4d12s18
         call enough('updateg. 12',bc,ibc)
         jbcoff=ibufs                                                   3d13s18
         ns1=0
         ns2=0
         if(c3.ne.0d0)then                                              4d10s18
          do is3=1,nsdlk1                                               4d10s18
           call ilimts(noc(isblk1(3,is3)),noc(isblk1(4,is3)),mynprocg,  4d10s18
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        4d10s18
           nhere0=ih+1-il                                               4d10s18
           call ilimts(noc(isblk1(3,is3)),nvirtc(isblk1(4,is3)),        4d10s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               4d10s18
           nhere3=ih+1-il                                               4d10s18
           i3xa(1,is3)=jbcoff                                           4d10s18
           if(noc(isblk1(2,is3)).ne.0.and.noc(isblk1(4,is3)).ne.0)then  4d10s18
c     (vm|op)
            need=nhere0*nvirtc(isblk1(1,is3))*noc(isblk1(2,is3))        4d10s18
            jbcoff=jbcoff+need                                          4d10s18
           end if                                                       4d10s18
           i3xa(2,is3)=jbcoff                                           4d10s18
           if(noc(isblk1(1,is3)).ne.0.and.noc(isblk1(4,is3)).ne.0)then  4d10s18
c     (nv|op)
            need=nhere0*nvirtc(isblk1(2,is3))*noc(isblk1(1,is3))        4d10s18
            jbcoff=jbcoff+need                                          4d10s18
           end if                                                       4d10s18
           i3xa(3,is3)=jbcoff                                           4d10s18
           if(noc(isblk1(4,is3)).ne.0)then                              4d10s18
c     (vv'|op)
            need=nhere0*nvirtc(isblk1(1,is3))*nvirtc(isblk1(2,is3))     4d10s18
            jbcoff=jbcoff+need                                          4d10s18
           end if                                                       4d10s18
           i3xa(4,is3)=jbcoff                                           4d13s18
           if(noc(isblk1(4,is3)).ne.0)then                              4d13s18
            do isj=1,nsdlk                                              4d13s18
             if(isblk(3,isj).eq.isblk1(3,is3).and.                      4d13s18
     $            isblk(4,isj).eq.isblk1(4,is3))then                      4d13s18
              iok=1                                                     4d13s18
             else if(isblk(3,isj).eq.isblk1(4,is3).and.                 4d13s18
     $             isblk(4,isj).eq.isblk1(3,is3).and.                   4d13s18
     $             (isblk(1,isj).eq.isblk1(1,is3).or.                   4d13s18
     $             isblk(2,isj).eq.isblk1(1,is3)))then                  4d13s18
              need=nhere0*noc(isblk1(1,is3))*noc(isblk1(2,is3))         4d13s18
              jbcoff=jbcoff+need                                        4d13s18
              go to 12121                                               4d13s18
             end if                                                     4d13s18
            end do                                                      4d13s18
12121       continue                                                    4d13s18
           end if                                                       4d13s18
           i3xb(1,is3)=jbcoff                                           4d10s18
           if(noc(isblk1(2,is3)).ne.0)then
c     (v'm|ov)
            need=nhere3*nvirtc(isblk1(1,is3))*noc(isblk1(2,is3))        4d10s18
            jbcoff=jbcoff+need                                          4d10s18
           end if
           i3xb(2,is3)=jbcoff                                           4d10s18
           if(noc(isblk1(1,is3)).ne.0)then                              4d10s18
c     (nv'|ov)
            need=nhere3*nvirtc(isblk1(2,is3))*noc(isblk1(1,is3))        4d10s18
            jbcoff=jbcoff+need                                          4d10s18
           end if                                                       4d10s18
           call enough('updateg. 13',bc,ibc)
          end do                                                        4d10s18
         end if                                                         4d10s18
         do is=1,nsdlkk                                                 3d13s18
          i1xa(is)=jbcoff                                               3d13s18
          n12=noc(isblkk(1,is))*noc(isblkk(2,is))                       3d13s18
          if(n12.gt.0)then                                              3d15s18
          ncl=0
          if(noc(isblkk(4,is)).gt.0)then                                3d13s18
           call ilimts(noc(isblkk(4,is)),nvirtc(isblkk(3,is)),mynprocg,  3d13s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d13s18
           jbcoff=jbcoff+n12*(ih-il+1)                                   3d13s18
           ncl=ncl+ih-il+1
          end if                                                        3d13s18
          i1xb(is)=jbcoff                                               3d13s18
          i4ok(is)=jbcoff                                               3d14s18
          if(noc(isblkk(3,is)).gt.0)then                                3d13s18
           call ilimts(noc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,  3d13s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         3d13s18
           jbcoff=jbcoff+n12*(ih-il+1)                                   3d13s18
           ncl=ncl+ih-il+1
           i4ok(is)=jbcoff                                              3d14s18
           if(noc(isblkk(4,is)).gt.0)then                               3d14s18
            call ilimts(noc(isblkk(3,is)),noc(isblkk(4,is)),mynprocg,   3d14s18
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        3d14s18
            jbcoff=jbcoff+n12*(ih+1-il)                                 3d14s18
            ncl=ncl+ih+1-il                                             3d14s18
           end if                                                       3d14s18
          end if                                                        3d13s18
          ntt=n12*ncl
          ns=0
          do j=1,mynprocg
           ns=ns+ior(j,is)
          end do
          fact1=dfloat(ntt)/dfloat(n12)
          fact2=dfloat(ns)/dfloat(n12)
          ns1=ns1+ns
          ns2=ns2+ntt
          end if
         end do                                                         3d13s18
         jbufs=ibufs                                                    3d13s18
         do ip=0,mynprocg-1                                             3d13s18
          ibc(isend+ip)=jbufs                                           3d13s18
          jbufs=jbufs+ibc(nsend+ip)                                     3d13s18
         end do                                                         3d13s18
        else                                                            3d13s18
         ibcoff=ibufr                                                   3d13s18
        end if                                                          3d13s18
       end do                                                           2d22s18
       ibcoff=ibufr                                                     3d16s18
      end if                                                            2d22s18
c
c     ipass = 1 means count, ipass = 2 means do
c
      isend=ibcoff                                                      1d10s18
      irecv=isend+mynprocg                                              1d10s18
      nsend=irecv+mynprocg                                              1d10s18
      nrecv=nsend+mynprocg                                              1d10s18
      ibcoff=nrecv+mynprocg                                             1d10s18
      call enough('updateg. 14',bc,ibc)
      do i=0,mynprocg-1                                                 1d10s18
       bc(nsend+i)=0                                                    1d10s18
       bc(nrecv+i)=0                                                    1d10s18
       do j=1,idbk
        ios(i+1,j)=0
        ior(i+1,j)=0
       end do                                                           1d10s18
      end do                                                            1d10s18
cccccccccccccccccccccccccccccccccccccccccc
      do ipass=1,2                                                      1d10s18
c
c     transformation of first two indices and sending
c
       if(c3.ne.0d0)then                                                4d10s18
        me2mea=0
        me2meb=0
        me2mec=0
        do is3=1,nsdlk1                                                 4d11s18
         if(nvirtc(isblk1(1,is3)).ne.0.and.nvirtc(isblk1(2,is3)).ne.0.  4d10s18
     $     and.noc(isblk1(3,is3)).ne.0.and.nvirtc(isblk1(4,is3)).ne.0)  4d10s18
     $        then                                                      4d10s18
          ncol1=noc(isblk1(1,is3))*noc(isblk1(2,is3))
          nhere=ncol1/mynprocg                                           4d10s18
          nleft=ncol1-nhere*mynprocg                                     4d10s18
          nhere0=nhere                                                  4d12s18
          nhere=nhere+1                                                 4d10s18
          ncuta=nleft*nhere                                             4d10s18
          nherea=nhere                                                  4d10s18
          nlefta=nleft                                                  4d10s18
          nhere0a=nhere0                                                4d10s18
          ncol2=nvirtc(isblk1(1,is3))*noc(isblk1(2,is3))
          nhere=ncol2/mynprocg                                           4d10s18
          nleft=ncol2-nhere*mynprocg                                     4d10s18
          nhere0=nhere                                                  4d12s18
          nhere=nhere+1                                                 4d10s18
          ncutb=nleft*nhere                                             4d10s18
          nhereb=nhere                                                  4d10s18
          nleftb=nleft                                                  4d10s18
          nhere0b=nhere0                                                4d10s18
          if(isblk1(1,is3).ne.isblk1(2,is3))then                        4d10s18
           ncol3=noc(isblk1(1,is3))*nvirtc(isblk1(2,is3))               4d10s18
           nhere=ncol3/mynprocg                                           4d10s18
           nleft=ncol3-nhere*mynprocg                                     4d10s18
           nhere0=nhere                                                 4d12s18
           nhere=nhere+1                                                 4d10s18
           ncutc=nleft*nhere                                             4d10s18
           nherec=nhere                                                  4d10s18
           nleftc=nleft                                                  4d10s18
           nhere0c=nhere0                                                4d10s18
          end if
          if(ipass.eq.1)then                                            4d10s18
           call ilimts(noc(isblk1(3,is3)),noc(isblk1(4,is3)),mynprocg,  4d10s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         4d10s18
           noo=ih+1-il                                                  4d10s18
           call ilimts(noc(isblk1(3,is3)),nvirtc(isblk1(4,is3)),        4d10s18
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                4d10s18
           nov=ih+1-il                                                  4d10s18
           nsum=nov+noo                                                 4d10s18
           n1=nherea*nsum                                                4d10s18
           do ip=0,nlefta-1                                             4d12s18
            ibc(nsend+ip)=ibc(nsend+ip)+n1                              4d10s18
            if(ip.eq.mynowprog)then
             me2mea=me2mea+n1
            end if
           end do                                                       4d10s18
           n2=nhere0a*nsum                                               4d10s18
           do ip=nlefta,mynprocg-1                                      4d12s18
            ibc(nsend+ip)=ibc(nsend+ip)+n2                              4d10s18
            if(ip.eq.mynowprog)then
             me2mea=me2mea+n2
            end if
           end do                                                       4d10s18
           n1=nhereb*nsum                                                4d10s18
           do ip=0,nleftb-1                                             4d12s18
            ibc(nsend+ip)=ibc(nsend+ip)+n1                              4d10s18
            if(ip.eq.mynowprog)then
             me2meb=me2meb+n1
            end if
           end do                                                       4d10s18
           n2=nhere0b*nsum                                               4d10s18
           do ip=nleftb,mynprocg-1                                      4d12s18
            ibc(nsend+ip)=ibc(nsend+ip)+n2                              4d10s18
            if(ip.eq.mynowprog)then
             me2meb=me2meb+n2
            end if
           end do                                                       4d10s18
           if(isblk1(1,is3).ne.isblk1(2,is3))then                       4d10s18
            n1=nherec*nsum                                                4d10s18
            do ip=0,nleftc-1                                            4d12s18
             ibc(nsend+ip)=ibc(nsend+ip)+n1                              4d10s18
            if(ip.eq.mynowprog)then
             me2mec=me2mec+n1
            end if
            end do                                                       4d10s18
            n2=nhere0c*nsum                                               4d10s18
            do ip=nleftc,mynprocg-1                                     4d12s18
             ibc(nsend+ip)=ibc(nsend+ip)+n2                              4d10s18
            if(ip.eq.mynowprog)then
             me2mec=me2mec+n2
            end if
            end do                                                       4d10s18
           end if                                                       4d10s18
          else                                                          4d10s18
           if(noc(isblk1(4,is3)).ne.0)then                               4d10s18
            call ilimts(noc(isblk1(3,is3)),noc(isblk1(4,is3)),mynprocg,  4d10s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         4d10s18
            nhere=ih+1-il
            ifoo=ibcoff                                                  4d10s18
            ibcoff=ifoo+nbasdwsc(isblk1(1,is3))*nbasdwsc(isblk1(2,is3))  4d10s18
     $           *nhere                                                  4d10s18
            do i=0,nbasdwsc(isblk1(1,is3))*nbasdwsc(isblk1(2,is3))*
     $           nhere-1
             bc(ifoo+i)=xnan
            end do
            do i4=0,noc(isblk1(2,is3))-1
             do i3=0,nvirtc(isblk1(1,is3))-1
              iad1=ifoo+nhere*(i3+noc(isblk1(1,is3))
     $             +nbasdwsc(isblk1(1,is3))*i4)
              iad2=i3xa(1,is3)+nhere*(i3+nvirtc(isblk1(1,is3))*i4)      4d12s18
              do i12=0,nhere-1
               bc(iad1+i12)=bc(iad2+i12)
              end do
             end do
            end do
            do i4=0,nvirtc(isblk1(2,is3))-1                             4d12s18
             i4p=i4+noc(isblk1(2,is3))                                  4d12s18
             do i3=0,noc(isblk1(1,is3))-1                               4d12s18
              iad1=ifoo+nhere*(i3                                       4d12s18
     $             +nbasdwsc(isblk1(1,is3))*i4p)                        4d12s18
              iad2=i3xa(2,is3)+nhere*(i3+noc(isblk1(1,is3))*i4)         4d12s18
              do i12=0,nhere-1                                          4d12s18
               bc(iad1+i12)=bc(iad2+i12)                                4d12s18
              end do                                                    4d12s18
             end do                                                     4d12s18
            end do                                                      4d12s18
            do i4=0,nvirtc(isblk1(2,is3))-1                             4d12s18
             i4p=i4+noc(isblk1(2,is3))                                  4d12s18
             do i3=0,nvirtc(isblk1(1,is3))-1                            4d12s18
              i3p=i3+noc(isblk1(1,is3))                                 4d12s18
              iad1=ifoo+nhere*(i3p                                       4d12s18
     $             +nbasdwsc(isblk1(1,is3))*i4p)                        4d12s18
              iad2=i3xa(3,is3)+nhere*(i3+nvirtc(isblk1(1,is3))*i4)      4d12s18
              do i12=0,nhere-1                                          4d12s18
               bc(iad1+i12)=bc(iad2+i12)                                4d12s18
              end do                                                    4d12s18
             end do                                                     4d12s18
            end do                                                      4d12s18
            do is=1,nsdlk                                                4d10s18
             if(isblk(1,is).eq.isblk1(1,is3).and.                        4d10s18
     $            isblk(2,is).eq.isblk1(2,is3).and.                      4d10s18
     $            isblk(3,is).eq.isblk1(3,is3))then                      4d10s18
              i10=i1s                                                    4d10s18
              i1n=noc(isblk1(3,is3))                                    4d10s18
              i4o=ioooo(is)                                              4d10s18
              if(isblk(1,is).eq.isblk(2,is))then                         4d10s18
               nrowo=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2           4d10s18
               iswitch=0                                                 4d10s18
              else
               nrowo=noc(isblk(1,is))*noc(isblk(2,is))                   4d10s18
               iswitch=1                                                 4d10s18
              end if                                                     4d10s18
              do i2=i2s,i2e                                              4d10s18
               i2m=i2-1                                                 4d12s18
               if(i2.eq.i2e)i1n=i1e                                      4d10s18
               do i1=i10,i1n
                icol=i1+noc(isblk1(3,is3))*i2m-il                       4d12s18
                do i4=0,noc(isblk1(2,is3))-1                             4d10s18
                 do i3=0,noc(isblk1(1,is3))-1                            4d10s18
                  ix=max(i3,i4)                                          4d10s18
                  in=min(i3,i4)                                          4d10s18
                  ieq=((ix*(ix+1))/2)+in                                 4d10s18
                  inot=i3+noc(isblk(1,is))*i4                            4d10s18
                  iad1=i4o+(inot-ieq)*iswitch+ieq                        4d10s18
                  iad2=ifoo+icol+nhere*(i3+nbasdwsc(isblk1(1,is3))*i4)  4d12s18
                  bc(iad2)=bc(iad1)                                     4d12s18
                 end do                                                  4d10s18
                end do                                                   4d10s18
                i4o=i4o+nrowo                                            4d10s18
               end do
               i10=1
              end do
              go to 2121                                                4d13s18
             else if(isblk(2,is).eq.isblk1(1,is3).and.                   4d10s18
     $           isblk(1,is).eq.isblk1(2,is3).and.                      4d10s18
     $           isblk(3,is).eq.isblk1(3,is3))then                      4d10s18
              i10=i1s                                                    4d10s18
              i1n=noc(isblk1(3,is3))                                    4d10s18
              i4o=ioooo(is)                                              4d10s18
              nrowo=noc(isblk(1,is))*noc(isblk(2,is))                    4d10s18
              do i2=i2s,i2e                                              4d10s18
               i2m=i2-1                                                 4d12s18
               if(i2.eq.i2e)i1n=i1e                                      4d10s18
               do i1=i10,i1n
                icol=i1+noc(isblk1(3,is3))*i2m-il                       4d12s18
                do i3=0,noc(isblk1(1,is3))-1                            4d12s18
                 iad1=i4o+i3*noc(isblk1(2,is3))                         4d12s18
                 do i4=0,noc(isblk1(2,is3))-1                             4d10s18
                  iad2=ifoo+icol+nhere*(i4+nbasdwsc(isblk1(2,is3))*i3)  4d12s18
                  bc(iad2)=bc(iad1+i4)                                  4d12s18
                 end do                                                  4d10s18
                end do                                                   4d10s18
                i4o=i4o+nrowo                                            4d10s18
               end do
               i10=1
              end do
              go to 2121                                                4d13s18
             end if
            end do                                                       4d10s18
            if(noc(isblk1(1,is3))*noc(isblk1(2,is3)).ne.0)then
             nhereq=i3xb(1,is3)-i3xa(4,is3)                             4d13s18
             nherew=nhere*noc(isblk1(1,is3))*noc(isblk1(2,is3))         4d13s18
             if(nherew.gt.nhereq)then                                   4d13s18
             write(6,*)('could not find 4o for '),(isblk1(j,is3),j=1,4)
             write(6,*)('so it should be under i3ax(4 '),i3xa(4,is3)
             write(6,*)('number of words here = '),nhereq
             write(6,*)('number of words wanted = '),nherew
              call dws_sync
              call dws_finalize
              stop
             end if
             do i4=0,noc(isblk1(2,is3))-1
              do i3=0,noc(isblk1(1,is3))-1
               iad1=i3xa(4,is3)+nhere*(i3+noc(isblk1(1,is3))*i4)        4d13s18
               iad2=ifoo+nhere*(i3+nbasdwsc(isblk1(1,is3))*i4)          4d13s18
               do i12=0,nhere-1                                         4d13s18
                bc(iad2+i12)=bc(iad1+i12)                               4d13s18
               end do                                                   4d13s18
              end do                                                    4d13s18
             end do                                                     4d13s18
            end if
 2121       continue                                                    4d13s18
            nn=nbasdwsc(isblk1(1,is3))*nbasdwsc(isblk1(2,is3))
            itmp1=ibcoff                                                4d10s18
            ibcoff=itmp1+nn*nhere                                       4d10s18
            call enough('updateg. 15',bc,ibc)
            nr=nhere*nbasdwsc(isblk1(1,is3))
            if(min(nr,noc(isblk1(2,is3)),nbasdwsc(isblk1(2,is3)))        5d10s19
     $           .gt.0)then                                             5d10s19
            call dgemm('n','n',nr,noc(isblk1(2,is3)),                   4d10s18
     $           nbasdwsc(isblk1(2,is3)),1d0,bc(ifoo),nr,               4d13s18
     $           bc(itotm(isblk1(2,is3))),nbasdwsc(isblk1(2,is3)),0d0,  4d10s18
     $           bc(itmp1),nr,                                           4d10s18
     d' updateg.  5')
            end if                                                      2d25s19
            itmp2=ibcoff                                                4d10s18
            ibcoff=itmp2+nr*noc(isblk1(2,is3))
            call enough('updateg. 16',bc,ibc)
            do i2=0,noc(isblk1(2,is3))-1                                4d10s18
             do i1=0,nbasdwsc(isblk1(1,is3))-1                          4d10s18
              iad1=itmp1+nhere*(i1+nbasdwsc(isblk1(1,is3))*i2)          4d13s18
              iad2=itmp2+nhere*(i2+noc(isblk1(2,is3))*i1)               4d10s18
              do i34=0,nhere-1                                          4d10s18
               bc(iad2+i34)=bc(iad1+i34)                                4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
            nr=nhere*noc(isblk1(2,is3))
            if(nr.gt.0.and.nbasdwsc(isblk1(1,is3)).gt.0)then            2d25s19
            call dgemm('n','n',nr,nbasdwsc(isblk1(1,is3)),              4d10s18
     $           nbasdwsc(isblk1(1,is3)),1d0,bc(itmp2),nr,              4d10s18
     $           bc(iumat(isblk1(1,is3))),nbasdwsc(isblk1(1,is3)),0d0,  4d10s18
     $           bc(itmp1),nr,                                          4d13s18
     d' updateg.  6')
            end if                                                      2d25s19
            do i2=0,noc(isblk1(2,is3))-1                                4d10s18
             do i1=0,nvirtc(isblk1(1,is3))-1                            4d10s18
              i1p=i1+noc(isblk1(1,is3))                                 4d10s18
              iad1=itmp1+nhere*(i2+noc(isblk1(2,is3))*i1p)              4d13s18
              icol=i1+nvirtc(isblk1(1,is3))*i2                          4d10s18
              if(icol.lt.ncutb)then                                     4d10s18
               ip=icol/nhereb                                           4d10s18
              else                                                      4d10s18
               ip=nleftb+(icol-ncutb)/nhere0b                           4d10s18
              end if                                                    4d10s18
              if(ip.lt.0.or.ip.ge.mynprocg)then
               write(6,*)('ip out of range: '),ip,i1,i2,icol,ncol2
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ilooper.ge.1)then
               call dws_sync
               call dws_finalize
               stop
              end if
              if(ibc(isend+ip)-ioffsx.lt.0.or.                          11d1s22
     $             ibc(isend+ip)-ioffsx.gt.maxbc)then                   11d1s22
               write(6,*)('isend out of range '),ibc(isend+ip)-ioffsx   11d1s22
               write(6,*)('for ip = '),ip,isend
               call dws_sync
               call dws_finalize
               stop
              end if
              do i34=0,nhere-1                                          4d10s18
               bc(ibc(isend+ip))=bc(iad1+i34)                           4d10s18
               ibc(isend+ip)=ibc(isend+ip)+1                            4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
            itmp3=ibcoff                                                4d10s18
            ibcoff=itmp3+nr*noc(isblk1(1,is3))                          4d10s18
            call enough('updateg. 17',bc,ibc)
            if(nr.gt.0.and.noc(isblk1(1,is3)).gt.0)then                 2d25s19
            call dgemm('n','n',nr,noc(isblk1(1,is3)),noc(isblk1(1,is3)),4d10s18
     $           1d0,bc(itmp1),nr,bc(itmat(isblk1(1,is3))),             4d13s18
     $           noc(isblk1(1,is3)),0d0,bc(itmp3),nr,                   4d10s18
     d' updateg.  7')
            end if                                                      2d25s19
            do i2=0,noc(isblk1(2,is3))-1                                4d10s18
             do i1=0,noc(isblk1(1,is3))-1                               4d10s18
              iad1=itmp3+nhere*(i2+noc(isblk1(2,is3))*i1)               4d10s18
              icol=i1+noc(isblk1(1,is3))*i2                             4d10s18
              if(icol.lt.ncuta)then                                     4d10s18
               ip=icol/nherea                                           4d10s18
              else                                                      4d10s18
               ip=nlefta+(icol-ncuta)/nhere0a                           4d10s18
              end if                                                    4d10s18
              do i34=0,nhere-1                                          4d10s18
               bc(ibc(isend+ip))=bc(iad1+i34)                           4d10s18
               ibc(isend+ip)=ibc(isend+ip)+1                            4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
            ibcoff=itmp2                                                4d10s18
            if(isblk1(1,is3).ne.isblk1(2,is3))then                      4d10s18
             do i2=0,nbasdwsc(isblk1(2,is3))-1                          4d10s18
              do i1=0,nbasdwsc(isblk1(1,is3))-1                         4d10s18
               iad1=ifoo+nhere*(i1+nbasdwsc(isblk1(1,is3))*i2)          4d13s18
               iad2=itmp1+nhere*(i2+nbasdwsc(isblk1(2,is3))*i1)         4d13s18
               do i34=0,nhere-1                                         4d10s18
                bc(iad2+i34)=bc(iad1+i34)                               4d10s18
               end do                                                   4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
             nr=nhere*nbasdwsc(isblk1(2,is3))                           4d10s18
             if(min(nr,noc(isblk1(1,is3)),nbasdwsc(isblk1(1,is3)))       5d10s19
     $            .gt.0)then                                            5d10s19
             call dgemm('n','n',nr,noc(isblk1(1,is3)),                  4d10s18
     $            nbasdwsc(isblk1(1,is3)),1d0,bc(itmp1),nr,             4d13s18
     $            bc(itotm(isblk1(1,is3))),nbasdwsc(isblk1(1,is3)),     4d10s18
     $            0d0,bc(ifoo),nr,                                      4d13s18
     d' updateg.  8')
             end if                                                     2d25s19
             do i2=0,nbasdwsc(isblk1(2,is3))-1                          4d10s18
              do i1=0,noc(isblk1(1,is3))-1                              4d10s18
               iad1=ifoo+nhere*(i2+nbasdwsc(isblk1(2,is3))*i1)          4d13s18
               iad2=itmp1+nhere*(i1+noc(isblk1(1,is3))*i2)              4d13s18
               do i34=0,nhere-1                                         4d10s18
                bc(iad2+i34)=bc(iad1+i34)                               4d10s18
               end do                                                   4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
             nr=nhere*noc(isblk1(1,is3))                                4d10s18
             jumat=iumat(isblk1(2,is3))+noc(isblk1(2,is3))              4d10s18
     $            *nbasdwsc(isblk1(2,is3))                              4d10s18
             if(nr.gt.0.and.nbasdwsc(isblk1(2,is3)).gt.0)then           2d25s19
             call dgemm('n','n',nr,nvirtc(isblk1(2,is3)),               4d10s18
     $            nbasdwsc(isblk1(2,is3)),1d0,bc(itmp1),nr,             4d13s18
     $            bc(jumat),nbasdwsc(isblk1(2,is3)),0d0,bc(ifoo),nr,    4d13s18
     d' updateg.  9')
             end if                                                     2d25s19
             do i2=0,nvirtc(isblk1(2,is3))-1                            4d10s18
              do i1=0,noc(isblk1(1,is3))-1                              4d10s18
               iad1=ifoo+nhere*(i1+noc(isblk1(1,is3))*i2)               4d13s18
               icol=i1+noc(isblk1(1,is3))*i2                            4d10s18
               if(icol.lt.ncutc)then                                    4d10s18
                ip=icol/nherec                                          4d10s18
               else                                                     4d10s18
                ip=nleftc+(icol-ncutc)/nhere0c                          4d10s18
               end if                                                   4d10s18
               do i34=0,nhere-1                                         4d10s18
                bc(ibc(isend+ip))=bc(iad1+i34)                          4d10s18
                ibc(isend+ip)=ibc(isend+ip)+1                           4d10s18
               end do                                                   4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
            end if                                                      4d10s18
            ibcoff=ifoo                                                 4d10s18
           end if                                                       4d10s18
           call ilimts(noc(isblk1(3,is3)),nvirtc(isblk1(4,is3)),        4d10s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               4d10s18
           nhere=ih+1-il                                                4d10s18
           ifov=ibcoff                                                  4d10s18
           nrowf=nbasdwsc(isblk1(1,is3))*nbasdwsc(isblk1(2,is3))        4d10s18
           need=nhere*nrowf                                             4d10s18
           ibcoff=ifov+need                                             4d10s18
           call enough('updateg. 18',bc,ibc)
           do i=0,need-1                                                4d10s18
            bc(ifov+i)=xnan                                             4d10s18
           end do                                                       4d10s18
           do i4=0,noc(isblk1(2,is3))-1                                 4d13s18
            do i3=0,nvirtc(isblk1(1,is3))-1                             4d13s18
             i3p=i3+noc(isblk1(1,is3))                                  4d13s18
             iad1=i3xb(1,is3)+nhere*(i3+nvirtc(isblk1(1,is3))*i4)       4d13s18
             iad2=ifov+nhere*(i3p+nbasdwsc(isblk1(1,is3))*i4)           4d13s18
             do i12=0,nhere-1                                           4d13s18
              bc(iad2+i12)=bc(iad1+i12)                                 4d13s18
             end do                                                     4d13s18
            end do                                                      4d13s18
           end do                                                       4d13s18
           do i4=0,nvirtc(isblk1(2,is3))-1                              4d13s18
            i4p=i4+noc(isblk1(2,is3))                                   4d13s18
            do i3=0,noc(isblk1(1,is3))-1                                4d13s18
             iad1=i3xb(2,is3)+nhere*(i3+noc(isblk1(1,is3))*i4)          4d13s18
             iad2=ifov+nhere*(i3+nbasdwsc(isblk1(1,is3))*i4p)           4d13s18
             do i12=0,nhere-1                                           4d13s18
              bc(iad2+i12)=bc(iad1+i12)                                 4d13s18
             end do                                                     4d13s18
            end do                                                      4d13s18
           end do                                                       4d13s18
           i10=i1s                                                      4d10s18
           i1n=noc(isblk1(3,is3))                                       4d13s18
           jfov=ifov                                                    4d10s18
           if(isblk1(1,is3).eq.isblk1(2,is3))then                       4d10s18
            nrow3=(nvirtc(isblk1(1,is3))*(nvirtc(isblk1(1,is3))+1))/2   4d10s18
            iswitch=0                                                   4d10s18
            nrowx=(noc(isblk1(1,is3))*(noc(isblk1(1,is3))+1))/2         4d10s18
           else                                                         4d10s18
            nrow3=nvirtc(isblk1(1,is3))*nvirtc(isblk1(2,is3))           4d10s18
            nrowx=noc(isblk1(1,is3))*noc(isblk1(2,is3))                 4d10s18
            iswitch=1                                                   4d10s18
           end if                                                       4d10s18
           j3=i3x(is3)
           jonex=ionex(is3)                                             4d10s18
           do i2=i2s,i2e                                                4d10s18
            if(i2.eq.i2e)i1n=i1e                                        4d10s18
            do i1=i10,i1n                                               4d10s18
             do i4=0,noc(isblk1(2,is3))-1                               4d10s18
              do i3=0,noc(isblk1(1,is3))-1                              4d10s18
               ix=max(i3,i4)                                            4d10s18
               in=min(i3,i4)                                            4d10s18
               ieq=((ix*(ix+1))/2)+in                                   4d10s18
               inot=i3+noc(isblk1(1,is3))*i4                            4d10s18
               iad1=jonex+(inot-ieq)*iswitch+ieq                        4d10s18
               iad2=jfov+nhere*(i3+nbasdwsc(isblk1(1,is3))*i4)          4d13s18
               bc(iad2)=bc(iad1)                                        4d13s18
              end do                                                    4d10s18
             end do                                                     4d10s18
             do i4=0,nvirtc(isblk1(2,is3))-1                            4d10s18
              i4p=i4+noc(isblk1(2,is3))                                 4d13s18
              do i3=0,nvirtc(isblk1(1,is3))-1                           4d10s18
               i3p=i3+noc(isblk1(1,is3))                                4d13s18
               ix=max(i3,i4)                                            4d10s18
               in=min(i3,i4)                                            4d10s18
               ieq=((ix*(ix+1))/2)+in                                   4d10s18
               inot=i3+nvirtc(isblk1(1,is3))*i4                         4d10s18
               iad1=j3+(inot-ieq)*iswitch+ieq                           4d10s18
               iad2=jfov+nhere*(i3p+nbasdwsc(isblk1(1,is3))*i4p)        4d13s18
               bc(iad2)=bc(iad1)                                        4d13s18
              end do                                                    4d10s18
             end do                                                     4d10s18
             jfov=jfov+1                                                4d13s18
             j3=j3+nrow3                                                4d10s18
             jonex=jonex+nrowx                                          4d10s18
            end do                                                      4d10s18
            i10=1                                                       4d10s18
           end do                                                       4d10s18
           nn=nbasdwsc(isblk1(2,is3))*nbasdwsc(isblk1(1,is3))           4d13s18
           itmp1=ibcoff                                                 4d10s18
           ibcoff=itmp1+nn*nhere                                        4d10s18
           call enough('updateg. 19',bc,ibc)
           nr=nhere*nbasdwsc(isblk1(1,is3))
           if(min(nr,noc(isblk1(2,is3)),nbasdwsc(isblk1(2,is3)))        5d10s19
     $          .gt.0)then                                              5d10s19
           call dgemm('n','n',nr,noc(isblk1(2,is3)),                    4d10s18
     $           nbasdwsc(isblk1(2,is3)),1d0,bc(ifov),nr,               4d13s18
     $           bc(itotm(isblk1(2,is3))),nbasdwsc(isblk1(2,is3)),0d0,  4d10s18
     $           bc(itmp1),nr,                                          4d13s18
     d' updateg. 10')
           end if                                                       2d25s19
           itmp2=ibcoff                                                 4d10s18
           ibcoff=itmp2+nr*noc(isblk1(2,is3))
           call enough('updateg. 20',bc,ibc)
           do i2=0,noc(isblk1(2,is3))-1                                 4d10s18
            do i1=0,nbasdwsc(isblk1(1,is3))-1                           4d10s18
             iad1=itmp1+nhere*(i1+nbasdwsc(isblk1(1,is3))*i2)           4d13s18
             iad2=itmp2+nhere*(i2+noc(isblk1(2,is3))*i1)                4d10s18
             do i34=0,nhere-1                                           4d10s18
              bc(iad2+i34)=bc(iad1+i34)                                 4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           end do                                                       4d10s18
           nr=nhere*noc(isblk1(2,is3))
           if(nr.gt.0.and.nbasdwsc(isblk1(1,is3)).gt.0)then             2d25s19
           call dgemm('n','n',nr,nbasdwsc(isblk1(1,is3)),               4d10s18
     $           nbasdwsc(isblk1(1,is3)),1d0,bc(itmp2),nr,              4d10s18
     $           bc(iumat(isblk1(1,is3))),nbasdwsc(isblk1(1,is3)),0d0,  4d10s18
     $           bc(itmp1),nr,                                          4d13s18
     d' updateg. 11')
           end if                                                       2d25s19
           do i2=0,noc(isblk1(2,is3))-1                                 4d10s18
            do i1=0,nvirtc(isblk1(1,is3))-1                             4d10s18
             i1p=i1+noc(isblk1(1,is3))                                  4d10s18
             iad1=itmp1+nhere*(i2+noc(isblk1(2,is3))*i1p)               4d13s18
             icol=i1+nvirtc(isblk1(1,is3))*i2                           4d10s18
             if(icol.lt.ncutb)then                                      4d10s18
              ip=icol/nhereb                                            4d10s18
             else                                                       4d10s18
              ip=nleftb+(icol-ncutb)/nhere0b                            4d10s18
             end if                                                     4d10s18
             do i34=0,nhere-1                                           4d10s18
              bc(ibc(isend+ip))=bc(iad1+i34)                            4d10s18
              ibc(isend+ip)=ibc(isend+ip)+1                             4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           end do                                                       4d10s18
           itmp3=ibcoff                                                 4d10s18
           ibcoff=itmp3+nr*noc(isblk1(1,is3))                           4d10s18
           call enough('updateg. 21',bc,ibc)
           if(nr.gt.0.and.noc(isblk1(1,is3)).gt.0)then                  2d25s19
           call dgemm('n','n',nr,noc(isblk1(1,is3)),noc(isblk1(1,is3)), 4d10s18
     $           1d0,bc(itmp1),nr,bc(itmat(isblk1(1,is3))),              4d10s18
     $           noc(isblk1(1,is3)),0d0,bc(itmp3),nr,                   4d10s18
     d' updateg. 12')
           end if                                                       2d25s19
           do i2=0,noc(isblk1(2,is3))-1                                 4d10s18
            do i1=0,noc(isblk1(1,is3))-1                                4d10s18
             iad1=itmp3+nhere*(i2+noc(isblk1(2,is3))*i1)                4d10s18
             icol=i1+noc(isblk1(1,is3))*i2                              4d10s18
             if(icol.lt.ncuta)then                                      4d10s18
              ip=icol/nherea                                            4d10s18
             else                                                       4d10s18
              ip=nlefta+(icol-ncuta)/nhere0a                            4d10s18
             end if                                                     4d10s18
             do i34=0,nhere-1                                           4d10s18
              bc(ibc(isend+ip))=bc(iad1+i34)                            4d10s18
              ibc(isend+ip)=ibc(isend+ip)+1                             4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           end do                                                       4d10s18
           ibcoff=itmp2                                                 4d10s18
           if(isblk1(1,is3).ne.isblk1(2,is3))then                       4d10s18
            do i2=0,nbasdwsc(isblk1(2,is3))-1                           4d10s18
             do i1=0,nbasdwsc(isblk1(1,is3))-1                          4d10s18
              iad1=ifov+nhere*(i1+nbasdwsc(isblk1(1,is3))*i2)           4d13s18
              iad2=itmp1+nhere*(i2+nbasdwsc(isblk1(2,is3))*i1)          4d13s18
              do i34=0,nhere-1                                          4d10s18
               bc(iad2+i34)=bc(iad1+i34)                                4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
            nr=nhere*nbasdwsc(isblk1(2,is3))                            4d10s18
            if(min(nr,noc(isblk1(1,is3)),nbasdwsc(isblk1(1,is3)))       5d10s19
     $           .gt.0)then                                             5d10s19
            call dgemm('n','n',nr,noc(isblk1(1,is3)),                   4d10s18
     $            nbasdwsc(isblk1(1,is3)),1d0,bc(itmp1),nr,             4d13s18
     $            bc(itotm(isblk1(1,is3))),nbasdwsc(isblk1(1,is3)),     4d10s18
     $            0d0,bc(ifov),nr,                                      4d13s18
     d' updateg. 13')
            end if                                                      2d25s19
            do i2=0,nbasdwsc(isblk1(2,is3))-1                           4d10s18
             do i1=0,noc(isblk1(1,is3))-1                               4d10s18
              iad1=ifov+nhere*(i2+nbasdwsc(isblk1(2,is3))*i1)           4d13s18
              iad2=itmp1+nhere*(i1+noc(isblk1(1,is3))*i2)               4d13s18
              do i34=0,nhere-1                                          4d10s18
               bc(iad2+i34)=bc(iad1+i34)                                4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
            nr=nhere*noc(isblk1(1,is3))                                 4d10s18
            jumat=iumat(isblk1(2,is3))+noc(isblk1(2,is3))               4d10s18
     $            *nbasdwsc(isblk1(2,is3))                              4d10s18
            if(nr.gt.0.and.nbasdwsc(isblk1(2,is3)).gt.0)then            2d25s19
            call dgemm('n','n',nr,nvirtc(isblk1(2,is3)),                4d10s18
     $            nbasdwsc(isblk1(2,is3)),1d0,bc(itmp1),nr,             4d13s18
     $            bc(jumat),nbasdwsc(isblk1(2,is3)),0d0,bc(ifov),nr,    4d13s18
     d' updateg. 14')
            end if                                                      2d25s19
            do i2=0,nvirtc(isblk1(2,is3))-1                             4d10s18
             do i1=0,noc(isblk1(1,is3))-1                               4d10s18
              iad1=ifov+nhere*(i1+noc(isblk1(1,is3))*i2)                4d13s18
              icol=i1+noc(isblk1(1,is3))*i2                             4d10s18
              if(icol.lt.ncutc)then                                     4d10s18
               ip=icol/nherec                                           4d10s18
              else                                                      4d10s18
               ip=nleftc+(icol-ncutc)/nhere0c                           4d10s18
              end if                                                    4d10s18
              do i34=0,nhere-1                                          4d10s18
               bc(ibc(isend+ip))=bc(iad1+i34)                           4d10s18
               ibc(isend+ip)=ibc(isend+ip)+1                            4d10s18
              end do                                                    4d10s18
             end do                                                     4d10s18
            end do                                                      4d10s18
           end if                                                       4d10s18
           ibcoff=ifov                                                  4d10s18
          end if                                                        4d10s18
         end if                                                         4d10s18
        end do                                                          4d10s18
       end if                                                           4d10s18
       if(c2.ne.0d0)then                                                3d15s18
c
c     for j
c
        do is=1,nsdlk                                                   3d15s18
         if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.        3d15s18
     $      nvirtc(isblk(3,is)).ne.0.and.nvirtc(isblk(4,is)).ne.0)then  3d15s18
          call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg, 3d15s18
     $        mynowprog,il,ih,i1s,i2e,i2s,i2e)                          3d15s18
          nhere=ih+1-il                                                 3d15s18
          nrow=noc(isblk(1,is))*noc(isblk(2,is))                        3d15s18
          il0=il
          if(ipass.eq.1)then                                            3d15s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           do ip=0,mynprocg-1                                           1d26s18
            if(ip.ge.nleft)then                                         1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngot*nhere                     1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngot*nhere                       1d26s18
            else                                                        1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngotp*nhere                    1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngotp*nhere                      1d26s18
            end if                                                      1d26s18
           end do                                                       1d26s18
          else                                                          3d15s18
           itmp1=ibcoff                                                  3d15s18
           itmp2=itmp1+nhere*nrow                                        3d15s18
           ibcoff=itmp2+nhere*nrow                                       3d15s18
           call enough('updateg. 22',bc,ibc)
           if(isblk(1,is).eq.isblk(2,is))then
            ii=jmats(is)                                                 3d15s18
            do i12=0,nhere-1
             do i4=0,noc(isblk(1,is))-1                                  3d15s18
              do i3=0,i4                                                 3d15s18
               iad1=itmp1+i12+nhere*(i3+noc(isblk(1,is))*i4)             3d15s18
               bc(iad1)=bc(ii)                                           3d15s18
               iad1=itmp1+i12+nhere*(i4+noc(isblk(1,is))*i3)             3d15s18
               bc(iad1)=bc(ii)                                           3d15s18
               ii=ii+1                                                   3d15s18
              end do                                                     3d15s18
             end do                                                      3d15s18
            end do                                                       3d15s18
           else                                                          3d15s18
            ii=jmats(is)                                                3d16s18
            do i12=0,nhere-1                                             3d15s18
             do i34=0,nrow-1                                             3d15s18
              iad1=itmp1+i12+nhere*i34                                         3d15s18
              bc(iad1)=bc(ii)                                           3d16s18
              ii=ii+1                                                   3d16s18
             end do                                                      3d15s18
            end do                                                       3d15s18
           end if                                                        3d15s18
           nx=nhere*noc(isblk(1,is))                                     3d15s18
           if(nx.gt.0.and.noc(isblk(2,is)).gt.0)then                    2d25s19
           call dgemm('n','n',nx,noc(isblk(2,is)),noc(isblk(2,is)),1d0,   3d15s18
     $         bc(itmp1),nx,bc(itmat(isblk(2,is))),noc(isblk(2,is)),0d0,3d15s18
     $         bc(itmp2),nx,                                            3d15s18
     d' updateg. 15')
           end if                                                       2d25s19
           do i12=0,nhere-1                                              3d15s18
            do i4=0,noc(isblk(2,is))-1                                   3d15s18
             do i3=0,noc(isblk(1,is))-1                                  3d15s18
              iad1=itmp1+i12+nhere*(i4+noc(isblk(2,is))*i3)               3d15s18
              iad2=itmp2+i12+nhere*(i3+noc(isblk(1,is))*i4)               3d15s18
              bc(iad1)=bc(iad2)                                         3d16s18
             end do                                                     3d16s18
            end do                                                      3d16s18
           end do                                                       3d16s18
           nx=nhere*noc(isblk(2,is))                                    3d16s18
           if(nx.gt.0.and.noc(isblk(1,is)).gt.0)then                    2d25s19
           call dgemm('n','n',nx,noc(isblk(1,is)),noc(isblk(1,is)),1d0, 3d16s18
     $          bc(itmp1),nx,bc(itmat(isblk(1,is))),noc(isblk(1,is)),   3d16s18
     $          0d0,bc(itmp2),nx,                                       3d16s18
     d' updateg. 16')
           end if                                                       2d25s19
           icol=1                                                       1d26s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           ips=0                                                        1d26s18
           mytop=ngot                                                   1d26s18
           if(nleft.gt.0)mytop=ngotp                                    1d26s18
           do i4=0,noc(isblk(2,is))-1                                   3d16s18
            do i3=0,noc(isblk(1,is))-1                                  3d16s18
             if(icol.gt.mytop)then                                      1d26s18
              ips=ips+1                                                 1d26s18
              if(ips.ge.nleft)then                                      1d26s18
               mytop=mytop+ngot                                         1d26s18
              else                                                      1d26s18
               mytop=mytop+ngotp                                        1d26s18
              end if                                                    1d26s18
             end if                                                     1d26s18
             do i12=0,nhere-1                                           1d26s18
              i12p=i12+il0-1                                            1d26s18
              iad1=itmp2+i12+nhere*(i4+noc(isblk(2,is))*i3)             3d16s18
              bc(ibc(isend+ips))=bc(iad1)                               1d26s18
              ibc(isend+ips)=ibc(isend+ips)+1                           1d26s18
             end do                                                     1d26s18
             icol=icol+1                                                1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           ibcoff=itmp1                                                 1d26s18
          end if                                                        3d16s18
         end if                                                         3d15s18
        end do                                                          3d15s18
c
c     for k and entourge
c
        do is=1,nsdlkk                                                  3d16s18
         if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.      3d16s18
     $      nvirtc(isblkk(3,is)).ne.0.and.nvirtc(isblkk(4,is)).ne.0)then3d16s18
          call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),        3d16s18
     $        mynprocg,mynowprog,il,ih,i1s,i2e,i2s,i2e)                 3d16s18
          nhere=ih+1-il                                                 3d15s18
          nrow=noc(isblkk(1,is))*noc(isblkk(2,is))                        3d15s18
          il0=il
          if(ipass.eq.1)then                                            3d15s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           do ip=0,mynprocg-1                                           1d26s18
            if(ip.ge.nleft)then                                         1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngot*nhere                     1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngot*nhere                       1d26s18
            else                                                        1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngotp*nhere                    1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngotp*nhere                      1d26s18
            end if                                                      1d26s18
           end do                                                       1d26s18
          else                                                          3d15s18
           itmp1=ibcoff                                                  3d15s18
           itmp2=itmp1+nhere*nrow                                        3d15s18
           ibcoff=itmp2+nhere*nrow                                       3d15s18
           call enough('updateg. 23',bc,ibc)
           ii=kmats(is)                                                 3d16s18
           do i12=0,nhere-1                                             3d15s18
            do i34=0,nrow-1                                             3d15s18
             iad1=itmp1+i12+nhere*i34                                         3d15s18
             bc(iad1)=bc(ii)                                            3d16s18
             ii=ii+1                                                    3d16s18
            end do                                                      3d15s18
           end do                                                       3d15s18
           nx=nhere*noc(isblkk(1,is))                                     3d15s18
           if(nx.gt.0.and.noc(isblkk(2,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(2,is)),noc(isblkk(2,is)),
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(2,is))),               3d16s18
     $          noc(isblkk(2,is)),0d0,bc(itmp2),nx,                     3d16s18
     d' updateg. 17')
           end if                                                       2d25s19
           do i12=0,nhere-1                                              3d15s18
            do i4=0,noc(isblkk(2,is))-1                                   3d15s18
             do i3=0,noc(isblkk(1,is))-1                                  3d15s18
              iad1=itmp1+i12+nhere*(i4+noc(isblkk(2,is))*i3)               3d15s18
              iad2=itmp2+i12+nhere*(i3+noc(isblkk(1,is))*i4)               3d15s18
              bc(iad1)=bc(iad2)                                         3d16s18
             end do                                                     3d16s18
            end do                                                      3d16s18
           end do                                                       3d16s18
           nx=nhere*noc(isblkk(2,is))                                    3d16s18
           if(nx.gt.0.and.noc(isblkk(1,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(1,is)),noc(isblkk(1,is)),   3d16s18
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(1,is))),               3d16s18
     $          noc(isblkk(1,is)),0d0,bc(itmp2),nx,                     3d16s18
     d' updateg. 18')
           end if                                                       2d25s19
           icol=1                                                       1d26s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           ips=0                                                        1d26s18
           mytop=ngot                                                   1d26s18
           if(nleft.gt.0)mytop=ngotp                                    1d26s18
           do i4=0,noc(isblkk(2,is))-1                                   3d16s18
            do i3=0,noc(isblkk(1,is))-1                                  3d16s18
             if(icol.gt.mytop)then                                      1d26s18
              ips=ips+1                                                 1d26s18
              if(ips.ge.nleft)then                                      1d26s18
               mytop=mytop+ngot                                         1d26s18
              else                                                      1d26s18
               mytop=mytop+ngotp                                        1d26s18
              end if                                                    1d26s18
             end if                                                     1d26s18
             do i12=0,nhere-1                                           1d26s18
              i12p=i12+il0-1                                            1d26s18
              iad1=itmp2+i12+nhere*(i4+noc(isblkk(2,is))*i3)             3d16s18
              bc(ibc(isend+ips))=bc(iad1)                               1d26s18
              ibc(isend+ips)=ibc(isend+ips)+1                           1d26s18
             end do                                                     1d26s18
             icol=icol+1                                                1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           ibcoff=itmp1                                                 1d26s18
          end if                                                        3d16s18
         end if                                                         3d15s18
         if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.      3d16s18
     $      nvirtc(isblkk(3,is)).ne.0.and.noc(isblkk(4,is)).ne.0)then   3d16s18
          call ilimts(noc(isblkk(4,is)),nvirtc(isblkk(3,is)),           3d16s18
     $        mynprocg,mynowprog,il,ih,i1s,i2e,i2s,i2e)                 3d16s18
          nhere=ih+1-il                                                 3d15s18
          nrow=noc(isblkk(1,is))*noc(isblkk(2,is))                        3d15s18
          il0=il
          if(ipass.eq.1)then                                            3d15s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           do ip=0,mynprocg-1                                           1d26s18
            if(ip.ge.nleft)then                                         1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngot*nhere                     1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngot*nhere                       1d26s18
            else                                                        1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngotp*nhere                    1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngotp*nhere                      1d26s18
            end if                                                      1d26s18
           end do                                                       1d26s18
          else                                                          3d15s18
           itmp1=ibcoff                                                  3d15s18
           itmp2=itmp1+nhere*nrow                                        3d15s18
           ibcoff=itmp2+nhere*nrow                                       3d15s18
           call enough('updateg. 24',bc,ibc)
           ii=i1xa(is)                                                  3d16s18
           do i12=0,nhere-1                                             3d15s18
            do i34=0,nrow-1                                             3d15s18
             iad1=itmp1+i12+nhere*i34                                         3d15s18
             bc(iad1)=bc(ii)                                            3d16s18
             ii=ii+1                                                    3d16s18
            end do                                                      3d15s18
           end do                                                       3d15s18
           nx=nhere*noc(isblkk(1,is))                                   3d23s18
           if(nx.gt.0.and.noc(isblkk(2,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(2,is)),noc(isblkk(2,is)),   3d23s18
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(2,is))),               3d23s18
     $          noc(isblkk(2,is)),0d0,bc(itmp2),nx,                     3d28s18
     d' updateg. 19')
           end if                                                       2d25s19
           do i4=0,noc(isblkk(2,is))-1                                   3d15s18
            do i3=0,noc(isblkk(1,is))-1                                  3d15s18
             iad1=itmp1+nhere*(i4+noc(isblkk(2,is))*i3)                 3d28s18
             iad2=itmp2+nhere*(i3+noc(isblkk(1,is))*i4)                 3d28s18
             do i12=0,nhere-1                                              3d15s18
              bc(iad1+i12)=bc(iad2+i12)                                 3d28s18
             end do                                                     3d16s18
            end do                                                      3d16s18
           end do                                                       3d16s18
           nx=nhere*noc(isblkk(2,is))                                   3d23s18
           if(nx.gt.0.and.noc(isblkk(1,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(1,is)),noc(isblkk(1,is)),   3d23s18
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(1,is))),               3d23s18
     $          noc(isblkk(1,is)),0d0,bc(itmp2),nx,                     3d23s18
     d' updateg. 20')
           end if                                                       2d25s19
           icol=1                                                       1d26s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           ips=0                                                        1d26s18
           mytop=ngot                                                   1d26s18
           if(nleft.gt.0)mytop=ngotp                                    1d26s18
           do i4=0,noc(isblkk(2,is))-1                                   3d16s18
            do i3=0,noc(isblkk(1,is))-1                                  3d16s18
             if(icol.gt.mytop)then                                      1d26s18
              ips=ips+1                                                 1d26s18
              if(ips.ge.nleft)then                                      1d26s18
               mytop=mytop+ngot                                         1d26s18
              else                                                      1d26s18
               mytop=mytop+ngotp                                        1d26s18
              end if                                                    1d26s18
             end if                                                     1d26s18
             do i12=0,nhere-1                                           1d26s18
              i12p=i12+il0-1                                            1d26s18
              iad1=itmp2+i12+nhere*(i4+noc(isblkk(2,is))*i3)            3d23s18
              bc(ibc(isend+ips))=bc(iad1)                               1d26s18
              ibc(isend+ips)=ibc(isend+ips)+1                           1d26s18
             end do                                                     1d26s18
             icol=icol+1                                                1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           ibcoff=itmp1                                                 1d26s18
          end if                                                        3d16s18
         end if                                                         3d16s18
         if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.      3d16s18
     $      nvirtc(isblkk(4,is)).ne.0.and.noc(isblkk(3,is)).ne.0)then   3d16s18
          call ilimts(noc(isblkk(3,is)),nvirtc(isblkk(4,is)),           3d16s18
     $        mynprocg,mynowprog,il,ih,i1s,i2e,i2s,i2e)                 3d16s18
          nhere=ih+1-il                                                 3d15s18
          nrow=noc(isblkk(1,is))*noc(isblkk(2,is))                        3d15s18
          il0=il
          if(ipass.eq.1)then                                            3d15s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           do ip=0,mynprocg-1                                           1d26s18
            if(ip.ge.nleft)then                                         1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngot*nhere                     1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngot*nhere                       1d26s18
            else                                                        1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngotp*nhere                    1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngotp*nhere                      1d26s18
            end if                                                      1d26s18
           end do                                                       1d26s18
          else                                                          3d15s18
           itmp1=ibcoff                                                  3d15s18
           itmp2=itmp1+nhere*nrow                                        3d15s18
           ibcoff=itmp2+nhere*nrow                                       3d15s18
           call enough('updateg. 25',bc,ibc)
           ii=i1xb(is)                                                  3d16s18
           do i12=0,nhere-1                                             3d15s18
            do i34=0,nrow-1                                             3d15s18
             iad1=itmp1+i12+nhere*i34                                         3d15s18
             bc(iad1)=bc(ii)                                            3d16s18
             ii=ii+1                                                    3d16s18
            end do                                                      3d15s18
           end do                                                       3d15s18
           nx=nhere*noc(isblkk(1,is))                                     3d15s18
           if(nx.gt.0.and.noc(isblkk(2,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(2,is)),noc(isblkk(2,is)),   3d16s18
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(2,is))),               3d16s18
     $          noc(isblkk(2,is)),0d0,bc(itmp2),nx,                     3d16s18
     d' updateg. 21')
           end if                                                       2d25s19
           do i12=0,nhere-1                                              3d15s18
            do i4=0,noc(isblkk(2,is))-1                                   3d15s18
             do i3=0,noc(isblkk(1,is))-1                                  3d15s18
              iad1=itmp1+i12+nhere*(i4+noc(isblkk(2,is))*i3)               3d15s18
              iad2=itmp2+i12+nhere*(i3+noc(isblkk(1,is))*i4)               3d15s18
              bc(iad1)=bc(iad2)                                         3d16s18
             end do                                                     3d16s18
            end do                                                      3d16s18
           end do                                                       3d16s18
           nx=nhere*noc(isblkk(2,is))                                    3d16s18
           if(nx.gt.0.and.noc(isblkk(1,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(1,is)),noc(isblkk(1,is)),   3d16s18
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(1,is))),               3d16s18
     $          noc(isblkk(1,is)),0d0,bc(itmp2),nx,                     3d16s18
     d' updateg. 22')
           end if                                                       2d25s19
           icol=1                                                       1d26s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           ips=0                                                        1d26s18
           mytop=ngot                                                   1d26s18
           if(nleft.gt.0)mytop=ngotp                                    1d26s18
           do i4=0,noc(isblkk(2,is))-1                                   3d16s18
            do i3=0,noc(isblkk(1,is))-1                                  3d16s18
             if(icol.gt.mytop)then                                      1d26s18
              ips=ips+1                                                 1d26s18
              if(ips.ge.nleft)then                                      1d26s18
               mytop=mytop+ngot                                         1d26s18
              else                                                      1d26s18
               mytop=mytop+ngotp                                        1d26s18
              end if                                                    1d26s18
             end if                                                     1d26s18
             do i12=0,nhere-1                                           1d26s18
              i12p=i12+il0-1                                            1d26s18
              iad1=itmp2+i12+nhere*(i4+noc(isblkk(2,is))*i3)             3d16s18
              bc(ibc(isend+ips))=bc(iad1)                               1d26s18
              ibc(isend+ips)=ibc(isend+ips)+1                           1d26s18
             end do                                                     1d26s18
             icol=icol+1                                                1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           ibcoff=itmp1                                                 1d26s18
          end if                                                        3d16s18
         end if                                                         3d16s18
         if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.      3d16s18
     $      noc(isblkk(4,is)).ne.0.and.noc(isblkk(3,is)).ne.0)then      3d16s18
          call ilimts(noc(isblkk(3,is)),noc(isblkk(4,is)),              3d16s18
     $        mynprocg,mynowprog,il,ih,i1s,i2e,i2s,i2e)                 3d16s18
          nhere=ih+1-il                                                 3d15s18
          nrow=noc(isblkk(1,is))*noc(isblkk(2,is))                        3d15s18
          il0=il
          if(ipass.eq.1)then                                            3d15s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           do ip=0,mynprocg-1                                           1d26s18
            if(ip.ge.nleft)then                                         1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngot*nhere                     1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngot*nhere                       1d26s18
            else                                                        1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngotp*nhere                    1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngotp*nhere                      1d26s18
            end if                                                      1d26s18
           end do                                                       1d26s18
          else                                                          3d15s18
           itmp1=ibcoff                                                  3d15s18
           itmp2=itmp1+nhere*nrow                                        3d15s18
           ibcoff=itmp2+nhere*nrow                                       3d15s18
           call enough('updateg. 26',bc,ibc)
           ii=i4ok(is)                                                  3d16s18
           do i12=0,nhere-1                                             3d15s18
            do i34=0,nrow-1                                             3d15s18
             iad1=itmp1+i12+nhere*i34                                         3d15s18
             bc(iad1)=bc(ii)                                            3d16s18
             ii=ii+1                                                    3d16s18
            end do                                                      3d15s18
           end do                                                       3d15s18
           nx=nhere*noc(isblkk(1,is))                                     3d15s18
           if(nx.gt.0.and.noc(isblkk(2,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(2,is)),noc(isblkk(2,is)),   3d16s18
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(2,is))),               3d16s18
     $          noc(isblkk(2,is)),0d0,bc(itmp2),nx,                     3d16s18
     d' updateg. 23')
           end if                                                       2d25s19
           do i12=0,nhere-1                                              3d15s18
            do i4=0,noc(isblkk(2,is))-1                                   3d15s18
             do i3=0,noc(isblkk(1,is))-1                                  3d15s18
              iad1=itmp1+i12+nhere*(i4+noc(isblkk(2,is))*i3)               3d15s18
              iad2=itmp2+i12+nhere*(i3+noc(isblkk(1,is))*i4)               3d15s18
              bc(iad1)=bc(iad2)                                         3d16s18
             end do                                                     3d16s18
            end do                                                      3d16s18
           end do                                                       3d16s18
           nx=nhere*noc(isblkk(2,is))                                    3d16s18
           if(nx.gt.0.and.noc(isblkk(1,is)).gt.0)then                   2d25s19
           call dgemm('n','n',nx,noc(isblkk(1,is)),noc(isblkk(1,is)),   3d16s18
     $          1d0,bc(itmp1),nx,bc(itmat(isblkk(1,is))),               3d16s18
     $          noc(isblkk(1,is)),0d0,bc(itmp2),nx,                     3d16s18
     d' updateg. 24')
           end if                                                       2d25s19
           icol=1                                                       1d26s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           ips=0                                                        1d26s18
           mytop=ngot                                                   1d26s18
           if(nleft.gt.0)mytop=ngotp                                    1d26s18
           do i4=0,noc(isblkk(2,is))-1                                   3d16s18
            do i3=0,noc(isblkk(1,is))-1                                  3d16s18
             if(icol.gt.mytop)then                                      1d26s18
              ips=ips+1                                                 1d26s18
              if(ips.ge.nleft)then                                      1d26s18
               mytop=mytop+ngot                                         1d26s18
              else                                                      1d26s18
               mytop=mytop+ngotp                                        1d26s18
              end if                                                    1d26s18
             end if                                                     1d26s18
             do i12=0,nhere-1                                           1d26s18
              i12p=i12+il0-1                                            1d26s18
              iad1=itmp2+i12+nhere*(i4+noc(isblkk(2,is))*i3)             3d16s18
              bc(ibc(isend+ips))=bc(iad1)                               1d26s18
              ibc(isend+ips)=ibc(isend+ips)+1                           1d26s18
             end do                                                     1d26s18
             icol=icol+1                                                1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           ibcoff=itmp1                                                 1d26s18
          end if                                                        3d16s18
         end if                                                         3d15s18
        end do                                                          3d15s18
       end if                                                           3d15s18
c
c     ionex
c
       if(c1.ne.0d0)then                                                3d15s18
        do is=1,nsdlk1                                                  1d26s18
         if(noc(isblk1(1,is)).ne.0.and.noc(isblk1(2,is)).ne.0.and.      1d26s18
     $     noc(isblk1(3,is)).ne.0.and.nvirtc(isblk1(4,is)).ne.0)then    1d26s18
          call ilimts(noc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  1d26s18
     $       mynowprog,il,ih,i1s,i2e,i2s,i2e)                           1d26s18
          nrow=noc(isblk1(1,is))*noc(isblk1(2,is))                      1d26s18
          nhere=ih+1-il                                                 1d26s18
          il0=il                                                        1d26s18
          if(ipass.eq.1)then                                            1d26s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           do ip=0,mynprocg-1                                           1d26s18
            if(ip.ge.nleft)then                                         1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngot*nhere                     1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngot*nhere                       1d26s18
            else                                                        1d26s18
             ibc(nsend+ip)=ibc(nsend+ip)+ngotp*nhere                    1d26s18
             ios(ip+1,is)=ios(ip+1,is)+ngotp*nhere                      1d26s18
            end if                                                      1d26s18
           end do                                                       1d26s18
          else                                                          1d26s18
           itmp1=ibcoff                                                 1d26s18
           itmp2=itmp1+nrow*nhere                                       1d26s18
           ibcoff=itmp2+nrow*nhere                                      1d26s18
           call enough('updateg. 27',bc,ibc)
           nr1=noc(isblk1(1,is))*nhere                                  1d26s18
           if(isblk1(1,is).eq.isblk1(2,is))then                         1d26s18
            nn=(noc(isblk1(1,is))*(noc(isblk1(1,is))+1))/2              1d26s18
            i1x=ionex(is)                                               1d26s18
            do i12=0,nhere-1                                            1d26s18
             do i4=0,noc(isblk1(1,is))-1                                1d26s18
              do i3=0,i4                                                1d26s18
               iad34=itmp1+i12+nhere*(i3+noc(isblk1(1,is))*i4)          1d26s18
               iad43=itmp1+i12+nhere*(i4+noc(isblk1(1,is))*i3)          1d26s18
               bc(iad34)=bc(i1x)                                        1d26s18
               bc(iad43)=bc(i1x)                                        1d26s18
               i1x=i1x+1                                                1d26s18
              end do                                                    1d26s18
             end do                                                     1d26s18
            end do                                                      1d26s18
           else                                                         1d26s18
            nn=noc(isblk1(1,is))*noc(isblk1(2,is))                      1d26s18
            i1x=ionex(is)                                               1d26s18
            do i12=0,nhere-1                                            1d26s18
             do i4=0,noc(isblk1(2,is))-1                                1d26s18
              do i3=0,noc(isblk1(1,is))-1                               1d26s18
               iad=itmp1+i12+nhere*(i3+noc(isblk1(1,is))*i4)            1d26s18
               bc(iad)=bc(i1x)                                          1d26s18
               i1x=i1x+1                                                1d26s18
              end do                                                    1d26s18
             end do                                                     1d26s18
            end do                                                      1d26s18
           end if                                                       1d26s18
           if(nr1.gt.0.and.noc(isblk1(2,is)).gt.0)then                  2d25s19
           call dgemm('n','n',nr1,noc(isblk1(2,is)),noc(isblk1(2,is)),  1d26s18
     $          1d0,bc(itmp1),nr1,bc(itmat(isblk1(2,is))),              1d26s18
     $          noc(isblk1(2,is)),0d0,bc(itmp2),nr1,                    1d26s18
     d' updateg. 25')
           end if                                                       2d25s19
           do i12=0,nhere-1                                             1d26s18
            do i4=0,noc(isblk1(2,is))-1                                 1d26s18
             do i3=0,noc(isblk1(1,is))-1                                1d26s18
              iad1=itmp2+i12+nhere*(i3+noc(isblk1(1,is))*i4)            1d26s18
              iad2=itmp1+i12+nhere*(i4+noc(isblk1(2,is))*i3)            1d26s18
              bc(iad2)=bc(iad1)                                         1d26s18
             end do                                                     1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           nr2=nhere*noc(isblk1(2,is))                                  1d26s18
           if(nr2.gt.0.and.noc(isblk1(1,is)).gt.0)then                  2d25s19
           call dgemm('n','n',nr2,noc(isblk1(1,is)),noc(isblk1(1,is)),  1d26s18
     $          1d0,bc(itmp1),nr2,bc(itmat(isblk1(1,is))),              1d26s18
     $          noc(isblk1(1,is)),0d0,bc(itmp2),nr2,                    1d26s18
     d' updateg. 26')
           end if                                                       2d25s19
           do i12=0,nhere-1                                             1d26s18
            do i4=0,noc(isblk1(2,is))-1                                 1d26s18
             do i3=0,noc(isblk1(1,is))-1                                1d26s18
              iad1=itmp1+i12+nhere*(i3+noc(isblk1(1,is))*i4)            1d26s18
              iad2=itmp2+i12+nhere*(i4+noc(isblk1(2,is))*i3)            1d26s18
              bc(iad1)=bc(iad2)                                         1d26s18
             end do                                                     1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           icol=1                                                       1d26s18
           ngot=nrow/mynprocg                                           1d26s18
           ngotp=ngot+1                                                 1d26s18
           nleft=nrow-ngot*mynprocg                                     1d26s18
           ips=0                                                        1d26s18
           mytop=ngot                                                   1d26s18
           if(nleft.gt.0)mytop=ngotp                                    1d26s18
           do i4=0,noc(isblk1(2,is))-1                                  1d26s18
            do i3=0,noc(isblk1(1,is))-1                                 1d26s18
             if(icol.gt.mytop)then                                      1d26s18
              ips=ips+1                                                 1d26s18
              if(ips.ge.nleft)then                                      1d26s18
               mytop=mytop+ngot                                         1d26s18
              else                                                      1d26s18
               mytop=mytop+ngotp                                        1d26s18
              end if                                                    1d26s18
             end if                                                     1d26s18
             do i12=0,nhere-1                                           1d26s18
              i12p=i12+il0-1                                            1d26s18
              i4x=i12p/noc(isblk1(3,is))                                1d26s18
              i3u=i12p-i4x*noc(isblk1(3,is))                            1d26s18
              iad1=itmp2+i12+nhere*(i4+noc(isblk1(2,is))*i3)            1d26s18
              bc(ibc(isend+ips))=bc(iad1)                               1d26s18
              ibc(isend+ips)=ibc(isend+ips)+1                           1d26s18
             end do                                                     1d26s18
             icol=icol+1                                                1d26s18
            end do                                                      1d26s18
           end do                                                       1d26s18
           ibcoff=itmp2                                                 1d26s18
          end if                                                        1d26s18
         end if                                                         1d26s18
        end do                                                          1d26s18
       end if                                                           1d26s18
c
c     4o
c
       do is=1,nsdlk
        if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.         1d10s18
     $     noc(isblk(3,is)).ne.0.and.noc(isblk(4,is)).ne.0)then         1d10s18
         call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,        1d10s18
     $       mynowprog,il,ih,i1s,i2e,i2s,i2e)                           1d10s18
         nrow=noc(isblk(1,is))*noc(isblk(2,is))                         1d10s18
         nhere=ih+1-il                                                  1d10s18
         il0=il                                                         1d10s18
         if(ipass.eq.1)then
          ngot=nrow/mynprocg                                            1d10s18
          ngotp=ngot+1                                                  1d10s18
          nleft=nrow-ngot*mynprocg                                      1d10s18
          do ip=0,mynprocg-1                                            1d10s18
           if(ip.ge.nleft)then                                          1d10s18
            ibc(nsend+ip)=ibc(nsend+ip)+ngot*nhere                      1d10s18
            ios(ip+1,is)=ios(ip+1,is)+ngot*nhere
           else                                                         1d10s18
            ibc(nsend+ip)=ibc(nsend+ip)+ngotp*nhere                     1d10s18
            ios(ip+1,is)=ios(ip+1,is)+ngotp*nhere
           end if                                                       1d10s18
          end do                                                        1d10s18
         else
          itmp1=ibcoff                                                  1d10s18
          itmp2=itmp1+nrow*nhere                                        1d10s18
          ibcoff=itmp2+nrow*nhere                                       1d10s18
          call enough('updateg. 28',bc,ibc)
          nr1=noc(isblk(1,is))*nhere                                    1d10s18
          if(isblk(1,is).eq.isblk(2,is))then                            1d10s18
           nn=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 1d10s18
           i4o=ioooo(is)                                                1d10s18
           do i12=0,nhere-1                                             1d10s18
            do i4=0,noc(isblk(1,is))-1                                  1d10s18
             do i3=0,i4                                                 1d10s18
              iad34=itmp1+i12+nhere*(i3+noc(isblk(1,is))*i4)             1d10s18
              iad43=itmp1+i12+nhere*(i4+noc(isblk(1,is))*i3)             1d10s18
              bc(iad34)=bc(i4o)                                         1d10s18
              bc(iad43)=bc(i4o)                                         1d10s18
              i4o=i4o+1                                                 1d10s18
             end do                                                     1d10s18
            end do                                                      1d10s18
           end do                                                       1d10s18
          else                                                          1d10s18
           nn=noc(isblk(1,is))*noc(isblk(2,is))                         1d10s18
           i4o=ioooo(is)                                                1d10s18
           do i12=0,nhere-1                                             1d10s18
            do i4=0,noc(isblk(2,is))-1                                  1d10s18
             do i3=0,noc(isblk(1,is))-1                                 1d10s18
              iad=itmp1+i12+nhere*(i3+noc(isblk(1,is))*i4)              1d10s18
              bc(iad)=bc(i4o)                                           1d10s18
              i4o=i4o+1                                                 1d10s18
             end do                                                     1d10s18
            end do                                                      1d10s18
           end do                                                       1d10s18
          end if                                                        1d10s18
          if(nr1.gt.0.and.noc(isblk(2,is)).gt.0)then                    2d25s19
          call dgemm('n','n',nr1,noc(isblk(2,is)),noc(isblk(2,is)),1d0, 1d10s18
     $        bc(itmp1),nr1,bc(itmat(isblk(2,is))),noc(isblk(2,is)),0d0,1d10s18
     $         bc(itmp2),nr1,                                           1d10s18
     d' updateg. 27')
          end if                                                        2d25s19
          do i12=0,nhere-1                                              1d10s18
           do i4=0,noc(isblk(2,is))-1                                   1d10s18
            do i3=0,noc(isblk(1,is))-1                                  1d10s18
             iad1=itmp2+i12+nhere*(i3+noc(isblk(1,is))*i4)              1d10s18
             iad2=itmp1+i12+nhere*(i4+noc(isblk(2,is))*i3)              1d10s18
             bc(iad2)=bc(iad1)                                          1d10s18
            end do                                                      1d10s18
           end do                                                       1d10s18
          end do                                                        1d10s18
          nr2=nhere*noc(isblk(2,is))                                    1d10s18
          if(nr2.gt.0.and.noc(isblk(1,is)).gt.0)then                    2d25s19
          call dgemm('n','n',nr2,noc(isblk(1,is)),noc(isblk(1,is)),1d0, 1d10s18
     $         bc(itmp1),nr2,bc(itmat(isblk(1,is))),noc(isblk(1,is)),   1d10s18
     $         0d0,bc(itmp2),nr2,                                       1d10s18
     d' updateg. 28')
          end if                                                        2d25s19
          icol=1                                                        1d10s18
          ngot=nrow/mynprocg                                            1d10s18
          ngotp=ngot+1                                                  1d10s18
          nleft=nrow-ngot*mynprocg                                      1d10s18
          ips=0                                                         1d10s18
          mytop=ngot                                                    1d10s18
          if(nleft.gt.0)mytop=ngotp                                     1d10s18
          do i4=0,noc(isblk(2,is))-1                                    1d10s18
           do i3=0,noc(isblk(1,is))-1                                   1d10s18
            if(icol.gt.mytop)then                                       1d10s18
             ips=ips+1                                                  1d10s18
             if(ips.ge.nleft)then                                       1d10s18
              mytop=mytop+ngot                                          1d10s18
             else                                                       1d10s18
              mytop=mytop+ngotp                                         1d10s18
             end if                                                     1d10s18
            end if                                                      1d10s18
            do i12=0,nhere-1                                            1d10s18
             i12p=i12+il0-1
             i4x=i12p/noc(isblk(3,is))
             i3u=i12p-i4x*noc(isblk(3,is))
             iad1=itmp2+i12+nhere*(i4+noc(isblk(2,is))*i3)              1d10s18
             bc(ibc(isend+ips))=bc(iad1)                                1d10s18
             ibc(isend+ips)=ibc(isend+ips)+1                            1d10s18
            end do                                                      1d10s18
            icol=icol+1                                                 1d11s18
           end do                                                       1d10s18
          end do                                                        1d10s18
          ibcoff=itmp2                                                  1d10s18
         end if
        end if                                                          1d10s18
       end do
       if(ipass.eq.2)then                                               1d10s18
        jbufs=0                                                         1d10s18
        do ip=0,mynprocg-1                                              1d10s18
         ibc(isend+ip)=jbufs                                            1d10s18
         jbufs=jbufs+ibc(nsend+ip)                                      1d10s18
        end do                                                          1d10s18
        call dws_all2allvb8(bc(ibufs),ibc(nsend),ibc(isend),bc(ibufr),   1d10s18
     $       ibc(nrecv),ibc(irecv))                                     1d10s18
        jbufr=ibufr                                                     1d10s18
        do ip=0,mynprocg-1                                              1d10s18
         ibc(irecv+ip)=jbufr                                            1d10s18
         jbufr=jbufr+ibc(nrecv+ip)                                      1d10s18
        end do                                                          1d10s18
       end if                                                           1d10s18
c
c     receiving
c
       if(ipass.eq.1)then                                               1d11s18
        if(c3.ne.0d0)then                                               4d11s18
         me2mea=0
         me2meb=0
         me2mec=0
         do is3=1,nsdlk1                                                4d11s18
          if(nvirtc(isblk1(1,is3)).ne.0.and.nvirtc(isblk1(2,is3)).ne.0. 4d11s18
     $         and.noc(isblk1(3,is3)).ne.0.and.                         4d11s18
     $         nvirtc(isblk1(4,is3)).ne.0)then                          4d11s18
           ncol1=noc(isblk1(1,is3))*noc(isblk1(2,is3))                  4d11s18
           nhere=ncol1/mynprocg                                          4d10s18
           nleft=ncol1-nhere*mynprocg                                    4d10s18
           nhere0=max(1,nhere)                                          4d12s18
           nhere=nhere+1                                                4d10s18
           nherea=nhere                                                 4d10s18
           nlefta=nleft                                                 4d10s18
           nhere0a=nhere0                                               4d10s18
           ncol2=nvirtc(isblk1(1,is3))*noc(isblk1(2,is3))
           nhere=ncol2/mynprocg                                          4d10s18
           nleft=ncol2-nhere*mynprocg                                    4d10s18
           nhere0=max(1,nhere)                                          4d12s18
           nhere=nhere+1                                                4d10s18
           nhereb=nhere                                                 4d10s18
           nleftb=nleft                                                 4d10s18
           nhere0b=nhere0                                               4d10s18
           if(isblk1(1,is3).ne.isblk1(2,is3))then                       4d10s18
            ncol3=noc(isblk1(1,is3))*nvirtc(isblk1(2,is3))              4d10s18
            nhere=ncol3/mynprocg                                          4d10s18
            nleft=ncol3-nhere*mynprocg                                    4d10s18
            nhere0=max(1,nhere)                                         4d12s18
            nhere=nhere+1                                                4d10s18
            nherec=nhere                                                 4d10s18
            nleftc=nleft                                                 4d10s18
            nhere0c=nhere0                                               4d10s18
            if(mynowprog.ge.nleftc)nherec=nherec-1                      4d12s18
           end if                                                       4d11s18
           if(mynowprog.ge.nlefta)nherea=nherea-1                       4d12s18
           if(mynowprog.ge.nleftb)nhereb=nhereb-1                       4d12s18
           do ip=0,mynprocg-1                                           4d11s18
            call ilimts(noc(isblk1(3,is3)),noc(isblk1(4,is3)),mynprocg, 4d11s18
     $           ip,il,ih,i1s,i1e,i2s,i2e)                              4d11s18
            noo=ih+1-il                                                 4d11s18
            call ilimts(noc(isblk1(3,is3)),nvirtc(isblk1(4,is3)),       4d11s18
     $           mynprocg,ip,il,ih,i1s,i1e,i2s,i2e)                     4d11s18
            nov=ih+1-il                                                 4d11s18
            nsum=noo+nov                                                4d11s18
            if(ip.eq.mynowprog)then
             me2mea=me2mea+nherea*nsum
             me2meb=me2meb+nhereb*nsum
            end if
            ibc(nrecv+ip)=ibc(nrecv+ip)+nsum*nherea                     4d11s18
            ibc(nrecv+ip)=ibc(nrecv+ip)+nsum*nhereb                     4d11s18
            if(isblk1(1,is3).ne.isblk1(2,is3))then                      4d12s18
             ibc(nrecv+ip)=ibc(nrecv+ip)+nsum*nherec                    4d11s18
             if(ip.eq.mynowprog)then
              me2mec=me2mec+nherec*nsum
             end if
            end if                                                      4d11s18
           end do                                                       4d11s18
          end if                                                        4d11s18
         end do                                                         4d11s18
        end if                                                          4d11s18
        if(c2.ne.0d0)then                                               3d16s18
c
c     jmats
c
         do is=1,nsdlk                                                   3d16s18
          if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.       3d16s18
     $       nvirtc(isblk(3,is)).ne.0.and.nvirtc(isblk(4,is)).ne.0)then 3d16s18
           call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,      3d16s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d27s18
           ncol=ih+1-il                                                 1d27s18
           il0=il                                                       1d27s18
           do ip=0,mynprocg-1                                           1d27s18
            call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),        3d16s18
     $           mynprocg,ip,il,ih,i1s,i1e,i2s,i2e)                     3d16s18
            nhere34=ih+1-il                                             1d27s18
            ibc(nrecv+ip)=ibc(nrecv+ip)+nhere34*ncol                    1d27s18
            ior(ip+1,is)=ior(ip+1,is)+nhere34*ncol                      1d27s18
           end do                                                       1d27s18
          end if                                                        1d27s18
         end do                                                          3d16s18
c
c     kmats and entourage
c
         do is=1,nsdlkk                                                 3d16s18
          if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.       3d16s18
     $      nvirtc(isblkk(3,is)).ne.0.and.nvirtc(isblkk(4,is)).ne.0)then 3d16s18
           call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,      3d16s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d27s18
           ncol=ih+1-il                                                 1d27s18
           il0=il                                                       1d27s18
           do ip=0,mynprocg-1                                           1d27s18
            call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),        3d16s18
     $           mynprocg,ip,il,ih,i1s,i1e,i2s,i2e)                     3d16s18
            nhere34=ih+1-il                                             1d27s18
            ibc(nrecv+ip)=ibc(nrecv+ip)+nhere34*ncol                    1d27s18
            ior(ip+1,is)=ior(ip+1,is)+nhere34*ncol                      1d27s18
           end do                                                       1d27s18
          end if                                                        1d27s18
          if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.       3d16s18
     $       nvirtc(isblkk(3,is)).ne.0.and.noc(isblkk(4,is)).ne.0)then  3d16s18
           call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,      3d16s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d27s18
           ncol=ih+1-il                                                 1d27s18
           il0=il                                                       1d27s18
           do ip=0,mynprocg-1                                           1d27s18
            call ilimts(noc(isblkk(4,is)),nvirtc(isblkk(3,is)),         3d16s18
     $           mynprocg,ip,il,ih,i1s,i1e,i2s,i2e)                     3d16s18
            nhere34=ih+1-il                                             1d27s18
            ibc(nrecv+ip)=ibc(nrecv+ip)+nhere34*ncol                    1d27s18
            ior(ip+1,is)=ior(ip+1,is)+nhere34*ncol                      1d27s18
           end do                                                       1d27s18
          end if                                                        1d27s18
          if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.       3d16s18
     $       nvirtc(isblkk(4,is)).ne.0.and.noc(isblkk(3,is)).ne.0)then  3d16s18
           call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,      3d16s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d27s18
           ncol=ih+1-il                                                 1d27s18
           il0=il                                                       1d27s18
           do ip=0,mynprocg-1                                           1d27s18
            call ilimts(noc(isblkk(3,is)),nvirtc(isblkk(4,is)),         3d16s18
     $           mynprocg,ip,il,ih,i1s,i1e,i2s,i2e)                     3d16s18
            nhere34=ih+1-il                                             1d27s18
            ibc(nrecv+ip)=ibc(nrecv+ip)+nhere34*ncol                    1d27s18
            ior(ip+1,is)=ior(ip+1,is)+nhere34*ncol                      1d27s18
           end do                                                       1d27s18
          end if                                                        1d27s18
          if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.       3d16s18
     $       noc(isblkk(4,is)).ne.0.and.noc(isblkk(3,is)).ne.0)then     3d16s18
           call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,      3d16s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d27s18
           ncol=ih+1-il                                                 1d27s18
           il0=il                                                       1d27s18
           do ip=0,mynprocg-1                                           1d27s18
            call ilimts(noc(isblkk(3,is)),noc(isblkk(4,is)),            3d16s18
     $           mynprocg,ip,il,ih,i1s,i1e,i2s,i2e)                     3d16s18
            nhere34=ih+1-il                                             1d27s18
            ibc(nrecv+ip)=ibc(nrecv+ip)+nhere34*ncol                    1d27s18
            ior(ip+1,is)=ior(ip+1,is)+nhere34*ncol                      1d27s18
           end do                                                       1d27s18
          end if                                                        1d27s18
         end do                                                          3d16s18
        end if                                                          3d16s18
c
c     onex
c
        if(c1.ne.0d0)then                                               1d27s18
         do is=1,nsdlk1                                                 1d27s18
          if(noc(isblk1(1,is)).ne.0.and.noc(isblk1(2,is)).ne.0.and.     1d27s18
     $       noc(isblk1(3,is)).ne.0.and.nvirtc(isblk1(4,is)).ne.0)then  1d27s18
           call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,    1d27s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d27s18
           ncol=ih+1-il                                                 1d27s18
           il0=il                                                       1d27s18
           do ip=0,mynprocg-1                                           1d27s18
            call ilimts(noc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,1d27s18
     $          ip,il,ih,i1s,i1e,i2s,i2e)                               1d27s18
            nhere34=ih+1-il                                             1d27s18
            ibc(nrecv+ip)=ibc(nrecv+ip)+nhere34*ncol                    1d27s18
            ior(ip+1,is)=ior(ip+1,is)+nhere34*ncol                      1d27s18
           end do                                                       1d27s18
          end if                                                        1d27s18
         end do                                                         1d27s18
        end if                                                          1d27s18
c
c     4o
c
        do is=1,nsdlk                                                    1d10s18
         if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.         1d10s18
     $       noc(isblk(3,is)).ne.0.and.noc(isblk(4,is)).ne.0)then       1d10s18
          call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,        1d10s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d10s18
          ncol=ih+1-il                                                   1d10s18
          il0=il
          do ip=0,mynprocg-1                                            1d10s18
           call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,      1d10s18
     $          ip,il,ih,i1s,i1e,i2s,i2e)                               1d10s18
           nhere34=ih+1-il                                                1d10s18
           ibc(nrecv+ip)=ibc(nrecv+ip)+nhere34*ncol                     1d10s18
           ior(ip+1,is)=ior(ip+1,is)+nhere34*ncol
          end do                                                        1d10s18
         end if                                                         1d11s18
        end do                                                          1d11s18
       else                                                             1d11s18
        do ip=0,mynprocg-1                                               1d11s18
         if(c3.ne.0d0)then                                              4d11s18
          do is3=1,nsdlk1                                               4d16s18
           if(nvirtc(isblk1(1,is3)).ne.0.and.nvirtc(isblk1(2,is3)).ne.0.  4d10s18
     $     and.noc(isblk1(3,is3)).ne.0.and.nvirtc(isblk1(4,is3)).ne.0)  4d10s18
     $        then                                                      4d10s18
            if(noc(isblk1(4,is3)).ne.0)then                             4d13s18
             call ilimts(nvirtc(isblk1(1,is3)),noc(isblk1(2,is3)),      4d13s18
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)              4d13s18
             nhere=ih+1-il                                                4d13s18
             call ilimts(noc(isblk1(3,is3)),noc(isblk1(4,is3)),           4d13s18
     $          mynprocg,ip,jl,jh,j1s,j1e,j2s,j2e)
             jl=jl-1                                                    4d13s18
             mhere=jh-jl                                                4d13s18
             do i12=0,nhere-1                                           4d13s18
              j10=j1s                                                   4d13s18
              j1n=noc(isblk1(3,is3))                                    4d16s18
              do j2=j2s,j2e                                             4d16s18
               if(j2.eq.j2e)j1n=j1e                                     4d16s18
               j2m=j2-1                                                 4d13s18
               do j1=j10,j1n                                            4d13s18
                j1m=j1-1                                                4d13s18
                iad1=iok3x(1,is3)+i12+nhere*(j1m                        4d16s18
     $               +noc(isblk1(3,is3))*j2m)                           4d13s18
                bc(iad1)=bc(ibc(irecv+ip))                               4d13s18
                ibc(irecv+ip)=ibc(irecv+ip)+1                            4d13s18
               end do                                                    4d13s18
               j10=1                                                    4d13s18
              end do                                                    4d13s18
             end do                                                     4d13s18
             call ilimts(noc(isblk1(1,is3)),noc(isblk1(2,is3)),         4d13s18
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)              4d13s18
             nhere=ih+1-il                                                4d13s18
             do i12=0,nhere-1                                           4d13s18
              j10=j1s                                                   4d13s18
              j1n=noc(isblk1(3,is3))                                    4d13s18
              do j2=j2s,j2e                                             4d13s18
               if(j2.eq.j2e)j1n=j1e                                     4d13s18
               j2m=j2-1                                                 4d13s18
               do j1=j10,j1n                                            4d13s18
                j1m=j1-1                                                4d16s18
                iad1=iok3x(2,is3)+i12+nhere*(j1m+noc(isblk1(3,is3))*j2m)4d16s18
                bc(iad1)=bc(ibc(irecv+ip))                              4d13s18
                ibc(irecv+ip)=ibc(irecv+ip)+1                             4d13s18
               end do                                                   4d13s18
               j10=1                                                    4d13s18
              end do                                                    4d13s18
             end do
             if(isblk1(1,is3).ne.isblk1(2,is3))then                     4d13s18
              call ilimts(noc(isblk1(1,is3)),nvirtc(isblk1(2,is3)),     4d13s18
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)              4d13s18
              nhere=ih+1-il                                                4d13s18
              do i12=0,nhere-1                                           4d13s18
               j10=j1s                                                   4d13s18
               j1n=noc(isblk1(3,is3))                                    4d13s18
               do j2=j2s,j2e                                             4d13s18
                j2m=j2-1                                                4d13s18
                if(j2.eq.j2e)j1n=j1e                                     4d13s18
                do j1=j10,j1n                                            4d13s18
                 j1m=j1-1                                               4d16s18
                 iad1=iok3x(3,is3)+i12+nhere*(j1m+noc(isblk1(3,is3))    4d16s18
     $                *j2m)                                             4d16s18
                 bc(iad1)=bc(ibc(irecv+ip))                              4d13s18
                 ibc(irecv+ip)=ibc(irecv+ip)+1                             4d13s18
                end do                                                   4d13s18
                j10=1                                                    4d13s18
               end do                                                    4d13s18
              end do
             end if                                                     4d13s18
            end if                                                      4d13s18
            call ilimts(nvirtc(isblk1(1,is3)),noc(isblk1(2,is3)),       4d13s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               4d13s18
            nhere=ih+1-il                                                 4d13s18
            call ilimts(noc(isblk1(3,is3)),nvirtc(isblk1(4,is3)),       4d13s18
     $         mynprocg,ip,jl,jh,j1s,j1e,j2s,j2e)
            jl=jl-1                                                     4d13s18
            mhere=jh-jl                                                 4d13s18
            do i12=0,nhere-1                                            4d13s18
             j10=j1s                                                    4d13s18
             j1n=noc(isblk1(3,is3))                                     4d13s18
             do j2=j2s,j2e                                              4d13s18
              if(j2.eq.j2e)j1n=j1e                                       4d13s18
              j2m=j2-1                                                  4d13s18
              j2p=j2m+noc(isblk1(4,is3))
              do j1=j10,j1n                                             4d13s18
               j1m=j1-1                                                 4d13s18
               iad1=iok3x(1,is3)+i12+nhere*(j1m                         4d16s18
     $              +noc(isblk1(3,is3))*j2p)                            4d13s18
               bc(iad1)=bc(ibc(irecv+ip))                                4d13s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                             4d13s18
              end do                                                     4d13s18
              j10=1                                                     4d13s18
             end do                                                     4d13s18
            end do                                                      4d13s18
            call ilimts(noc(isblk1(1,is3)),noc(isblk1(2,is3)),          4d13s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               4d13s18
            nhere=ih+1-il                                                 4d13s18
            do i12=0,nhere-1                                            4d13s18
             j10=j1s                                                    4d13s18
             j1n=noc(isblk1(3,is3))                                     4d13s18
             do j2=j2s,j2e                                              4d13s18
              if(j2.eq.j2e)j1n=j1e                                      4d13s18
              j2m=j2-1                                                  4d13s18
              j2p=j2m+noc(isblk1(4,is3))                                4d13s18
              do j1=j10,j1n                                             4d13s18
               j1m=j1-1                                                 4d16s18
               iad1=iok3x(2,is3)+i12+nhere*(j1m+noc(isblk1(3,is3))*j2p) 4d16s18
               bc(iad1)=bc(ibc(irecv+ip))                               4d13s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                              4d13s18
              end do                                                    4d13s18
              j10=1                                                     4d13s18
             end do                                                     4d13s18
            end do
            if(isblk1(1,is3).ne.isblk1(2,is3))then                      4d13s18
             call ilimts(noc(isblk1(1,is3)),nvirtc(isblk1(2,is3)),      4d13s18
     $          mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               4d13s18
             nhere=ih+1-il                                                 4d13s18
             do i12=0,nhere-1                                            4d13s18
              j10=j1s                                                    4d13s18
              j1n=noc(isblk1(3,is3))                                     4d13s18
              do j2=j2s,j2e                                              4d13s18
               if(j2.eq.j2e)j1n=j1e                                      4d13s18
               j2m=j2-1                                                 4d13s18
               j2p=j2m+noc(isblk1(4,is3))                               4d13s18
               do j1=j10,j1n                                             4d13s18
                j1m=j1-1                                                4d16s18
                iad1=iok3x(3,is3)+i12+nhere*(j1m+noc(isblk1(3,is3))*j2p)4d16s18
                bc(iad1)=bc(ibc(irecv+ip))                               4d13s18
                ibc(irecv+ip)=ibc(irecv+ip)+1                              4d13s18
               end do                                                    4d13s18
               j10=1                                                     4d13s18
              end do                                                     4d13s18
             end do
            end if                                                      4d13s18
           end if                                                       4d11s18
          end do                                                        4d11s18
         end if                                                         4d11s18
         if(c2.ne.0d0)then                                              3d16s18
c
c     jmats
c
          do is=1,nsdlk                                                 3d16s18
           if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.      3d16s18
     $        nvirtc(isblk(3,is)).ne.0.and.nvirtc(isblk(4,is)).ne.0)then3d16s18
            call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),        3d16s18
     $          mynprocg,ip,ih,il,i1s,i1e,i2s,i2e)                      3d16s18
            call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,     3d16s18
     $         mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                     1d27s18
            ncol=ih+1-il                                                1d27s18
            do i34=0,ncol-1                                             1d27s18
             i10=i1s                                                    1d27s18
             i1n=nvirtc(isblk(3,is))                                    3d16s18
             do i2=i2s,i2e                                              1d27s18
              i2m=i2-1+noc(isblk(4,is))                                 3d16s18
              if(i2.eq.i2e)i1n=i1e                                      1d27s18
              do i1=i10,i1n                                             1d27s18
               i1m=i1-1+noc(isblk(3,is))                                3d16s18
               iad=ijnn(is)+i34+ncol*(i1m+nbasdwsc(isblk(3,is))*i2m)    3d16s18
               bc(iad)=bc(ibc(irecv+ip))                                1d27s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                            1d27s18
              end do                                                    1d27s18
              i10=1                                                     1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
           end if                                                       1d27s18
          end do                                                        1d27s18
c
c     kmats and entourage
c
          do is=1,nsdlkk                                                 3d16s18
           if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.    3d16s18
     $      nvirtc(isblkk(3,is)).ne.0.and.nvirtc(isblkk(4,is)).ne.0)then3d16s18
            call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),        3d16s18
     $          mynprocg,ip,ih,il,i1s,i1e,i2s,i2e)                      3d16s18
            call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,     3d16s18
     $         mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                     1d27s18
            ncol=ih+1-il                                                1d27s18
            do i34=0,ncol-1                                             1d27s18
             i10=i1s                                                    1d27s18
             i1n=nvirtc(isblkk(3,is))                                    3d16s18
             do i2=i2s,i2e                                              1d27s18
              i2m=i2-1+noc(isblkk(4,is))                                3d16s18
              if(i2.eq.i2e)i1n=i1e                                      1d27s18
              do i1=i10,i1n                                             1d27s18
               i1m=i1-1+noc(isblkk(3,is))                               3d16s18
               iad=iokx(is)+i34+ncol*(i1m+nbasdwsc(isblkk(3,is))*i2m)    3d16s18
               bc(iad)=bc(ibc(irecv+ip))                                1d27s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                            1d27s18
              end do                                                    1d27s18
              i10=1                                                     1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
           end if                                                       1d27s18
           if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.    3d16s18
     $      nvirtc(isblkk(3,is)).ne.0.and.noc(isblkk(4,is)).ne.0)then   3d16s18
            call ilimts(noc(isblkk(4,is)),nvirtc(isblkk(3,is)),         3d16s18
     $          mynprocg,ip,ih,il,i1s,i1e,i2s,i2e)                      3d16s18
            call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,     3d16s18
     $         mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                     1d27s18
            ncol=ih+1-il                                                1d27s18
            do i34=0,ncol-1                                             1d27s18
             i10=i1s                                                    1d27s18
             i1n=noc(isblkk(4,is))                                      3d16s18
             do i2=i2s,i2e                                              1d27s18
              i2m=i2-1+noc(isblkk(3,is))                                3d16s18
              if(i2.eq.i2e)i1n=i1e                                      1d27s18
              do i1=i10,i1n                                             1d27s18
               i1m=i1-1                                                 3d16s18
               iad=iokx(is)+i34+ncol*(i2m+nbasdwsc(isblkk(3,is))*i1m)    3d16s18
               bc(iad)=bc(ibc(irecv+ip))                                1d27s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                            1d27s18
              end do                                                    1d27s18
              i10=1                                                     1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
           end if                                                       1d27s18
           if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.    3d16s18
     $      nvirtc(isblkk(4,is)).ne.0.and.noc(isblkk(3,is)).ne.0)then   3d16s18
            call ilimts(noc(isblkk(3,is)),nvirtc(isblkk(4,is)),         3d16s18
     $          mynprocg,ip,ih,il,i1s,i1e,i2s,i2e)                      3d16s18
            call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,     3d16s18
     $         mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                     1d27s18
            ncol=ih+1-il                                                1d27s18
            do i34=0,ncol-1                                             1d27s18
             i10=i1s                                                    1d27s18
             i1n=noc(isblkk(3,is))                                      3d16s18
             do i2=i2s,i2e                                              1d27s18
              i2m=i2-1+noc(isblkk(4,is))                                3d16s18
              if(i2.eq.i2e)i1n=i1e                                      1d27s18
              do i1=i10,i1n                                             1d27s18
               i1m=i1-1                                                 3d16s18
               iad=iokx(is)+i34+ncol*(i1m+nbasdwsc(isblkk(3,is))*i2m)    3d16s18
               bc(iad)=bc(ibc(irecv+ip))                                1d27s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                            1d27s18
              end do                                                    1d27s18
              i10=1                                                     1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
           end if                                                       1d27s18
           if(noc(isblkk(1,is)).ne.0.and.noc(isblkk(2,is)).ne.0.and.    3d16s18
     $      noc(isblkk(4,is)).ne.0.and.noc(isblkk(3,is)).ne.0)then      3d16s18
            call ilimts(noc(isblkk(3,is)),noc(isblkk(4,is)),            3d16s18
     $          mynprocg,ip,ih,il,i1s,i1e,i2s,i2e)                      3d16s18
            call ilimts(noc(isblkk(1,is)),noc(isblkk(2,is)),mynprocg,     3d16s18
     $         mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                     1d27s18
            ncol=ih+1-il                                                1d27s18
            do i34=0,ncol-1                                             1d27s18
             i10=i1s                                                    1d27s18
             i1n=noc(isblkk(3,is))                                      3d16s18
             do i2=i2s,i2e                                              1d27s18
              i2m=i2-1                                                  3d16s18
              if(i2.eq.i2e)i1n=i1e                                      1d27s18
              do i1=i10,i1n                                             1d27s18
               i1m=i1-1                                                 3d16s18
               iad=iokx(is)+i34+ncol*(i1m+nbasdwsc(isblkk(3,is))*i2m)    3d16s18
               bc(iad)=bc(ibc(irecv+ip))                                1d27s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                            1d27s18
              end do                                                    1d27s18
              i10=1                                                     1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
           end if                                                       1d27s18
          end do                                                        1d27s18
         end if                                                         3d16s18
c
c     onex
c
         if(c1.ne.0d0)then                                              1d27s18
          do is=1,nsdlk1                                                1d27s18
           if(noc(isblk1(1,is)).ne.0.and.noc(isblk1(2,is)).ne.0.and.    1d27s18
     $        noc(isblk1(3,is)).ne.0.and.nvirtc(isblk1(4,is)).ne.0)then 1d27s18
            call ilimts(noc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,1d27s18
     $        ip,ih,il,i1s,i1e,i2s,i2e)                                 1d27s18
            call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,   1d27s18
     $         mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                     1d27s18
            ncol=ih+1-il                                                1d27s18
            do i34=0,ncol-1                                             1d27s18
             i34p=i34+il-1                                              1d27s18
             i2p=i34p/noc(isblk1(1,is))                                 1d27s18
             i1p=i34p-i2p*noc(isblk1(1,is))                             1d27s18
             i10=i1s                                                    1d27s18
             i1n=noc(isblk1(3,is))                                      1d27s18
             do i2=i2s,i2e                                              1d27s18
              i2m=i2-1                                                  1d27s18
              if(i2.eq.i2e)i1n=i1e                                      1d27s18
              do i1=i10,i1n                                             1d27s18
               i1m=i1-1                                                 1d27s18
               iad=iovnn(is)+i34+ncol*(i1m+noc(isblk1(3,is))*i2m)       1d27s18
               bc(iad)=bc(ibc(irecv+ip))                                1d27s18
               ibc(irecv+ip)=ibc(irecv+ip)+1                            1d27s18
              end do                                                    1d27s18
              i10=1                                                     1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
           end if                                                       1d27s18
          end do                                                        1d27s18
         end if                                                         1d27s18
c
c     4o
c
         do is=1,nsdlk                                                    1d10s18
          if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.         1d10s18
     $       noc(isblk(3,is)).ne.0.and.noc(isblk(4,is)).ne.0)then       1d10s18
           call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,         1d11s18
     $       ip,ih,il,i1s,i1e,i2s,i2e)                                  1d11s18
           call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,        1d10s18
     $        mynowprog,il,ih,i1s0,i1e0,i2s0,i2e0)                      1d10s18
           ncol=ih+1-il                                                   1d10s18
           do i34=0,ncol-1                                              1d10s18
            i34p=i34+il-1
            i2p=i34p/noc(isblk(1,is))
            i1p=i34p-i2p*noc(isblk(1,is))
            i10=i1s                                                      1d10s18
            i1n=noc(isblk(3,is))                                         1d10s18
            do i2=i2s,i2e                                                1d10s18
             i2m=i2-1                                                    1d10s18
             if(i2.eq.i2e)i1n=i1e                                        1d10s18
             do i1=i10,i1n                                               1d10s18
              i1m=i1-1                                                   1d10s18
              iad=ioonn(is)+i34+ncol*(i1m+noc(isblk(3,is))*i2m)             1d10s18
              bc(iad)=bc(ibc(irecv+ip))                                 1d10s18
              ibc(irecv+ip)=ibc(irecv+ip)+1                             1d10s18
             end do                                                     1d10s18
             i10=1                                                       1d10s18
            end do                                                      1d10s18
           end do                                                       1d10s18
          end if                                                        1d11s18
         end do                                                         1d11s18
        end do
c
       end if                                                           1d11s18
       if(ipass.eq.1)then                                               1d10s18
        ntots=0                                                         1d10s18
        ntotr=0                                                         1d10s18
        do ip=0,mynprocg-1                                              1d10s18
         ntots=ntots+ibc(nsend+ip)                                      1d10s18
         ntotr=ntotr+ibc(nrecv+ip)                                      1d10s18
        end do                                                          1d10s18
        ibufs=ibcoff                                                    1d10s18
        ibufr=ibufs+ntots                                               1d10s18
        ibcoff=ibufr+ntotr                                              1d10s18
        call enough('updateg. 29',bc,ibc)
        jbufs=ibufs                                                     1d10s18
        jbufr=0                                                         1d10s18
        do ip=0,mynprocg-1                                              1d10s18
         ibc(isend+ip)=jbufs                                             1d10s18
         jbufs=jbufs+ibc(nsend+ip)                                      1d10s18
         ibc(irecv+ip)=jbufr                                            1d10s18
         jbufr=jbufr+ibc(nrecv+ip)                                      1d10s18
        end do                                                          1d10s18
       end if                                                           1d10s18
      end do                                                            1d10s18
cccccccccccccccccccccccccccccccccccccccccc
      ibcoff=isend                                                      1d10s18
      if(c3.ne.0d0)then                                                 4d11s18
c
c     we have everything we need for this contribution
c
       do is3=1,nsdlk1                                                  4d13s18
        call ilimts(nvirtc(isblk1(1,is3)),noc(isblk1(2,is3)),mynprocg,  4d16s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           4d13s18
        nhere=ih+1-il                                                   4d13s18
        if(nhere.gt.0.and.noc(isblk1(3,is3)).ne.0.and.                  4d13s18
     $       nbasdwsc(isblk1(4,is3)).ne.0)then                          4d13s18
         itmp=ibcoff                                                    4d13s18
         ncol=noc(isblk1(3,is3))*noc(isblk1(4,is3))                     4d16s18
         ibcoff=itmp+nhere*ncol                                         4d13s18
         call enough('updateg. 30',bc,ibc)
         nr=nhere*noc(isblk1(3,is3))                                    4d16s18
         if(min(nr,noc(isblk1(4,is3)),nbasdwsc(isblk1(4,is3))).gt.0)then5d10s19
         call dgemm('n','n',nr,noc(isblk1(4,is3)),
     $        nbasdwsc(isblk1(4,is3)),1d0,bc(iok3x(1,is3)),nr,          4d16s18
     $        bc(itotm(isblk1(4,is3))),nbasdwsc(isblk1(4,is3)),0d0,     4d16s18
     $        bc(itmp),nr,                                              4d16s18
     d' updateg. 29')
         end if                                                         2d25s19
         do i4=0,noc(isblk1(4,is3))-1                                   4d16s18
          do i3=0,noc(isblk1(3,is3))-1                                  4d16s18
           iad1=itmp+nhere*(i3+noc(isblk1(3,is3))*i4)                   4d16s18
           iad2=iok3x(1,is3)+nhere*(i4+noc(isblk1(4,is3))*i3)           4d16s18
           do i12=0,nhere-1                                             4d16s18
            bc(iad2+i12)=bc(iad1+i12)                                   4d16s18
           end do                                                       4d16s18
          end do                                                        4d16s18
         end do                                                         4d16s18
         nr=nhere*noc(isblk1(4,is3))                                    4d16s18
         if(nr.gt.0.and.noc(isblk1(3,is3)).gt.0)then                    2d25s19
         call dgemm('n','n',nr,noc(isblk1(3,is3)),noc(isblk1(3,is3)),   4d16s18
     $        c3,bc(iok3x(1,is3)),nr,bc(itmat(isblk1(3,is3))),          4d16s18
     $        noc(isblk1(3,is3)),0d0,bc(itmp),nr,                       4d16s18
     d' updateg. 30')
         end if                                                         2d25s19
         do is1=1,nsdlk1                                                4d16s18
          if(isblk1(1,is3).eq.isblk1(4,is1).and.                        4d16s18
     $       isblk1(2,is3).eq.isblk1(3,is1))then                        4d16s18
           if(isblk1(3,is3).eq.isblk1(2,is1))then                       4d16s18
            i10=i1s                                                     4d16s18
            i1n=nvirtc(isblk1(1,is3))                                   4d16s18
            jtmp=itmp                                                   4d16s18
            do i2=i2s,i2e                                               4d16s18
             i2m=i2-1                                                   4d16s18
             if(i2.eq.i2e)i1n=i1e                                       4d16s18
             do i1=i10,i1n                                              4d16s18
              i1m=i1-1                                                  4d16s18
              do i3=0,noc(isblk1(3,is3))-1                              4d16s18
               iad1=jtmp+nhere*noc(isblk1(4,is3))*i3                    4d16s18
               iad2=ionex3(is1)+noc(isblk1(4,is3))*(i3                  4d16s18
     $               +noc(isblk1(3,is3))*(i2m+noc(isblk1(3,is1))*i1m))  4d16s18
               do i4=0,noc(isblk1(4,is3))-1                             4d16s18
                kad1=iad1+nhere*i4                                      4d16s18
                bc(iad2+i4)=bc(iad2+i4)+bc(kad1)                        4d16s18
               end do                                                   4d16s18
              end do                                                    4d16s18
              jtmp=jtmp+1                                               4d16s18
             end do                                                     4d16s18
             i10=1                                                      4d16s18
            end do                                                      4d16s18
           end if                                                       4d16s18
           if(isblk1(3,is3).eq.isblk1(1,is1))then                       4d16s18
            i10=i1s                                                     4d16s18
            i1n=nvirtc(isblk1(1,is3))                                   4d16s18
            jtmp=itmp                                                   4d16s18
            do i2=i2s,i2e                                               4d16s18
             i2m=i2-1                                                   4d16s18
             if(i2.eq.i2e)i1n=i1e                                       4d16s18
             do i1=i10,i1n                                              4d16s18
              i1m=i1-1                                                  4d16s18
              do i3=0,noc(isblk1(3,is3))-1                              4d16s18
               iad1=jtmp+nhere*noc(isblk1(4,is3))*i3                    4d16s18
               iad2=ionex3(is1)+i3+noc(isblk1(3,is3))*(                 4d16s18
     $               +noc(isblk1(4,is3))*(i2m+noc(isblk1(3,is1))*i1m))  4d16s18
               do i4=0,noc(isblk1(4,is3))-1                             4d16s18
                jad2=iad2+i4*noc(isblk1(1,is1))                         4d16s18
                kad1=iad1+nhere*i4                                      4d16s18
                bc(jad2)=bc(jad2)+bc(kad1)                              4d16s18
               end do                                                   4d16s18
              end do                                                    4d16s18
              jtmp=jtmp+1                                               4d16s18
             end do                                                     4d16s18
             i10=1                                                      4d16s18
            end do                                                      4d16s18
           end if                                                       4d16s18
          end if                                                        4d16s18
         end do                                                         4d16s18
         ibcoff=itmp                                                    4d16s18
        end if                                                          4d13s18
        call ilimts(noc(isblk1(1,is3)),noc(isblk1(2,is3)),mynprocg,     4d16s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           4d13s18
        nhere=ih+1-il                                                   4d13s18
        if(nhere.gt.0.and.noc(isblk1(3,is3)).ne.0.and.                  4d13s18
     $       nbasdwsc(isblk1(4,is3)).ne.0)then                          4d13s18
         itmp=ibcoff                                                    4d13s18
         ncol=noc(isblk1(3,is3))*nbasdwsc(isblk1(4,is3))                4d16s18
         ibcoff=itmp+nhere*ncol                                         4d13s18
         call enough('updateg. 31',bc,ibc)
         nr=nhere*noc(isblk1(3,is3))                                    4d16s18
         if(nr.gt.0.and.nbasdwsc(isblk1(4,is3)).gt.0)then               2d25s19
         call dgemm('n','n',nr,nbasdwsc(isblk1(4,is3)),                 4d16s18
     $        nbasdwsc(isblk1(4,is3)),c3,bc(iok3x(2,is3)),nr,           4d16s18
     $        bc(iumat(isblk1(4,is3))),nbasdwsc(isblk1(4,is3)),         4d16s18
     $        0d0,bc(itmp),nr,                                          4d16s18
     d' updateg. 31')
         end if                                                         2d25s19
         do i4=0,nbasdwsc(isblk1(4,is3))-1                              4d16s18
          do i3=0,noc(isblk1(3,is3))-1                                  4d16s18
           iad1=itmp+nhere*(i3+noc(isblk1(3,is3))*i4)                   4d16s18
           iad2=iok3x(2,is3)+nhere*(i4+nbasdwsc(isblk1(4,is3))*i3)      4d16s18
           do i12=0,nhere-1                                             4d16s18
            bc(iad2+i12)=bc(iad1+i12)                                   4d16s18
           end do                                                       4d16s18
          end do                                                        4d16s18
         end do                                                         4d16s18
         nr=nhere*nbasdwsc(isblk1(4,is3))                               4d16s18
         if(nr.gt.0.and.noc(isblk1(3,is3)).gt.0)then                    2d25s19
         call dgemm('n','n',nr,noc(isblk1(3,is3)),noc(isblk1(3,is3)),   4d16s18
     $        1d0,bc(iok3x(2,is3)),nr,bc(itmat(isblk1(3,is3))),         4d16s18
     $        noc(isblk1(3,is3)),0d0,bc(itmp),nr,                       4d16s18
     d' updateg. 32')
         end if                                                         2d25s19
         do is1=1,nsdlk1                                                4d16s18
          if(isblk1(4,is1).eq.isblk1(4,is3).and.                        4d16s18
     $         isblk1(3,is1).eq.isblk1(3,is3))then                      4d16s18
           if(isblk1(1,is1).eq.isblk1(1,is3))then                       4d16s18
            i10=i1s                                                     4d16s18
            i1n=noc(isblk1(1,is3))                                      4d16s18
            nnn=nhere*nbasdwsc(isblk1(4,is3))                           4d17s18
            nnnn=noc(isblk1(1,is1))*noc(isblk1(2,is1))                  4d16s18
            jtmp=itmp                                                   4d17s18
            do i2=i2s,i2e                                               4d16s18
             i2m=i2-1                                                   4d16s18
             if(i2.eq.i2e)i1n=i1e                                       4d16s18
             do i1=i10,i1n                                              4d16s18
              i1m=i1-1                                                  4d16s18
              do i4=0,nvirtc(isblk1(4,is3))-1                           4d16s18
               i4p=i4+noc(isblk1(4,is3))                                4d17s18
               iad1=jtmp+nhere*i4p                                      4d17s18
               iad2=ionex3(is1)+i1m+noc(isblk1(1,is1))*(i2m             4d16s18
     $              +noc(isblk1(2,is1))*noc(isblk1(3,is1))*i4)          4d16s18
               do i3=0,noc(isblk1(3,is3))-1                              4d16s18
                jad1=iad1+nnn*i3                                        4d16s18
                jad2=iad2+nnnn*i3                                       4d16s18
                bc(jad2)=bc(jad2)+bc(jad1)                              4d17s18
               end do                                                   4d16s18
              end do                                                    4d16s18
              jtmp=jtmp+1                                               4d16s18
             end do                                                     4d16s18
             i10=1                                                      4d16s18
            end do                                                      4d16s18
           else if(isblk1(2,is1).eq.isblk1(1,is3))then                       4d16s18
            i10=i1s                                                     4d16s18
            i1n=noc(isblk1(1,is3))                                      4d16s18
            nnn=nhere*nbasdwsc(isblk1(4,is3))                           4d17s18
            nnnn=noc(isblk1(1,is1))*noc(isblk1(2,is1))                  4d16s18
            jtmp=itmp                                                   4d17s18
            do i2=i2s,i2e                                               4d16s18
             i2m=i2-1                                                   4d16s18
             if(i2.eq.i2e)i1n=i1e                                       4d16s18
             do i1=i10,i1n                                              4d16s18
              i1m=i1-1                                                  4d16s18
              do i4=0,nvirtc(isblk1(4,is3))-1                           4d16s18
               i4p=i4+noc(isblk1(4,is3))                                4d17s18
               iad1=jtmp+nhere*i4p                                      4d17s18
               iad2=ionex3(is1)+i2m+noc(isblk1(1,is1))*(i1m             4d16s18
     $              +noc(isblk1(2,is1))*noc(isblk1(3,is1))*i4)          4d16s18
               do i3=0,noc(isblk1(3,is3))-1                              4d16s18
                jad1=iad1+nnn*i3                                        4d16s18
                jad2=iad2+nnnn*i3                                       4d16s18
                bc(jad2)=bc(jad2)+bc(jad1)                              4d17s18
               end do                                                   4d16s18
              end do                                                    4d16s18
              jtmp=jtmp+1                                               4d16s18
             end do                                                     4d16s18
             i10=1                                                      4d16s18
            end do                                                      4d16s18
           end if                                                       4d16s18
          end if                                                        4d16s18
         end do                                                         4d16s18
         do i4=0,noc(isblk1(4,is3))-1                                   4d16s18
          do i3=0,noc(isblk1(3,is3))-1                                  4d16s18
           iad1=itmp+nhere*(i4+nbasdwsc(isblk1(4,is3))*i3)              4d16s18
           iad2=iok3x(2,is3)+nhere*(i3+noc(isblk1(3,is3))*i4)           4d16s18
           do i12=0,nhere-1                                             4d16s18
            bc(iad2+i12)=bc(iad1+i12)                                   4d16s18
           end do                                                       4d16s18
          end do                                                        4d16s18
         end do                                                         4d16s18
         nr=nhere*noc(isblk1(3,is3))                                    4d16s18
         if(nr.gt.0.and.noc(isblk1(4,is3)).gt.0)then                    2d25s19
         call dgemm('n','n',nr,noc(isblk1(4,is3)),noc(isblk1(4,is3)),   4d16s18
     $        1d0,bc(iok3x(2,is3)),nr,bc(itmat(isblk1(4,is3))),                 4d16s18
     $        noc(isblk1(4,is3)),0d0,bc(itmp),nr,                       4d16s18
     d' updateg. 33')
         end if                                                         2d25s19
         do isj=1,nsdlk                                                 4d16s18
          if(isblk(1,isj).eq.isblk1(1,is3).and.                         4d16s18
     $         isblk(2,isj).eq.isblk1(2,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(3,is3))then                       4d16s18
           if(isblk(1,isj).eq.isblk(2,isj))then                         4d16s18
            nrowj=(noc(isblk(1,isj))*(noc(isblk(1,isj))+1))/2           4d16s18
            iswitch=0                                                   4d16s18
           else                                                         4d16s18
            nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                   4d16s18
            iswitch=1                                                   4d16s18
           end if                                                       4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             if(i1.le.i2.or.isblk1(1,is3).ne.isblk1(2,is3))then         4d16s18
              i1m=i1-1                                                  4d16s18
              ieq=((i2*i2m)/2)+i1m                                      4d16s18
              inot=i1m+noc(isblk(1,isj))*i2m                            4d16s18
              iad2=ioooo2(isj)+(inot-ieq)*iswitch+ieq                   4d16s18
              do i3=0,noc(isblk1(3,is3))-1                              4d16s18
               iad1=jtmp+nhere*i3                                       4d16s18
               jad2=iad2+nrowj*i3                                        4d16s18
               do i4=0,noc(isblk1(4,is3))-1                             4d16s18
                kad2=jad2+nrowj*noc(isblk(3,isj))*i4                    4d16s18
                kad1=iad1+nhere*noc(isblk1(3,is3))*i4                   4d16s18
                bc(kad2)=bc(kad2)+bc(kad1)                              4d16s18
               end do                                                   4d16s18
              end do                                                    4d16s18
             end if                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          else if(isblk(1,isj).eq.isblk1(2,is3).and.                    4d16s18
     $         isblk(2,isj).eq.isblk1(1,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(3,is3))then                       4d16s18
           nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                    4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             i1m=i1-1                                                   4d16s18
             iad2=ioooo2(isj)+i2m+noc(isblk(1,isj))*i1m                 4d16s18
             do i3=0,noc(isblk1(3,is3))-1                               4d16s18
              iad1=jtmp+nhere*i3                                        4d16s18
              jad2=iad2+nrowj*i3                                        4d16s18
              do i4=0,noc(isblk1(4,is3))-1                              4d16s18
               kad2=jad2+nrowj*noc(isblk(3,isj))*i4                     4d16s18
               kad1=iad1+nhere*noc(isblk1(3,is3))*i4                    4d16s18
               bc(kad2)=bc(kad2)+bc(kad1)                               4d16s18
              end do                                                    4d16s18
             end do                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          end if                                                        4d16s18
          if(isblk(1,isj).eq.isblk1(1,is3).and.                         4d16s18
     $         isblk(2,isj).eq.isblk1(2,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(4,is3))then                       4d16s18
           if(isblk(1,isj).eq.isblk(2,isj))then                         4d16s18
            nrowj=(noc(isblk(1,isj))*(noc(isblk(1,isj))+1))/2           4d16s18
            iswitch=0                                                   4d16s18
           else                                                         4d16s18
            nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                   4d16s18
            iswitch=1                                                   4d16s18
           end if                                                       4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             if(i1.le.i2.or.isblk1(1,is3).ne.isblk1(2,is3))then         4d16s18
              i1m=i1-1                                                  4d16s18
              ieq=((i2*i2m)/2)+i1m                                      4d16s18
              inot=i1m+noc(isblk(1,isj))*i2m                            4d16s18
              iad2=ioooo2(isj)+(inot-ieq)*iswitch+ieq                   4d16s18
              do i3=0,noc(isblk1(3,is3))-1                              4d16s18
               iad1=jtmp+nhere*i3                                       4d16s18
               jad2=iad2+nrowj*noc(isblk(3,isj))*i3                     4d16s18
               do i4=0,noc(isblk1(4,is3))-1                             4d16s18
                kad2=jad2+nrowj*i4                                      4d16s18
                kad1=iad1+nhere*noc(isblk1(3,is3))*i4                   4d16s18
                bc(kad2)=bc(kad2)+bc(kad1)                              4d16s18
               end do                                                   4d16s18
              end do                                                    4d16s18
             end if                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          else if(isblk(1,isj).eq.isblk1(2,is3).and.                    4d16s18
     $         isblk(2,isj).eq.isblk1(1,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(4,is3))then                       4d16s18
           nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                    4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             i1m=i1-1                                                   4d16s18
             iad2=ioooo2(isj)+i2m+noc(isblk(1,isj))*i1m                 4d16s18
             do i3=0,noc(isblk1(3,is3))-1                               4d16s18
              iad1=jtmp+nhere*i3                                        4d16s18
              jad2=iad2+nrowj*noc(isblk(3,isj))*i3                      4d16s18
              do i4=0,noc(isblk1(4,is3))-1                              4d16s18
               kad2=jad2+nrowj*i4                                       4d16s18
               kad1=iad1+nhere*noc(isblk1(3,is3))*i4                    4d16s18
               bc(kad2)=bc(kad2)+bc(kad1)                               4d16s18
              end do                                                    4d16s18
             end do                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          end if                                                        4d16s18
          if(isblk(1,isj).eq.isblk1(3,is3).and.                         4d16s18
     $         isblk(2,isj).eq.isblk1(4,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(1,is3))then                       4d16s18
           if(isblk(1,isj).eq.isblk(2,isj))then                         4d16s18
            nrowj=(noc(isblk(1,isj))*(noc(isblk(1,isj))+1))/2           4d16s18
            iswitch=0                                                   4d16s18
           else                                                         4d16s18
            nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                   4d16s18
            iswitch=1                                                   4d16s18
           end if                                                       4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             i1m=i1-1                                                   4d16s18
             iad2=ioooo2(isj)+nrowj*(i1m+noc(isblk(3,isj))*i2m)         4d16s18
             do i4=0,noc(isblk1(4,is3))-1                               4d16s18
              i3top=(1-iswitch)*i4+iswitch*(noc(isblk1(3,is3))-1)       4d17s18
              iad1=jtmp+nhere*noc(isblk1(3,is3))*i4                     4d17s18
              ixx=(i4*(i4+1))/2                                         4d16s18
              jad2=iad2+(noc(isblk(1,isj))*i4-ixx)*iswitch+ixx          4d17s18
              do i3=0,i3top                                             4d17s18
               kad1=iad1+nhere*i3                                       4d17s18
               bc(jad2+i3)=bc(jad2+i3)+bc(kad1)                         4d16s18
              end do                                                    4d16s18
             end do                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          else if(isblk(1,isj).eq.isblk1(3,is3).and.                    4d16s18
     $         isblk(2,isj).eq.isblk1(4,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(2,is3))then                       4d16s18
           nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                    4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             i1m=i1-1                                                   4d16s18
             iad2=ioooo2(isj)+nrowj*(i2m+noc(isblk(3,isj))*i1m)         4d16s18
             do i4=0,noc(isblk1(4,is3))-1                               4d16s18
              iad1=jtmp+nhere*noc(isblk1(3,is3))*i4                     4d17s18
              jad2=iad2+noc(isblk(1,isj))*i3                            4d16s18
              do i3=0,noc(isblk1(3,is3))-1                              4d16s18
               kad1=iad1+nhere*i3                                       4d17s18
               bc(jad2+i3)=bc(jad2+i3)+bc(kad1)                         4d17s18
              end do                                                    4d16s18
             end do                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          end if                                                        4d16s18
          if(isblk(1,isj).eq.isblk1(4,is3).and.                         4d16s18
     $         isblk(2,isj).eq.isblk1(3,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(1,is3))then                       4d16s18
           if(isblk(1,isj).eq.isblk(2,isj))then                         4d16s18
            nrowj=(noc(isblk(1,isj))*(noc(isblk(1,isj))+1))/2           4d16s18
            iswitch=0                                                   4d16s18
           else                                                         4d16s18
            nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                   4d16s18
            iswitch=1                                                   4d16s18
           end if                                                       4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             i1m=i1-1                                                   4d16s18
             iad2=ioooo2(isj)+nrowj*(i1m+noc(isblk(3,isj))*i2m)         4d16s18
             do i3=0,noc(isblk1(3,is3))-1                               4d17s18
              i4top=(1-iswitch)*i3+iswitch*(noc(isblk1(4,is3))-1)       4d17s18
              iad1=jtmp+nhere*i3                                        4d17s18
              ixx=(i3*(i3+1))/2                                         4d16s18
              jad2=iad2+(noc(isblk(1,isj))*i3-ixx)*iswitch+ixx          4d17s18
              do i4=0,i4top                                             4d17s18
               kad1=iad1+nhere*noc(isblk1(3,is3))*i4                    4d17s18
               bc(jad2+i4)=bc(jad2+i4)+bc(kad1)                         4d16s18
              end do                                                    4d16s18
             end do                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          else if(isblk(1,isj).eq.isblk1(4,is3).and.                    4d16s18
     $         isblk(2,isj).eq.isblk1(3,is3).and.                       4d16s18
     $         isblk(3,isj).eq.isblk1(2,is3))then                       4d16s18
           nrowj=noc(isblk(1,isj))*noc(isblk(2,isj))                    4d16s18
           i10=i1s                                                      4d16s18
           i1n=noc(isblk1(1,is3))                                       4d16s18
           jtmp=itmp                                                    4d16s18
           do i2=i2s,i2e                                                4d16s18
            if(i2.eq.i2e)i1n=i1e                                        4d16s18
            i2m=i2-1                                                    4d16s18
            do i1=i10,i1n                                               4d16s18
             i1m=i1-1                                                   4d16s18
             iad2=ioooo2(isj)+nrowj*(i2m+noc(isblk(3,isj))*i1m)         4d16s18
             do i3=0,noc(isblk1(3,is3))-1                               4d16s18
              iad1=jtmp+nhere*i3                                        4d17s18
              jad2=iad2+noc(isblk(1,isj))*i3                            4d16s18
              do i4=0,noc(isblk1(4,is3))-1                              4d17s18
               kad1=iad1+nhere*noc(isblk1(3,is3))*i4                     4d16s18
               bc(jad2+i4)=bc(jad2+i4)+bc(kad1)                         4d17s18
              end do                                                    4d16s18
             end do                                                     4d16s18
             jtmp=jtmp+1                                                4d16s18
            end do                                                      4d16s18
            i10=1                                                       4d16s18
           end do                                                       4d16s18
          end if                                                        4d16s18
         end do                                                         4d16s18
         ibcoff=itmp                                                    4d16s18
        end if                                                          4d17s18
        if(isblk1(1,is3).ne.isblk1(2,is3))then                          4d16s18
          call ilimts(noc(isblk1(1,is3)),nvirtc(isblk1(2,is3)),mynprocg,4d16s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           4d13s18
          nhere=ih+1-il                                                   4d13s18
          if(nhere.gt.0.and.noc(isblk1(3,is3)).ne.0.and.                  4d13s18
     $         nbasdwsc(isblk1(4,is3)).ne.0)then                          4d13s18
           itmp=ibcoff                                                    4d13s18
           ncol=noc(isblk1(3,is3))*noc(isblk1(4,is3))                     4d16s18
           ibcoff=itmp+nhere*ncol                                         4d13s18
           call enough('updateg. 32',bc,ibc)
           nr=nhere*noc(isblk1(3,is3))                                    4d16s18
           if(min(nr,noc(isblk1(4,is3)),nbasdwsc(isblk1(4,is3)))        5d10s19
     $          .gt.0)then                                              5d10s19
           call dgemm('n','n',nr,noc(isblk1(4,is3)),
     $        nbasdwsc(isblk1(4,is3)),1d0,bc(iok3x(3,is3)),nr,          4d16s18
     $        bc(itotm(isblk1(4,is3))),nbasdwsc(isblk1(4,is3)),0d0,     4d16s18
     $        bc(itmp),nr,                                              4d16s18
     d' updateg. 34')
           end if                                                       2d25s19
           do i4=0,noc(isblk1(4,is3))-1                                   4d16s18
            do i3=0,noc(isblk1(3,is3))-1                                  4d16s18
             iad1=itmp+nhere*(i3+noc(isblk1(3,is3))*i4)                   4d16s18
             iad2=iok3x(3,is3)+nhere*(i4+noc(isblk1(4,is3))*i3)           4d16s18
             do i12=0,nhere-1                                             4d16s18
              bc(iad2+i12)=bc(iad1+i12)                                   4d16s18
             end do                                                       4d16s18
            end do                                                        4d16s18
           end do                                                         4d16s18
           nr=nhere*noc(isblk1(4,is3))                                    4d16s18
           if(nr.gt.0.and.noc(isblk1(3,is3)).gt.0)then                  2d25s19
           call dgemm('n','n',nr,noc(isblk1(3,is3)),noc(isblk1(3,is3)),   4d16s18
     $          c3,bc(iok3x(3,is3)),nr,bc(itmat(isblk1(3,is3))),          4d16s18
     $          noc(isblk1(3,is3)),0d0,bc(itmp),nr,                       4d16s18
     d' updateg. 35')
           end if                                                       2d25s19
           do is1=1,nsdlk1                                              4d16s18
            if(isblk1(4,is1).eq.isblk1(2,is3).and.                      4d16s18
     $         isblk1(3,is1).eq.isblk1(1,is3))then                      4d16s18
             if(isblk1(1,is1).eq.isblk1(3,is3))then
              i10=i1s                                                   4d16s18
              i1n=noc(isblk1(1,is3))                                    4d16s18
              jtmp=itmp                                                 4d16s18
              do i2=i2s,i2e                                             4d16s18
               i2m=i2-1                                                 4d16s18
               if(i2.eq.i2e)i1n=i1e                                     4d16s18
               do i1=i10,i1n                                            4d16s18
                i1m=i1-1                                                4d17s18
                do i3=0,noc(isblk1(3,is3))-1                            4d16s18
                 iad1=jtmp+nhere*noc(isblk1(4,is3))*i3                  4d16s18
                 iad2=ionex3(is1)+i3+noc(isblk1(1,is1))                 4d16s18
     $                *noc(isblk1(2,is1))*(i1m+noc(isblk1(3,is1))*i2m)  4d17s18
                 do i4=0,noc(isblk1(4,is3))-1                           4d16s18
                  kad1=iad1+nhere*i4                                    4d16s18
                  kad2=iad2+noc(isblk1(1,is1))*i4                       4d16s18
                  bc(kad2)=bc(kad2)+bc(kad1)                            4d16s18
                 end do                                                 4d16s18
                end do                                                  4d16s18
                jtmp=jtmp+1                                             4d16s18
               end do                                                   4d16s18
               i10=1                                                    4d16s18
              end do                                                    4d16s18
             end if                                                     4d16s18
             if(isblk1(2,is1).eq.isblk1(3,is3))then
              i10=i1s                                                   4d16s18
              i1n=noc(isblk1(1,is3))                                    4d16s18
              jtmp=itmp                                                 4d16s18
              do i2=i2s,i2e                                             4d16s18
               i2m=i2-1                                                 4d16s18
               if(i2.eq.i2e)i1n=i1e                                     4d16s18
               do i1=i10,i1n                                            4d16s18
                i1m=i1-1                                                4d17s18
                do i3=0,noc(isblk1(3,is3))-1                            4d16s18
                 iad1=jtmp+nhere*noc(isblk1(4,is3))*i3                  4d16s18
                 iad2=ionex3(is1)+noc(isblk1(1,is1))*(i3                4d16s18
     $                +noc(isblk1(2,is1))*(i1m+noc(isblk1(3,is1))*i2m)) 4d17s18
                 do i4=0,noc(isblk1(4,is3))-1                           4d16s18
                  kad1=iad1+nhere*i4                                    4d16s18
                  bc(iad2+i4)=bc(iad2+i4)+bc(kad1)                         4d16s18
                 end do                                                 4d16s18
                end do                                                  4d16s18
                jtmp=jtmp+1                                             4d16s18
               end do                                                   4d16s18
               i10=1                                                    4d16s18
              end do                                                    4d16s18
             end if                                                     4d16s18
            end if                                                      4d16s18
           end do                                                       4d16s18
           ibcoff=itmp                                                  4d16s18
          end if                                                        4d16s18
        end if                                                          4d13s18
       end do                                                           4d13s18
      end if                                                            4d11s18
      if(c2.ne.0d0)then                                                 3d16s18
c
c     for jmats. need to grab parts of onex and 4o first.
c
       do is=1,nsdlk                                                    3d16s18
        if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0              3d16s18
     $   .and.nbasdwsc(isblk(3,is)).ne.0.and.nbasdwsc(isblk(4,is)).ne.0)3d16s18
     $       then                                                       3d16s18
         call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,        3d16s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           3d16s18
         nhere=ih+1-il                                                  3d16s18
         if(noc(isblk(3,is)).ne.0.and.nvirtc(isblk(4,is)).ne.0)then     3d16s18
          do isx=1,nsdlk1                                               3d16s18
           if(isblk(1,is).eq.isblk1(1,isx).and.                         3d16s18
     $          isblk(2,is).eq.isblk1(2,isx).and.                       3d16s18
     $          isblk(3,is).eq.isblk1(3,isx))then                       3d16s18
            do i2=0,nvirtc(isblk1(4,isx))-1                               3d16s18
             i2p=i2+noc(isblk1(4,isx))                                  3d16s18
             do i1=0,noc(isblk1(3,isx))-1                                3d16s18
              iad1=iovnn(isx)+nhere*(i1+noc(isblk1(3,isx))*i2)          3d16s18
              iad2=ijnn(is)+nhere*(i1+nbasdwsc(isblk(3,is))*i2p)        3d16s18
              do i34=0,nhere-1                                          3d16s18
               bc(iad2+i34)=bc(iad1+i34)                                3d16s18
              end do                                                    3d16s18
             end do                                                     3d16s18
            end do                                                      3d16s18
           end if                                                       3d16s18
          end do                                                        3d16s18
         end if                                                         3d16s18
         if(noc(isblk(4,is)).ne.0.and.nvirtc(isblk(3,is)).ne.0)then     3d16s18
          do isx=1,nsdlk1                                               3d16s18
           if(isblk(1,is).eq.isblk1(1,isx).and.                         3d16s18
     $          isblk(2,is).eq.isblk1(2,isx).and.                       3d16s18
     $          isblk(3,is).eq.isblk1(4,isx))then                       3d16s18
            do i2=0,nvirtc(isblk1(4,isx))-1                               3d16s18
             i2p=i2+noc(isblk1(4,isx))                                  3d16s18
             do i1=0,noc(isblk1(3,isx))-1                                3d16s18
              iad1=iovnn(isx)+nhere*(i1+noc(isblk1(3,isx))*i2)          3d16s18
              iad2=ijnn(is)+nhere*(i2p+nbasdwsc(isblk(3,is))*i1)        3d16s18
              do i34=0,nhere-1                                          3d16s18
               bc(iad2+i34)=bc(iad1+i34)                                3d16s18
              end do                                                    3d16s18
             end do                                                     3d16s18
            end do                                                      3d16s18
           end if                                                       3d16s18
          end do                                                        3d16s18
         end if                                                         3d16s18
         if(noc(isblk(4,is)).ne.0.and.noc(isblk(3,is)).ne.0)then        3d16s18
          do i2=0,noc(isblk(4,is))-1                                    3d16s18
           do i1=0,noc(isblk(3,is))-1                                   3d16s18
            iad1=ioonn(is)+nhere*(i1+noc(isblk(3,is))*i2)               3d16s18
            iad2=ijnn(is)+nhere*(i1+nbasdwsc(isblk(3,is))*i2)           3d16s18
            do i34=0,nhere-1                                            3d16s18
             bc(iad2+i34)=bc(iad1+i34)                                  3d16s18
            end do                                                      3d16s18
           end do                                                       3d16s18
          end do                                                        3d16s18
         end if                                                         3d16s18
         ncol=nbasdwsc(isblk(3,is))*nbasdwsc(isblk(4,is))               3d16s18
         itmp1=ibcoff                                                   3d16s18
         ibcoff=itmp1+nhere*ncol                                        3d16s18
         call enough('updateg. 33',bc,ibc)
         nr=nhere*nbasdwsc(isblk(3,is))                                 3d16s18
         if(nr.gt.0.and.nbasdwsc(isblk(4,is)).gt.0)then                 2d25s19
         call dgemm('n','n',nr,nbasdwsc(isblk(4,is)),                   3d16s18
     $        nbasdwsc(isblk(4,is)),1d0,bc(ijnn(is)),nr,                3d16s18
     $    bc(iumat(isblk(4,is))),nbasdwsc(isblk(4,is)),0d0,bc(itmp1),nr,3d16s18
     d' updateg. 36')
         end if                                                         2d25s19
         do i4=0,nbasdwsc(isblk(4,is))-1                                3d16s18
          do i3=0,nbasdwsc(isblk(3,is))-1                               3d16s18
           iad1=itmp1+nhere*(i3+nbasdwsc(isblk(3,is))*i4)               3d16s18
           iad2=ijnn(is)+nhere*(i4+nbasdwsc(isblk(4,is))*i3)            3d16s18
           do i12=0,nhere-1                                             3d16s18
            bc(iad2+i12)=bc(iad1+i12)                                   3d16s18
           end do                                                       3d16s18
          end do                                                        3d16s18
         end do                                                         3d16s18
         nr=nhere*nbasdwsc(isblk(4,is))                                 3d16s18
         if(nr.gt.0.and.nbasdwsc(isblk(3,is)).gt.0)then                 2d25s19
         call dgemm('n','n',nr,nbasdwsc(isblk(3,is)),                   3d16s18
     $        nbasdwsc(isblk(3,is)),1d0,bc(ijnn(is)),nr,                3d16s18
     $    bc(iumat(isblk(3,is))),nbasdwsc(isblk(3,is)),0d0,bc(itmp1),nr, 3d16s18
     d' updateg. 37')
         end if                                                         2d25s19
         nr=nhere*nbasdwsc(isblk(4,is))                                 3d16s18
         if(nr.gt.0.and.noc(isblk(3,is)).gt.0)then                      2d25s19
         call dgemm('n','n',nr,noc(isblk(3,is)),noc(isblk(3,is)),c2,    3d16s18
     $        bc(itmp1),nr,bc(itmat(isblk(3,is))),noc(isblk(3,is)),0d0, 3d16s18
     $        bc(ijnn(is)),nr,                                          3d16s18
     d' updateg. 38')
         end if                                                         2d25s19
         do is1=1,nsdlk1                                                3d16s18
          if(isblk1(1,is1).eq.isblk(1,is).and.                          3d16s18
     $         isblk1(2,is1).eq.isblk(2,is).and.isblk1(3,is1).eq.       3d16s18
     $         isblk(3,is))then                                         3d16s18
           do i4=0,nvirtc(isblk(4,is))-1                                 3d16s18
            i4p=i4+noc(isblk(4,is))                                      3d16s18
            do i3=0,noc(isblk(3,is))-1                                   3d16s18
             iad1=ijnn(is)+nhere*(i4p+nbasdwsc(isblk(4,is))*i3)          3d16s18
             iad2=ionex2(is1)+nhere*(i3+noc(isblk(3,is))*i4)            3d16s18
             do irow=0,nhere-1                                          3d16s18
              bc(iad2+irow)=bc(iad2+irow)+bc(iad1+irow)                 3d16s18
             end do                                                      3d16s18
            end do                                                       3d16s18
           end do                                                        3d16s18
          else if(isblk1(1,is1).eq.isblk(2,is).and.                     3d16s18
     $         isblk1(2,is1).eq.isblk(1,is).and.isblk1(3,is1).eq.       3d16s18
     $         isblk(3,is))then                                         3d16s18
           write(6,*)('i''m in the second block for ionex2 ...')
           if(bc(132).ne.-132d0)then
            call dws_sync
            call dws_finalize
            stop
           end if
           do i4=0,nvirtc(isblk(4,is))-1                                 3d16s18
            i4p=i4+noc(isblk(4,is))                                      3d16s18
            do i3=0,noc(isblk(3,is))-1                                   3d16s18
             iad1=ijnn(is)+nhere*(i4p+nbasdwsc(isblk(4,is))*i3)          3d16s18
             iad2=ionex2(is1)+nhere*(i3+noc(isblk(3,is))*i4)            3d16s18
             do irow=0,nhere-1                                          3d16s18
              bc(iad2+irow)=bc(iad2+irow)+bc(iad1+irow)                 3d16s18
             end do                                                      3d16s18
            end do                                                       3d16s18
           end do                                                        3d16s18
          end if                                                        3d16s18
         end do                                                         3d16s18
         itmp2=ibcoff                                                   3d16s18
         ibcoff=itmp2+nhere*noc(isblk(3,is))*noc(isblk(4,is))           3d16s18
         call enough('updateg. 34',bc,ibc)
         do i4=0,noc(isblk(4,is))-1                                     3d16s18
          do i3=0,noc(isblk(3,is))-1                                    3d16s18
           iad1=ijnn(is)+nhere*(i4+nbasdwsc(isblk(4,is))*i3)            3d16s18
           iad2=itmp2+nhere*(i3+noc(isblk(3,is))*i4)                    3d16s18
           do i12=0,nhere-1                                             3d16s18
            bc(iad2+i12)=bc(iad1+i12)                                   3d16s18
           end do                                                       3d16s18
          end do                                                        3d16s18
         end do                                                         3d16s18
         nr=nhere*noc(isblk(3,is))                                      3d16s18
         if(nr.gt.0.and.noc(isblk(4,is)).gt.0)then                      2d25s19
         call dgemm('n','n',nr,noc(isblk(4,is)),noc(isblk(4,is)),1d0,   3d16s18
     $        bc(itmp2),nr,bc(itmat(isblk(4,is))),noc(isblk(4,is)),0d0, 3d16s18
     $        bc(ijnn(is)),nr,                                          3d16s18
     d' updateg. 39')
         end if                                                         2d25s19
         nrow1=noc(isblk(1,is))*noc(isblk(2,is))                        3d16s18
         ncol=noc(isblk(3,is))*noc(isblk(4,is))
         if(isblk(1,is).eq.isblk(2,is))then
          nrow1=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2               3d16s18
         else                                                           3d16s18
          nrow1=noc(isblk(1,is))*noc(isblk(2,is))                       3d16s18
         end if                                                         3d16s18
         do i4=0,noc(isblk(4,is))-1                                     3d16s18
          do i3=0,noc(isblk(3,is))-1                                    3d16s18
           iad2=ijnn(is)+nhere*(i3+noc(isblk(3,is))*i4)                 3d16s18
           iad1=ioooo2(is)+nrow1*(i3+noc(isblk(3,is))*i4)               3d16s18
           i10=i1s                                                      3d16s18
           i1n=noc(isblk(1,is))                                         3d19s18
           do i2=i2s,i2e                                                3d16s18
            i2m=i2-1                                                    3d16s18
            if(i2.eq.i2e)i1n=i1e                                        3d16s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d16s18
             do i1=i10,min(i2,i1n)                                      3d19s18
              i1m=i1-1                                                  3d16s18
              irow=((i2*i2m)/2)+i1m                                     3d16s18
              jrow=i1+noc(isblk(1,is))*i2m-il                           3d16s18
              orig=bc(iad1+irow)
              bc(iad1+irow)=bc(iad1+irow)+bc(iad2+jrow)                 3d16s18
             end do                                                     3d16s18
            else                                                        3d16s18
             do i1=i10,i1n                                               3d16s18
              i1m=i1-1                                                   3d16s18
              irow=i1m+noc(isblk(1,is))*i2m                              3d16s18
              jrow=i1+noc(isblk(1,is))*i2m-il                            3d16s18
              orig=bc(iad1+irow)
              bc(iad1+irow)=bc(iad1+irow)+bc(iad2+jrow)                  3d16s18
             end do                                                      3d16s18
            end if                                                      3d16s18
            i10=1                                                       3d16s18
           end do                                                       3d16s18
          end do                                                        3d16s18
         end do                                                         3d16s18
         do iso=1,nsdlk                                                 3d16s18
          if(isblk(1,is).eq.isblk(3,iso).and.isblk(2,is).eq.isblk(4,iso)3d16s18
     $         .and.isblk(3,is).eq.isblk(1,iso))then                    3d16s18
           if(isblk(1,iso).eq.isblk(2,iso))then
            nrow1=(noc(isblk(1,iso))*(noc(isblk(1,iso))+1))/2           3d16s18
            iswitch=0
           else                                                         3d16s18
            nrow1=noc(isblk(1,iso))*noc(isblk(2,iso))                       3d16s18
            iswitch=1                                                   3d16s18
           end if                                                       3d16s18
           do i4=0,noc(isblk(4,is))-1                                     3d16s18
            do i3=0,iswitch*(noc(isblk(3,is))-1)+i4*(1-iswitch)         3d19s18
             iad2=ijnn(is)+nhere*(i3+noc(isblk(3,is))*i4)                 3d16s18
             iad3=ijnn(is)+nhere*(i4+noc(isblk(3,is))*i3)                 3d16s18
             ieq=((i4*(i4+1))/2)+i3                                     3d16s18
             inot=i3+noc(isblk(1,iso))*i4                               3d16s18
             iad1=ioooo2(iso)+(inot-ieq)*iswitch+ieq                    3d16s18
             i10=i1s                                                      3d16s18
             i1n=noc(isblk(1,is))                                       3d19s18
             do i2=i2s,i2e                                                3d16s18
              i2m=i2-1                                                    3d16s18
              if(i2.eq.i2e)i1n=i1e                                        3d16s18
              do i1=i10,i1n                                              3d16s18
               i1m=i1-1                                                  3d16s18
               irow=nrow1*(i1m+noc(isblk(1,is))*i2m)                    3d16s18
               jrow=i1+noc(isblk(1,is))*i2m-il                           3d16s18
               orig=bc(iad1+irow)
               bc(iad1+irow)=bc(iad1+irow)+bc(iad2+jrow)                 3d16s18
              end do                                                     3d16s18
              i10=1                                                       3d16s18
             end do                                                       3d16s18
            end do                                                        3d16s18
           end do                                                         3d16s18
          else if(isblk(2,is).eq.isblk(3,iso).and.                      3d16s18
     $         isblk(1,is).eq.isblk(4,iso)                              3d16s18
     $         .and.isblk(3,is).eq.isblk(1,iso))then                    3d16s18
           nrow1=noc(isblk(1,iso))*noc(isblk(2,iso))                       3d16s18
           do i4=0,noc(isblk(4,is))-1                                     3d16s18
            do i3=0,noc(isblk(3,is))-1                                  3d16s18
             iad2=ijnn(is)+nhere*(i3+noc(isblk(3,is))*i4)                 3d16s18
             iad1=ioooo2(iso)+i3+noc(isblk(3,is))*i4                    3d16s18
             i10=i1s                                                      3d16s18
             i1n=noc(isblk(1,is))                                       3d19s18
             do i2=i2s,i2e                                                3d16s18
              i2m=i2-1                                                    3d16s18
              if(i2.eq.i2e)i1n=i1e                                        3d16s18
              do i1=i10,i1n                                              3d16s18
               i1m=i1-1                                                  3d16s18
               irow=nrow1*(i2m+noc(isblk(2,is))*i1m)                    3d16s18
               jrow=i1+noc(isblk(1,is))*i2m-il                           3d16s18
               orig=bc(iad1+irow)
               bc(iad1+irow)=bc(iad1+irow)+bc(iad2+jrow)                 3d16s18
              end do                                                     3d16s18
              i10=1                                                       3d16s18
             end do                                                       3d16s18
            end do                                                        3d16s18
           end do                                                         3d16s18
          else if(isblk(2,is).eq.isblk(3,iso).and.                      3d16s18
     $         isblk(1,is).eq.isblk(4,iso)                              3d16s18
     $         .and.isblk(4,is).eq.isblk(1,iso))then                    3d16s18
           nrow1=noc(isblk(1,iso))*noc(isblk(2,iso))                       3d16s18
           do i4=0,noc(isblk(4,is))-1                                     3d16s18
            do i3=0,noc(isblk(3,is))-1                                  3d16s18
             iad2=ijnn(is)+nhere*(i3+noc(isblk(3,is))*i4)                 3d16s18
             iad1=ioooo2(iso)+i4+noc(isblk(4,is))*i3                    3d16s18
             i10=i1s                                                      3d16s18
             i1n=noc(isblk(1,is))                                       3d19s18
             do i2=i2s,i2e                                                3d16s18
              i2m=i2-1                                                    3d16s18
              if(i2.eq.i2e)i1n=i1e                                        3d16s18
              do i1=i10,i1n                                              3d16s18
               i1m=i1-1                                                  3d16s18
               irow=nrow1*(i2m+noc(isblk(2,is))*i1m)                    3d16s18
               jrow=i1+noc(isblk(1,is))*i2m-il                           3d16s18
               orig=bc(iad1+irow)
               bc(iad1+irow)=bc(iad1+irow)+bc(iad2+jrow)                 3d16s18
              end do                                                     3d16s18
              i10=1                                                       3d16s18
             end do                                                       3d16s18
            end do                                                        3d16s18
           end do                                                         3d16s18
          else if(isblk(1,is).eq.isblk(3,iso).and.                      3d16s18
     $         isblk(2,is).eq.isblk(4,iso)                              3d16s18
     $         .and.isblk(4,is).eq.isblk(1,iso))then                    3d16s18
           nrow1=noc(isblk(1,iso))*noc(isblk(2,iso))                       3d16s18
           do i4=0,noc(isblk(4,is))-1                                     3d16s18
            do i3=0,noc(isblk(3,is))-1                                  3d16s18
             iad2=ijnn(is)+nhere*(i3+noc(isblk(3,is))*i4)                 3d16s18
             iad1=ioooo2(iso)+i4+noc(isblk(4,is))*i3                    3d16s18
             i10=i1s                                                      3d16s18
             i1n=noc(isblk(1,is))                                       3d19s18
             do i2=i2s,i2e                                                3d16s18
              i2m=i2-1                                                    3d16s18
              if(i2.eq.i2e)i1n=i1e                                        3d16s18
              do i1=i10,i1n                                              3d16s18
               i1m=i1-1                                                  3d16s18
               irow=nrow1*(i1m+noc(isblk(1,is))*i2m)                    3d16s18
               jrow=i1+noc(isblk(1,is))*i2m-il                           3d16s18
               orig=bc(iad1+irow)
               bc(iad1+irow)=bc(iad1+irow)+bc(iad2+jrow)                 3d16s18
              end do                                                     3d16s18
              i10=1                                                       3d16s18
             end do                                                       3d16s18
            end do                                                        3d16s18
           end do                                                         3d16s18
          end if
         end do                                                         3d16s18
         if(isblk(3,is).ne.isblk(4,is))then                             3d19s18
          do i3=0,nvirtc(isblk(3,is))-1                                  3d16s18
           i3p=i3+noc(isblk(3,is))                                       3d16s18
           do i4=0,noc(isblk(4,is))-1                                    3d16s18
            iad1=itmp1+nhere*(i4+nbasdwsc(isblk(4,is))*i3p)              3d16s18
            iad2=ijnn(is)+nhere*(i3+nvirtc(isblk(3,is))*i4)              3d16s18
            do i12=0,nhere-1                                             3d16s18
             bc(iad2+i12)=bc(iad1+i12)                                   3d16s18
            end do                                                       3d16s18
           end do                                                        3d16s18
          end do                                                         3d16s18
          nr=nhere*nvirtc(isblk(3,is))                                   3d16s18
          if(nr.gt.0.and.noc(isblk(4,is)).gt.0)then                     2d25s19
          call dgemm('n','n',nr,noc(isblk(4,is)),noc(isblk(4,is)),c2,    3d16s18
     $        bc(ijnn(is)),nr,bc(itmat(isblk(4,is))),noc(isblk(4,is)),  3d16s18
     $        0d0,bc(itmp1),nr,                                         3d16s18
     d' updateg. 40')
          end if                                                        2d25s19
          do is1=1,nsdlk1                                                3d16s18
           if(isblk1(1,is1).eq.isblk(1,is).and.                          3d16s18
     $         isblk1(2,is1).eq.isblk(2,is).and.isblk1(3,is1).eq.       3d16s18
     $         isblk(4,is))then                                         3d16s18
            do i4=0,noc(isblk(4,is))-1                                   3d16s18
             do i3=0,nvirtc(isblk(3,is))-1                               3d16s18
              iad1=itmp1+nhere*(i3+nvirtc(isblk(3,is))*i4)               3d19s18
              iad2=ionex2(is1)+nhere*(i4+noc(isblk(4,is))*i3)            3d16s18
              do irow=0,nhere-1                                          3d16s18
               bc(iad2+irow)=bc(iad2+irow)+bc(iad1+irow)                 3d16s18
              end do                                                      3d16s18
             end do                                                       3d16s18
            end do                                                        3d16s18
           else if(isblk1(1,is1).eq.isblk(2,is).and.                     3d16s18
     $         isblk1(2,is1).eq.isblk(1,is).and.isblk1(3,is1).eq.       3d16s18
     $         isblk(4,is))then                                         3d16s18
            write(6,*)('we are in other ijnn onex2 block ...')
            if(bc(132).ne.-132d0)then
             call dws_sync
             call dws_finalize
             stop
            end if
            nrow1=noc(isblk1(1,is1))*noc(isblk1(2,is1))                  3d16s18
            do i4=0,noc(isblk(4,is))-1                                   3d16s18
             do i3=0,nvirtc(isblk(3,is))-1                               3d16s18
              iad1=tmp1+nhere*(i4+noc(isblk(4,is))*i3)                   3d16s18
              iad2=ionex2(is1)+nhere*(i4+noc(isblk(4,is))*i3)            3d16s18
              do irow=0,nhere-1                                          3d16s18
               bc(iad2+irow)=bc(iad2+irow)+bc(iad1+irow)                 3d16s18
              end do
             end do                                                       3d16s18
            end do                                                        3d16s18
           end if                                                        3d16s18
          end do                                                         3d16s18
         end if
         ibcoff=itmp1                                                   3d16s18
        end if                                                          3d16s18
       end do                                                           3d16s18
c
c     for kmats
c     recall K_{nm}^{AB}=(nB|mA)
c
       do isk=1,nsdlkk                                                  3d20s18
        call ilimts(noc(isblkk(1,isk)),noc(isblkk(2,isk)),mynprocg,     3d20s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           3d20s18
        nhere=ih+1-il                                                   3d20s18
        nrowb=nbasdwsc(isblkk(3,isk))*nbasdwsc(isblkk(4,isk))           3d20s18
        itmp=ibcoff                                                     3d20s18
        ibcoff=itmp+nrowb*nhere                                         3d20s18
        nr=nhere*nbasdwsc(isblkk(3,isk))
        if(nr.gt.0.and.nbasdwsc(isblkk(4,isk)).gt.0)then                2d25s19
        call dgemm('n','n',nr,nbasdwsc(isblkk(4,isk)),
     $       nbasdwsc(isblkk(4,isk)),1d0,bc(iokx(isk)),nr,              3d20s18
     $       bc(iumat(isblkk(4,isk))),nbasdwsc(isblkk(4,isk)),0d0,      3d20s18
     $       bc(itmp),nr,                                               3d20s18
     d' updateg. 41')
        end if                                                          2d25s19
        do i4=0,nbasdwsc(isblkk(4,isk))-1
         do i3=0,nbasdwsc(isblkk(3,isk))-1
          iad1=itmp+nhere*(i3+nbasdwsc(isblkk(3,isk))*i4)
          iad2=iokx(isk)+nhere*(i4+nbasdwsc(isblkk(4,isk))*i3)
          do i12=0,nhere-1
           bc(iad2+i12)=bc(iad1+i12)
          end do
         end do
        end do
        nr=nhere*nbasdwsc(isblkk(4,isk))
        if(nr.gt.0.and.nbasdwsc(isblkk(3,isk)).gt.0)then                2d25s19
        call dgemm('n','n',nr,nbasdwsc(isblkk(3,isk)),
     $       nbasdwsc(isblkk(3,isk)),c2,bc(iokx(isk)),nr,               3d20s18
     $       bc(iumat(isblkk(3,isk))),nbasdwsc(isblkk(3,isk)),0d0,      3d20s18
     $       bc(itmp),nr,                                               3d20s18
     d' updateg. 42')
        end if                                                          2d25s19
        itmpt=ibcoff                                                    3d20s18
        ibcoff=itmpt+nhere*nbasdwsc(isblkk(4,isk))*noc(isblkk(3,isk))   3d20s18
        call enough('updateg. 35',bc,ibc)
        if(nr.gt.0.and.noc(isblkk(3,isk)).gt.0)then                     2d25s19
        call dgemm('n','n',nr,noc(isblkk(3,isk)),noc(isblkk(3,isk)),1d0,3d20s18
     $       bc(itmp),nr,bc(itmat(isblkk(3,isk))),noc(isblkk(3,isk)),   3d20s18
     $       0d0,bc(itmpt),nr,                                          3d20s18
     d' updateg. 43')
        end if                                                          2d25s19
c
c     onex2 contribution
c     recall we have K_{nm}^{AB}=(nB|mA), so currently under itmpt we
c     have (nB|ma)=(ma|nB) stored with aB indicies reversed.
c     and we want to add this into (ma|nB) and (am|nB).
c     but first 2 indices of
c     onex2 are distributed while nm are distributed for K. So we need
c     full onex3.
c
        do isx=1,nsdlk1                                                 3d29s18
         if(isblk1(1,isx).eq.isblkk(2,isk).and.                         3d30s18
     $        isblk1(2,isx).eq.isblkk(3,isk).and.                       3d30s18
     $        isblk1(3,isx).eq.isblkk(1,isk))then                       3d30s18
          do i3=0,noc(isblkk(3,isk))-1                                  3d29s18
           do i4=0,nvirtc(isblkk(4,isk))-1                              3d29s18
            i10=i1s                                                     3d30s18
            i1n=noc(isblkk(1,isk))                                      3d30s18
            do i2=i2s,i2e                                               3d30s18
             i2m=i2-1                                                   3d30s18
             if(i2.eq.i2e)i1n=i1e                                       3d30s18
             do i1=i10,i1n                                              3d30s18
              i1m=i1-1                                                  3d30s18
              i4p=i4+noc(isblkk(4,isk))                                 3d30s18
              iad1=ionex3(isx)+i2m+noc(isblk1(1,isx))*(i3+              3d30s18
     $             noc(isblk1(2,isx))*(i1m+noc(isblk1(3,isx))*i4))      3d30s18
              iad2=itmpt+i1+noc(isblkk(1,isk))*i2m-il                   3d30s18
     $             +nhere*(i4p+nbasdwsc(isblkk(4,isk))*i3)              3d30s18
              bc(iad1)=bc(iad1)+bc(iad2)                                3d30s18
             end do                                                     3d30s18
             i10=1                                                      3d30s18
            end do                                                      3d29s18
           end do                                                       3d29s18
          end do                                                        3d29s18
         end if                                                         3d29s18
         if(isblk1(1,isx).eq.isblkk(3,isk).and.                         3d30s18
     $      isblk1(2,isx).eq.isblkk(2,isk).and.                         3d30s18
     $        isblk1(3,isx).eq.isblkk(1,isk))then                       3d30s18
          do i3=0,noc(isblkk(3,isk))-1                                  3d29s18
           do i4=0,nvirtc(isblkk(4,isk))-1                              3d29s18
            i4p=i4+noc(isblkk(4,isk))                                   3d29s18
            i10=i1s                                                     3d29s18
            i1n=noc(isblkk(1,isk))                                      3d29s18
            do i2=i2s,i2e                                               3d29s18
             i2m=i2-1                                                   3d29s18
             if(i2.eq.i2e)i1n=i1e                                       3d29s18
             do i1=i10,i1n                                              3d29s18
              i1m=i1-1                                                  3d29s18
              iad1=itmpt+i1+noc(isblkk(1,isk))*i2m-il                   3d29s18
     $             +nhere*(i4p+nbasdwsc(isblkk(4,isk))*i3)              3d29s18
              iad2=ionex3(isx)+i3+noc(isblk1(1,isx))*(i2m               3d30s18
     $             +noc(isblk1(2,isx))*(i1m+noc(isblk1(3,isx))*i4))     3d30s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d29s18
             end do                                                     3d29s18
             i10=1                                                      3d29s18
            end do                                                      3d29s18
           end do                                                       3d29s18
          end do                                                        3d29s18
         end if                                                         3d29s18
        end do                                                          3d29s18
        do i3=0,noc(isblkk(3,isk))-1
         do i4=0,noc(isblkk(4,isk))-1
          iad1=itmpt+nhere*(i4+nbasdwsc(isblkk(4,isk))*i3)              3d20s18
          iad2=iokx(isk)+nhere*(i3+noc(isblkk(3,isk))*i4)               3d20s18
          do i12=0,nhere-1                                              3d20s18
           bc(iad2+i12)=bc(iad1+i12)                                    3d20s18
          end do                                                        3d20s18
         end do                                                         3d20s18
        end do                                                          3d20s18
        nr=nhere*noc(isblkk(3,isk))
        if(nr.gt.0.and.noc(isblkk(4,isk)).gt.0)then                     2d25s19
        call dgemm('n','n',nr,noc(isblkk(4,isk)),noc(isblkk(4,isk)),1d0,3d20s18
     $       bc(iokx(isk)),nr,bc(itmat(isblkk(4,isk))),
     $       noc(isblkk(4,isk)),0d0,bc(itmpt),nr,                       3d20s18
     d' updateg. 44')
        end if                                                          2d25s19
        nrows=noc(isblkk(3,isk))*noc(isblkk(4,isk))                     3d20s18
c
c     we have K_{ij}^{kl}=(il|jk)
c
        do is=1,nsdlk                                                   3d20s18
         if(isblk(1,is).eq.isblkk(1,isk).and.                           3d20s18
     $        isblk(2,is).eq.isblkk(4,isk).and.                         3d26s18
     $        isblk(3,is).eq.isblkk(2,isk))then                         3d26s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d20s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d20s18
           iswitch=0                                                    3d20s18
          else                                                          3d20s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                          3d20s18
           iswitch=1                                                    3d20s18
          end if                                                        3d20s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d26s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             itop=i1m                                                   3d26s18
            else                                                        3d26s18
             itop=noc(isblkk(4,isk))-1                                  3d26s18
            end if                                                      3d26s18
            do i4=0,itop                                                3d26s18
             ieq=((i1m*(i1m+1))/2)+i4                                   3d26s18
             inot=i1m+noc(isblk(1,is))*i4                               3d26s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i3=0,noc(isblkk(3,isk))-1                               3d20s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i2m+noc(isblk(3,is))*i3                              3d26s18
              iad2=jrow+nr*jcol                                         3d20s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         else if(isblk(1,is).eq.isblkk(2,isk).and.                           3d20s18
     $        isblk(2,is).eq.isblkk(3,isk).and.                         3d26s18
     $        isblk(3,is).eq.isblkk(1,isk))then                         3d26s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d26s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d26s18
           iswitch=0                                                    3d26s18
          else                                                          3d26s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                         3d26s18
           iswitch=1                                                    3d26s18
          end if                                                        3d26s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d28s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             itop=i2m                                                   3d26s18
            else                                                        3d26s18
             itop=noc(isblkk(3,isk))-1                                  3d26s18
            end if                                                      3d26s18
            do i3=0,itop                                                3d26s18
             ieq=((i2m*(i2m+1))/2)+i3                                   3d26s18
             inot=i2m+noc(isblk(1,is))*i3                               3d26s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i4=0,noc(isblkk(4,isk))-1                               3d26s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i1m+noc(isblk(3,is))*i4                              3d26s18
              iad2=jrow+nr*jcol                                         3d26s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         end if
         if(isblk(1,is).eq.isblkk(4,isk).and.                           3d26s18
     $        isblk(2,is).eq.isblkk(1,isk).and.                         3d26s18
     $        isblk(3,is).eq.isblkk(2,isk))then                         3d26s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d20s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d20s18
           iswitch=0                                                    3d20s18
          else                                                          3d20s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                          3d20s18
           iswitch=1                                                    3d20s18
          end if                                                        3d20s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d26s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             ibot=i1m                                                   3d26s18
            else                                                        3d26s18
             ibot=0                                                     3d26s18
            end if                                                      3d26s18
            do i4=ibot,noc(isblkk(4,isk))-1                             3d26s18
             ix=max(i4,i1m)
             in=min(i4,i1m)
             ieq=((ix*(ix+1))/2)+in                                     3d26s18
             inot=i4+noc(isblk(1,is))*i1m                               3d26s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i3=0,noc(isblkk(3,isk))-1                               3d20s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i2m+noc(isblk(3,is))*i3                              3d26s18
              iad2=jrow+nr*jcol                                         3d26s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         else if(isblk(1,is).eq.isblkk(3,isk).and.                      3d26s18
     $        isblk(2,is).eq.isblkk(2,isk).and.                         3d27s18
     $        isblk(3,is).eq.isblkk(1,isk))then                         3d27s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d26s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d26s18
           iswitch=0                                                    3d26s18
          else                                                          3d26s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                          3d20s18
           iswitch=1                                                    3d26s18
          end if                                                        3d26s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d20s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             ibot=i2m                                                   3d27s18
            else                                                        3d26s18
             ibot=0                                                     3d27s18
            end if                                                      3d26s18
            do i3=ibot,noc(isblkk(3,isk))-1                             3d27s18
             ieq=((i3*(i3+1))/2)+i2m                                    3d27s18
             inot=i3+noc(isblk(1,is))*i2m                               3d26s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i4=0,noc(isblkk(4,isk))-1                                3d20s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i1m+noc(isblk(3,is))*i4                              3d27s18
              iad2=jrow+nr*jcol                                         3d26s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         end if
         if(isblk(1,is).eq.isblkk(1,isk).and.                           3d26s18
     $        isblk(2,is).eq.isblkk(4,isk).and.                         3d20s18
     $        isblk(3,is).eq.isblkk(3,isk))then                         3d26s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d20s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d20s18
           iswitch=0                                                    3d20s18
          else                                                          3d20s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                          3d20s18
           iswitch=1                                                    3d20s18
          end if                                                        3d20s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d20s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             itop=i1m                                                   3d26s18
            else                                                        3d26s18
             itop=noc(isblkk(4,isk))-1                                  3d26s18
            end if                                                      3d26s18
            do i4=0,itop                                                3d26s18
             ieq=((i1m*i1)/2)+i4                                        3d26s18
             inot=i1m+noc(isblk(1,is))*i4                               3d26s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i3=0,noc(isblkk(3,isk))-1                               3d26s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i3+noc(isblk(3,is))*i2m                              3d26s18
              iad2=jrow+nr*jcol                                         3d26s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         else if(isblk(1,is).eq.isblkk(2,isk).and.                      3d26s18
     $        isblk(2,is).eq.isblkk(3,isk).and.                         3d20s18
     $        isblk(3,is).eq.isblkk(4,isk))then                         3d26s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d26s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d26s18
           iswitch=0                                                    3d26s18
          else                                                          3d26s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                          3d20s18
           iswitch=1                                                    3d26s18
          end if                                                        3d26s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d20s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             itop=i2m                                                   3d26s18
            else                                                        3d26s18
             itop=noc(isblkk(3,isk))-1                                  3d26s18
            end if                                                      3d26s18
            do i3=0,itop                                                3d26s18
             ieq=((i2m*i2)/2)+i3                                        3d26s18
             inot=i2m+noc(isblk(1,is))*i3                               3d26s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i4=0,noc(isblkk(4,isk))-1                                3d20s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i4+noc(isblk(3,is))*i1m                              3d26s18
              iad2=jrow+nr*jcol                                         3d26s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         end if
         if(isblk(1,is).eq.isblkk(4,isk).and.                           3d26s18
     $        isblk(2,is).eq.isblkk(1,isk).and.                         3d26s18
     $        isblk(3,is).eq.isblkk(3,isk))then                         3d26s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d20s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d20s18
           iswitch=0                                                    3d20s18
          else                                                          3d20s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                          3d20s18
           iswitch=1                                                    3d20s18
          end if                                                        3d20s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d20s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             ibot=i1m                                                   3d26s18
            else                                                        3d26s18
             ibot=0                                                     3d26s18
            end if                                                      3d26s18
            do i4=ibot,noc(isblkk(4,isk))-1                             3d26s18
             ieq=((i4*(i4+1))/2)+i1m                                    3d26s18
             inot=i4+noc(isblk(1,is))*i1m                               3d26s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i3=0,noc(isblkk(3,isk))-1                               3d26s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i3+noc(isblk(3,is))*i2m                              3d26s18
              iad2=jrow+nr*jcol                                         3d26s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         else if(isblk(1,is).eq.isblkk(3,isk).and.                      3d26s18
     $        isblk(2,is).eq.isblkk(2,isk).and.                         3d27s18
     $        isblk(3,is).eq.isblkk(4,isk))then                         3d26s18
          i10=i1s                                                       3d20s18
          if(isblk(1,is).eq.isblk(2,is))then                            3d26s18
           nr=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 3d26s18
           iswitch=0                                                    3d26s18
          else                                                          3d26s18
           nr=noc(isblk(1,is))*noc(isblk(2,is))                          3d20s18
           iswitch=1                                                    3d26s18
          end if                                                        3d26s18
          i1n=noc(isblkk(1,isk))                                        3d20s18
          do i2=i2s,i2e                                                 3d20s18
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e                                         3d20s18
           do i1=i10,i1n                                                3d20s18
            i1m=i1-1                                                    3d20s18
            irow=itmpt+i1+noc(isblkk(1,isk))*i2m-il                     3d20s18
            if(isblk(1,is).eq.isblk(2,is))then                          3d26s18
             ibot=i2m                                                   3d27s18
            else                                                        3d26s18
             ibot=0                                                     3d27s18
            end if                                                      3d26s18
            do i3=ibot,noc(isblkk(3,isk))-1                             3d27s18
             ieq=((i3*(i3+1))/2)+i2m                                    3d27s18
             inot=i3+noc(isblk(1,is))*i2m                               3d27s18
             jrow=ioooo2(is)+(inot-ieq)*iswitch+ieq                     3d26s18
             do i4=0,noc(isblkk(4,isk))-1                                3d20s18
              icol=i3+noc(isblkk(3,isk))*i4                             3d20s18
              iad1=irow+nhere*icol                                      3d20s18
              jcol=i4+noc(isblk(3,is))*i1m                              3d27s18
              iad2=jrow+nr*jcol                                         3d26s18
              bc(iad2)=bc(iad2)+bc(iad1)                                3d20s18
             end do                                                     3d20s18
            end do                                                      3d20s18
           end do                                                       3d20s18
           i10=1                                                        3d20s18
          end do                                                        3d20s18
         end if
        end do                                                          3d20s18
        if(isblkk(3,isk).ne.isblkk(4,isk))then                          3d29s18
c
c     other onex2 contribution.                                         3d29s18
c     recall itmp has full indices in umat basis, with last two indices 3d29s18
c     swapped. Now unswap them and transform last index into tmat basis
c
         do i3=0,nvirtc(isblkk(3,isk))-1                                 3d29s18
          i3p=i3+noc(isblkk(3,isk))                                      3d29s18
          do i4=0,noc(isblkk(4,isk))-1                                   3d29s18
           iad1=itmp+nhere*(i4+nbasdwsc(isblkk(4,isk))*i3p)              3d29s18
           iad2=iokx(isk)+nhere*(i3+nvirtc(isblkk(3,isk))*i4)            3d29s18
           do i12=0,nhere-1                                              3d29s18
            bc(iad2+i12)=bc(iad1+i12)                                    3d29s18
           end do                                                        3d29s18
          end do                                                         3d29s18
         end do                                                          3d29s18
         nr=nhere*nvirtc(isblkk(3,isk))                                  3d29s18
         if(nr.gt.0.and.noc(isblkk(4,isk)).gt.0)then                    2d25s19
         call dgemm('n','n',nr,noc(isblkk(4,isk)),noc(isblkk(4,isk)),1d03d29s18
     $        ,bc(iokx(isk)),nr,                                        3d29s18
     $        bc(itmat(isblkk(4,isk))),noc(isblkk(4,isk)),0d0,           3d29s18
     $        bc(itmp),nr,                                               3d29s18
     d' updateg. 45')
         end if                                                         2d25s19
c
c     currently under itmp we K_{nm}^{Ab}=(nb|mA)
c     and we want to add this into (nb|mA) and (bn|mA),
c     but first 2 indices of
c     onex2 are distributed while nm are distributed for K. So we need
c     full onex3
c
         do isx=1,nsdlk1                                                 3d29s18
          if(isblk1(1,isx).eq.isblkk(1,isk).and.
     $         isblk1(2,isx).eq.isblkk(4,isk).and.                      3d30s18
     $         isblk1(3,isx).eq.isblkk(2,isk))then                      3d30s18
           do i4=0,noc(isblkk(4,isk))-1                                  3d29s18
            do i3=0,nvirtc(isblkk(3,isk))-1                              3d29s18
             i10=i1s                                                    3d30s18
             i1n=noc(isblkk(1,isk))                                     3d30s18
             do i2=i2s,i2e                                              3d30s18
              i2m=i2-1                                                  3d30s18
              if(i2.eq.i2e)i1n=i1e                                      3d30s18
              do i1=i10,i1n                                             3d30s18
               i1m=i1-1                                                 3d30s18
               iad1=itmp+i1+noc(isblkk(1,isk))*i2m-il                   3d30s18
     $              +nhere*(i3+nvirtc(isblkk(3,isk))*i4)                3d30s18
               iad2=ionex3(isx)+i1m+noc(isblk1(1,isx))*(i4              3d30s18
     $              +noc(isblk1(2,isx))*(i2m+noc(isblk1(3,isx))*i3))    3d30s18
               bc(iad2)=bc(iad2)+bc(iad1)                               3d30s18
              end do                                                    3d30s18
              i10=1                                                     3d30s18
             end do                                                      3d29s18
            end do                                                       3d29s18
           end do                                                        3d29s18
          end if                                                         3d29s18
          if(isblk1(1,isx).eq.isblkk(4,isk).and.                        3d30s18
     $       isblk1(2,isx).eq.isblkk(1,isk).and.                        3d29s18
     $       isblk1(3,isx).eq.isblkk(2,isk))then                        3d30s18
           do i4=0,noc(isblkk(4,isk))-1                                 3d29s18
            do i3=0,nvirtc(isblkk(3,isk))-1                             3d29s18
             i10=i1s                                                    3d29s18
             i1n=noc(isblkk(1,isk))                                     3d29s18
             do i2=i2s,i2e                                              3d29s18
              i2m=i2-1
              if(i2.eq.i2e)i1n=i1e                                      3d29s18
              do i1=i10,i1n                                             3d29s18
               i1m=i1-1
               iad1=itmp+i1+noc(isblkk(1,isk))*i2m-il                   3d29s18
     $              +nhere*(i3+nvirtc(isblkk(3,isk))*i4)                3d29s18
               iad2=ionex3(isx)+i4+noc(isblk1(1,isx))*(                 3d30s18
     $              i1m+noc(isblk1(2,isx))*(i2m+noc(isblk1(3,isx))*i3)) 3d30s18
               bc(iad2)=bc(iad2)+bc(iad1)                               3d29s18
              end do                                                    3d29s18
              i10=1                                                     3d29s18
             end do                                                     3d29s18
            end do                                                      3d29s18
           end do                                                       3d29s18
          end if                                                        3d29s18
         end do                                                          3d29s18
        end if                                                          3d29s18
       end do                                                           3d20s18
      end if                                                            3d16s18
c
c     for onex. need to grab parts of 4o first.                         1d27s18
c
      if(c1.ne.0d0)then                                                 1d27s18
       do is=1,nsdlk1                                                   1d27s18
        if(noc(isblk1(1,is)).ne.0.and.noc(isblk1(2,is)).ne.0.and.       1d28s18
     $      noc(isblk1(3,is)).ne.0.and.nbasdwsc(isblk1(4,is)).ne.0)then 1d27s18
         call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,      1d27s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           1d27s18
         ncol=noc(isblk1(3,is))*nbasdwsc(isblk1(4,is))                  1d27s18
         nhere=ih+1-il                                                  1d27s18
         ncx=noc(isblk1(3,is))*nvirtc(isblk1(4,is))
         ncn=noc(isblk1(1,is))*noc(isblk1(2,is))
         itmp1=ibcoff                                                   1d27s18
         itmp2=itmp1+nhere*ncol                                         1d27s18
         ibcoff=itmp2+nhere*ncol                                        1d27s18
         do i=0,nhere*ncol-1
          bc(itmp1+i)=0d0
         end do
         call enough('updateg. 36',bc,ibc)
         jj=iovnn(is)                                                   1d27s18
         do i2=0,nvirtc(isblk1(4,is))-1                                 1d27s18
          i2p=i2+noc(isblk1(4,is))                                      1d27s18
          do i1=0,noc(isblk1(3,is))-1                                   1d27s18
           iad1=itmp1+nhere*(i1+noc(isblk1(3,is))*i2p)                  1d27s18
           do i34=0,nhere-1                                             1d27s18
            bc(iad1+i34)=bc(jj)                                         1d27s18
            jj=jj+1                                                     1d27s18
           end do                                                       1d27s18
          end do                                                        1d27s18
         end do                                                         1d27s18
         ncx=noc(isblk1(3,is))*nbasdwsc(isblk1(4,is))
         do is4=1,nsdlk                                                 1d27s18
          if(isblk1(1,is).eq.isblk(1,is4).and.                          1d27s18
     $         isblk1(2,is).eq.isblk(2,is4))then                        1d28s18
           if(isblk1(3,is).eq.isblk(3,is4).and.                         1d27s18
     $         isblk1(4,is).eq.isblk(4,is4))then                        1d27s18
            do i2=0,noc(isblk1(4,is))-1                                 1d27s18
             do i1=0,noc(isblk1(3,is))-1                                1d27s18
              iad1=itmp1+nhere*(i1+noc(isblk1(3,is))*i2)                1d27s18
              iad2=ioonn(is4)+nhere*(i1+noc(isblk1(3,is))*i2)           1d27s18
              do i34=0,nhere-1                                          1d27s18
               bc(iad1+i34)=bc(iad2+i34)                                1d27s18
              end do                                                    1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
            go to 27                                                    1d27s18
           else if(isblk1(3,is).eq.isblk(4,is4).and.                    1d27s18
     $         isblk1(4,is).eq.isblk(3,is4))then                        1d27s18
            do i2=0,noc(isblk1(4,is))-1                                 1d27s18
             do i1=0,noc(isblk1(3,is))-1                                1d27s18
              iad1=itmp1+nhere*(i1+noc(isblk1(3,is))*i2)                1d27s18
              iad2=ioonn(is4)+nhere*(i2+noc(isblk1(4,is))*i1)           2d7s18
              do i34=0,nhere-1                                          1d27s18
               bc(iad1+i34)=bc(iad2+i34)                                1d27s18
              end do                                                    1d27s18
             end do                                                     1d27s18
            end do                                                      1d27s18
            go to 27                                                    1d27s18
           end if                                                       1d27s18
          end if                                                        1d27s18
         end do                                                         1d27s18
   27    continue                                                       1d27s18
c                                                                       1d27s18
c     we will mult by umat rather than itotm for we need the v part     1d27s18
c     with umat for occ-virt derivative.                                1d27s18
c
         nr=nhere*noc(isblk1(3,is))                                     1d27s18
         nb=nbasdwsc(isblk1(4,is))                                      1d27s18
         no=noc(isblk1(3,is))                                           1d27s18
         if(nr.gt.0.and.nb.gt.0)then                                    2d25s19
         call dgemm('n','n',nr,nb,nb,c1,bc(itmp1),nr,                   1d27s18
     $        bc(iumat(isblk1(4,is))),nb,0d0,bc(itmp2),nr,              1d27s18
     d' updateg. 46')
         end if                                                         2d25s19
         do i2=0,nb-1                                                   1d27s18
          do i1=0,no-1                                                  1d27s18
           iad1=itmp1+nhere*(i2+nb*i1)                                  1d27s18
           iad2=itmp2+nhere*(i1+noc(isblk1(3,is))*i2)                   1d27s18
           do i34=0,nhere-1                                             1d27s18
            bc(iad1+i34)=bc(iad2+i34)                                   1d27s18
           end do                                                       1d27s18
          end do                                                        1d27s18
         end do                                                         1d27s18
         nr=nhere*nb                                                    1d27s18
         if(nr.gt.0.and.no.gt.0)then                                    2d25s19
         call dgemm('n','n',nr,no,no,1d0,bc(itmp1),nr,                  1d27s18
     $        bc(itmat(isblk1(3,is))),no,0d0,bc(itmp2),nr,              1d27s18
     d' updateg. 47')
         end if                                                         2d25s19
         nv=nvirtc(isblk1(4,is))                                        1d27s18
         do i2=0,nv-1                                                   1d27s18
          i2p=i2+noc(isblk1(4,is))                                      1d27s18
          do i1=0,no-1                                                  1d27s18
           iad1=iovnn(is)+nhere*(i1+no*i2)                              1d27s18
           iad2=itmp2+nhere*(i2p+nb*i1)                                 1d27s18
           do i34=0,nhere-1                                             1d27s18
            bc(iad1+i34)=bc(iad2+i34)                                   1d27s18
           end do                                                       1d27s18
          end do                                                        1d27s18
         end do                                                         1d27s18
c
c     higher order contributions
c
         if(c2.ne.0d0)then                                              3d16s18
          ntot=nhere*no*nv                                              3d16s18
          nrow=noc(isblk1(1,is))*noc(isblk1(2,is))                      3d29s18
          call dws_gsumf(bc(ionex3(is)),nrow*no*nv)                     3d29s18
          do i34=0,no*nv-1                                              3d29s18
           iad1=ionex2(is)+nhere*i34                                    3d29s18
           iad2=ionex3(is)+nrow*i34+il-1                                3d29s18
           do i12=0,nhere-1                                             3d29s18
            bc(iad1+i12)=bc(iad1+i12)+bc(iad2+i12)                      3d29s18
           end do                                                       3d29s18
          end do                                                        3d29s18
          do i=0,ntot-1                                                 3d16s18
           bc(iovnn(is)+i)=bc(iovnn(is)+i)+bc(ionex2(is)+i)             3d16s18
          end do                                                        3d16s18
         end if                                                         3d16s18
         ncx=nvirtc(isblk1(3,is))*noc(isblk1(4,is))
         do i2=0,noc(isblk1(4,is))-1                                    1d27s18
          do i1=0,no-1                                                  1d27s18
           iad1=itmp1+nhere*(i1+no*i2)                                  1d27s18
           iad2=itmp2+nhere*(i2+nbasdwsc(isblk1(4,is))*i1)                   1d27s18
           do i34=0,nhere-1                                             1d27s18
            bc(iad1+i34)=bc(iad2+i34)                                   1d27s18
           end do                                                       1d27s18
          end do                                                        1d27s18
         end do                                                         1d27s18
         nr=nhere*no                                                    1d27s18
         no4=noc(isblk1(4,is))                                          1d27s18
         ncx=noc(isblk1(3,is))*noc(isblk1(4,is))
         if(nr.gt.0.and.no4.gt.0)then                                   2d25s19
         call dgemm('n','n',nr,no4,no4,1d0,bc(itmp1),nr,                1d27s18
     $        bc(itmat(isblk1(4,is))),no4,0d0,bc(itmp2),nr,             1d27s18
     d' updateg. 48')
         end if                                                         2d25s19
         if(isblk1(3,is).eq.isblk1(4,is))then                           2d7s18
          do i4=0,no-1                                                  2d7s18
           do i3=0,no-1                                                 2d7s18
            iad1=itmp2+nhere*(i3+no*i4)                                 2d7s18
            iad2=itmp2+nhere*(i4+no*i3)                                 2d7s18
            do i12=0,nhere-1                                            2d7s18
             sum=0.5d0*(bc(iad1+i12)+bc(iad2+i12))                      2d7s18
             bc(iad1+i12)=sum                                           2d7s18
             bc(iad2+i12)=sum                                           2d7s18
            end do                                                      2d7s18
           end do                                                       2d7s18
          end do                                                        2d7s18
         end if                                                         2d7s18
         do is4=1,nsdlk                                                 1d27s18
          if(isblk(1,is4).eq.isblk(2,is4))then                          2d6s18
           iswitch=0                                                    2d6s18
           nrowo=(noc(isblk(1,is4))*(noc(isblk(1,is4))+1))/2            2d6s18
          else                                                          2d6s18
           nrowo=noc(isblk(1,is4))*noc(isblk(2,is4))                    2d6s18
           iswitch=1                                                    2d6s18
          end if                                                        2d6s18
          if(isblk(1,is4).eq.isblk1(4,is).and.                          1d27s18
     $       isblk(2,is4).eq.isblk1(3,is))then                          1d28s18
           if(isblk(3,is4).eq.isblk1(1,is))then                         1d27s18
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            na1=na1+1
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              do i4=0,noc(isblk1(4,is))-1                               1d28s18
               i3top=(noc(isblk1(3,is))-1)*iswitch+i4*(1-iswitch)       2d6s18
               do i3=0,i3top                                            2d6s18
                ix=max(i3,i4)                                           2d6s18
                in=min(i3,i4)                                           2d6s18
                ieq=((ix*(ix+1))/2)+in                                  2d6s18
                inot=i4+noc(isblk(1,is4))*i3                            2d6s18
                iad1=ioooo2(is4)+(inot-ieq)*iswitch+ieq                 2d6s18
     $               +nrowo*(i1m+noc(isblk(3,is4))*i2m)                 2d6s18
                iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                bc(iad1)=bc(iad1)+bc(iad2)                              2d7s18
               end do                                                   1d28s18
              end do                                                    1d28s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           else if(isblk(4,is4).eq.isblk1(1,is))then                    1d27s18
            na2=na2+1
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              do i4=0,noc(isblk1(4,is))-1                               1d28s18
               i3top=(noc(isblk1(3,is))-1)*iswitch+i4*(1-iswitch)       2d6s18
               do i3=0,i3top                                            2d6s18
                ix=max(i4,i3)                                           2d6s18
                in=min(i4,i3)                                           2d6s18
                ieq=((ix*(ix+1))/2)+in                                  2d6s18
                inot=i4+noc(isblk(1,is4))*i3                            2d6s18
                iad1=ioooo2(is4)+(inot-ieq)*iswitch+ieq                 2d6s18
     $               +nrowo*(i2m+noc(isblk(3,is4))*i1m)                 2d6s18
                iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                bc(iad1)=bc(iad1)+bc(iad2)                              2d7s18
               end do                                                   1d28s18
              end do                                                    1d28s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           end if                                                       1d27s18
          end if
          if(isblk(2,is4).eq.isblk1(4,is).and.                          1d27s18
     $       isblk(1,is4).eq.isblk1(3,is))then                          1d28s18
           if(isblk(3,is4).eq.isblk1(1,is))then                         1d28s18
            na3=na3+1
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              do i4=0,noc(isblk1(4,is))-1                               1d28s18
               i3top=(noc(isblk1(3,is))-1)*iswitch+i4*(1-iswitch)       2d6s18
               do i3=0,i3top                                            2d6s18
                ix=max(i3,i4)                                           2d6s18
                in=min(i3,i4)                                           2d6s18
                ieq=((ix*(ix+1))/2)+in                                  2d6s18
                inot=i3+noc(isblk(1,is4))*i4                            2d6s18
                iad1=ioooo2(is4)+(inot-ieq)*iswitch+ieq                 2d6s18
     $               +nrowo*(i1m+noc(isblk(3,is4))*i2m)                 2d6s18
                iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                bc(iad1)=bc(iad1)+bc(iad2)                              2d7s18
               end do                                                   1d28s18
              end do                                                    1d28s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           else if(isblk(3,is4).eq.isblk1(2,is))then                    1d28s18
            na4=na4+1
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              do i4=0,noc(isblk1(4,is))-1                               1d28s18
               i3top=(noc(isblk1(3,is))-1)*iswitch+i4*(1-iswitch)       2d6s18
               do i3=0,i3top                                            2d6s18
                ix=max(i3,i4)                                           2d6s18
                in=min(i3,i4)                                           2d6s18
                ieq=((ix*(ix+1))/2)+in                                  2d6s18
                inot=i3+noc(isblk(1,is4))*i4                            2d6s18
                iad1=ioooo2(is4)+(inot-ieq)*iswitch+ieq                 2d6s18
     $               +nrowo*(i2m+noc(isblk(3,is4))*i1m)                 2d6s18
                iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                bc(iad1)=bc(iad1)+bc(iad2)                              2d7s18
               end do                                                   1d28s18
              end do                                                    1d28s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           end if                                                       1d28s18
          end if
          if(isblk(3,is4).eq.isblk1(4,is).and.                          1d27s18
     $       isblk(4,is4).eq.isblk1(3,is))then                          1d28s18
           if(isblk(1,is4).eq.isblk1(1,is))then                         1d28s18
            nb1=nb1+1
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              if((iswitch.eq.0.and.i1m.ge.i2m).or.iswitch.eq.1)then     2d6s18
               ieq=((i1m*(i1m+1))/2)+i2m                                2d6s18
               inot=i1m+noc(isblk(1,is4))*i2m                           2d6s18
               iii=(inot-ieq)*iswitch+ieq                               2d6s18
               do i4=0,noc(isblk1(4,is))-1                               1d28s18
                do i3=0,noc(isblk1(3,is))-1                              1d28s18
                 iad1=ioooo2(is4)+iii+nrowo*(i4+noc(isblk(3,is4))*i3)   2d6s18
                 iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                 bc(iad1)=bc(iad1)+bc(iad2)                             2d7s18
                end do                                                   1d28s18
               end do                                                    1d28s18
              end if                                                    2d6s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           else if(isblk(1,is4).eq.isblk1(2,is))then                    1d28s18
            nb2=nb2+1
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              if((iswitch.eq.0.and.i1m.ge.i2m).or.iswitch.eq.1)then     2d6s18
               ieq=((i1m*(i1m+1))/2)+i2m                                2d6s18
               inot=i2m+noc(isblk(1,is4))*i1m                           2d6s18
               iii=(inot-ieq)*switch+ieq                                2d6s18
               do i4=0,noc(isblk1(4,is))-1                               1d28s18
                do i3=0,noc(isblk1(3,is))-1                              1d28s18
                 iad1=ioooo2(is4)+iii+nrowo*(i4+noc(isblk(3,is4))*i3)   2d6s18
                 iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                 bc(iad1)=bc(iad1)+bc(iad2)                             2d7s18
                end do                                                   1d28s18
               end do                                                    1d28s18
              end if                                                    2d6s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           end if                                                       1d28s18
          end if
          if(isblk(4,is4).eq.isblk1(4,is).and.                          1d27s18
     $       isblk(3,is4).eq.isblk1(3,is))then                           1d27s18
           if(isblk(1,is4).eq.isblk1(1,is))then                         1d28s18
            nb3=nb3+1
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              if((iswitch.eq.0.and.i1m.ge.i2m).or.iswitch.eq.1)then     2d6s18
               ieq=((i1m*(i1m+1))/2)+i2m                                2d6s18
               inot=i1m+noc(isblk(1,is4))*i2m                           2d6s18
               iii=(inot-ieq)*iswitch+ieq                               2d6s18
               do i4=0,noc(isblk1(4,is))-1                               1d28s18
                do i3=0,noc(isblk1(3,is))-1                              1d28s18
                 iad1=ioooo2(is4)+iii+nrowo*(i3+noc(isblk(3,is4))*i4)   2d6s18
                 iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                 bc(iad1)=bc(iad1)+bc(iad2)                             2d7s18
                end do                                                   1d28s18
               end do                                                    1d28s18
              end if                                                    2d6s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           else if(isblk(1,is4).eq.isblk1(2,is))then                    1d28s18
            nb4=nb4+1
            i10=i1s                                                     1d28s18
            i1n=noc(isblk1(1,is))                                       1d28s18
            irow=0                                                      1d28s18
            do i2=i2s,i2e                                               1d28s18
             i2m=i2-1                                                   1d28s18
             if(i2.eq.i2e)i1n=i1e                                       1d28s18
             do i1=i10,i1n                                              1d28s18
              i1m=i1-1                                                  1d28s18
              if((iswitch.eq.0.and.i1m.ge.i2m).or.iswitch.eq.1)then     2d6s18
               ieq=((i1m*(i1m+1))/2)+i2m                                2d6s18
               inot=i2m+noc(isblk(1,is4))*i1m                           2d6s18
               iii=(inot-ieq)*iswitch+ieq                               2d6s18
               do i4=0,noc(isblk1(4,is))-1                               1d28s18
                do i3=0,noc(isblk1(3,is))-1                              1d28s18
                 iad1=ioooo2(is4)+iii+nrowo*(i3+noc(isblk(3,is4))*i4)   2d6s18
                 iad2=itmp2+irow+nhere*(i3+noc(isblk1(3,is))*i4)         1d28s18
                 bc(iad1)=bc(iad1)+bc(iad2)                             2d7s18
                end do                                                   1d28s18
               end do                                                    1d28s18
              end if                                                    2d6s18
              irow=irow+1                                               1d28s18
             end do                                                     1d28s18
             i10=1                                                      1d28s18
            end do                                                      1d28s18
           end if                                                       1d28s18
          end if
         end do                                                         1d27s18
         ibcoff=itmp1                                                   1d28s18
        end if                                                          1d27s18
       end do                                                           1d27s18
      end if                                                            1d27s18
c
c     for 4o.                                                           1d27s18
c
      do is=1,nsdlk                                                     1d11s18
       if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.           1d10s18
     $     noc(isblk(3,is)).ne.0.and.noc(isblk(4,is)).ne.0)then         1d10s18
        call ilimts(noc(isblk(1,is)),noc(isblk(2,is)),mynprocg,          1d10s18
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            1d11s18
        ncol=ih+1-il                                                     1d10s18
        nrow=noc(isblk(3,is))*noc(isblk(4,is))                          1d10s18
        itmp=ibcoff                                                     1d11s18
        ibcoff=itmp+nrow*ncol                                           1d11s18
        call enough('updateg. 37',bc,ibc)
        nr4=ncol*noc(isblk(3,is))                                       1d11s18
        if(nr4.gt.0.and.noc(isblk(4,is)).gt.0)then                      2d25s19
        call dgemm('n','n',nr4,noc(isblk(4,is)),noc(isblk(4,is)),1d0,   1d11s18
     $       bc(ioonn(is)),nr4,bc(itmat(isblk(4,is))),noc(isblk(4,is))  1d11s18
     $       ,0d0,bc(itmp),nr4,                                         1d11s18
     d' updateg. 49')
        end if                                                          2d25s19
        do i4=0,noc(isblk(4,is))-1                                      1d11s18
         do i3=0,noc(isblk(3,is))-1                                     1d11s18
          iad1=itmp+ncol*(i3+noc(isblk(3,is))*i4)                       1d11s18
          iad2=ioonn(is)+ncol*(i4+noc(isblk(4,is))*i3)                  1d11s18
          do i12=0,ncol-1                                               1d11s18
           bc(iad2+i12)=bc(iad1+i12)                                    1d11s18
          end do                                                        1d11s18
         end do                                                         1d11s18
        end do                                                          1d11s18
        nr3=ncol*noc(isblk(4,is))                                       1d11s18
        if(nr3.gt.0.and.noc(isblk(3,is)).gt.0)then                      2d25s19
        call dgemm('n','n',nr3,noc(isblk(3,is)),noc(isblk(3,is)),1d0,   1d11s18
     $       bc(ioonn(is)),nr3,bc(itmat(isblk(3,is))),noc(isblk(3,is))  1d11s18
     $       ,0d0,bc(itmp),nr3,                                         1d11s18
     d' updateg. 50')
        end if                                                          2d25s19
        i10=i1s                                                         1d11s18
        i1n=noc(isblk(1,is))                                            1d11s18
        if(isblk(1,is).eq.isblk(2,is))then                              1d11s18
         nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 1d11s18
         iswitch=0                                                      1d11s18
        else                                                            1d11s18
         nrow=noc(isblk(1,is))*noc(isblk(2,is))                         1d11s18
         iswitch=1                                                      1d11s18
        end if                                                          1d11s18
        ii=itmp                                                         1d11s18
        do i2=i2s,i2e                                                   1d11s18
         if(i2.eq.i2e)i1n=i1e                                           1d11s18
         do i1=i10,i1n                                                  1d11s18
          if(isblk(1,is).ne.isblk(2,is).or.i1.le.i2)then                1d11s18
           ieq=((i2*(i2-1))/2)+i1-1                                     1d11s18
           inot=i1-1+noc(isblk(1,is))*(i2-1)                            1d11s18
           irow=(inot-ieq)*iswitch+ieq                                  1d11s18
           do i4=0,noc(isblk(4,is))-1                                   1d11s18
            do i3=0,noc(isblk(3,is))-1                                  1d11s18
             iad1=ii+ncol*(i4+noc(isblk(4,is))*i3)                      1d11s18
             iad2=ioooo2(is)+irow+nrow*(i3+noc(isblk(3,is))*i4)         1d11s18
             bc(iad2)=bc(iad2)+bc(iad1)*c0                              1d11s18
            end do                                                      1d11s18
           end do                                                       1d11s18
          end if                                                        1d11s18
          ii=ii+1                                                       1d11s18
         end do                                                         1d11s18
         i10=1                                                          1d11s18
        end do                                                          1d11s18
        ibcoff=itmp                                                     1d11s18
       end if
      end do                                                            1d11s18
      call dws_gsumf(bc(ioooo2(1)),ntot2)                               1d10s18
      if(icall.eq.-1)write(6,*)('oooo2 after gs: ')
      do is=1,nsdlk
       if(noc(isblk(1,is)).ne.0.and.noc(isblk(2,is)).ne.0.and.          1d10s18
     $    noc(isblk(3,is)).ne.0.and.noc(isblk(4,is)).ne.0)then          1d10s18
        if(icall.eq.-1)
     $      write(6,*)('integral type '),is,(isblk(j,is),j=1,4)
        if(isblk(1,is).ne.isblk(2,is))then
         nrow=noc(isblk(1,is))*noc(isblk(2,is))
        else
         nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 1d10s18
        end if
        ncol=noc(isblk(3,is))*noc(isblk(4,is))
        if(icall.eq.-1)call prntm2(bc(ioooo2(is)),nrow,ncol,nrow)
        if(isblk(1,is).eq.isblk(2,is).and.icall.eq.-1)then
         nrow2=noc(isblk(1,is))*noc(isblk(1,is))
         itmp=ibcoff
         ibcoff=itmp+noc(isblk(1,is))*noc(isblk(2,is))*ncol
         call enough('updateg. 38',bc,ibc)
         do i=0,ncol-1
          do j=0,noc(isblk(1,is))-1
           do k=0,j
            iad1=ioooo2(is)+((j*(j+1))/2)+k+nrow*i
            iad2=itmp+k+noc(isblk(1,is))*(j+noc(isblk(1,is))*i)
            iad3=itmp+j+noc(isblk(1,is))*(k+noc(isblk(1,is))*i)
            bc(iad2)=bc(iad1)
            bc(iad3)=bc(iad1)
           end do
          end do
         end do
         call prntm2(bc(itmp),nrow2,ncol,nrow2)
         if(nsymb.eq.1.and.call.ne.-1)then
         rms1=0d0
         rms2=0d0
         rms3=0d0
         do i1=0,noc(1)-1
          do i2=0,noc(1)-1
           do i3=0,noc(1)-1
            do i4=0,noc(1)-1
             iad1=itmp+i4+noc(1)*(i3+noc(1)*(i2+noc(1)*i1))
             iad2=itmp+i3+noc(1)*(i4+noc(1)*(i2+noc(1)*i1))
             iad3=itmp+i4+noc(1)*(i3+noc(1)*(i1+noc(1)*i2))
             iad4=itmp+i2+noc(1)*(i1+noc(1)*(i4+noc(1)*i3))
             rms1=rms1+(bc(iad1)-bc(iad2))**2
             rms2=rms2+(bc(iad1)-bc(iad3))**2
             rms3=rms3+(bc(iad1)-bc(iad4))**2
            end do
           end do
          end do
         end do
         rms1=sqrt(rms1)
         rms2=sqrt(rms2)
         rms3=sqrt(rms3)
         write(6,*)('symmetry tests: '),rms1,rms2,rms3
         call printa(bc(itmp),noc,0,noc,0,noc,0,noc,0,bc(ibcoff))
         end if
         ibcoff=itmp
        end if
        call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,         1d10s18
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           1d10s18
        nhere=ih+1-il                                                   1d10s18
        if(nhere.gt.0)then                                              1d11s18
         ioooo2(is)=ioooo2(is)+nrow*(il-1)                               1d10s18
        end if                                                          1d11s18
       end if                                                           1d10s18
      end do                                                            1d10s18
      irms=0                                                            2d16s18
      jh02=ih02                                                         1d8s18
      jh02h=ih02h                                                       1d16s17
      do isb=1,nsymb
         icasa1=0
         icasa2=0
         icasa3=0
         icasa4=0
         icasb1=0
         icasb2=0
         icasb3=0
         icasb4=0
         icasc1=0
         icasc2=0
         icasc3=0
         icasc4=0
         icasd1=0
         icasd2=0
         icasd3=0
         icasd4=0
       nad=iacto(isb)*idoub(isb)
       nvo=nvirtc(isb)*noc(isb)
       ntot=nad+nvo
       ig=ignew(isb)
       ilook=ig+9+3
       if(nad.gt.0)then
        do i=0,nad-1
         bc(ig+i)=0d0
        end do
       end if                                                           1d5s18
       igp=ig+nad                                                       1d5s18
       do i=0,nvo-1
        bc(igp+i)=0d0
       end do
       if(mynowprog.eq.0.and.idoit(1).ne.0)then                         1d8s18
c
c     dd h0 part
c
        do id=0,idoub(isb)-1
         do ia=0,iacto(isb)-1
          iad=jh02+ia+idoub(isb)+nbasdwsc(isb)*id                        1d8s18
          iad2=ig+ia+iacto(isb)*id                                      1d8s18
          bc(iad2)=4d0*bc(iad)                                          1d8s18
         end do
        end do
        do idt=0,idoub(isb)-1                                           4d4s18
         iad1=igp+nvirtc(isb)*idt                                       4d4s18
         do id=0,idoub(isb)-1                                           4d4s18
          iad2=jh02h+noc(isb)+nbasdwsc(isb)*id                          4d4s18
          iad3=itmat(isb)+idt+noc(isb)*id                               4d4s18
          fact=4d0*bc(iad3)                                             4d4s18
          do iv=0,nvirtc(isb)-1                                         4d4s18
           bc(iad1+iv)=bc(iad1+iv)+fact*bc(iad2+iv)                     4d4s18
          end do                                                        4d4s18
         end do                                                         4d4s18
        end do                                                          4d4s18
        do ia=0,iacto(isb)-1                                             1d16s17
         iap=ia+idoub(isb)                                              1d16s17
         iad1=igp+nvirtc(isb)*iap                                       1d16s17
         do id=0,idoub(isb)-1                                           1d16s17
          iad2=jh02h+noc(isb)+nbasdwsc(isb)*id                          1d16s17
          iad3=itmat(isb)+iap+noc(isb)*id                               1d16s17
          fact=4d0*bc(iad3)                                             1d16s17
          do iv=0,nvirtc(isb)-1                                         1d16s17
           bc(iad1+iv)=bc(iad1+iv)+fact*bc(iad2+iv)                     1d16s17
          end do                                                        1d16s17
         end do                                                         1d16s17
        end do                                                          1d16s17
       end if
       if(idoit(2).ne.0)then
c
c     dddd part
c
        ilook=ig
        if(isb.eq.1)ilook=ilook+iacto(1)
        do is=1,nsdlk
         if(isblk(1,is).eq.isb.and.isblk(2,is).eq.isb)then
          is34=isblk(3,is)
          call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,ih,     1d5s18
     $         i1s,i1e,i2s,i2e)                                         1d5s18
          nrow=(noc(isb)*(noc(isb)+1))/2
          i4o=ioooo2(is)
          i10=i1s
          i1n=noc(is34)
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            if(i1.eq.i2.and.i1.le.idoub(is34))then
             do id=0,idoub(isb)-1
              do ia=0,iacto(isb)-1
               iap=ia+idoub(isb)
               iad2=i4o+((iap*(iap+1))/2)+id
               iad1=ig+ia+iacto(isb)*id
               orig=bc(iad1)
               bc(iad1)=bc(iad1)+bc(iad2)*8d0                           1d26s18
              end do
             end do
            end if
            i4o=i4o+nrow
           end do
           i10=1
          end do
         end if
         if(isblk(1,is).eq.isb.and.isblk(3,is).eq.isb)then
          is24=isblk(2,is)
          call ilimts(noc(isb),noc(is24),mynprocg,mynowprog,il,ih,      1d5s18
     $         i1s,i1e,i2s,i2e)                                         1d5s18
          if(isb.eq.is24)then
           nrow=(noc(isb)*(noc(isb)+1))/2
           iswitch=0
          else
           nrow=noc(isb)*noc(is24)
           iswitch=1
          end if
          i4o=ioooo2(is)
          i10=i1s
          i1n=noc(isb)
          do i2=i2s,i2e
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            if(i1.le.idoub(isb).and.i2.le.idoub(is24))then
             i1m=i1-1
             do ia=0,iacto(isb)-1
              iap=ia+idoub(isb)
              ix=max(iap,i2m)
              in=min(iap,i2m)
              ieq=((ix*(ix+1))/2)+in
              inot=iap+noc(isb)*i2m                                     1d8s18
              iad2=i4o+(inot-ieq)*iswitch+ieq
              iad1=ig+ia+iacto(isb)*i1m
              bc(iad1)=bc(iad1)-bc(iad2)*4d0                            1d26s18
             end do
            end if
            i4o=i4o+nrow
           end do
           i10=1
          end do
         else if(isblk(2,is).eq.isb.and.isblk(3,is).eq.isb)then
          is14=isblk(1,is)
          call ilimts(noc(isb),noc(is24),mynprocg,mynowprog,il,ih,      1d5s18
     $         i1s,i1e,i2s,i2e)                                         1d5s18
          nrow=noc(isb)*noc(is14)
          i4o=ioooo2(is)
          i10=i1s
          i1n=noc(isb)
          do i2=i2s,i2e
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            if(i1.le.idoub(isb).and.i2.le.idoub(is14))then
             i1m=i1-1
             do ia=0,iacto(isb)-1
              iad2=i4o+i2m+noc(is14)*(ia+idoub(isb))                    1d8s18
              iad1=ig+ia+iacto(isb)*i1m
              bc(iad1)=bc(iad1)-bc(iad2)*4d0                            1d26s18
             end do
            end if
            i4o=i4o+nrow
           end do
           i10=1
          end do
         else if(isblk(2,is).eq.isb.and.isblk(4,is).eq.isb)then
          is13=isblk(1,is)
          call ilimts(noc(is13),noc(isb),mynprocg,mynowprog,il,ih,      1d5s18
     $         i1s,i1e,i2s,i2e)                                         1d5s18
          nrow=noc(isb)*noc(is13)
          i4o=ioooo2(is)
          i10=i1s
          i1n=noc(is13)
          do i2=i2s,i2e
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            if(i1.le.idoub(is13).and.i2.le.idoub(isb))then
             i1m=i1-1
             do ia=0,iacto(isb)-1
              iad2=i4o+i1m+noc(is13)*(ia+idoub(isb))                    1d8s18
              iad1=ig+ia+iacto(isb)*i2m
              bc(iad1)=bc(iad1)-bc(iad2)*4d0                            1d26s18
             end do
            end if
            i4o=i4o+nrow
           end do
           i10=1
          end do
         else if(isblk(1,is).eq.isb.and.isblk(4,is).eq.isb)then
          is23=isblk(2,is)
          call ilimts(noc(is23),noc(isb),mynprocg,mynowprog,il,ih,      1d5s18
     $         i1s,i1e,i2s,i2e)                                         1d5s18
          nrow=noc(isb)*noc(is23)
          i4o=ioooo2(is)
          i10=i1s
          i1n=noc(is23)
          do i2=i2s,i2e
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            if(i1.le.idoub(is23).and.i2.le.idoub(isb))then
             i1m=i1-1
             do ia=0,iacto(isb)-1
              iad2=i4o+ia+idoub(isb)+noc(isb)*i1m                       1d8s18
              iad1=ig+ia+iacto(isb)*i2m
              bc(iad1)=bc(iad1)-bc(iad2)*4d0                            1d26s18
             end do
            end if
            i4o=i4o+nrow
           end do
           i10=1
          end do
         end if
        end do
c
c     for occ - virt rotations
c
        if(c1.ne.0d0)then                                               2d13s18
         do is=1,nsdlk1                                                 2d13s18
          if(isblk1(4,is).eq.isb.and.isblk1(3,is).eq.isb)then           2d13s18
           is12=isblk1(1,is)                                            2d13s18
           call ilimts(noc(is12),noc(is12),mynprocg,mynowprog,il,ih,    2d13s18
     $          i1s,i1e,i2s,i2e)                                        2d13s18
           nrow=ih+1-il                                                 2d13s18
           i10=i1s                                                      2d13s18
           i1n=noc(is12)                                                2d13s18
           i1x=iovnn(is)                                                2d13s18
           do i2=i2s,i2e                                                2d13s18
            if(i2.eq.i2e)i1n=i1e                                        2d13s18
            do i1=i10,i1n                                               2d13s18
             if(i1.eq.i2.and.i1.le.idoub(is12))then                     2d13s18
              do iv=0,nvirtc(isb)-1                                     2d13s18
               do id=0,idoub(isb)-1                                     2d13s18
                iad2=i1x+nrow*(id+noc(isb)*iv)                          2d13s18
                do idt=0,idoub(isb)-1                                   2d13s18
                 iad1=igp+iv+nvirtc(isb)*idt                            2d13s18
                 iad3=itmat(isb)+idt+noc(isb)*id                        2d13s18
                 bc(iad1)=bc(iad1)+bc(iad2)*bc(iad3)*8d0                2d13s18
                end do                                                  2d13s18
                do iat=0,iacto(isb)-1                                   2d13s18
                 iatp=iat+idoub(isb)                                    2d13s18
                 iad1=igp+iv+nvirtc(isb)*iatp                           2d13s18
                 iad3=itmat(isb)+iatp+noc(isb)*id                       2d13s18
                 bc(iad1)=bc(iad1)+bc(iad2)*bc(iad3)*8d0                2d13s18
                end do                                                  2d13s18
               end do                                                   2d13s18
              end do                                                    2d13s18
             end if                                                     2d13s18
             i1x=i1x+1                                                  2d13s18
            end do                                                      2d13s18
            i10=1                                                       2d13s18
           end do                                                       2d13s18
          end if                                                        2d12s18
          if(isblk1(4,is).eq.isb)then                                   2d13s18
           if(isblk1(2,is).eq.isb)then                                  2d13s18
            call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,   2d13s18
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        2d13s18
            nrow=ih+1-il                                                 2d13s18
            i10=i1s                                                      2d13s18
            i1n=noc(isblk1(1,is))                                       2d13s18
            i1x=iovnn(is)                                                2d13s18
            do i2=i2s,i2e                                                2d13s18
             i2m=i2-1                                                   2d13s18
             if(i2.eq.i2e)i1n=i1e                                        2d13s18
             do i1=i10,i1n                                               2d13s18
              if(i2.le.idoub(isblk1(2,is)).and                          2d13s18
     $             .i1.le.idoub(isblk1(3,is)))then                      2d13s18
               i1m=i1-1                                                  2d13s18
               do iv=0,nvirtc(isb)-1                                     2d13s18
                do idt=0,idoub(isb)-1                                    2d13s18
                 iad1=igp+iv+nvirtc(isb)*idt                             2d13s18
                 iad2=i1x+nrow*(i1m+noc(isblk1(3,is))*iv)                2d13s18
                 iad3=itmat(isb)+idt+noc(isb)*i2m                        2d13s18
                 bc(iad1)=bc(iad1)-4d0*bc(iad3)*bc(iad2)                 2d13s18
                end do                                                   2d13s18
                do iat=0,iacto(isb)-1                                   2d13s18
                 iatp=iat+idoub(isb)                                    2d13s18
                 iad1=igp+iv+nvirtc(isb)*iatp                           2d13s18
                 iad2=i1x+nrow*(i1m+noc(isblk1(3,is))*iv)                2d13s18
                 iad3=itmat(isb)+iatp+noc(isb)*i2m                      2d13s18
                 bc(iad1)=bc(iad1)-4d0*bc(iad3)*bc(iad2)                 2d13s18
                end do                                                  2d13s18
               end do                                                    2d13s18
              end if                                                    2d13s18
              i1x=i1x+1                                                 2d13s18
             end do                                                     2d13s18
             i10=1                                                      2d13s18
            end do                                                      2d13s18
           else if(isblk1(1,is).eq.isb)then                             2d13s18
            call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,   2d13s18
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        2d13s18
            nrow=ih+1-il                                                 2d13s18
            i10=i1s                                                      2d13s18
            i1n=noc(isblk1(1,is))                                       2d13s18
            i1x=iovnn(is)                                                2d13s18
            do i2=i2s,i2e                                                2d13s18
             i2m=i2-1                                                   2d13s18
             if(i2.eq.i2e)i1n=i1e                                        2d13s18
             do i1=i10,i1n                                               2d13s18
              i1m=i1-1                                                  2d13s18
              if(i1.le.idoub(isb).and                                   2d13s18
     $             .i2.le.idoub(isblk1(3,is)))then                      2d13s18
               do iv=0,nvirtc(isb)-1                                     2d13s18
                do idt=0,idoub(isb)-1                                   2d13s18
                 iad1=igp+iv+nvirtc(isb)*idt                             2d13s18
                 iad2=i1x+nrow*(i2m+noc(isblk1(3,is))*iv)                2d13s18
                 iad3=itmat(isb)+idt+noc(isb)*i1m                        2d13s18
                 bc(iad1)=bc(iad1)-4d0*bc(iad3)*bc(iad2)                 2d13s18
                end do                                                   2d13s18
                do iat=0,iacto(isb)-1                                   2d13s18
                 iatp=iat+idoub(isb)                                    2d13s18
                 iad1=igp+iv+nvirtc(isb)*iatp                           2d13s18
                 iad2=i1x+nrow*(i2m+noc(isblk1(3,is))*iv)                2d13s18
                 iad3=itmat(isb)+iatp+noc(isb)*i1m                      2d13s18
                 bc(iad1)=bc(iad1)-4d0*bc(iad3)*bc(iad2)                 2d13s18
                end do                                                   2d13s18
               end do                                                    2d13s18
              end if                                                    2d13s18
              i1x=i1x+1                                                 2d13s18
             end do                                                     2d13s18
             i10=1                                                      2d13s18
            end do                                                      2d13s18
           end if                                                       2d13s18
          end if                                                        2d13s18
         end do                                                         2d12s18
        end if                                                          2d13s18
       end if
       if(idoit(3).ne.0.and.mynowprog.eq.0)then                         1d8s18
c
c     aa part
c
        do id=0,idoub(isb)-1                                              1d8s18
         do ia=0,iacto(isb)-1                                             1d8s18
          iap=ia+idoub(isb)                                             1d8s18
          iad=ig+ia+iacto(isb)*id                                       1d8s18
          iad2=iden(isb)+iacto(isb)*ia                                  1d8s18
          iad3=jh02+idoub(isb)+nbasdwsc(isb)*id                         1d8s18
          do ja=0,iacto(isb)-1                                          1d8s18
           bc(iad)=bc(iad)-2d0*bc(iad3+ja)*bc(iad2+ja)                  1d8s18
          end do                                                        1d8s18
         end do                                                         1d8s18
        end do                                                          1d8s18
        do ia2=0,iacto(isb)-1                                           1d26s18
         iad4=iden(isb)+iacto(isb)*ia2                                  1d26s18
         iad2=itmat(isb)+noc(isb)*(ia2+idoub(isb))                      1d26s18
         do ia3=0,iacto(isb)-1                                          1d26s18
          fact1=2d0*bc(iad4+ia3)                                        1d26s18
          iad3=jh02h+noc(isb)+nbasdwsc(isb)*(ia3+idoub(isb))            1d26s18
          do id=0,noc(isb)-1                                            4d5s18
           iad1=igp+nvirtc(isb)*id                                      1d26s18
           fact2=fact1*bc(iad2+id)                                      1d26s18
           do iv=0,nvirtc(isb)-1                                        1d26s18
            orig=bc(iad1+iv)
            bc(iad1+iv)=bc(iad1+iv)+fact2*bc(iad3+iv)                   1d26s18
           end do                                                       1d26s18
          end do                                                        1d26s18
         end do                                                         1d26s18
        end do                                                          1d26s18
       end if                                                           1d8s18
       if(idoit(4).ne.0)then                                            1d8s18
        if(isb.eq.1)then
         ilook=ig+1
        else
         ilook=-1
        end if
c
c     ddaa part
c
        do is=1,nsdlk                                                   1d8s18
         if(isblk(1,is).eq.isb.and.isblk(2,is).eq.isb)then              1d8s18
          is34=isblk(3,is)                                              1d8s18
          call ilimts(noc(is34),noc(is34),mynprocg,mynowprog,il,ih,     1d8s18
     $         i1s,i1e,i2s,i2e)                                         1d8s18
          i10=i1s                                                       1d8s18
          i1n=noc(is34)                                                 1d8s18
          i4o=ioooo2(is)                                                1d8s18
          nrow=(noc(isb)*(noc(isb)+1))/2                                1d8s18
          do i2=i2s,i2e                                                 1d8s18
           i2m=i2-1-idoub(is34)                                         1d8s18
           if(i2.eq.i2e)i1n=i1e                                         1d8s18
           do i1=i10,i1n                                                1d8s18
            if(i1.gt.idoub(is34).and.i2.gt.idoub(is34))then             1d8s18
             i1m=i1-1-idoub(is34)                                       1d8s18
             iad2=iden(is34)+i1m+iacto(is34)*i2m                        1d8s18
             do id=0,idoub(isb)-1                                       1d8s18
              do ia=0,iacto(isb)-1                                      1d8s18
               iap=ia+idoub(isb)                                        1d8s18
               iad1=ig+ia+iacto(isb)*id                                 1d8s18
               iad3=i4o+((iap*(iap+1))/2)+id                            1d8s18
               bc(iad1)=bc(iad1)+4d0*bc(iad2)*bc(iad3)                  1d26s18
              end do                                                    1d8s18
             end do                                                     1d8s18
            else if(i1.le.idoub(is34).and.i2.eq.i1)then                 1d11s18
             do id=0,idoub(isb)-1                                       1d11s18
              do ia=0,iacto(isb)-1                                      1d11s18
               iap=ia+idoub(isb)                                        1d11s18
               iad3=i4o+((iap*(iap+1))/2)+id                              1d11s18
               fact=-4d0*bc(iad3)                                       1d26s18
               iad2=iden(isb)+iacto(isb)*ia                             1d11s18
               iad1=ig+iacto(isb)*id                                    1d16s17
               do iat=0,iacto(isb)-1                                    1d11s18
                bc(iad1+iat)=bc(iad1+iat)+fact*bc(iad2+iat)             1d16s17
               end do                                                   1d11s18
              end do                                                    1d11s18
             end do                                                     1d11s18
            end if                                                      1d8s18
            i4o=i4o+nrow                                                1d8s18
           end do                                                       1d8s18
           i10=1                                                        1d8s18
          end do                                                        1d8s18
         end if                                                         1d11s18
         if(isblk(1,is).eq.isb.and.isblk(3,is).eq.isb)then              1d8s18
          is24=isblk(2,is)                                              1d8s18
          call ilimts(noc(isb),noc(is24),mynprocg,mynowprog,il,ih,      1d8s18
     $         i1s,i1e,i2s,i2e)                                         1d8s18
          if(is24.eq.isb)then                                           1d8s18
           nrow=(noc(isb)*(noc(isb)+1))/2                               1d8s18
           iswitch=0                                                    1d8s18
          else                                                          1d8s18
           nrow=noc(isb)*noc(is24)                                      1d8s18
           iswitch=1                                                    1d8s18
          end if                                                        1d8s18
          i10=i1s                                                       1d8s18
          i1n=noc(isb)                                                  1d8s18
          i4o=ioooo2(is)                                                1d8s18
          do i2=i2s,i2e                                                 1d8s18
           if(i2.eq.i2e)i1n=i1e                                         1d8s18
           i2n=i2-1                                                     1d11s18
           i2m=i2n-idoub(is24)                                          1d11s18
           do i1=i10,i1n                                                1d8s18
            if(i1.le.idoub(isb).and.i2.gt.idoub(is24))then              1d8s18
             i1m=i1-1                                                   1d8s18
             do i4=0,iacto(is24)-1                                      1d8s18
              i4p=i4+idoub(is24)                                        1d8s18
              iad2=iden(is24)+i4+iacto(is24)*i2m                        1d8s18
              do i3=0,iacto(isb)-1                                      1d8s18
               i3p=i3+idoub(isb)                                        1d8s18
               ix=max(i3p,i4p)                                          1d8s18
               in=min(i3p,i4p)                                          1d8s18
               ieq=((ix*(ix+1))/2)+in                                   1d8s18
               inot=i3p+noc(isb)*i4p                                    1d8s18
               iad1=i4o+(inot-ieq)*iswitch+ieq                          1d8s18
               iad3=ig+i3+iacto(isb)*i1m                                1d8s18
               bc(iad3)=bc(iad3)-2d0*bc(iad2)*bc(iad1)                  1d26s18
              end do                                                    1d8s18
             end do                                                     1d8s18
            else if(i1.gt.idoub(isb).and.i2.le.idoub(is24))then         1d11s18
             i1m=i1-1-idoub(isb)                                        1d11s18
             iad2=iden(isb)+iacto(isb)*i1m                              1d11s18
             do id=0,idoub(isb)-1                                       1d11s18
              ix=max(id,i2n)                                            1d11s18
              in=min(id,i2n)                                            1d11s18
              ieq=((ix*(ix+1))/2)+in                                    1d11s18
              inot=id+noc(isb)*i2n                                      1d11s18
              iad1=i4o+(inot-ieq)*iswitch+ieq                           1d11s18
              fact=2d0*bc(iad1)                                         1d26s18
              iad3=ig+iacto(isb)*id                                     1d11s18
              do iat=0,iacto(isb)-1                                      1d11s18
               bc(iad3+iat)=bc(iad3+iat)+fact*bc(iad2+iat)              1d11s18
              end do                                                    1d11s18
             end do                                                     1d11s18
            end if                                                      1d8s18
            i4o=i4o+nrow                                                1d8s18
           end do                                                       1d8s18
           i10=1                                                        1d8s18
          end do                                                        1d8s18
         else if(isblk(2,is).eq.isb.and.isblk(3,is).eq.isb)then         1d8s18
          is14=isblk(1,is)                                              1d8s18
          call ilimts(noc(isb),noc(is14),mynprocg,mynowprog,il,ih,      1d11s18
     $         i1s,i1e,i2s,i2e)                                         1d8s18
          nrow=noc(isb)*noc(is14)                                       1d8s18
          i10=i1s                                                       1d8s18
          i1n=noc(isb)                                                  1d8s18
          i4o=ioooo2(is)                                                1d8s18
          do i2=i2s,i2e                                                 1d8s18
           if(i2.eq.i2e)i1n=i1e                                         1d8s18
           i2n=i2-1                                                     1d11s18
           i2m=i2n-idoub(is14)                                          1d11s18
           do i1=i10,i1n                                                1d8s18
            if(i1.le.idoub(isb).and.i2.gt.idoub(is14))then              1d8s18
             i1m=i1-1                                                   1d8s18
             do i4=0,iacto(is14)-1                                      1d8s18
              i4p=i4+idoub(is14)                                        1d8s18
              iad2=iden(is14)+i4+iacto(is14)*i2m                        1d8s18
              do i3=0,iacto(isb)-1                                      1d8s18
               iad1=i4o+i4p+noc(is14)*i3p                               1d8s18
               iad3=ig+i3+iacto(isb)*i1m                                1d8s18
               bc(iad3)=bc(iad3)-2d0*bc(iad2)*bc(iad1)                  1d26s18
              end do                                                    1d8s18
             end do                                                     1d8s18
            else if(i1.gt.idoub(isb).and.i2.le.idoub(is14))then         1d11s18
             i1m=i1-1-idoub(isb)                                        1d11s18
             iad2=iden(isb)+iacto(isb)*i1m                              1d11s18
             do id=0,idoub(isb)-1                                       1d11s18
              iad1=i4o+i2n+noc(is14)*id                                 1d11s18
              fact=2d0*bc(iad1)                                         1d26s18
              iad3=ig+iacto(isb)*id                                     1d11s18
              do iat=0,iacto(isb)-1                                      1d11s18
               bc(iad3+iat)=bc(iad3+iat)+fact*bc(iad2+iat)              1d11s18
              end do                                                    1d11s18
             end do                                                     1d11s18
            end if                                                      1d8s18
            i4o=i4o+nrow                                                1d8s18
           end do                                                       1d8s18
           i10=1                                                        1d8s18
          end do                                                        1d8s18
         else if(isblk(2,is).eq.isb.and.isblk(4,is).eq.isb)then         1d8s18
          is13=isblk(1,is)                                              1d8s18
          call ilimts(noc(is13),noc(isb),mynprocg,mynowprog,il,ih,      1d8s18
     $         i1s,i1e,i2s,i2e)                                         1d8s18
          nrow=noc(isb)*noc(is13)                                       1d8s18
          i10=i1s                                                       1d8s18
          i1n=noc(is13)                                                 1d8s18
          i4o=ioooo2(is)                                                1d8s18
          do i2=i2s,i2e                                                 1d8s18
           if(i2.eq.i2e)i1n=i1e                                         1d8s18
           i2m=i2-1                                                     1d8s18
           do i1=i10,i1n                                                1d8s18
            if(i2.le.idoub(isb).and.i1.gt.idoub(is13))then              1d8s18
             i1m=i1-1-idoub(is13)                                       1d8s18
             do i4=0,iacto(is13)-1                                      1d8s18
              i4p=i4+idoub(is13)                                        1d8s18
              iad2=iden(is13)+i4+iacto(is13)*i1m                        1d8s18
              do i3=0,iacto(isb)-1                                      1d8s18
               i3p=i3+idoub(isb)                                        1d11s18
               iad1=i4o+i4p+noc(is13)*i3p                               1d8s18
               iad3=ig+i3+iacto(isb)*i2m                                1d8s18
               bc(iad3)=bc(iad3)-2d0*bc(iad2)*bc(iad1)                  1d26s18
              end do                                                    1d8s18
             end do                                                     1d8s18
            else if(i2.gt.idoub(isb).and.i1.le.idoub(is13))then         1d11s18
             i2n=i2m-idoub(isb)                                         1d11s18
             i1m=i1-1                                                   1d11s18
             iad2=iden(isb)+iacto(isb)*i2n                              1d11s18
             do id=0,idoub(isb)-1                                       1d11s18
              iad1=i4o+i1m+noc(is13)*id                                 1d11s18
              fact=2d0*bc(iad1)                                         1d26s18
              iad3=ig+iacto(isb)*id                                     1d11s18
              do iat=0,iacto(isb)-1                                      1d11s18
               bc(iad3+iat)=bc(iad3+iat)+fact*bc(iad2+iat)              1d11s18
              end do                                                    1d11s18
             end do                                                     1d11s18
            end if                                                      1d8s18
            i4o=i4o+nrow                                                1d8s18
           end do                                                       1d8s18
           i10=1                                                        1d8s18
          end do                                                        1d8s18
         else if(isblk(1,is).eq.isb.and.isblk(4,is).eq.isb)then         1d8s18
          is23=isblk(2,is)                                              1d8s18
          call ilimts(noc(is23),noc(isb),mynprocg,mynowprog,il,ih,      1d8s18
     $         i1s,i1e,i2s,i2e)                                         1d8s18
          nrow=noc(isb)*noc(is23)                                       1d8s18
          i10=i1s                                                       1d8s18
          i1n=noc(is23)                                                 1d8s18
          i4o=ioooo2(is)                                                1d8s18
          do i2=i2s,i2e                                                 1d8s18
           if(i2.eq.i2e)i1n=i1e                                         1d8s18
           i2m=i2-1                                                     1d8s18
           do i1=i10,i1n                                                1d8s18
            if(i2.le.idoub(isb).and.i1.gt.idoub(is23))then              1d8s18
             i1m=i1-1-idoub(is23)                                       1d8s18
             do i4=0,iacto(is23)-1                                      1d8s18
              i4p=i4+idoub(is23)                                        1d8s18
              iad2=iden(is23)+i4+iacto(is23)*i1m                        1d8s18
              do i3=0,iacto(isb)-1                                      1d8s18
               iad1=i4o+i3p+noc(isb)*i4p                                1d8s18
               iad3=ig+i3+iacto(isb)*i2m                                1d8s18
               bc(iad3)=bc(iad3)-2d0*bc(iad2)*bc(iad1)                  1d26s18
              end do                                                    1d8s18
             end do                                                     1d8s18
            else if(i2.gt.idoub(isb).and.i1.le.idoub(is23))then         1d11s18
             i2n=i2m-idoub(isb)                                         1d11s18
             i1m=i1-1                                                   1d11s18
             iad2=iden(isb)+iacto(isb)*i2n                              1d11s18
             do id=0,idoub(isb)-1                                       1d11s18
              iad1=i4o+id+noc(isb)*i1m                                  1d11s18
              fact=2d0*bc(iad1)                                         1d26s18
              iad3=ig+iacto(isb)*id                                     1d11s18
              do iat=0,iacto(isb)-1                                      1d11s18
               bc(iad3+iat)=bc(iad3+iat)+fact*bc(iad2+iat)              1d11s18
              end do                                                    1d11s18
             end do                                                     1d11s18
            end if                                                      1d8s18
            i4o=i4o+nrow                                                1d8s18
           end do                                                       1d8s18
           i10=1                                                        1d8s18
          end do                                                        1d8s18
         end if                                                         1d8s18
        end do                                                           1d8s18
c                                                                       2d13s18
c     for noc-virt rotations                                            2d13s18
c                                                                       2d13s18
        if(c1.ne.0d0)then                                               2d13s18
         ilook=-1
         if(isb.eq.2)ilook=igp+1-1+15*(3-1)
         do is=1,nsdlk1                                                 2d13s18
          if(isblk1(4,is).eq.isb)then                                   2d13s18
           if(isblk1(3,is).eq.isb)then                                  2d13s18
            is12=isblk1(1,is)                                           2d13s18
            call ilimts(noc(is12),noc(is12),mynprocg,mynowprog,il,ih,   2d13s18
     $           i1s,i1e,i2s,i2e)                                       2d13s18
            nrow=ih+1-il                                                2d13s18
            i10=i1s                                                     2d13s18
            i1n=noc(is12)                                               2d13s18
            i1x=iovnn(is)                                               2d13s18
            do i2=i2s,i2e                                               2d13s18
             if(i2.eq.i2e)i1n=i1e                                       2d13s18
             i2m=i2-idoub(is12)-1                                       2d13s18
             do i1=i10,i1n                                              2d13s18
              if(i1.eq.i2.and.i1.le.idoub(is12))then                    2d13s18
               do iv=0,nvirtc(isb)-1                                    2d13s18
                do ia=0,iacto(isb)-1                                    2d13s18
                 do iap=0,iacto(isb)-1                                  2d13s18
                  iad1=iden(isb)+iap+iacto(isb)*ia                      2d13s18
                  iapp=iap+idoub(isb)                                   2d13s18
                  do idt=0,idoub(isb)-1                                 2d13s18
                   iad2=itmat(isb)+idt+noc(isb)*iapp                    2d13s18
                   iad3=igp+iv+nvirtc(isb)*idt                          2d13s18
                   iad4=i1x+nrow*(ia+idoub(isb)+noc(isb)*iv)            2d13s18
                   bc(iad3)=bc(iad3)+4d0*bc(iad1)*bc(iad2)*bc(iad4)     2d13s18
                  end do                                                2d13s18
                  do iat=0,iacto(isb)-1                                 2d13s18
                   iatp=iat+idoub(isb)                                  2d13s18
                   iad2=itmat(isb)+iatp+noc(isb)*iapp                   2d13s18
                   iad3=igp+iv+nvirtc(isb)*iatp                         2d13s18
                   iad4=i1x+nrow*(ia+idoub(isb)+noc(isb)*iv)            2d13s18
                   bc(iad3)=bc(iad3)+4d0*bc(iad1)*bc(iad2)*bc(iad4)     2d13s18
                  end do                                                2d13s18
                 end do                                                 2d13s18
                end do                                                  2d13s18
               end do                                                   2d13s18
              end if                                                    2d13s18
              if(i1.gt.idoub(is12).and.i2.gt.idoub(is12))then           2d13s18
               i1m=i1-idoub(is12)-1                                     2d13s18
               iad4=iden(is12)+i1m+iacto(is12)*i2m                      2d13s18
               fact=4d0*bc(iad4)                                        2d13s18
               do iv=0,nvirtc(isb)-1                                    2d13s18
                do id=0,idoub(isb)-1                                    2d13s18
                 do idt=0,idoub(isb)-1                                  2d13s18
                  iad1=itmat(isb)+idt+noc(isb)*id                       2d13s18
                  iad2=i1x+nrow*(id+noc(isb)*iv)                        2d13s18
                  iad3=igp+iv+nvirtc(isb)*idt                           2d13s18
                  bc(iad3)=bc(iad3)+fact*bc(iad1)*bc(iad2)              2d13s18
                 end do                                                 2d13s18
                 do iat=0,iacto(isb)-1                                  2d13s18
                  iatp=iat+idoub(isb)                                   2d13s18
                  iad1=itmat(isb)+iatp+noc(isb)*id                      2d13s18
                  iad2=i1x+nrow*(id+noc(isb)*iv)                        2d13s18
                  iad3=igp+iv+nvirtc(isb)*iatp                          2d13s18
                  bc(iad3)=bc(iad3)+fact*bc(iad1)*bc(iad2)              2d13s18
                 end do                                                 2d13s18
                end do                                                  2d13s18
               end do                                                   2d13s18
              end if                                                    2d13s18
              i1x=i1x+1                                                 2d13s18
             end do                                                     2d13s18
             i10=1                                                      2d13s18
            end do                                                      2d13s18
           end if                                                       2d13s18
           if(isblk1(2,is).eq.isb)then                                  2d13s18
            call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,   2d13s18
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        2d13s18
            nrow=ih+1-il                                                2d13s18
            i10=i1s                                                     2d13s18
            i1n=noc(isblk1(1,is))                                       2d13s18
            i1x=iovnn(is)                                               2d13s18
            do i2=i2s,i2e                                               2d13s18
             if(i2.eq.i2e)i1n=i1e                                       2d13s18
             i2n=i2-1                                                   2d13s18
             i2m=i2n-idoub(isb)                                         2d13s18
             do i1=i10,i1n                                              2d13s18
              if(i2.gt.idoub(isb).and.i1.le.idoub(isblk1(1,is)))then    2d13s18
               i1m=i1-1                                                 2d13s18
               do iv=0,nvirtc(isb)-1                                    2d13s18
                iad2=i1x+nrow*(i1m+noc(isblk1(1,is))*iv)                2d13s18
                fmul=-2d0*bc(iad2)                                      2d13s18
                do iap=0,iacto(isb)-1                                   2d13s18
                 iad1=iden(isb)+iap+iacto(isb)*i2m                      2d13s18
                 fact=fmul*bc(iad1)                                     2d13s18
                 iad3=itmat(isb)+noc(isb)*(idoub(isb)+iap)              2d13s18
                 do idt=0,idoub(isb)-1                                  2d13s18
                  iad4=igp+iv+nvirtc(isb)*idt                           2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+idt)                   2d13s18
                 end do                                                 2d13s18
                 do iat=0,iacto(isb)-1                                  2d13s18
                  iatp=iat+idoub(isb)                                   2d13s18
                  iad4=igp+iv+nvirtc(isb)*iatp                          2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+iatp)                  2d13s18
                 end do                                                 2d13s18
                end do                                                  2d13s18
               end do                                                   2d13s18
              end if                                                    2d13s18
              if(i2.le.idoub(isb).and.i1.gt.idoub(isblk1(1,is)))then    2d13s18
               i1m=i1-1-idoub(isblk1(1,is))                             2d13s18
               iad3=itmat(isb)+noc(isb)*i2n                             2d13s18
               do iap=0,iacto(isblk1(3,is))-1                           2d13s18
                iad1=iden(isblk1(3,is))+i1m+iacto(isblk1(3,is))*iap     2d13s18
                fmul=-2d0*bc(iad1)                                      2d13s18
                do iv=0,nvirtc(isb)-1                                    2d13s18
                 iad2=i1x+nrow*(                                        2d13s18
     $                iap+idoub(isblk1(3,is))+noc(isblk1(3,is))*iv)     2d13s18
                 fact=fmul*bc(iad2)                                     2d13s18
                 do idt=0,idoub(isb)-1                                  2d13s18
                  iad4=igp+iv+nvirtc(isb)*idt                           2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+idt)                   2d13s18
                 end do                                                 2d13s18
                 do iat=0,iacto(isb)-1                                  2d13s18
                  iatp=iat+idoub(isb)                                   2d13s18
                  iad4=igp+iv+nvirtc(isb)*iatp                          2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+iatp)                  2d13s18
                 end do                                                 2d13s18
                end do                                                  2d13s18
               end do                                                   2d13s18
              end if                                                    2d13s18
              i1x=i1x+1                                                 2d13s18
             end do                                                     2d13s18
             i10=1                                                      2d13s18
            end do                                                      2d13s18
           else if(isblk1(1,is).eq.isb)then                             2d13s18
            call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,   2d13s18
     $          mynowprog,il,ih,i1s,i1e,i2s,i2e)                        2d13s18
            nrow=ih+1-il                                                2d13s18
            i10=i1s                                                     2d13s18
            i1n=noc(isblk1(1,is))                                       2d13s18
            i1x=iovnn(is)                                               2d13s18
            do i2=i2s,i2e                                               2d13s18
             if(i2.eq.i2e)i1n=i1e                                       2d13s18
             i2m=i2-1                                                   2d13s18
             i2n=i2m-idoub(isblk1(2,is))                                2d13s18
             do i1=i10,i1n                                              2d13s18
              if(i1.gt.idoub(isb).and.i2.le.idoub(isblk1(2,is)))then    2d14s18
               i1m=i1-1-idoub(isb)                                      2d13s18
               do iv=0,nvirtc(isb)-1                                    2d13s18
                iad2=i1x+nrow*(i2m+noc(isblk1(2,is))*iv)                2d13s18
                fmul=-2d0*bc(iad2)                                      2d13s18
                do iap=0,iacto(isb)-1                                   2d13s18
                 iad1=iden(isb)+iap+iacto(isb)*i1m                      2d13s18
                 fact=fmul*bc(iad1)                                     2d13s18
                 iad3=itmat(isb)+noc(isb)*(idoub(isb)+iap)              2d13s18
                 do idt=0,idoub(isb)-1                                  2d13s18
                  iad4=igp+iv+nvirtc(isb)*idt                           2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+idt)                   2d13s18
                 end do                                                 2d13s18
                 do iat=0,iacto(isb)-1                                  2d13s18
                  iatp=iat+idoub(isb)                                   2d13s18
                  iad4=igp+iv+nvirtc(isb)*iatp                          2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+iatp)                  2d13s18
                 end do                                                 2d13s18
                end do                                                  2d13s18
               end do                                                   2d13s18
              end if                                                    2d13s18
              if(i1.le.idoub(isb).and.i2.gt.idoub(isblk1(2,is)))then    2d13s18
               i1m=i1-1                                                 2d13s18
               iad3=itmat(isb)+noc(isb)*i1m                             2d13s18
               do iap=0,iacto(isblk1(3,is))-1                           2d13s18
                iad1=iden(isblk1(3,is))+iap+iacto(isblk1(3,is))*i2n     2d13s18
                fmul=-2d0*bc(iad1)                                      2d13s18
                do iv=0,nvirtc(isb)-1                                    2d13s18
                 iad2=i1x+nrow*(                                        2d13s18
     $                iap+idoub(isblk1(3,is))+noc(isblk1(3,is))*iv)     2d13s18
                 fact=fmul*bc(iad2)                                     2d13s18
                 do idt=0,idoub(isb)-1                                  2d13s18
                  iad4=igp+iv+nvirtc(isb)*idt                           2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+idt)                   2d13s18
                 end do                                                 2d13s18
                 do iat=0,iacto(isb)-1                                  2d13s18
                  iatp=iat+idoub(isb)                                   2d13s18
                  iad4=igp+iv+nvirtc(isb)*iatp                          2d13s18
                  bc(iad4)=bc(iad4)+fact*bc(iad3+iatp)                  2d13s18
                 end do                                                 2d13s18
                end do                                                  2d13s18
               end do                                                   2d13s18
              end if                                                    2d13s18
              i1x=i1x+1                                                 2d13s18
             end do                                                     2d13s18
             i10=1                                                      2d13s18
            end do                                                      2d13s18
           end if                                                       2d13s18
          end if                                                        2d13s18
         end do                                                         2d13s18
        end if                                                          2d13s18
       end if                                                           1d8s18
       if(idoit(5).ne.0)then                                            1d11s18
c
c     aaaa part
c
        if(isb.eq.1)then                                                1d11s18
         ilook=ig+1
        else
         ilook=-1
        end if
        do is=1,nsdlk                                                   1d11s18
         if(isblk(1,is).eq.isblk(2,is))then                             1d11s18
          fmul=-4d0                                                     1d11s18
         else                                                           1d11s18
          fmul=-2d0                                                     1d11s18
         end if                                                         1d11s18
         if(isblk(1,is).eq.isb.and.j2den(is).gt.0)then                  1d11s18
          call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,       1d11s18
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          1d11s18
          i10=i1s                                                       1d11s18
          i1n=noc(isblk(3,is))                                          1d11s18
          if(isblk(1,is).eq.isblk(2,is))then                            1d11s18
           nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2               1d11s18
           nrowd=(iacto(isblk(1,is))*(iacto(isblk(1,is))+1))/2          1d11s18
           iswitch=0                                                    1d11s18
          else                                                          1d11s18
           nrow=noc(isblk(1,is))*noc(isblk(2,is))                       1d11s18
           nrowd=iacto(isblk(1,is))*iacto(isblk(2,is))                  1d11s18
           iswitch=1                                                    1d11s18
          end if                                                        1d11s18
          i4o=ioooo2(is)                                                1d11s18
          do i2=i2s,i2e                                                 1d11s18
           i2m=i2-1-idoub(isblk(4,is))                                  1d11s18
           if(i2.eq.i2e)i1n=i1e                                         1d11s18
           do i1=i10,i1n                                                1d11s18
            if(i1.gt.idoub(isblk(3,is)).and.i2.gt.idoub(isblk(4,is)))
     $           then                                                   1d11s18
             i1m=i1-1-idoub(isblk(3,is))                                1d11s18
             ix=max(i1m,i2m)                                            1d11s18
             in=min(i1m,i2m)                                            1d11s18
             ieq=((ix*(ix+1))/2)+in                                     1d11s18
             inot=i1m+iacto(isblk(3,is))*i2m                            1d11s18
             j2d=j2den(is)+nrowd*((inot-ieq)*iswitch+ieq)               1d11s18
             do ia=0,iacto(isblk(2,is))-1                               1d11s18
              iap=ia+idoub(isblk(2,is))                                 1d11s18
              do id=0,idoub(isb)-1                                      1d11s18
               ieq=((iap*(iap+1))/2)+id                                 1d11s18
               inot=id+noc(isb)*iap                                     1d11s18
               iad1=i4o+(inot-ieq)*iswitch+ieq                          1d11s18
               fact=fmul*bc(iad1)                                       1d26s18
               iad3=ig+iacto(isb)*id                                    1d11s18
               do iat=0,iacto(isb)-1                                    1d11s18
                ix=max(ia,iat)                                          1d11s18
                in=min(ia,iat)                                          1d11s18
                ieq=((ix*(ix+1))/2)+in                                  1d11s18
                inot=iat+iacto(isb)*ia                                  1d11s18
                iad2=j2d+(inot-ieq)*iswitch+ieq                         1d11s18
                bc(iad3+iat)=bc(iad3+iat)+fact*bc(iad2)                 1d11s18
               end do                                                   1d11s18
              end do                                                    1d11s18
             end do                                                     1d11s18
            end if                                                      1d11s18
            i4o=i4o+nrow                                                1d11s18
           end do                                                       1d11s18
           i10=1                                                        1d11s18
          end do                                                        1d11s18
         else if(isblk(2,is).eq.isb.and.j2den(is).gt.0)then             1d11s18
          call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,       1d11s18
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          1d11s18
          i10=i1s                                                       1d11s18
          i1n=noc(isblk(3,is))                                          1d11s18
          nrow=noc(isblk(1,is))*noc(isblk(2,is))                        1d11s18
          nrowd=iacto(isblk(1,is))*iacto(isblk(2,is))                   1d11s18
          i4o=ioooo2(is)                                                1d11s18
          do i2=i2s,i2e                                                 1d11s18
           i2m=i2-1-idoub(isblk(4,is))                                  1d11s18
           if(i2.eq.i2e)i1n=i1e                                         1d11s18
           do i1=i10,i1n                                                1d11s18
            if(i1.gt.idoub(isblk(3,is)).and.i2.gt.idoub(isblk(4,is)))
     $           then                                                   1d11s18
             i1m=i1-1-idoub(isblk(3,is))                                1d11s18
             j2d=j2den(is)+nrowd*(i1m+iacto(isblk(3,is))*i2m)           1d11s18
             do ia=0,iacto(isblk(1,is))-1                               1d11s18
              iap=ia+idoub(isblk(1,is))                                 1d11s18
              do id=0,idoub(isb)-1                                      1d11s18
               iad1=i4o+iap+noc(isblk(1,is))*id                         1d11s18
               fact=fmul*bc(iad1)                                       1d26s18
               iad3=ig+iacto(isb)*id                                    1d11s18
               do iat=0,iacto(isb)-1                                    1d11s18
                iad2=j2d+ia+iacto(isblk(1,is))*iat                      1d16s17
                bc(iad3+iat)=bc(iad3+iat)+fact*bc(iad2)                 1d11s18
               end do                                                   1d11s18
              end do                                                    1d11s18
             end do                                                     1d11s18
            end if                                                      1d11s18
            i4o=i4o+nrow                                                1d11s18
           end do                                                       1d11s18
           i10=1                                                        1d11s18
          end do                                                        1d11s18
         end if                                                         1d11s18
        end do                                                          1d11s18
c
c     occ-virt part                                                     2d15s18
c
        if(c1.ne.0d0)then                                               2d15s18
         ilook=0
         if(nsymb.eq.1)ilook=igp+22-1+66*(3-1)
         if(isb.eq.2)ilook=igp
         do is=1,nsdlk1                                                 2d15s18
          if(isblk1(4,is).eq.isb)then                                   2d15s18
           call ilimts(noc(isblk1(1,is)),noc(isblk1(2,is)),mynprocg,     2d15s18
     $         mynowprog,il,ih,i1s,i1e,i2s,i2e)                         2d15s18
           nrow=ih+1-il                                                 2d15s18
           do isd=1,nsdlk                                               2d15s18
            if(j2den(isd).gt.0)then                                     2d15s18
             if(isblk(1,isd).eq.isblk(2,isd))then                         2d15s18
              nrowd=(iacto(isblk(1,isd))*(iacto(isblk(1,isd))+1))/2        2d15s18
              iswitch=0                                                 2d15s18
              fmul=1d0                                                  2d16s18
             else                                                       2d15s18
              nrowd=iacto(isblk(1,isd))*iacto(isblk(2,isd))               2d15s18
              iswitch=1                                                 2d15s18
              fmul=0.5d0                                                2d16s18
             end if                                                     2d15s18
             fmul=fmul*4d0                                              2d16s18
             if(isblk(3,isd).eq.isblk(4,isd))then                         2d15s18
              jswitch=0                                                 2d15s18
             else                                                       2d15s18
              jswitch=1                                                 2d15s18
             end if                                                     2d15s18
             if(isblk1(1,is).eq.isblk(1,isd).and.                       2d15s18
     $           isblk1(2,is).eq.isblk(2,isd).and.                      2d15s18
     $           isblk1(3,is).eq.isblk(3,isd))then                      2d15s18
              i10=i1s                                                   2d15s18
              i1n=noc(isblk1(1,is))                                     2d15s18
              i1x=iovnn(is)                                             2d15s18
              do i2=i2s,i2e                                             2d15s18
               i2m=i2-1-idoub(isblk1(2,is))                             2d15s18
               if(i2.eq.i2e)i1n=i1e                                     2d15s18
               do i1=i10,i1n                                            2d15s18
                if(i1.gt.idoub(isblk1(1,is)).and.                       2d15s18
     $               i2.gt.idoub(isblk1(2,is)))then                     2d15s18
                 i1m=i1-1-idoub(isblk1(1,is))                           2d15s18
                 ix=max(i1m,i2m)                                        2d15s18
                 in=min(i1m,i2m)                                        2d15s18
                 ieq=((ix*(ix+1))/2)+in                                 2d15s18
                 inot=i1m+iacto(isblk(1,isd))*i2m                       2d15s18
                 i2den=j2den(isd)+(inot-ieq)*iswitch+ieq                2d15s18
                 do ia3=0,iacto(isblk1(3,is))-1                         2d15s18
                  ia3p=ia3+idoub(isblk1(3,is))                          2d15s18
                  do ia4=0,iacto(isb)-1                                 2d15s18
                   ia4p=ia4+idoub(isb)                                  2d15s18
                   ix=max(ia3,ia4)                                      2d15s18
                   in=min(ia3,ia4)                                      2d15s18
                   ieq=((ix*(ix+1))/2)+in                               2d15s18
                   inot=ia3+iacto(isblk(3,isd))*ia4                     2d15s18
                   k2den=i2den+nrowd*((inot-ieq)*jswitch+ieq)           2d15s18
                   fact=bc(k2den)*fmul                                  2d16s18
                   do iv=0,nvirtc(isb)-1                                 2d15s18
                    iad3=i1x+nrow*(ia3p+noc(isblk1(3,is))*iv)           2d16s18
                    fact2=fact*bc(iad3)                                 2d16s18
                    do idt=0,idoub(isb)-1                                2d15s18
                     iad1=igp+iv+nvirtc(isb)*idt                         2d15s18
                     iad2=itmat(isb)+idt+noc(isb)*ia4p                  2d15s18
                     bc(iad1)=bc(iad1)+bc(iad2)*fact2                   2d16s18
  121                format(a2,6i3,3f16.8,2x,4i1,2x,4i1)
                    end do                                              2d15s18
                    do iat=0,iacto(isb)-1                               2d16s18
                     iatp=iat+idoub(isb)                                2d16s18
                     iad1=igp+iv+nvirtc(isb)*iatp                       2d16s18
                     iad2=itmat(isb)+iatp+noc(isb)*ia4p                 2d16s18
                     bc(iad1)=bc(iad1)+bc(iad2)*fact2                   2d16s18
                    end do                                              2d15s18
                   end do                                               2d15s18
                  end do                                                2d15s18
                 end do                                                 2d15s18
                end if                                                  2d15s18
                i1x=i1x+1                                               2d15s18
               end do                                                   2d15s18
               i10=1                                                    2d15s18
              end do                                                    2d15s18
              icasa1=icasa1+1
             else if(isblk1(1,is).eq.isblk(1,isd).and.                       2d15s18
     $           isblk1(2,is).eq.isblk(2,isd).and.                      2d15s18
     $           isblk1(3,is).eq.isblk(4,isd))then                      2d15s18
              i10=i1s                                                   2d15s18
              i1n=noc(isblk1(1,is))                                     2d15s18
              i1x=iovnn(is)                                             2d15s18
              do i2=i2s,i2e                                             2d15s18
               i2m=i2-1-idoub(isblk1(2,is))                             2d15s18
               if(i2.eq.i2e)i1n=i1e                                     2d15s18
               do i1=i10,i1n                                            2d15s18
                if(i1.gt.idoub(isblk1(1,is)).and.                       2d15s18
     $               i2.gt.idoub(isblk1(2,is)))then                     2d15s18
                 i1m=i1-1-idoub(isblk1(1,is))                           2d15s18
                 ix=max(i1m,i2m)                                        2d15s18
                 in=min(i1m,i2m)                                        2d15s18
                 ieq=((ix*(ix+1))/2)+in                                 2d15s18
                 inot=i1m+iacto(isblk(1,isd))*i2m                       2d15s18
                 i2den=j2den(isd)+(inot-ieq)*iswitch+ieq                2d15s18
                 do ia3=0,iacto(isblk1(3,is))-1                         2d15s18
                  ia3p=ia3+idoub(isblk1(3,is))                          2d15s18
                  do ia4=0,iacto(isb)-1                                 2d15s18
                   ia4p=ia4+idoub(isb)                                  2d15s18
                   inot=ia4+iacto(isblk(3,isd))*ia3                     2d16s18
                   k2den=i2den+nrowd*inot                               2d16s18
                   fact=bc(k2den)*fmul                                  2d16s18
                   do iv=0,nvirtc(isb)-1                                 2d15s18
                    iad3=i1x+nrow*(ia3p+noc(isblk1(3,is))*iv)           2d16s18
                    fact2=fact*bc(iad3)                                 2d16s18
                    do idt=0,idoub(isb)-1                                2d15s18
                     iad1=igp+iv+nvirtc(isb)*idt                         2d15s18
                     iad2=itmat(isb)+idt+noc(isb)*ia4p                  2d15s18
                     bc(iad1)=bc(iad1)+bc(iad2)*fact2                   2d16s18
                    end do                                              2d15s18
                    do iat=0,iacto(isb)-1                               2d16s18
                     iatp=iat+idoub(isb)                                2d16s18
                     iad1=igp+iv+nvirtc(isb)*iatp                       2d16s18
                     iad2=itmat(isb)+iatp+noc(isb)*ia4p                 2d16s18
                     bc(iad1)=bc(iad1)+bc(iad2)*fact2                   2d16s18
                    end do                                              2d16s18
                   end do                                               2d15s18
                  end do                                                2d15s18
                 end do                                                 2d15s18
                end if                                                  2d15s18
                i1x=i1x+1                                               2d15s18
               end do                                                   2d15s18
               i10=1                                                    2d15s18
              end do                                                    2d15s18
              icasa4=icasa4+1
             else if(isblk1(1,is).eq.isblk(2,isd).and.                       2d15s18
     $           isblk1(2,is).eq.isblk(1,isd).and.                      2d15s18
     $           isblk1(3,is).eq.isblk(3,isd))then                      2d15s18
              icasa2=icasa2+1
             else if(isblk1(1,is).eq.isblk(4,isd).and.                       2d15s18
     $           isblk1(2,is).eq.isblk(3,isd).and.                      2d15s18
     $           isblk1(3,is).eq.isblk(1,isd))then                      2d15s18
              icasb2=icasb2+1
             else if(isblk1(1,is).eq.isblk(2,isd).and.                       2d15s18
     $           isblk1(2,is).eq.isblk(1,isd).and.                      2d15s18
     $           isblk1(3,is).eq.isblk(4,isd))then                      2d15s18
              icasa3=icasa3+1
             else if(isblk1(1,is).eq.isblk(4,isd).and.                       2d15s18
     $           isblk1(2,is).eq.isblk(3,isd).and.                      2d15s18
     $           isblk1(3,is).eq.isblk(2,isd))then                      2d15s18
              icasb3=icasb3+1
             else if(isblk1(1,is).eq.isblk(3,isd).and.                       2d15s18
     $           isblk1(2,is).eq.isblk(4,isd).and.                      2d15s18
     $           isblk1(3,is).eq.isblk(2,isd))then                      2d15s18
              icasb4=icasb4+1
             end if                                                     2d15s18
            end if
           end do                                                       2d15s18
          end if                                                        2d15s18
         end do                                                         2d15s18
        end if                                                          2d15s18
       end if                                                           1d11s18
       call dws_gsumf(bc(ignew(isb)),ntot)
       rmssz(isb)=0d0                                                   4d2s18
       do i=0,ntot-1                                                    2d16s18
        rmssz(isb)=rmssz(isb)+bc(ignew(isb)+i)**2                       4d2s18
       end do                                                           2d16s18
       if(ntot.ne.0)then                                                4d2s18
        rmssz(isb)=sqrt(rmssz(isb)/dfloat(ntot))                        4d2s18
       end if                                                           4d2s18
       irms=irms+ntot                                                   2d16s18
       if(nad.gt.0.and.lprint)then                                      4d5s18
        write(6,*)('act-doub part: '),loc(bc(ig)),ig
        call prntm2(bc(ig),iacto(isb),idoub(isb),iacto(isb))
        if(nsymb.eq.1)call printa(bc(ig),iacto,idoub,1,0,idoub,0,1,0,
     $       bc(ibcoff))
       end if
       if(nvo.gt.0.and.lprint)then                                      4d5s18
        write(6,*)('virt-occ part: '),loc(bc(igp)),igp
        call prntm2(bc(igp),nvirtc(isb),noc(isb),nvirtc(isb))
        if(nsymb.eq.1)call printa(bc(igp),nvirtc,noc,1,0,noc,0,1,0,
     $       bc(ibcoff))
       end if
       jh02=jh02+nbasdwsc(isb)*nbasdwsc(isb)                            1d8s18
       jh02h=jh02h+nbasdwsc(isb)*noc(isb)                               1d16s17
      end do
      ibcoff=ih02h                                                      1d16s17
      if(lprint)write(6,*)('rmssize of updated gradient: '),            4d5s18
     $     (rmssz(i),i=1,nsymb)                                         4d5s18
      if(lprint)then                                                    2d14s20
       write(6,*)('ending h0: ')                                        2d14s20
       jh0mo=ih02                                                       2d14s20
       do isb=1,nsymb                                                   2d14s20
        if(noc(isb).gt.0)then                                           2d14s20
         write(6,*)('for symmetry block '),isb                          2d14s20
         call prntm2(bc(jh0mo),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         jh0mo=jh0mo+nbasdwsc(isb)**2                                        2d14s20
        end if                                                          2d14s20
       end do                                                           2d14s20
       write(6,*)('ending oooo: ')
       do is=1,nsdlk                                                    2d14s20
        if(isblk(1,is).eq.isblk(2,is))then                              2d14s20
         nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 2d14s20
        else                                                            2d14s20
         nrow=noc(isblk(1,is))*noc(isblk(2,is))                         2d14s20
        end if                                                          2d14s20
        if(isblk(3,is).eq.isblk(4,is))then                              2d14s20
         ncol=(noc(isblk(3,is))*(noc(isblk(3,is))+1))/2                 2d14s20
        else                                                            2d14s20
         ncol=noc(isblk(3,is))*noc(isblk(4,is))                         2d14s20
        end if                                                          2d14s20
        if(min(nrow,ncol).gt.0)then                                     2d14s20
         write(6,*)('integral type '),(isblk(j,is),j=1,4)
         call prntm2(bc(ioooo2(is)),nrow,ncol,nrow)                      2d14s20
        end if                                                          2d14s20
       end do                                                           2d14s20
      end if                                                            2d14s20
      return
      end
