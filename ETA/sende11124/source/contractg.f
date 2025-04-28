c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine contractg(x,ngaus,ibdat,ibstor,isstor,idorel,iflag,
     $   nbasisp,nbasisc,iapair,nlzz,iorbsym,iorbsymz,mrow,icode,bc,ibc)11d9s22
      implicit real*8 (a-h,o-z)
      integer*8 ibstor,isstor                                           4d27s21
      logical ldeb                                                      4d30s18
c
c     take matrix x in primitive basis and contract accourding to
c     basis set contraction coefficients.
c     if iflag=0, we do similarity transformation for all symmetries,
c     and nbasdws is set to the number of contracted functions.
c     we will also update the orbital qns if nlzz ne 0.                 4d19s21
c     if iflag gt 0, we are transforming the orthogonalizing
c     transformation so the bra refers to the uncontracted basis.
c     that is, input rows are contracted and then uncontracted.         4d27s21
c     nbasdws on input will be the number of contracted functions.
c     this is done only for symmetry block iflag.                       4d27s21
c     if iflag=-1, we just do column contraction for all symmetries.    4d27s21
c     the number of rows in this case is given by mrow.                 4d27s21
c
      dimension x(*),ipt(8),ibstor(*),isstor(*),joff(8),                4d25s18
     $     nbasisp(8),nbasisc(8),ntmp(8),iapair(3,*),                   4d27s21
     $     nbasiscu(8),iorbsym(*),iorbsymz(*),ipt1(8),mrow(*)           4d27s21
      include "common.store"
      include "common.hf"
      include "common.print"                                            1d22s20
      common/singcm/iuse,nff
      if(iprtr(13).eq.0)then                                            1d22s20
       ldeb=.false.                                                       4d30s18
      else                                                              1d22s20
       ldeb=.true.                                                       4d30s18
      end if                                                            1d22s20
      if(ldeb)write(6,*)('in contractg with flag '),iflag,              10d11s22
     $     (' code '),icode                                             10d11s22
c
c     icoff is the offset in bstor and sstor for contracted fcns.
c
      xnan=-4d0
      icoff=0                                                           4s23s18
      if(ldeb)then                                                      11d15s22
       write(6,*)('ibstor,isstor for prims ')
       write(6,*)('idorel = '),idorel                                   11d15s22
       write(6,*)('nsymb = '),nsymb
       write(6,*)('nbasisc: '),(nbasisc(isb),isb=1,nsymb)
       write(6,*)('nlorcont: '),(nlorcont(isb),isb=1,nsymb)
      end if                                                            11d15s22
      nbct=0                                                            4d25s18
      do isb=1,nsymb                                                    4s23s18
       ntmp(isb)=nbasisp(isb)                                           4d25s18
       if(idorel.eq.0)then                                              3d3s20
        nbasiscu(isb)=nbasisc(isb)                                       3d3s20
       else                                                             3d3s20
        nbasiscu(isb)=nlorcont(isb)                                     3d29s21
       end if                                                           3d3s20
       nbct=nbct+nbasiscu(isb)                                          3d3s20
       do j=1,nbasisp(isb)
        jp=icoff+j
        if(ldeb)write(6,*)jp,ibstor(jp),isstor(jp),nbasisp(isb),
     $       nbasiscu(isb)                                              3d3s20
       end do
       icoff=icoff+nbasisp(isb)                                         4s23s18
      end do                                                            4s23s18
      if(ldeb)write(6,*)('for contracts '),icoff
      do j=1,nbct                                                       4d25s18
       jm=icoff+j                                                       4d24s18
       if(ldeb)write(6,*)j,jm,ibstor(jm),isstor(jm)
      end do
      ioff=0
      if(iflag.le.0)then                                                4d27s21
       ioff1=0                                                          4d19s21
       do isb=1,nsymb
        joff(isb)=0                                                      9d26s17
        nhere=nbasisp(isb)                                              9d26s17
        if(idorel.ne.0)nhere=nhere*2
        ipt(isb)=ioff
        ipt1(isb)=ioff1                                                 4d19s21
        ioff1=ioff1+nhere                                               4d19s21
        if(iflag.eq.0)then                                              4d27s21
         ioff=ioff+nhere*nhere
        else                                                            4d27s21
         ioff=ioff+nhere*mrow(isb)                                      4d27s21
        end if                                                          4d27s21
       end do
       itmp=ibcoff
       itmp2=itmp+ioff
       ibcoff=itmp2+ioff
       if(nlzz.ne.0.and.iflag.eq.0)then                                 4d27s21
        itmplz=ibcoff                                                   4d19s21
        ibcoff=itmplz+ioff1                                             4d19s21
        if(nlzz.ne.0)then                                               4d19s21
         itmpll=ibcoff                                                  4d19s21
         ibcoff=itmpll+ioff1                                            4d19s21
        end if                                                          4d19s21
       end if                                                           4d19s21
       call enough('contractg.  1',bc,ibc)
       do i=0,ioff-1
        bc(itmp+i)=xnan
       end do
       if(iflag.eq.0)then                                               4d27s21
        npass=2                                                          9d26s17
       else                                                             4d27s21
        npass=1                                                         4d27s21
       end if                                                           4d27s21
      else
       ioff=0                                                           9d26s17
       npass=1                                                          9d26s17
       itmp=ibcoff                                                      4s23s18
       nrow=nbasisc(iflag)                                              1d30s18
       ncol=nbasisp(iflag)                                              1d30s18
       if(idorel.ne.0)then                                              1d30s18
        nrow=nrow*2                                                     1d30s18
        ncol=ncol*2                                                     1d30s18
       end if                                                           1d30s18
       ibcoff=itmp+nrow*ncol                                            1d30s18
       do i=0,nrow*ncol-1                                               1d30s18
        bc(itmp+i)=xnan
       end do
      end if                                                            9d26s17
      call enough('contractg.  2',bc,ibc)
      if(iflag.le.0)then                                                4d27s21
       do isb=1,nsymb                                                    9d26s17
        nhere=nbasisp(isb)                                              9d26s17
        if(idorel.ne.0)nhere=nhere*2                                     9d26s17
        ifrom=ipt(isb)+1                                                 9d26s17
        ito=ipt(isb)+itmp2                                               9d26s17
        if(iflag.eq.0)then                                              4d27s21
         do i=1,nhere                                                     9d26s17
          do j=1,i                                                        9d26s17
           ji=ito+i-1+nhere*(j-1)                                         9d26s17
           ij=ito+j-1+nhere*(i-1)                                         9d26s17
           bc(ji)=x(ifrom)
           bc(ij)=x(ifrom)
           ifrom=ifrom+1
          end do
         end do
        else                                                            4d27s21
         do i=0,nhere*mrow(isb)-1                                       4d27s21
          bc(ito+i)=x(ifrom+i)                                          4d27s21
         end do                                                         4d27s21
        end if                                                          4d27s21
        if(ldeb)then                                                    4d30s18
         write(6,*)('for symmetry block '),isb
         write(6,*)('starting matrix = '),ito
         if(iflag.eq.0)then                                             4d27s21
          call prntm2(bc(ito),nhere,nhere,nhere)
         else                                                           4d27s21
          call prntm2(bc(ito),mrow(isb),nhere,mrow(isb))                4d27s21
         end if                                                         4d27s21
        end if                                                          4d30s18
       end do                                                            9d26s17
      end if
      ng7=ngaus*7                                                       5d2s18
      ng8=ngaus*8
      do ipass=1,npass
       if(ldeb)write(6,*)('for pass no. '),ipass                        4d30s18
       do isb=1,nsymb                                                   3d29s21
        nbasiscu(isb)=nlorcont(isb)                                     3d29s21
       end do                                                           3d29s21
       jbdat=ibdat
       ioff=0
       koff=0                                                            4s23s18
       ignext=0
       do ig=1,ngaus
        l=ibc(jbdat)
        nl=2*l+1
        kbdat=jbdat+ng8
        icx=ibc(kbdat)
        nhere=ibc(jbdat+ng7)                                            5d2s18
        if(ldeb)write(6,*)('shell '),ig,(' has l '),l,(' and icx = '),  4d30s18
     $       icx,(' iapair '),iapair(1,nhere),kbdat                           5d2s18
        if(iapair(1,nhere).eq.0)then                                    5d2s18
         nseq=1                                                         5d2s18
        else                                                            5d2s18
         nseq=2                                                         5d2s18
        end if                                                          5d2s18
        if(icx.gt.0)then
         nz=0
         ngxx=0                                                         5d8s18
         do jg=1,ngaus-ig
          jcx=ibc(kbdat+jg)
          ngxx=max(ngxx,-jcx)                                           5d8s18
          if(ldeb)write(6,*)jg,jcx
          if(jcx.eq.0)nz=nz+1
          if(jcx.ne.0.and.nz.gt.0.or.jcx.gt.0.or.jcx.eq.-1)then         1d22s20
           jjg=jg                                                       1d22s20
           go to 1                                                      1d22s20
          end if                                                        1d22s20
         end do
         if(ldeb)then
         write(6,*)('could not find end of contraction ')
         write(6,*)('ngxx = '),ngxx
         write(6,*)('set it to last+1 ')                                  1d22s20
         end if
         jjg=ngaus-ig+1                                                 1d22s20
    1    continue
         jg=jjg                                                         1d22s20
         if(ldeb)then                                                   4d30s18
          write(6,*)('jg = '),jg
          write(6,*)('nz = '),nz
         end if                                                         4d30s18
         numt=jg                                                        1d22s20
         numl=jg                                                        1d28s19
         ignext=ig+numt                                                 1d28s19
         if(ldeb)write(6,*)('ignext = '),ignext                         4d30s18
         jg=numt                                                        1d28s19
         jgr=jg                                                         2d14s19
         if(idorel.ne.0)jgr=jgr*2                                       2d14s19
         nleft=ngxx                                                     1d22s20
         if(nleft.eq.0)nleft=1                                          10d27s20
         if(ldeb)then                                                   4d30s18
         write(6,*)('we will contract '),jgr,(' gaussians down to '),   2d14s19
     $       nleft,(' fcns ')
         write(6,*)('with the coefficients '),ibdat+icx
         call prntm2(bc(ibdat+icx),jgr,nleft,jgr)                       2d14s19
         end if                                                         4d30s18
c
c     recall ibstor, isstor refer to primitive fcns and
c     ibstor(+icoff, isstor refer to contracted fcns ...
         do iseq=1,nseq                                                 5d2s18
          do j=1,nl                                                       4s23s18
           do i=0,jg-1                                                  5d2s18
            iad=ioff+j+nl*(iseq-1+nseq*i)                               5d2s18
            jad=koff+j+nl*(iseq-1+nseq*i)+icoff                         5d2s18
            if(ldeb)write(6,2)iad,j,i+1,iseq,ibstor(iad),                5d2s18
     $         isstor(iad),ibstor(jad),isstor(jad),jad                      4s23s18
    2       format('fcn ',4i5,' is at ',i5,' in symmetry block ',i1,3i5)  4s23s18
           end do
           if(iflag.le.0)then                                           4d27s21
            if(ldeb)write(6,*)('contract ')                               4d30s18
            if(ipass.eq.1.and.nlzz.ne.0.and.iflag.eq.0)then             4d27s21
             iad=ioff+j+nl*(iseq-1)                                     4d19s21
             isb=isstor(iad)
             icoll=ibstor(iad)
             jad=iorbsym(isb)+icoll-1                                   4d19s21
             lzqn=ibc(jad)                                               4d19s21
             if(nlzz.ne.2)then                                          4d19s21
              jad=iorbsymz(isb)+icoll-1                                 4d19s21
              llqn=ibc(jad)                                             4d19s21
              if(ldeb)write(6,*)('l,lz for this contraction: '),llqn,   4d19s21
     $             lzqn                                                 4d19s21
             else                                                       4d19s21
              if(ldeb)write(6,*)('lz for this contraction: '),lzqn
             end if
            end if                                                      4d19s21
            do i=0,jg-1                                                   11d19s17
             iad=ioff+j+nl*(iseq-1+nseq*i)                              5d2s18
             isb=isstor(iad)
             icoll=ibstor(iad)
             if(iflag.eq.0)then                                         4d27s21
              nhere=ntmp(isb)                                             4d25s18
              if(idorel.ne.0.and.ipass.eq.1)nhere=nhere*2                2d14s19
             else                                                       4d27s21
              nhere=mrow(isb)                                           4d27s21
             end if                                                     4d27s21
             if(i.eq.0)then                                               11d19s17
              icpy=ibcoff                                                   11d19s17
              ibcoff=icpy+jgr*nhere                                     2d14s19
              call enough('contractg.  3',bc,ibc)
              jcpy=icpy                                                   11d19s17
             end if                                                       11d19s17
             locat=itmp2+ipt(isb)+nhere*(icoll-1)                         11d19s17
             do k=0,nhere-1                                               11d19s17
              bc(jcpy+k)=bc(locat+k)                                      11d19s17
             end do                                                       11d19s17
             if(idorel.ne.0)then                                        2d14s19
              locat=itmp2+ipt(isb)+nhere*(icoll+nbasisp(isb)-1)         2d14s19
              kcpy=jcpy+nhere*jg                                        2d14s19
              do k=0,nhere-1                                            2d14s19
               bc(kcpy+k)=bc(locat+k)                                   2d14s19
              end do                                                    2d14s19
             end if                                                     2d14s19
             jcpy=jcpy+nhere                                              11d19s17
            end do                                                        11d19s17
            itmpo=ibcoff                                                  11d19s17
            ibcoff=itmpo+nhere*nleft                                      11d19s17
            call enough('contractg.  4',bc,ibc)
            if(nhere.gt.0.and.jgr.gt.0)then                             2d22s19
             call dgemm('n','n',nhere,nleft,jgr,1d0,bc(icpy),nhere,      2d14s19
     $         bc(ibdat+icx),jgr,0d0,bc(itmpo),nhere,                   2d14s19
     d' contractg.  1')
            end if                                                      2d22s19
            if(ldeb)then
             write(6,*)('to yield ')                                      4d30s18
             call prntm2(bc(itmpo),nhere,nleft,nhere)                      11d19s17
             write(6,*)('koff = '),koff
            end if                                                       4d30s18
            jtmpo=itmpo                                                   11d19s17
            icolx=0                                                      4d24s18
            do i=0,nleft-1                                                11d19s17
             iad=koff+j+nl*(iseq-1+nseq*i)+icoff                        5d2s18
             icoll=ibstor(iad)                                            4s23s18
             if(ldeb)write(6,*)i,iad-icoff,icoll,iad,isstor(iad)            12d13s22
             icolx=max(icolx,icoll)                                      4d24s18
             isb=isstor(iad)                                              4s23s18
             joff(isb)=joff(isb)+1                                        4s23s18
             jtmp=itmp+ipt(isb)+nhere*(icoll-1)                           4s23s18
             if(nlzz.gt.0.and.ipass.eq.1.and.iflag.eq.0)then            4d27s21
              jlz=itmplz+ipt1(isb)+icoll-1                              4d19s21
              ibc(jlz)=lzqn                                             4d19s21
              if(nlzz.ne.2)then                                         4d19s21
               jllz=itmpll+ipt1(isb)+icoll-1                            4d19s21
               ibc(jllz)=llqn                                           4d19s21
              end if
             end if                                                     4d19s21
             do k=0,nhere-1                                               11d19s17
              bc(jtmp+k)=bc(jtmpo+k)                                      11d19s17
             end do                                                       11d19s17
             jtmpo=jtmpo+nhere                                            11d19s17
            end do                                                        11d19s17
            if(ldeb)then                                                 4d30s18
             write(6,*)('tmp so fara '),ipt(isb),itmp,itmp+ipt(isb)
             call prntm2(bc(itmp+ipt(isb)),nhere,icolx,nhere)             4d24s18
            end if                                                       4d30s18
            ibcoff=icpy                                                  4d30s18
           else if(isstor(ioff+j+nl*(iseq-1)).eq.iflag)then             5d2s18
            if(ldeb)write(6,*)('uncontract ')                           4d27s21
            if(icode.eq.0)then                                          10d12s22
             itmpi=ibcoff                                                  4s23s18
             nhere=nbasisc(iflag)                                         4d25s18
             itmpo=itmpi+nleft*nhere                                      4d25s18
             ibcoff=itmpo+jgr*nhere                                      2d14s19
             call enough('contractg.  5',bc,ibc)
             nrow=nbasisc(iflag)                                         1d30s18
             do i=0,nleft-1                                               4d25s18
              iad=koff+j+nl*(iseq-1+nseq*i)+icoff                        5d2s18
              irow=ibstor(iad)                                            4d25s18
              do icol=0,nhere-1                                           4d25s18
               iad1=itmpi+i+nleft*icol                                    4d25s18
               iad2=irow+nrow*icol                                       1d30s18
               bc(iad1)=x(iad2)                                           4d25s18
              end do                                                      4d25s18
             end do                                                       4d25s18
             if(ldeb)then                                                 4d30s18
              write(6,*)('uncontracting '),nrow                          10d11s22
              call prntm2(bc(itmpi),nleft,nhere,nleft)                     4d25s18
             end if                                                       4d30s18
            else                                                        10d12s22
             nrowx=jg                                                   10d12s22
             nnrowx=nbasisp(iflag)                                      10d12s22
             if(idorel.ne.0)then                                        10d12s22
              nrowx=nrowx*2                                             10d12s22
              nnrowx=nnrowx*2                                           10d12s22
             end if                                                     10d12s22
             itmpi=ibcoff                                               10d12s22
             ibcoff=itmpi+nrowx*nbasisc(iflag)                          10d12s22
             call enough('contractg.  6',bc,ibc)
             do i=0,jg-1                                                10d12s22
              iad=ioff+j+nl*(iseq-1+nseq*i)                             10d12s22
              irow=ibstor(iad)                                            4d25s18
              jtmpi=itmpi+i                                             10d12s22
              do k=0,nbasisc(iflag)-1                                   10d12s22
               bc(jtmpi)=x(irow)                                        10d12s22
               jtmpi=jtmpi+nrowx                                        10d12s22
               irow=irow+nnrowx                                         10d12s22
              end do                                                    10d12s22
              if(idorel.ne.0)then                                       10d12s22
               irow=ibstor(iad)+nbasisp(iflag)                          10d12s22
               jtmpi=itmpi+i+jg                                         10d12s22
               do k=0,nbasisc(iflag)-1                                  10d12s22
                bc(jtmpi)=x(irow)                                        10d12s22
                jtmpi=jtmpi+nrowx                                        10d12s22
                irow=irow+nnrowx                                         10d12s22
               end do                                                    10d12s22
              end if                                                    10d12s22
             end do                                                     10d12s22
             if(ldeb)then                                                 4d30s18
              write(6,*)('recontracting ')
              call prntm2(bc(itmpi),nrowx,mrow(iflag),nrowx)            2d2s23
             end if                                                       4d30s18
            end if                                                      10d12s22
            if(jgr.gt.0.and.nleft.gt.0)then                             2d22s19
             if(icode.eq.0)then                                         10d11s22
              call dgemm('n','n',jgr,nhere,nleft,1d0,bc(ibdat+icx),jgr,   2d14s19
     $          bc(itmpi),nleft,0d0,bc(itmpo),jgr,                      2d14s19
     d' contractg.  2')
              if(ldeb)then                                                 4d30s18
               write(6,*)('to yield')
               call prntm2(bc(itmpo),jgr,nhere,jgr)                       2d14s19
              end if                                                       4d30s18
             else                                                       10d11s22
              if(ldeb)then                                              10d12s22
               write(6,*)('coefficients ')
               call prntm2(bc(ibdat+icx),jgr,nleft,jgr)
              end if                                                    10d12s22
              itmpo=ibcoff                                              10d12s22
              iwgt=itmpo+nleft*mrow(iflag)                              2d2s23
              iscr=iwgt+nrowx                                           10d12s22
              ibcoff=iscr+2*nleft+nleft*nleft+2*nleft*jgr+jgr           10d12s22
              call enough('contractg.  7',bc,ibc)
              do i=0,nrowx-1                                            10d12s22
               bc(iwgt+i)=1d0                                           10d12s22
              end do                                                    10d12s22
              iuse=0                                                    10d12s22
              call lsqfit2(bc(ibdat+icx),jgr,nleft,bc(itmpi),nrowx,     10d12s22
     $             mrow(iflag),jgr,bc(itmpo),nleft,bc(iscr),            2d6s23
     $             bc(iwgt),0,rmsovr,bc,ibc)                            3d27s23
              ibcoff=iwgt                                               10d12s22
              if(ldeb)then                                                 4d30s18
               write(6,*)('rmsovr '),rmsovr
               write(6,*)('to yield')
               call prntm2(bc(itmpo),nleft,mrow(iflag),nleft)           2d2s23
              end if                                                       4d30s18
             end if                                                     10d11s22
            end if                                                      2d22s19
            if(icode.ne.0)then                                          10d12s22
             do i=0,nleft-1                                             10d12s22
              iad=koff+j+nl*(iseq-1+nseq*i)+icoff                       10d12s22
              irow=itmp+ibstor(iad)-1                                   10d12s22
              jtmpo=itmpo+i                                             10d12s22
              do k=0,nbasisc(iflag)-1                                   10d12s22
               bc(irow)=bc(jtmpo)                                       10d12s22
               irow=irow+nbasisc(iflag)                                 10d12s22
               jtmpo=jtmpo+nleft                                        10d12s22
              end do                                                    10d12s22
             end do                                                      10d12s22
             if(ldeb)then                                               10d12s22
              write(6,*)('tmp so farr ')                                10d12s22
              call prntm2(bc(itmp),nbasisc(iflag),mrow(iflag),          2d2s23
     $             nbasisc(iflag))                                      10d12s22
             end if                                                     10d12s22
            else                                                        10d12s22
             irowx=0                                                      4d27s18
             ncol=nbasisp(iflag)                                         1d30s18
             if(idorel.ne.0)ncol=ncol*2                                  1d30s18
             do i=0,jg-1                                                 2d14s19
              iad=ioff+j+nl*(iseq-1+nseq*i)                              5d2s18
              irow=ibstor(iad)                                            4d25s18
              irowx=max(irowx,irow)
              do icol=0,nhere-1                                            4d25s18
               iad1=itmpo+i+jgr*icol                                     2d14s19
               iad2=itmp+irow-1+ncol*icol                                1d30s18
               bc(iad2)=bc(iad1)                                          4d25s18
              end do                                                       4d25s18
              if(idorel.ne.0)then                                        2d14s19
               ip=i+jg                                                   2d14s19
               irow=ibstor(iad)+nbasisp(iflag)                           2d14s19
               irowx=max(irowx,irow)
               do icol=0,nhere-1                                            4d25s18
                iad1=itmpo+ip+jgr*icol                                   2d14s19
                iad2=itmp+irow-1+ncol*icol                                1d30s18
                bc(iad2)=bc(iad1)                                          4d25s18
               end do                                                       4d25s18
              end if                                                     2d14s19
             end do                                                       4d25s18
             if(ldeb)then
              write(6,*)('tmp so farb ')
              call prntm2(bc(itmp),irowx,nrow,ncol)                      1d30s18
             end if
            end if                                                      10d12s22
            ibcoff=itmpi                                                  4s23s18
           end if                                                         9d26s17
          end do
         end do                                                         5d2s18
         koff=koff+nl*nleft*nseq                                        5d2s18
        else                                                             9d26s17
         if(ldeb)write(6,*)('ig vs. ignext '),ig,ignext
         if(ig.ge.ignext)then
          if(ldeb)write(6,*)('nl,nseq '),nl,nseq
          do iseq=1,nseq                                                5d2s18
           lu=1                                                         2d14s19
           nlr=nl                                                       2d14s19
           do l=1,nl
            iad=ioff+l+nl*(iseq-1)                                      5d2s18
            jad=koff+lu+nlr*(iseq-1)+icoff                              2d14s19
            lu=lu+1                                                     2d14s19
            if(iflag.le.0)then                                          4d27s21
             icol=ibstor(iad)
             icoln=ibstor(jad)
             if(ldeb)write(6,*)('copy column: '),icol,(' to '),icoln
             isb=isstor(iad)
             if(iflag.eq.0)then                                         4d27s21
              nhere=ntmp(isb)                                             4d25s18
              if(idorel.ne.0.and.ipass.eq.1)nhere=nhere*2                2d14s19
              if(nlzz.ne.0.and.ipass.eq.1)then                           4d19s21
               locat=iorbsym(isb)+icol-1                                 4d19s21
               locatn=itmplz+ipt1(isb)+icoln-1                           4d19s21
               ibc(locatn)=ibc(locat)                                    4d19s21
               if(nlzz.ne.2)then                                         4d19s21
                locat=iorbsymz(isb)+icol-1                               4d19s21
                locatn=itmpll+ipt1(isb)+icoln-1                          4d19s21
                ibc(locatn)=ibc(locat)                                   4d19s21
               end if                                                    4d19s21
              end if                                                     4d19s21
             else                                                       4d27s21
              nhere=mrow(isb)                                           4d27s21
             end if                                                     4d27s21
             locat=itmp2+ipt(isb)+nhere*(icol-1)                                 9d26s17
             locatn=itmp+ipt(isb)+nhere*(icoln-1)                         4s23s18
             do i=0,nhere-1
              bc(locatn+i)=bc(locat+i)
             end do
             if(idorel.ne.0)then                                        1d30s18
              icolr=icol+nbasisp(isb)                                   1d30s18
              nbasiscu(isb)=nbasiscu(isb)+1                             3d29s21
              icoln=nbasiscu(isb)                                       3d29s21
              if(ldeb)write(6,*)('copya column: '),icolr,(' to '),icoln  1d30s18
              if(nlzz.ne.0.and.ipass.eq.1.and.iflag.eq.0)then           4d27s21
               locat=iorbsym(isb)+icolr-1                               4d19s21
               locatn=itmplz+ipt1(isb)+icoln-1                           4d19s21
               ibc(locatn)=ibc(locat)                                    4d19s21
               if(nlzz.ne.2)then                                         4d19s21
                locat=iorbsymz(isb)+icolr-1                             4d19s21
                locatn=itmpll+ipt1(isb)+icoln-1                          4d19s21
                ibc(locatn)=ibc(locat)                                   4d19s21
               end if                                                    4d19s21
              end if                                                    4d19s21
              locat=itmp2+ipt(isb)+nhere*(icolr-1)                      1d30s18
              locatn=itmp+ipt(isb)+nhere*(icoln-1)                      2d14s19
              do i=0,nhere-1                                            1d30s18
               bc(locatn+i)=bc(locat+i)                                 1d30s18
              end do                                                    1d30s18
             end if                                                     1d30s18
            else if(isstor(iad).eq.iflag)then
             irowf=ibstor(iad)                                           4d25s18
             irowt=ibstor(jad)                                            4s23s18
             nhere=nbasisc(iflag)                                        4d25s18
             if(icode.ne.0)then                                         10d12s22
              if(ldeb)write(6,*)('copyb row '),irowf,(' to '),irowt     2d2s23
              irowf0=irowf                                              2d2s23
              jtmp=itmp+irowt-1                                         2d2s23
              nnrowx=nbasisp(iflag)                                      10d12s22
              if(idorel.ne.0)nnrowx=nnrowx*2                            12d15s22
              do i=0,mrow(iflag)-1                                      2d2s23
               bc(jtmp)=x(irowf)                                        2d2s23
               jtmp=jtmp+nbasisc(iflag)                                 10d12s22
               irowf=irowf+nnrowx                                       2d2s23
              end do                                                    10d12s22
              if(idorel.ne.0)then                                       2d2s23
               irowf=irowf0+nbasisp(iflag)                              2d2s23
               jtmp=itmp+nbasiscu(iflag)                                2d2s23
               nbasiscu(iflag)=nbasiscu(iflag)+1                        2d2s23
               if(ldeb)write(6,*)('copyb2 row '),irowf,                 2d2s23
     $             (' to '),nbasiscu(iflag)                             2d2s23
               do i=0,mrow(iflag)-1                                     2d2s23
                bc(jtmp)=x(irowf)                                       2d2s23
                jtmp=jtmp+nbasisc(iflag)                                2d2s23
                irowf=irowf+nnrowx                                      2d2s23
               end do                                                   2d2s23
              end if                                                    2d2s23
             else                                                       10d12s22
              if(ldeb)write(6,*)('copyb row '),irowt,(' to '),irowf
              nrow=nbasisp(iflag)                                        1d30s18
              if(idorel.ne.0)then                                        1d30s18
               nrow=nrow*2                                               1d30s18
              end if                                                     1d30s18
              do i=0,nhere-1                                              4d25s18
               bc(itmp+irowf-1+i*nrow)=                                  1d30s18
     $             x(irowt+i*nhere)                                     1d30s18
              end do                                                      4d25s18
              if(idorel.ne.0)then                                        1d30s18
               irowfr=irowf+nbasisp(iflag)                               1d30s18
               nbasiscu(iflag)=nbasiscu(iflag)+1                         3d29s21
               irowtr=nbasiscu(iflag)                                    3d29s21
               if(ldeb)write(6,*)('copyc row '),irowtr,(' to '),irowfr
               do i=0,nhere-1                                              4d25s18
                bc(itmp+irowfr-1+i*nrow)=                                1d30s18
     $             x(irowtr+i*nhere)                                    1d30s18
               end do                                                      4d25s18
              end if                                                     1d30s18
             end if                                                     10d12s22
            end if                                                       4d25s18
           end do
          end do
          koff=koff+nlr*nseq                                            2d14s19
         end if                                                          9d26s17
        end if
        jbdat=jbdat+1                                                    9d26s17
        ioff=ioff+nl*nseq                                               5d2s18
       end do
       if(iflag.le.0)then                                               4d27s21
        do isb=1,nsymb
         if(ldeb)write(6,*)('for symmetry '),isb                                4d25s18
         iad=itmp+ipt(isb)
         if(iflag.eq.0)then                                             4d27s21
          nrow=ntmp(isb)
         else                                                           4d27s21
          nrow=mrow(isb)                                                4d27s21
         end if                                                         4d27s21
         ncol=nbasisc(isb)                                              1d30s18
         if(idorel.ne.0.and.ipass.eq.1)then                             2d14s19
          nrow=nrow*2                                                   1d30s18
         end if                                                         1d30s18
         if(ldeb)call prntm2(bc(iad),nrow,ncol,nrow)                    1d30s18
         if(iflag.eq.0)then                                             4d27s21
          if(ipass.eq.1.and.nlzz.ne.0)then                               4d19s21
           do i=0,ncol-1                                                 4d19s21
            ibc(iorbsym(isb)+i)=ibc(itmplz+ipt1(isb)+i)                  4d19s21
           end do                                                        4d19s21
           if(nlzz.ne.2)then                                             4d19s21
            do i=0,ncol-1                                                 4d19s21
             ibc(iorbsymz(isb)+i)=ibc(itmpll+ipt1(isb)+i)                4d19s21
            end do                                                        4d19s21
           end if                                                        4d19s21
          end if                                                         4d19s21
          do i=0,ncol-1                                                  1d30s18
           do j=0,nrow-1                                                 1d30s18
            ji=iad+j+nrow*i                                              1d30s18
            ij=itmp2+ipt(isb)+i+ncol*j                                   1d30s18
            bc(ij)=bc(ji)                                                4d25s18
           end do                                                        4d25s18
          end do                                                         4d25s18
          if(ldeb)then
           write(6,*)('transposed ')                                      4d25s18
           call prntm2(bc(itmp2+ipt(isb)),ncol,nrow,ncol)                 1d30s18
          end if
         end if                                                         4d27s21
         ntmp(isb)=nbasisc(isb)                                         4d25s18
        end do
       end if                                                           4d24s18
      end do                                                            4d24s18
      if(iflag.le.0)then                                                4d27s21
       ioff0=1                                                          5d2s18
       do isb=1,nsymb
        ioff=ioff0                                                      5d2s18
        nhere=nbasisc(isb)                                              1d30s18
        if(iflag.eq.0)then                                              4d27s21
         do i=0,nhere-1                                                  1d30s18
          do j=0,i
           ji=itmp2+ipt(isb)+j+nhere*i                                   1d30s18
           x(ioff)=bc(ji)
           ioff=ioff+1
          end do
         end do
        else                                                            4d27s21
         jtmp2=itmp+ipt(isb)                                            4d28s21
         if(ldeb)then                                                   4d28s21
          write(6,*)('copy '),jtmp2                                     4d28s21
          call prntm2(bc(jtmp2),mrow(isb),nhere,mrow(isb))              4d28s21
         end if                                                         4d28s21
         do i=0,nhere*mrow(isb)-1                                       4d27s21
          x(ioff+i)=bc(jtmp2+i)                                         4d27s21
         end do                                                         4d27s21
         if(ldeb)then                                                   4d28s21
          write(6,*)('to '),ioff                                        4d28s21
          call prntm2(x(ioff),mrow(isb),nhere,mrow(isb))                4d28s21
         end if                                                         4d28s21
        end if                                                          4d27s21
        nspac=nbasisp(isb)                                              1d30s18
        if(idorel.ne.0)nspac=nspac*2                                    1d30s18
        if(iflag.eq.0)then                                              4d27s21
         ioff0=ioff0+nspac**2                                            1d30s18
        else                                                            4d27s21
         ioff0=ioff0+nspac*mrow(isb)                                    4d27s21
        end if                                                          4d27s21
       end do
      else                                                              4d25s18
       ncol=nbasisc(iflag)                                              1d30s18
       if(icode.ne.0)then                                               10d12s22
        nrow=nbasisc(iflag)                                             10d12s22
        ncol=mrow(iflag)                                                2d2s23
       else                                                             10d12s22
        nrow=nbasisp(iflag)                                              1d30s18
        if(idorel.ne.0)then                                              1d30s18
         nrow=nrow*2                                                     1d30s18
        end if                                                           1d30s18
       end if                                                           10d12s22
       if(ldeb)then
        write(6,*)('new matrix ')
        call prntm2(bc(itmp),nrow,ncol,nrow)                             1d30s18
       end if
       jtmp=itmp-1                                                      4d25s18
       do i=1,nrow*ncol                                                 1d30s18
        x(i)=bc(jtmp+i)                                                 4d25s18
       end do                                                           4d25s18
      end if                                                            4d24s18
      ibcoff=itmp
      return
      end
