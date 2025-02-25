c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine denmx(mdon,mdoo,ibasis,iptr,ncsf,nfcn,                 3d3s21
     $     vint,nctf,isymmrci,irel,ism,irefo,nvirt,nsymb,multh,ixw1,    3d3s21
     $     ixw2,nh0av,nroot,ihsdiag,nff0,iff0,nff1,iff1,nsing,ndoub,    3d9s21
     $     maxbx,norb,nff22,nfdat,vdnono,nec,ncsf2,ih0copy,nbasdws,     3d11s21
     $     nbasisp,ncomp,lwrite,sr2,srh,norbci,iorb,idoubo,bc,ibc)      11d10s22
      implicit real*8 (a-h,o-z)
c                                                                       3d3s21
c     compute one particle density matrix in mo basis.                  3d3s21
c
      logical lwrite                                                    2d24s21
      integer*8 ihsdiag(mdoo+1,nsymb,2)                                 3d3s21
      dimension ibasis(3,*),iptr(2,mdoo+1),ncsf(*),vint(nctf,*),        3d3s21
     $     irel(*),ism(*),irefo(*),nvirt(*),multh(8,8),nh0av(*),iden(8),3d3s21
     $     nff1(mdoo+1,*),iff1(*),nff0(*),iff0(*),vdnono(*),ncsf2(4,*), 3d11s21
     $     nbasdws(*),nbasisp(*),iorbno(8),iorb(*),idoubo(*)            11d2s22
      include "common.store"
      if(lwrite)write(6,*)('Hi, my name is denmx'),nroot,ibcoff         2d24s21
      ibcoffo=ibcoff                                                    3d3s21
      if(norbci.ne.0)then                                               1d24s22
       do isb=1,nsymb                                                   1d24s22
        iorbno(isb)=ibcoff                                              1d24s22
        ibcoff=ibcoff+nh0av(isb)*nh0av(isb)                             1d24s22
       end do                                                           1d24s22
      end if                                                            1d24s22
      do isb=1,nsymb                                                    3d3s21
       iden(isb)=ibcoff                                                 3d3s21
       ibcoff=iden(isb)+nh0av(isb)*nh0av(isb)*nroot                     3d3s21
      end do                                                            3d3s21
      call enough('denmx.  1',bc,ibc)
      do i=ibcoffo,ibcoff-1                                             3d3s21
       bc(i)=0d0                                                        3d3s21
      end do                                                            3d3s21
      nwds=ibcoff-ibcoffo                                               3d3s21
      call hccsfd1(vint,nctf,ibasis,ncsf,nfcn,iden,nroot,mdon,1,        3d3s21
     $     iptr,ixw1,ixw2,mdoo,nh0av,norb,ism,irel,bc,ibc,1)              11d10s22
      if(nsing.ne.0)then                                                3d4s21
       call hcisd1(ihsdiag,nff1,iff1,nff0,iff0,nctf,ncsf,mdon,mdoo,      3d3s21
     $      nsymb,multh,ixw1,ixw2,iden,nh0av,nvirt,vint,bc,ibc,1)       1d27s23
       call hcssd1(ihsdiag,nff1,iff1,ncsf,mdon,mdoo,nsymb,multh,        3d4s21
     $      ixw1,ixw2,iden,nh0av,nroot,ism,irel,irefo,nvirt,isymmrci,   3d4s21
     $      norb,maxbx,bc,ibc,1)                                        1d30s23
      end if                                                            3d4s21
      if(ndoub.ne.0)then                                                3d9s21
       call hcddjkd1(nff22,nfdat,vdnono,nsymb,mdon,mdoo,nec,multh,      3d9s21
     $     isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,ixw1,ixw2,norb,     3d9s21
     $     nroot,iden,nh0av,sr2,srh,ibc,bc,ibc,bc,1)                    1d27s23
       if(nsing.gt.0)then                                               3d11s21
        call hcdsd1(ihsdiag,nff1,iff1,nff22,nfdat,vdnono,nsymb,mdon,    3d11s21
     $      mdoo,nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,    3d11s21
     $      ixw1,ixw2,norb,nroot,iden,nh0av,maxbx,sr2,ibc,bc,ibc,bc,1)  1d27s23
       end if                                                           3d11s21
      end if                                                            3d9s21
      call dws_gsumf(bc(ibcoffo),nwds)                                  3d3s21
      itrace=ibcoff                                                     4d21s21
      ibcoff=itrace+nroot                                               4d21s21
      call enough('denmx.  2',bc,ibc)
      do i=itrace,ibcoff-1                                              4d21s21
       bc(i)=0d0                                                        4d21s21
      end do                                                            4d21s21
      jh0copy=ih0copy                                                   3d11s21
      hcore=0d0                                                         3d11s21
      ikin=ibcoff                                                       3d11s21
      ibcoff=ikin+nroot                                                 3d11s21
      call enough('denmx.  3',bc,ibc)
      do i=ikin,ibcoff-1                                                3d11s21
       bc(i)=0d0                                                        3d11s21
      end do                                                            3d11s21
      do isb=1,nsymb                                                    3d3s21
       if(nh0av(isb).gt.0)then                                          3d3s21
        if(lwrite)write(6,*)('for symmetry '),isb                       2d24s21
        ncor=nbasdws(isb)-nh0av(isb)                                    3d11s21
        do i=0,ncor-1                                                   3d11s21
         ii=jh0copy+i*(nbasdws(isb)+1)                                  3d11s21
         hcore=hcore+2d0*bc(ii)                                         3d11s21
        end do                                                          3d11s21
        nn=nh0av(isb)**2                                                3d3s21
        jden=iden(isb)                                                  3d3s21
        do ir=1,nroot                                                   3d3s21
         if(lwrite)write(6,*)('for root '),ir                           2d24s21
         jtrace=itrace+ir-1                                             4d21s21
         do i=0,nh0av(isb)-1                                            3d3s21
          ii=jden+i*(nh0av(isb)+1)                                      3d3s21
          bc(jtrace)=bc(jtrace)+bc(ii)                                  4d21s21
         end do                                                         3d3s21
         if(lwrite)write(6,*)('trace so far'),bc(jtrace)                2d24s21
         do i=0,nh0av(isb)-1                                            3d11s21
          do j=0,i-1                                                    3d11s21
           ji=jden+i+nh0av(isb)*j                                       3d11s21
           ij=jden+j+nh0av(isb)*i                                       3d11s21
           bc(ij)=bc(ji)                                                3d11s21
          end do                                                        3d11s21
         end do                                                         3d11s21
         jkin=ikin+ir-1                                                 3d11s21
         do i=0,nh0av(isb)-1                                            3d11s21
          ip=i+ncor                                                     3d11s21
          iad=jh0copy+ncor+nbasdws(isb)*ip
          jad=jden+nh0av(isb)*i                                         3d11s21
          do j=0,nh0av(isb)-1                                           3d11s21
           bc(jkin)=bc(jkin)+bc(iad+j)*bc(jad+j)                        3d11s21
          end do                                                        3d11s21
         end do                                                         3d11s21
         if(lwrite)then                                                 11d2s22
          idcopy=ibcoff                                                  3d12s21
          ieig=idcopy+nh0av(isb)*nh0av(isb)                              3d12s21
          ivec=ieig+nh0av(isb)                                           3d11s21
          isym=ivec+nh0av(isb)*nh0av(isb)                                3d11s21
          ibcoff=isym+nh0av(isb)                                         3d11s21
          call enough('denmx.  4',bc,ibc)
          do i=0,nh0av(isb)*nh0av(isb)-1                                 3d12s21
           bc(idcopy+i)=-bc(jden+i)                                       3d12s21
          end do                                                         3d12s21
          call diagx(nh0av(isb),bc(idcopy),bc(ieig),bc(ivec),ibc(isym), 11d14s22
     $         bc,ibc)                                                  11d14s22
          if(ir.eq.norbci)then                                           1d24s22
           do i=0,nh0av(isb)*nh0av(isb)-1                                1d24s22
            bc(iorbno(isb)+i)=bc(ivec+i)                                 1d24s22
           end do                                                        1d24s22
           write(6,*)('no ')
           call prntm2(bc(iorbno(isb)),nh0av(isb),nh0av(isb),nh0av(isb)) 10d26s22
           write(6,*)('iorb ')
           nbasdwsf=nbasisp(isb)*ncomp                                      3d11s21
           nrow=ncomp*(nh0av(isb)+idoubo(isb))                           10d26s22
           call prntm2(bc(iorb(isb)),nbasisp(isb),nrow,nbasisp(isb))     10d26s22
           itmpo=ibcoff                                                  10d26s22
           ibcoff=itmpo+nbasdwsf*nh0av(isb)                              10d26s22
           jorb=iorb(isb)+nbasdwsf*idoubo(isb)                           10d26s22
           call enough('denmx.4b',bc,ibc)
           call dgemm('n','n',nbasdwsf,nh0av(isb),nh0av(isb),1d0,        10d26s22
     $         bc(jorb),nbasdwsf,bc(iorbno(isb)),nh0av(isb),            10d26s22
     $         0d0,bc(itmpo),nbasdwsf)                                  10d26s22
           write(6,*)('no in ao basis ')                                 10d26s22
           call prntm2(bc(itmpo),nbasdwsf,nh0av(isb),nbasdwsf)           10d26s22
           ibcoff=itmpo                                                  10d26s22
          end if                                                        11d2s22
          write(6,*)('natural orbital occupations: ')                   2d24s21
          do i=0,nh0av(isb)-1                                            3d12s21
           bc(ieig+i)=-bc(ieig+i)                                        3d12s21
          end do                                                         3d12s21
          call prntm2(bc(ieig),1,nh0av(isb),1)                          2d24s21
          ibcoff=idcopy                                                  3d12s21
         end if                                                         1d24s22
         jden=jden+nn                                                   3d3s21
        end do                                                          3d3s21
       end if                                                           3d3s21
       nbasdwsf=nbasisp(isb)*ncomp                                      3d11s21
       jh0copy=jh0copy+nbasdwsf*nbasdwsf                                3d11s21
      end do                                                            3d3s21
      if(norbci.lt.0.and.lwrite)then                                    1d24s22
       write(6,*)('averaged over roots ')                               1d24s22
       do isb=1,nsymb                                                   1d24s22
        nn=nh0av(isb)*nh0av(isb)                                        1d24s22
        idcopy=ibcoff                                                   1d24s22
        ieig=idcopy+nn                                                  1d24s22
        isym=ieig+nh0av(isb)                                            1d24s22
        ibcoff=isym+nh0av(isb)                                          1d24s22
        call enough('denmx.  5',bc,ibc)
        nnm=nn-1                                                        1d24s22
        do i=0,nnm                                                      1d24s22
         bc(idcopy+i)=bc(iden(isb)+i)                                   1d24s22
        end do                                                          1d24s22
        jden=iden(isb)+nn                                               1d24s22
        do ir=2,nroot                                                   1d24s22
         do i=0,nnm                                                     1d24s22
          bc(idcopy+i)=bc(idcopy+i)+bc(jden+i)                          1d24s22
         end do                                                         1d24s22
         jden=jden+nn                                                   1d24s22
        end do                                                          1d24s22
        fact=-1d0/dfloat(nroot)                                         1d24s22
        do i=0,nnm                                                      1d24s22
         bc(idcopy+i)=bc(idcopy+i)*fact                                 1d24s22
        end do                                                          1d24s22
        call diagx(nh0av(isb),bc(idcopy),bc(ieig),bc(iorbno(isb)),      1d24s22
     $       ibc(isym),bc,ibc)                                          11d14s22
        write(6,*)('natural orbital occupations for symmetry  '),isb    1d24s22
        do i=0,nh0av(isb)-1                                             1d24s22
         bc(ieig+i)=-bc(ieig+i)                                         1d24s22
        end do                                                          1d24s22
        call prntm2(bc(ieig),1,nh0av(isb),1)                            1d24s22
        write(6,*)('no ')
        call prntm2(bc(iorbno(isb)),nh0av(isb),nh0av(isb),nh0av(isb))   10d26s22
        write(6,*)('iorb ')
        nbasdwsf=nbasisp(isb)*ncomp                                      3d11s21
        nrow=ncomp*(nh0av(isb)+idoubo(isb))                             10d26s22
        call prntm2(bc(iorb(isb)),nbasisp(isb),nrow,nbasisp(isb))       10d26s22
        itmpo=ibcoff                                                    10d26s22
        ibcoff=itmpo+nbasdwsf*nh0av(isb)                                10d26s22
        call enough('denmx.5b',bc,ibc)
        jorb=iorb(isb)+nbasdwsf*idoubo(isb)                             10d26s22
        call dgemm('n','n',nbasdwsf,nh0av(isb),nh0av(isb),1d0,          10d26s22
     $         bc(jorb),nbasdwsf,bc(iorbno(isb)),nh0av(isb),            10d26s22
     $         0d0,bc(itmpo),nbasdwsf)                                  10d26s22
        write(6,*)('no in ao basis ')                                   10d26s22
        call prntm2(bc(itmpo),nbasdwsf,nh0av(isb),nbasdwsf)             10d26s22
        ibcoff=itmpo                                                    10d26s22
        ibcoff=idcopy                                                   1d24s22
       end do                                                           1d24s22
      end if                                                            1d24s22
      if(lwrite)then                                                    11d9s22
       do ir=1,nroot                                                     3d11s21
        jkin=ikin+ir-1                                                   3d11s21
        write(6,*)('kinetic energy for root '),ir,hcore,                2d24s21
     $      bc(jkin),hcore+bc(jkin)                                     2d24s21
       end do                                                            3d11s21
      end if                                                            11d9s22
      ibcoff=ibcoffo                                                    3d3s21
      return
      end
