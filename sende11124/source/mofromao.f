c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mofromao(ivecs,ieigs,iorb,morb,idorel,ixinv,nbasisp,   2d15s19
     $     nbasisc,ismile,bc,ibc,nposs)                                 5d2s23
      implicit real*8 (a-h,o-z)
      include "common.hf"
      include "common.store"
      include "common.spher"
      include "common.print"                                            4d28s21
      logical lwrite                                                    1s1s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ivecs(8),ieigs(8),iorb(8),morb(8),ixinv(8),nbasisp(8),  2d15s19
     $     nbasisc(8),ismile(8),nposs(*)                                5d2s23
      ncomp=1
      if(iprtr(22).ne.0)then                                            4d28s21
       lwrite=.true.                                                     1s1s20
      else                                                              4d28s21
       lwrite=.false.                                                    1s1s20
      end if                                                            4d28s21
      if(idorel.ne.0)ncomp=2
      if(lwrite)write(6,*)('in mofromao ')
      do isb=1,nsymb
       if(lwrite)write(6,*)('for symmetry block '),isb
       nrow=nbasisp(isb)*ncomp                                          2d1s19
       nh=nbasdws(isb)                                                  2d1s19
       jvecs=ivecs(isb)+nrow*nbasisc(isb)                               2d15s19
       if(min(nrow,nh).gt.0)then                                        5d10s19
       if(lwrite)then                                                   1s1s20
        ihalf=ibcoff                                                    10d27s20
        ibcoff=ihalf+nrow*nh                                            10d27s20
        call enough('mofromao.  1',bc,ibc)
        write(6,*)('jvecs ')
        call prntm2(bc(jvecs),nrow,nrow,nrow)
        write(6,*)('ivecs ')
        call prntm2(bc(ivecs(isb)),nrow,nh,nrow)
        call dgemm('n','n',nrow,nh,nrow,1d0,bc(jvecs),nrow,             10d27s20
     $       bc(ivecs(isb)),nrow,0d0,bc(ihalf),nrow,                    10d27s20
     d' mofromao.  1')
        call dgemm('t','n',nh,nh,nrow,1d0,bc(ivecs(isb)),nrow,
     $       bc(ihalf),nrow,0d0,bc(ibcoff),nh,
     d' mofromao.  2')
       end if                                                           1s1s20
       morb(isb)=ibcoff
       ibcoff=morb(isb)+nh*nbasisc(isb)                                 2d9s20
       ibcoffo=ibcoff
       call enough('mofromao.  2',bc,ibc)
       itmp1=ibcoff
       ibcoff=itmp1+nrow*nbasisc(isb)                                   2d9s20
       call enough('mofromao.  3',bc,ibc)
        call dgemm('n','n',nrow,nh,nrow,1d0,bc(jvecs),nrow,              2d15s19
     $      bc(ivecs(isb)),nrow,0d0,bc(itmp1),nrow,                     2d15s19
     d' mofromao.  3')
       itmp2=ibcoff                                                     2d15s19
       ibcoff=itmp2+nrow*nbasisc(isb)                                   2d9s20
       call enough('mofromao.  4',bc,ibc)
       do i=0,nh-1                                                      2d15s19
        do j=0,nrow-1                                                   2d15s19
         ji=itmp1+j+nrow*i                                              2d15s19
         ij=itmp2+i+nh*j                                                2d15s19
         bc(ij)=bc(ji)                                                  2d15s19
        end do                                                          2d15s19
       end do                                                           2d15s19
        if(lwrite)then                                                  4d28s21
         write(6,*)('iorb ')
         call prntm2(bc(iorb(isb)),nrow,nh,nrow)
         ixi=ibcoff
         ibcoff=ixi+nrow*nh
         call dgemm('n','n',nrow,nh,nrow,1d0,bc(jvecs),nrow,
     $        bc(iorb(isb)),nrow,0d0,bc(ixi),nrow)
         write(6,*)('S*o')
         call prntm2(bc(ixi),nrow,nh,nrow)
         call dgemm('t','n',nh,nh,nrow,1d0,bc(iorb(isb)),nrow,
     $        bc(ixi),nrow,0d0,bc(ibcoff),nh)
         write(6,*)('ortho test? ')
         call prntm2(bc(ibcoff),nh,nh,nh)
         write(6,*)('trans from ao to mo: ')
         call prntm2(bc(itmp2),nh,nrow,nh)
        end if                                                          4d28s21
        call dgemm('n','n',nh,nh,nrow,1d0,bc(itmp2),nh,bc(iorb(isb)),    2d15s19
     $      nrow,0d0,bc(morb(isb)),nh,                                  2d15s19
     d' mofromao.  4')
        if(idorel.ne.0)then                                             2d1s23
         if(nposs(isb).ne.0)then                                        5d2s23
          igetfrom=ivecs(isb)                                           2d1s23
     $         +nbasisp(isb)*ncomp*(nbasisc(isb)*2+nbasisp(isb)*ncomp)  2d1s23
          itmp=ibcoff                                                   2d1s23
          ibcoff=itmp+nh*nh                                             2d1s23
          call enough('mofromao.tmp',bc,ibc)                            2d1s23
          do i=0,nh*nh-1                                                2d1s23
           bc(itmp+i)=bc(morb(isb)+i)                                   2d1s23
          end do                                                        2d1s23
          call dgemm('n','n',nbasisc(isb),nh,nh,1d0,bc(igetfrom),       2d1s23
     $         nbasisc(isb),bc(itmp),nh,0d0,bc(morb(isb)),nbasisc(isb)) 2d1s23
          ibcoff=itmp                                                   2d1s23
         end if                                                         2d1s23
        end if                                                          2d1s23
       end if                                                           5d10s19
       ibcoff=ibcoffo
      end do
      return
      end
