c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prephdd(nhdiag,nsymb,mdon,mdoo,nff2,ncsf,ncsf2,nvirt,  5d28s21
     $     isymmrci,multh,nroot,maxbxd,nlocald,bc,ibc)                  11d10s22
      implicit real*8 (a-h,o-z)                                         11d19s20
      integer*8 nhdiag,irow,icol                                        11d19s20
c
c     set up distributed storage of uncontracted doubles data
c     1: best vectors
c     2: H*best vectors
c
      dimension nhdiag(mdoo+1,nsymb,2),nff2(mdoo+1,nsymb),              6d8s21
     $     ncsf(*),ncsf2(4,*),nvirt(*),multh(8,8)                       5d28s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.store"
      nspace=0                                                          11d28s20
      irow=1                                                            11d19s20
      maxbxd=0                                                          6d7s21
      do isb=1,nsymb                                                    11d19s20
       isbv12=multh(isb,isymmrci)                                         11d19s20
       nvisv=0                                                          6d8s21
       nvnotv=0                                                         6d8s21
       do isbv1=1,nsymb                                                 6d8s21
        isbv2=multh(isbv12,isbv1)                                       6d8s21
        if(isbv2.ge.isbv1)then                                          6d8s21
         if(isbv1.eq.isbv2)then                                         6d8s21
          nvisv=nvisv+nvirt(isbv1)                                      6d8s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         6d8s21
         else                                                           6d8s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 6d8s21
         end if                                                         6d8s21
         nvnotv=nvnotv+nvv                                              6d8s21
        end if                                                          6d8s21
       end do                                                           6d8s21
       do ii=mdon+1,mdoo+1                                              5d23s23
        iarg=ii-mdon                                                    6d8s21
        nhdiag(ii,isb,1)=-1                                             6d8s21
        nhdiag(ii,isb,2)=-1                                             6d8s21
        if(nff2(ii,isb).gt.0)then                                       6d8s21
         icol=nff2(ii,isb)                                              6d8s21
         irow=nroot*(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)             6d8s21
         maxbxd=max(maxbxd,irow*icol)                                   6d8s21
         call ilimts(1,nff2(ii,isb),mynprocg,mynowprog,il,ih,i1s,i1e,   6d9s21
     $        i2s,i2e)                                                  6d9s21
         nlocald=nlocald+irow*(ih+1-il)                                 6d9s21
         call ddi_create(bc,ibc,irow,icol,nhdiag(ii,isb,1))             11d15s22
         call ddi_create(bc,ibc,irow,icol,nhdiag(ii,isb,2))             11d15s22
         call ddi_zero(bc,ibc,nhdiag(ii,isb,1))                         11d15s22
         call ddi_zero(bc,ibc,nhdiag(ii,isb,2))                         11d15s22
        end if                                                          11d19s20
       end do                                                           6d8s21
      end do                                                            11d19s20
      return                                                            11d19s20
      end
