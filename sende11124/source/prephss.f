c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prephss(nhdiag,nsymb,mdon,mdoo,nff1,ncsf,nvirt,        11d19s20
     $     isymmrci,multh,nroot,nlocal,maxbx,bc,ibc)                    11d10s22
      implicit real*8 (a-h,o-z)                                         11d19s20
      integer*8 nhdiag,irow,icol                                        11d19s20
c
c     set up distributed storage of singles data
c     1: best vectors
c     2: H*best vectors
c
      dimension nhdiag(mdoo+1,nsymb,2),nff1(mdoo+1,nsymb),ncsf(*),      11d19s20
     $     nvirt(*),multh(8,8)                                          11d10s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.store"                                            11d10s20
      maxbx=0                                                           1d29s21
      nlocal=0                                                          11d10s20
      do isb=1,nsymb                                                    11d19s20
       isbv=multh(isb,isymmrci)                                         11d19s20
       irow=nvirt(isbv)*nroot                                           11d10s20
       do ii=mdon+1,mdoo+1                                              5d23s23
        iarg=ii-mdon                                                    11d19s20
        ndblock=ncsf(iarg)*nff1(ii,isb)                                 11d19s20
        nhdiag(ii,isb,1)=-1                                             11d26s20
        nhdiag(ii,isb,2)=-1                                             11d26s20
        if(ndblock.gt.0)then                                            11d19s20
         maxbx=max(maxbx,irow*ndblock)                                  1d30s21
         icol=ndblock                                                   11d10s20
         icol4=icol                                                     11d10s20
         call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  11d10s20
         nhere=ih+1-il                                                  11d10s20
         nlocal=nlocal+nhere*nvirt(isbv)*nroot                          11d10s20
         call ddi_create(bc,ibc,irow,icol,nhdiag(ii,isb,1))             11d15s22
         call ddi_create(bc,ibc,irow,icol,nhdiag(ii,isb,2))             11d15s22
         call ddi_zero(bc,ibc,nhdiag(ii,isb,1))                         11d15s22
         call ddi_zero(bc,ibc,nhdiag(ii,isb,2))                         11d15s22
        end if                                                          11d19s20
       end do                                                           11d19s20
      end do                                                            11d19s20
      maxbx=maxbx+mynprocg*4                                            1d30s21
      return                                                            11d19s20
      end
