c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine makegs(igsdiag,nff1,mdon,mdoop,nsymb,nvirt,isymmrci,   8d21s21
     $     multh,ncsf,nroot,maxbx,npadddi,bc,ibc)                       11d10s22
      implicit real*8 (a-h,o-z)                                         8d21s21
      integer*8 igsdiag(mdoop,*),i18,i28                                8d21s21
      dimension nff1(mdoop,*),nvirt(*),multh(8,8),ncsf(*)               8d21s21
      include "common.store"
      do ii=mdon+1,mdoop                                                8d21s21
       iarg=ii-mdon                                                     8d21s21
       do isb=1,nsymb                                                   8d21s21
        igsdiag(ii,isb)=-1                                              8d21s21
        if(nff1(ii,isb).gt.0)then                                       8d21s21
         isbv=multh(isb,isymmrci)                                       8d21s21
         i28=nff1(ii,isb)*ncsf(iarg)                                    8d21s21
         i18=nvirt(isbv)*nroot                                          8d21s21
         nsz=i28*i18                                                    11d12s21
         maxbx=max(maxbx,nsz+npadddi)                                   11d15s21
         call ddi_create(bc,ibc,i18,i28,igsdiag(ii,isb))                11d15s22
        end if                                                          8d21s21
       end do                                                           8d21s21
      end do                                                            8d21s21
      call dws_synca                                                    8d21s21
      return                                                            8d21s21
      end                                                               8d21s21
