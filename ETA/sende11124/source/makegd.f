c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine makegd(igddiag,nff2,mdon,mdoop,nsymb,nvirt,isymmrci,   8d21s21
     $     multh,ncsf,ncsf2,nroot,maxbxd,npadddi,bc,ibc)                11d15s22
      implicit real*8 (a-h,o-z)                                         8d21s21
      integer*8 igddiag(mdoop,*),i18,i28                                8d21s21
      dimension nff2(mdoop,*),nvirt(*),multh(8,8),ncsf(*),ncsf2(*)      8d21s21
      do isb=1,nsymb                                                    8d21s21
       nvisv=0                                                          8d21s21
       nvnotv=0                                                         8d21s21
       isbv12=multh(isb,isymmrci)                                       8d21s21
       do isbv1=1,nsymb                                                 8d21s21
        isbv2=multh(isbv1,isbv12)                                       8d21s21
        if(isbv2.ge.isbv1)then                                          8d21s21
         if(isbv1.eq.isbv2)then                                         8d21s21
          nvisv=nvisv+nvirt(isbv1)                                      8d21s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d21s21
         else                                                           8d21s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d21s21
         end if                                                         8d21s21
         nvnotv=nvnotv+nvv                                              8d21s21
        end if                                                          8d21s21
       end do                                                           8d21s21
       do ii=mdon+1,mdoop                                                8d21s21
        iarg=ii-mdon                                                     8d21s21
        igddiag(ii,isb)=-1                                              8d21s21
        if(min(nvisv+nvnotv,nff2(ii,isb)).gt.0)then                     8d24s21
         i28=nff2(ii,isb)                                               8d21s21
         i18=nroot*(ncsf2(iarg)*nvisv+ncsf(iarg)*nvnotv)                8d21s21
         nsx=i18*i28                                                    11d15s21
         maxbxd=max(maxbxd,nsx+npadddi)                                 11d15s21
         call ddi_create(bc,ibc,i18,i28,igddiag(ii,isb))                11d15s22
        end if                                                          8d21s21
       end do                                                           8d21s21
      end do                                                            8d21s21
      call dws_synca                                                    8d21s21
      return                                                            8d21s21
      end                                                               8d21s21
