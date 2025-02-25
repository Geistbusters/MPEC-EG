c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine unmakegs(igsdiag,ihsdiag,nff1,mdon,mdoop,nsymb,nvirt,  8d21s21
     $     isymmrci,multh,ncsf,nroot,ioverwrite,bc,ibc)                 11d10s22
      implicit real*8 (a-h,o-z)                                         8d21s21
      integer*8 igsdiag(mdoop,*),ihsdiag(mdoop,*),i18,i28,i38,i48       8d21s21
      dimension nff1(mdoop,*),nvirt(*),multh(8,8),ncsf(*)               8d21s21
      include "common.store"                                            8d21s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      if(ioverwrite.ne.0)then                                           8d21s21
       i18=1                                                             8d21s21
       do ii=mdon+1,mdoop                                                8d21s21
        iarg=ii-mdon                                                     8d21s21
        do isb=1,nsymb                                                   8d21s21
         if(nff1(ii,isb).gt.0)then                                       8d21s21
          isbv=multh(isb,isymmrci)                                       8d21s21
          ncol=nff1(ii,isb)*ncsf(iarg)                                    8d21s21
          nrow=nvirt(isbv)*nroot                                          8d21s21
          call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  8d21s21
          nhere=ih+1-il                                                 8d21s21
          if(nhere.gt.0)then                                            8d21s21
           itmp=ibcoff                                                    8d21s21
           ibcoff=itmp+nhere*nrow                                       8d21s21
           i28=nrow                                                     8d21s21
           call enough('unmakegs.  1',bc,ibc)
           i38=il                                                       8d21s21
           i48=ih                                                       8d21s21
           call ddi_get(bc,ibc,igsdiag(ii,isb),i18,i28,i38,i48,bc(itmp))11d15s22
           call ddi_acc(bc,ibc,ihsdiag(ii,isb),i18,i28,i38,i48,bc(itmp))11d15s22
          end if                                                        8d21s21
         end if                                                          8d21s21
        end do                                                           8d21s21
       end do                                                            8d21s21
       call dws_synca                                                    8d21s21
      end if                                                            8d21s21
      do ii=mdoop,mdon+1,-1                                             8d21s21
       do isb=nsymb,1,-1                                                8d21s21
        if(igsdiag(ii,isb).gt.0)then                                    8d21s21
         call ddi_destroy(igsdiag(ii,isb))                              8d21s21
        end if                                                          8d21s21
       end do                                                           8d21s21
      end do                                                            8d21s21
      call dws_synca                                                    8d21s21
      return                                                            8d21s21
      end                                                               8d21s21
