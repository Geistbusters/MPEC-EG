c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine unmakegd(igddiag,ihddiag,nff2,mdon,mdoop,nsymb,nvirt,  8d21s21
     $     isymmrci,multh,ncsf,ncsf2,nroot,ioverwrite,bc,ibc)           11d10s22
      implicit real*8 (a-h,o-z)                                         8d21s21
      integer*8 igddiag(mdoop,*),ihddiag(mdoop,*),i18,i28,i38,i48       8d21s21
      dimension nff2(mdoop,*),nvirt(*),multh(8,8),ncsf(*),ncsf2(*)      8d21s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      include "common.store"                                            8d21s21
      if(ioverwrite.ne.0)then                                           8d21s21
       i18=1                                                            8d21s21
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
         if(min(nff2(ii,isb),nvisv+nvnotv).gt.0)then                    8d24s21
          ncol=nff2(ii,isb)                                             8d21s21
          nrow=nroot*(ncsf2(iarg)*nvisv+ncsf(iarg)*nvnotv)              8d21s21
          call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  8d21s21
          nhere=ih+1-il                                                 8d21s21
          if(nhere.gt.0)then                                            8d21s21
           itmp=ibcoff                                                  8d21s21
           ibcoff=itmp+nrow*nhere                                       8d21s21
           call enough('unmakegd.  1',bc,ibc)
           i28=nrow                                                     8d21s21
           i38=il                                                       8d21s21
           i48=ih                                                       8d21s21
           call ddi_get(bc,ibc,igddiag(ii,isb),i18,i28,i38,i48,bc(itmp))11d15s22
           call ddi_acc(bc,ibc,ihddiag(ii,isb),i18,i28,i38,i48,bc(itmp))11d15s22
           ibcoff=itmp                                                  8d21s21
          end if                                                        8d21s21
         end if                                                          8d21s21
        end do                                                           8d21s21
       end do                                                            8d21s21
       call dws_synca                                                    8d21s21
      end if                                                            8d21s21
      do isb=nsymb,1,-1                                                 8d21s21
       do ii=mdoop,mdon+1,-1                                            8d21s21
        if(igddiag(ii,isb).gt.0)then                                    8d21s21
         call ddi_destroy(igddiag(ii,isb))                              8d21s21
        end if                                                          8d21s21
       end do                                                           8d21s21
      end do                                                            8d21s21
      call dws_synca                                                    8d21s21
      return                                                            8d21s21
      end                                                               8d21s21
