c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumpdoubuc(iddi,nff2,iff2,ncsf,ncsf2,mdon,mdoo,nsymb,
     $     multh,nvirt,isymmrci,nroot,iflag,idcsf2,bc,ibc)              11d14s22
      implicit real*8 (a-h,o-z)                                         7d8s21
      integer*8 iddi(mdoo+1,nsymb),i18,i28,i38                          7d8s21
      dimension nff2(mdoo+1,nsymb,2),iff2(*),ncsf(*),multh(8,8),        7d8s21
     $     nvirt(*),ncsf2(idcsf2,*)                                     2d14s22
      include "common.store"                                            7d8s21
      if(iflag.eq.0)then                                                7d27s21
       mff2=0                                                            7d8s21
       do isb=1,nsymb                                                    7d8s21
        do ii=mdon+1,mdoo+1                                              7d8s21
         if(nff2(ii,isb,1).gt.0)mff2=mff2+nff2(ii,isb,1)                 7d8s21
        end do                                                           7d8s21
       end do                                                            7d8s21
       mff2=mff2*2                                                       7d8s21
       write(2)mff2                                                      7d8s21
       write(2)(iff2(i),i=1,mff2)                                        7d8s21
       write(2)(((nff2(ii,isb,j),ii=mdon+1,mdoo+1),isb=1,nsymb),j=1,2)   7d8s21
       write(2)(ncsf(iarg),iarg=1,mdoo+1-mdon)                           7d8s21
       write(2)(ncsf2(1,iarg),iarg=1,mdoo+1-mdon)                        7d21s21
      else                                                              7d27s21
       i18=1                                                             7d8s21
       do isb=1,nsymb                                                    7d8s21
        isbv12=multh(isb,isymmrci)                                       7d21s21
        nvisv=0                                                          7d21s21
        nvnotv=0                                                         7d21s21
        do isbv1=1,nsymb                                                 7d21s21
         isbv2=multh(isbv1,isbv12)                                       7d21s21
         if(isbv2.ge.isbv1)then                                          7d21s21
          if(isbv1.eq.isbv2)then                                         7d21s21
           nvisv=nvisv+nvirt(isbv1)                                      7d21s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         7d21s21
          else                                                           7d21s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                 7d21s21
          end if                                                         7d21s21
          nvnotv=nvnotv+nvv                                              7d21s21
         end if                                                          7d21s21
        end do                                                           7d21s21
        do ii=mdon+1,mdoo+1                                              7d8s21
         if(nff2(ii,isb,1).gt.0)then                                     7d8s21
          iarg=ii-mdon                                                   7d8s21
          nrow=nroot*(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)             7d21s21
          ncol=nff2(ii,isb,1)                                            7d21s21
          itmp=ibcoff                                                    7d8s21
          ibcoff=itmp+nrow*ncol                                          7d8s21
          call enough('dumpdoubuc.  1',bc,ibc)
          i28=nrow                                                       7d8s21
          i38=ncol                                                       7d8s21
          call ddi_get(bc,ibc,iddi(ii,isb),i18,i28,i18,i38,bc(itmp))    11d15s22
          nam=nrow*ncol-1                                                7d8s21
          write(2)(bc(itmp+i),i=0,nam)                                   7d8s21
          ibcoff=itmp                                                    7d8s21
         end if                                                          7d8s21
        end do                                                           7d8s21
       end do                                                            7d8s21
      end if                                                            7d27s21
      return                                                            7d8s21
      end                                                               7d8s21
