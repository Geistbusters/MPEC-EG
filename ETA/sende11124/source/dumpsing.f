c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumpsing(iddi,nff1,iff1,ncsf,mdon,mdoo,nsymb,multh,    7d8s21
     $     nvirt,isymmrci,nroot,vstry,iflag,bc,ibc)                     11d14s22
      implicit real*8 (a-h,o-z)                                         7d8s21
      integer*8 iddi(mdoo+1,nsymb),i18,i28,i38                          7d8s21
      dimension nff1(mdoo+1,nsymb,2),iff1(*),ncsf(*),multh(8,8),        7d8s21
     $     nvirt(*),vstry(*)                                            7d9s21
      include "common.store"                                            7d8s21
      if(iflag.eq.0)then                                                7d27s21
       mff1=0                                                            7d8s21
       do isb=1,nsymb                                                    7d8s21
        do ii=mdon+1,mdoo+1                                              7d8s21
         if(nff1(ii,isb,1).gt.0)mff1=mff1+nff1(ii,isb,1)                 7d8s21
        end do                                                           7d8s21
       end do                                                            7d8s21
       mff1=mff1*2                                                       7d8s21
       write(2)mff1                                                      7d8s21
       write(2)(iff1(i),i=1,mff1)                                        7d8s21
       write(2)(((nff1(ii,isb,j),ii=mdon+1,mdoo+1),isb=1,nsymb),j=1,2)   7d8s21
       write(2)(ncsf(iarg),iarg=1,mdoo+1-mdon)                           7d8s21
      else                                                              7d27s21
       i18=1                                                             7d8s21
       do isb=1,nsymb                                                    7d8s21
        isbv=multh(isb,isymmrci)                                         7d8s21
        do ii=mdon+1,mdoo+1                                              7d8s21
         if(nff1(ii,isb,1).gt.0)then                                     7d8s21
          iarg=ii-mdon                                                   7d8s21
          nrow=nroot*nvirt(isbv)                                         7d8s21
          ncol=ncsf(iarg)*nff1(ii,isb,1)                                 7d8s21
          itmp=ibcoff                                                    7d8s21
          ibcoff=itmp+nrow*ncol                                          7d8s21
          call enough('dumpsing.  1',bc,ibc)
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
