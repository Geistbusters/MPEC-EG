c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine countem2(nff2,nsymb,ncdx,mdon,mdoo,ncsf2,              11d19s20
     $     multh,nvirt,nvxv,bc,ibc)                                     11d14s22
      implicit real*8 (a-h,o-z)
      logical ldebug                                                    12d12s19
      dimension nff2(mdoo+1,*),ncdx(3,8,2),ncsf2(4,*)                      11d19s20
     $     ,multh(8,8),nvirt(8),nvxv(2,8),nkk(4)                        11d19s20
      include "common.mrci"                                             8d27s19
      include "common.store"
      ldebug=.false.                                                    1d25s21
      if(ldebug)write(6,*)('in countem2: ')                             12d12s19
      do isb=1,nsymb
       do l=1,4                                                         11d19s20
        nkk(l)=0                                                        11d19s20
       end do                                                           11d19s20
       do ii=mdon+1,mdoo+1                                              5d23s23
        iarg=ii-mdon                                                    11d19s20
        do l=1,4                                                        11d19s20
         nkk(l)=nkk(l)+nff2(ii,isb)*ncsf2(l,iarg)                       11d19s20
        end do                                                          11d19s20
       end do                                                           11d19s20
       isymvv=multh(isb,isymmrci)                                       8d27s19
       ncdx(1,isb,1)=nkk(1)                                             11d19s20
       do i=1,3
        ip=i+1                                                          11d19s20
        ncdx(i,isb,2)=nkk(ip)                                           11d19s20
       end do
       nvxv(1,isb)=0                                                    8d2s19
       nvxv(2,isb)=0                                                    8d2s19
       do jsb=1,nsymb                                                   8d2s19
        ksb=multh(jsb,isymvv)                                           8d27s19
        if(jsb.le.ksb)then                                              8d27s19
         if(jsb.ne.ksb)then                                              8d2s19
          nvvs=nvirt(jsb)*nvirt(ksb)                                     8d2s19
          nvvt=nvvs                                                      8d2s19
         else                                                            8d2s19
          nvvs=(nvirt(jsb)*(nvirt(jsb)+1))/2                             8d2s19
          nvvt=(nvirt(jsb)*(nvirt(jsb)-1))/2                             8d2s19
         end if                                                          8d2s19
         if(ldebug)write(6,*)isb,jsb,ksb,nvvs,nvvt                      12d12s19
         nvxv(1,isb)=nvxv(1,isb)+nvvs                                    8d2s19
         nvxv(2,isb)=nvxv(2,isb)+nvvt                                    8d2s19
        end if                                                          8d27s19
       end do                                                           8d2s19
       if(ldebug)write(6,*)nvxv(1,isb),nvxv(2,isb)                      12d12s19
      end do
      return
      end
