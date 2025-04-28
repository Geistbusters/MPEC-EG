c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diaghii1(tmpa,numa,nalpha,iaorb,ih0e,i2e,debug,nsymb,  8d8s06
     $     ndet,ism,irelo,bc,ibc)                                       11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 iaorb(nalpha,numa)
      dimension tmpa(numa)
      include "common.store"
      include "common.cas"
      integer*8 ism(1),irelo(1)                                         5d7s18
      dimension i2e(1),ih0e(1),ioff(8),ioffb(8),ndet(1)                 8d8s06
      do i=1,nsymb                                                      8d8s06
       if(i.eq.1)then                                                   8d8s06
        ioff(i)=ndet(i)                                                 8d8s06
        ioffb(i)=0                                                      8d8s06
       else                                                             8d8s06
        ioff(i)=ioff(i-1)+ndet(i)                                       8d8s06
        ioffb(i)=ioffb(i-1)+iacto(i-1)                                  8d8s06
       end if                                                           8d8s06
      end do                                                            8d8s06
      isym=1                                                            8d8s06
      do i=1,numa
       sum=0d0
       do j=1,nalpha
        jo=irelo(iaorb(j,i))                                            8d8s06
        jj=((jo*(jo+1))/2)-1
        nj=ism(iaorb(j,i))                                              8d8s06
        na=(iacto(nj)*(iacto(nj)+1))/2                                  8d8s06
        iadd=ih0e(nj)+(jo-1)*(iacto(nj)+1)                              8d8s06
        sum=sum+bc(iadd)                                                7d12s06
        do jp=1,nalpha
         jpo=irelo(iaorb(jp,i))                                         8d8s06
         jjp=((jpo*(jpo+1))/2)-1
         njp=ism(iaorb(jp,i))                                           8d8s06
         i2eu=ifind2(nj,nj,njp,njp,icase)                               8d15s06
   32    format(4i1,2i5)
         jeu=i2eu                                                       12d11s06
         iadd1=i2e(i2eu)+jj+jjp*na                                      8d8s06
         i2eu=ifind2(nj,njp,nj,njp,icase)                               8d15s06
         if(icase.eq.2.or.icase.eq.4)then                               7d5s18
          iadd2=i2e(i2eu)+jpo-1+iacto(njp)*(jo-1+iacto(nj)              8d15s06
     $         *(jpo-1+iacto(njp)*(jo-1)))
         else if(nj.eq.njp)then
          ix=max(jpo,jo)
          in=min(jpo,jo)
          kk=((ix*(ix-1))/2)+in-1
          iadd2=i2e(i2eu)+kk*(na+1)                                     8d14s06
         else                                                           8d8s06
          iadd2=i2e(i2eu)+jo-1+iacto(nj)*(jpo-1+iacto(njp)              8d14s06
     $         *(jo-1+iacto(nj)*(jpo-1)))                               8d8s06
         end if                                                         8d8s06
         orig=sum
         sum=sum+0.5d0*(bc(iadd1)-bc(iadd2))                            7d12s06
        end do
       end do
       tmpa(i)=sum
 3878  format(es21.14,10i3)
      end do
      return
      end
