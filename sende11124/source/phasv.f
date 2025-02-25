c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine phasv(v,id,nr,nc)
      implicit real*8 (a-h,o-z)
c
c     phase vectors so largest component is positive
c
      dimension v(id,nc)
      do 1 ir=1,nc
       x=abs(v(1,ir))                                                   10d17s23
       do il=2,nr                                                       10d17s23
        x=max(x,abs(v(il,ir)))                                          10d17s23
       end do                                                           10d17s23
       thrs=0.9999d0*x                                                  10d17s23
       do il=1,nr                                                       10d17s23
        if(abs(v(il,ir)).ge.thrs)then                                   10d17s23
         ix=il                                                          10d17s23
         go to 22                                                       10d17s23
        end if                                                          10d17s23
       end do                                                           10d17s23
   22  continue                                                         10d17s23
       if(v(ix,ir).lt.0d0)then
        do 3 il=1,nr
         v(il,ir)=-v(il,ir)
    3   continue
       end if
    1 continue
      return
      end
