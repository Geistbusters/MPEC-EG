c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine avgr(xin,navg,tovr)                                    3d1s21
      implicit real*8 (a-h,o-z)                                         3d1s21
      dimension xin(*)                                                  3d1s21
      tovr=0d0                                                          3d1s21
      do i=2,navg+1                                                     3d1s21
       im=i-1                                                           3d1s21
       tdelta=xin(i)-xin(im)                                            3d1s21
       tovr=tovr+tdelta                                                 3d1s21
      end do                                                            3d1s21
      tovr=tovr/dfloat(navg)                                            3d1s21
      return                                                            3d1s21
      end                                                               3d1s21
