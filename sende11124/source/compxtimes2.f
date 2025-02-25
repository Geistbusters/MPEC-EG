c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine compxtimes2(icmp1,icmp2,cmp3,xin,ndim,nmul,xout,ndim2, 2d22s21
     $     nrow,mcol,ff,icsf0,mcsf,bc,ibc)                              11d14s22
      implicit real*8 (a-h,o-z)
      dimension icmp1(2,*),icmp2(*),cmp3(*),xout(ndim2,*),              2d22s21
     $     xin(ndim,*)
      include "common.store"
      if(ff.eq.0d0)then                                                 12d5s20
       do i=1,mcol                                                      2d22s21
        do j=1,nrow                                                     2d22s21
         xout(j,i)=0d0                                                  12d5s20
        end do                                                          12d5s20
       end do                                                           12d5s20
      else                                                              12d5s20
       do i=1,mcol                                                      2d22s21
        do j=1,nrow                                                     2d22s21
         xout(j,i)=ff*xout(j,i)                                          12d5s20
        end do
       end do
      end if                                                            12d5s20
      if(mod(mcol,2).eq.0)then                                          2d22s21
       kbot=1                                                           12d6s20
      else                                                              12d6s20
       kbot=2                                                           12d6s20
      end if                                                            12d6s20
      kbotm=kbot-1                                                      12d6s20
      ii=1
      itop=icsf0+mcsf-1                                                 2d22s21
      ii=1
      do i=1,itop                                                       2d22s21
       if(i.ge.icsf0.and.i.le.itop)then                                 2d22s21
        im=i+1-icsf0                                                    2d22s21
        do k=1,kbotm                                                     12d6s20
         iiu=ii                                                          12d6s20
         do j=icmp1(1,i),icmp1(2,i)                                      12d6s20
          xout(icmp2(iiu),k)=xout(icmp2(iiu),k)+cmp3(iiu)*xin(im,k)     2d22s21
          iiu=iiu+1                                                      12d6s20
         end do                                                          12d6s20
        end do                                                           12d6s20
        do k=kbot,mcol,2                                                 2d22s21
         kp=k+1                                                          12d6s20
         iiu=ii                                                          12d6s20
         do j=icmp1(1,i),icmp1(2,i)                                      12d6s20
          xout(icmp2(iiu),k)=xout(icmp2(iiu),k)+cmp3(iiu)*xin(im,k)     2d22s21
          xout(icmp2(iiu),kp)=xout(icmp2(iiu),kp)+cmp3(iiu)*xin(im,kp)  2d22s21
          iiu=iiu+1                                                      12d6s20
         end do                                                          12d6s20
        end do                                                           12d6s20
        ii=iiu                                                          2d22s21
       else                                                             2d22s21
        iio=ii
        ii=ii+icmp1(2,i)+1-icmp1(1,i)                                   2d22s21
       end if                                                           2d22s21
      end do                                                            12d5s20
      return
      end
