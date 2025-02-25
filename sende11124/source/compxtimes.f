c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine compxtimes(icmp1,icmp2,cmp3,xin,ndim,nrow,xout,ndim2,  12d5s20
     $     ncsfb,ncsfk,ff,fff)                                          3d25s22
      implicit real*8 (a-h,o-z)
      dimension icmp1(2,*),icmp2(*),cmp3(*),xout(ndim2,*),              1d18s23
     $     xin(ndim,*)
      data icall/0/
      save icall
      icall=icall+1
      if(ff.eq.0d0)then                                                 12d5s20
       do i=1,ncsfk                                                     12d5s20
        do j=1,ncsfb                                                    12d5s20
         xout(j,i)=0d0                                                  12d5s20
        end do                                                          12d5s20
       end do                                                           12d5s20
      else                                                              12d5s20
       do i=1,ncsfk
        do j=1,ncsfb
         xout(j,i)=ff*xout(j,i)                                          12d5s20
        end do
       end do
      end if                                                            12d5s20
      kbot=mod(ncsfk,4)+1                                               3d11s21
      kbotm=kbot-1                                                      12d6s20
      ii=1
      do i=1,nrow                                                       12d6s20
       do k=1,kbotm                                                     12d6s20
        iiu=ii                                                          12d6s20
        fact=xin(i,k)*fff                                               3d25s22
        do j=icmp1(1,i),icmp1(2,i)                                      12d6s20
         xout(icmp2(iiu),k)=xout(icmp2(iiu),k)+cmp3(iiu)*fact           3d25s22
         iiu=iiu+1                                                      12d6s20
        end do                                                          12d6s20
       end do                                                           12d6s20
       do k=kbot,ncsfk,4                                                12d6s20
        kp=k+1                                                          12d6s20
        kq=kp+1                                                         3d11s21
        kr=kq+1                                                         3d11s21
        iiu=ii                                                          12d6s20
        fact=xin(i,k)*fff                                               3d25s22
        factp=xin(i,kp)*fff                                             3d25s22
        factq=xin(i,kq)*fff                                             3d25s22
        factr=xin(i,kr)*fff                                             3d25s22
        do j=icmp1(1,i),icmp1(2,i)                                      12d6s20
         xout(icmp2(iiu),k)=xout(icmp2(iiu),k)+cmp3(iiu)*fact           3d25s22
         xout(icmp2(iiu),kp)=xout(icmp2(iiu),kp)+cmp3(iiu)*factp        3d25s22
         xout(icmp2(iiu),kq)=xout(icmp2(iiu),kq)+cmp3(iiu)*factq        3d25s22
         xout(icmp2(iiu),kr)=xout(icmp2(iiu),kr)+cmp3(iiu)*factr        3d25s22
         iiu=iiu+1                                                      12d6s20
        end do                                                          12d6s20
       end do                                                           12d6s20
       ii=iiu                                                           12d6s20
      end do                                                            12d5s20
      return
      end
