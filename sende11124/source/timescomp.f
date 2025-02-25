c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine timescomp(xin,ndim,nrow,icmp1,icmp2,cmp3,              2d9s23
     $     xout,ndim2,ncol,ff,fff)                                      2d9s23
      implicit real*8 (a-h,o-z)
c
c     form xout=xout*ff+fff*xin*compressed
c
      dimension xin(ndim,*),icmp1(2,*),icmp2(*),cmp3(*),xout(ndim2,*)   2d9s23
      if(ff.eq.0d0)then                                                 2d9s23
       do i=1,ncol                                                      2d9s23
        do j=1,nrow                                                     2d9s23
         xout(j,i)=0d0                                                  2d9s23
        end do                                                          2d9s23
       end do                                                           2d9s23
      else if(abs(ff-1d0).gt.1d-14)then                                 2d9s23
       do i=1,ncol                                                      2d9s23
        do j=1,nrow                                                     2d9s23
         xout(j,i)=ff*xout(j,i)                                         2d9s23
        end do                                                          2d9s23
       end do                                                           2d9s23
      end if                                                            2d9s23
      iiu=1                                                             2d9s23
      do i=1,ncol                                                       2d9s23
       do j=icmp1(1,i),icmp1(2,i)
        fact=cmp3(iiu)*fff                                              2d9s23
        do k=1,nrow                                                     2d9s23
         xout(k,i)=xout(k,i)+xin(k,icmp2(iiu))*fact                     2d9s23
        end do                                                          2d9s23
        iiu=iiu+1                                                       2d9s23
       end do                                                           2d9s23
      end do                                                            2d9s23
      return                                                            2d9s23
      end                                                               2d9s23
