c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine compcomp(nrow,icmp1a,icmp2a,cmp3a,nmid,icmp1,icmp2,    2d9s23
     $     cmp3,xout,ndim,ncol,ff,fff,imap)                             2d9s23
      implicit real*8 (a-h,o-z)
c
c     form xout=xout*ff+fff*compresseda*compressed
c
      dimension icmp1a(2,*),icmp2a(*),cmp3a(*),icmp1(2,*),icmp2(*),     2d9s23
     $     cmp3(*),xout(ndim,*),imap(*)                                 2d9s23
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
      iiua=1                                                            2d9s23
      do i=1,nmid                                                       2d9s23
       imap(i)=iiua                                                     2d9s23
       iiua=iiua+(icmp1a(2,i)+1-icmp1a(1,i))                            2d9s23
      end do                                                            2d9s23
      iiu=1                                                             2d9s23
      do i=1,ncol                                                       2d9s23
       do j=icmp1(1,i),icmp1(2,i)
        fact=cmp3(iiu)*fff                                              2d9s23
        jju=imap(icmp2(iiu))                                            2d9s23
        do k=icmp1a(1,icmp2(iiu)),icmp1a(2,icmp2(iiu))                  2d9s23
         xout(icmp2a(jju),i)=xout(icmp2a(jju),i)+cmp3a(jju)*fact        2d9s23
         jju=jju+1                                                      2d9s23
        end do                                                          2d9s23
        iiu=iiu+1                                                       2d9s23
       end do                                                           2d9s23
      end do                                                            2d9s23
      return                                                            2d9s23
      end                                                               2d9s23
