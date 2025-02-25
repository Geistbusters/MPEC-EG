c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine uncompxu(icmp1,icmp2,cmp3,wp,ncsfb,ncsfk)
      implicit real*8 (a-h,o-z)
      dimension icmp1(2,*),icmp2(*),cmp3(*),wp(ncsfb,*)                 1d18s23
      do i=1,ncsfk
       do j=1,ncsfb
        wp(j,i)=0d0
       end do
      end do
      ii=1
      do i=1,ncsfk                                                      12d5s20
       do j=icmp1(1,i),icmp1(2,i)
        wp(icmp2(ii),i)=cmp3(ii)                                        12d5s20
        ii=ii+1
       end do
      end do
      return
      end
