c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine uncompxut(icmp1,icmp2,cmp3,wp,ncsfb,ncsfk)
      implicit real*8 (a-h,o-z)
      dimension icmp1(2,ncsfb),icmp2(*),cmp3(*),wp(ncsfb,*)
      do i=1,ncsfk
       do j=1,ncsfb
        wp(j,i)=0d0
       end do
      end do
      ii=1
      do i=1,ncsfb                                                      12d5s20
       do j=icmp1(1,i),icmp1(2,i)
        wp(i,icmp2(ii))=cmp3(ii)                                        12d5s20
        ii=ii+1
       end do
      end do
      return
      end
