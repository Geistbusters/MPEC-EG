c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine msl(lab,x,n1,n2)
      implicit real*8 (a-h,o-z)
      character*(*)lab
      dimension x(n1,n2),sums(4)
      do i=1,2
       sums(i)=0d0
      end do
      do i=1,n2
       do j=1,n1
        sums(1)=sums(1)+x(j,i)
        sums(2)=sums(2)+x(j,i)**2
       end do
      end do
      sums(3)=sums(1)
      sums(4)=sums(2)
      return
      end
