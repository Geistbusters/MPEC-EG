c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mxmt(a,na,b,nb,c)                                      11d23s20
      implicit real*8 (a-h,o-z)                                         11d23s20
c
c     a*b+c=c, c stored as triangle
c
      dimension a(na,*),b(nb,*),c(*)
      ij=0
      do j=1,na
       do k=1,nb
        do i=1,j
         c(ij+i)=c(ij+i)+a(i,k)*b(k,j)
        end do
       end do
       ij=ij+j
      end do
      return
      end
