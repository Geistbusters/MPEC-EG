c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine square(x,n)
      implicit real*8 (a-h,o-z)
c
c     take x stored as upper triangle by cols
c     and convert it into square storage.
c
      dimension x(1)
      do 1 ir=n,1,-1
       i1=(ir-1)*n
       i2=ir*(ir-1)/2
       do 2 il=ir,1,-1
        x(i1+il)=x(i2+il)
    2  continue
    1 continue
      do 3 ir=1,n
       do 4 il=1,ir
        i1=il+(ir-1)*n
        i2=ir+(il-1)*n
        x(i2)=x(i1)
    4  continue
    3 continue
      return
      end
