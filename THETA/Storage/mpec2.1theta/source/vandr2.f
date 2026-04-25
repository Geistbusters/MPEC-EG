c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine vandr2(x,q,n)
      implicit real  *8 (a-h,o-z)
c
c     solve the vandermode linear system using algorithm  5.6-2
c     from matrix computations.
c     input parameters:
c     x - "nodes"
c     q - "moments"
c     n - number of nodes or moments
c     output parameters:
c     x - unchanged
c     q - "weights"
c     n - unchanged.
c
      dimension x(n),q(n)
      do 1 k=1,n-1
       do 2 i=n,k+1,-1
        q(i)=q(i)-x(k)*q(i-1)
    2  continue
    1 continue
      do 3 k=n-1,1,-1
       do 4 i=k+1,n
        tmpdws=x(i)-x(i-k)                                              12d8s88
        swdpmt=1d0/tmpdws                                               12d8s88
        swdpmt=swdpmt*(2d0-tmpdws*swdpmt)                               12d8s88
        q(i)=q(i)*swdpmt                                                12d8s88
    4  continue
       do 5 i=k,n-1
        q(i)=q(i)-q(i+1)
    5  continue
    3 continue
      return
      end
