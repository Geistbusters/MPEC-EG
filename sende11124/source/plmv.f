c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine plmv(plm,l,m,g,nqd)                                    5d12s88
      implicit real  *8 (a-h,o-z)                                       5d12s88
c                                                                       5d12s88
c     this calculates the associated legendre functions                 5d12s88
c      m      m           m                                             5d12s88
c     p (g), p (g), ..., p (g)                                          5d12s88
c      m      m+1         l                                             5d12s88
c     and stores the result in plm. m must be positive.                 5d12s88
c     the definition used is taken from edmonds.                        5d12s88
c     note pmm is stored as the first column.                           5d12s88
c                                                                       5d12s88
      parameter (one=1.0d0)                                             5d12s88
      dimension plm(nqd,l+1-m),g(nqd)                                   5d12s88
c                                                                       5d12s88
c     calculate pmm                                                     5d12s88
c                                                                       5d12s88
      if(l.ge.m)then                                                    5d12s88
       xnorm=one                                                        5d12s88
       do 1 i=1,m                                                       5d12s88
        xnorm=xnorm* dfloat(2*i-1)                                      5d12s88
    1  continue                                                         5d12s88
       mh=m/2                                                           5d12s88
       if(mh*2.eq.m)then                                                5d12s88
        do 2 i=1,nqd                                                    5d12s88
         plm(i,1)=one-g(i)*g(i)                                         5d12s88
         if(plm(i,1).eq.0d0)then                                        2d12s92
          if(mh.eq.0)then                                               2d12s92
           plm(i,1)=xnorm                                               2d12s92
          else                                                          2d12s92
           plm(i,1)=0d0                                                 2d12s92
          end if                                                        2d12s92
         else
          plm(i,1)=xnorm*(plm(i,1)**mh)                                 2d12s92
         end if                                                         2d12s92
    2   continue                                                        5d12s88
       else                                                             5d12s88
        do 3 i=1,nqd                                                    5d12s88
         plm(i,1)=one-g(i)*g(i)                                         5d12s88
         if(plm(i,1).eq.0d0)then                                        2d12s92
          plm(i,1)=0d0                                                  2d12s92
         else                                                           2d12s92
          plm(i,1)=xnorm*(plm(i,1)**mh)*sqrt(plm(i,1))                  5d12s88
         end if                                                         2d12s92
    3   continue                                                        5d12s88
       end if                                                           5d12s88
      end if                                                            5d12s88
c                                                                       5d12s88
c     calculate pm+1m                                                   5d12s88
c                                                                       5d12s88
      if(l.ge.m+1)then                                                  5d12s88
       xnorm= dfloat(2*m+1)                                             5d12s88
       do 4 i=1,nqd                                                     5d12s88
        plm(i,2)=xnorm*g(i)*plm(i,1)                                    5d12s88
    4  continue                                                         5d12s88
      end if                                                            5d12s88
c                                                                       5d12s88
c     calculate the remaining functions by recursion                    5d12s88
c                                                                       5d12s88
      if(l.ge.m+2)then                                                  5d12s88
       do 5 j=3,l-m+1                                                   5d12s88
        tmpdws= dfloat(j-1)                                             12d8s88
        swdpmt=1d0/tmpdws                                               12d8s88
        swdpmt=swdpmt*(2d0-swdpmt*tmpdws)                               12d8s88
        xnorm1= dfloat(2*(j+m)-3)*swdpmt                                12d8s88
        xnorm2= dfloat(j+2*m-2)*swdpmt                                  12d8s88
        do 5 i=1,nqd                                                    5d12s88
         plm(i,j)=g(i)*xnorm1*plm(i,j-1)-xnorm2*plm(i,j-2)              5d12s88
    5  continue                                                         5d12s88
      end if                                                            5d12s88
      return                                                            5d12s88
      end                                                               1d27s88
