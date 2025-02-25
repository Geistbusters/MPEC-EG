c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      real*8 function cleb2(j1,m1,j2,m2,j3,m3)
      implicit real*8 (a-h,o-z)
c
c     the 2 means args are twice values so always integers.
c
      cleb2=f3j(j1,j2,j3,m1,m2,-m3,1)
      cleb2=cleb2*sqrt(dfloat(j3+1))
      ipp=j2-j1-m3
      ipp=ipp/2
      if(mod(ipp,2).ne.0)cleb2=-cleb2
      return
      end
