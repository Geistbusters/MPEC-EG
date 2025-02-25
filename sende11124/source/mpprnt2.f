c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mpprnt2(a,n)
      implicit real*8   (a-h,o-z)
c
c     print a matrix stored in packed form - the lower triangle by
c     rows or the upper triangle by cols.                               4d24s96
c
      dimension a(n*(n+1)/2)
      iunxp=6
      j1=1
c     numcol=9                                                          9d13s88
      numcol=5                                                          9d13s88
      numadd=numcol-1                                                   9d13s88
      j2=min0(numcol,n)                                                 9d13s88
      write(iunxp,20)n,n                                                9d30s91
   20 format(1x,2i5,'@')                                                2d11s89
    1 continue
       write(iunxp,2)(j,j=j1,j2)                                        9d30s91
    2  format(1x,3hm-n,9(4x,i5,5x))                                     2d11s89
       do 3 i=j1,n                                                      10d6s93
        write(iunxp,4)i,(a(((i*(i-1))/2)+j),j=j1,min0(i,j2))            4d24s96
    4   format(1x,i4,1x,1p9e14.6)                                       11d7s91
    3  continue
       if(j2.eq.n)then                                                  2d11s89
        write(iunxp,121)                                                9d30s91
  121   format(1x,'$')                                                  2d11s89
        return                                                          2d11s89
       end if                                                           2d11s89
       j1=j2+1
       j2=min0(j1+numadd,n)                                             9d13s88
      go to 1
      end
