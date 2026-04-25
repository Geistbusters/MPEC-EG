c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine step(npts,h,c,ider)                                    2d15s91
      implicit real  *8 (a-h,o-z)
c
c      calculate the coefficients in the f.d. formula for the ider'th   2d15s91
c      derivative.
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0)                         11d25s88
      dimension h(npts),c(npts)                                         11d25s88
      hmax=abs(h(1))                                                    11d25s88
      if(ider.gt.npts)then
       write(6,2)ider,npts                                              2d15s91
    2  format(/1x,'in step, order of derivative = ',i5,                 2d15s91
     $        /1x,'exceeds number of points = ',i5)                     2d15s91
       stop                                                             2d15s91
      end if                                                            2d15s91
      do 77 i=1,npts
       hmax=max(hmax,abs(h(i)))                                         6d2s98
       c(i)=zero                                                        11d25s88
   77 continue
      c(ider+1)=one                                                     2d15s91
      hscal=one/hmax                                                    7d31s88
      hscal=hscal*(2d0-hscal*hmax)                                      12d8s88
      do 79 i=1,npts
       h(i)=h(i)*hscal                                                  11d25s88
   79 continue
      call vandr2(h,c,npts)                                             11d25s88
      fact=one
      do 1 i=1,ider                                                     2d15s91
       fact=fact*hscal*dfloat(i)                                        2d15s91
    1 continue                                                          2d15s91
      do 83 i=1,npts
       c(i)=c(i)*fact                                                   11d25s88
       h(i)=h(i)*hmax                                                   5d2s90
   83 continue
      return
      end
