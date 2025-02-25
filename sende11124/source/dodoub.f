c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dodoub(idat2,ndoub,cvec,g,numb,nroot,nhere,i2e,iacto,  3d21s17
     $                  ism,irelo,bc,ibc)                               11d14s22
      implicit real*8 (a-h,o-z)
      integer*2 idat2(4,*),idwstmp2(2)                                  1d19s23
      integer*1 idwstmp1(4)                                             5d30s06
      equivalence (idwstmp2,idwstmp1)                                   5d31s06
      integer*8 ism,irelo                                               12d13s06
      include "common.store"                                            10d6s06
      dimension cvec(nroot*nhere,numb),g(nroot*nhere,numb),             3d14s17
     $     i2e(*),iacto(8),ism(*),irelo(*),phs(2)                       1d19s23
      data phs/1d0,-1d0/
      data icall/0/
      save icall
      do i=1,ndoub
       i1signed=idat2(1,i)                                              3d12s17
       i1=iabs(i1signed)                                                3d12s17
       iphdws=(((idat2(1,i)-i1)/(2*i1)))+2                              3d12s17
       if(iphdws.lt.1.or.iphdws.gt.2)then
        write(6,*)i1signed,iphdws
        stop 'dodoub'
       end if
       i2=idat2(2,i)
       idwstmp2(1)=idat2(3,i)                                           5d30s06
       idwstmp2(2)=idat2(4,i)                                           5d30s06
       n1=irelo(idwstmp1(1))                                            9d1s06
       n2=irelo(idwstmp1(2))                                            9d1s06
       n3=irelo(idwstmp1(3))                                            9d1s06
       n4=irelo(idwstmp1(4))                                            9d1s06
       isx1=ism(idwstmp1(1))                                             9d1s06
       isx2=ism(idwstmp1(2))                                             9d1s06
       isx3=ism(idwstmp1(3))                                             9d1s06
       isx4=ism(idwstmp1(4))                                             9d1s06
       i2eu=ifind2(isx1,isx3,isx2,isx4,icase)                           12d19s06
       if(isx1.eq.isx3)then                                             12d19s06
        in=min(n1,n3)                                                   12d19s06
        ix=max(n1,n3)                                                   12d19s06
        ii1=((ix*(ix-1))/2)+in-1                                        9d1s06
        in1save=in
        ix1save=ix
        in=min(n2,n4)                                                   12d19s06
        ix=max(n2,n4)                                                   12d19s06
        in2save=in
        ix2save=ix
        ii2=((ix*(ix-1))/2)+in-1                                        9d1s06
        na=(iacto(isx1)*(iacto(isx1)+1))/2                              9d1s06
        iad1=i2e(i2eu)+ii1+na*ii2                                       9d1s06
        ii1save=ii1
        ii2save=ii2
        nasave=na
       else                                                             9d1s06
        if(icase.eq.1)then                                              9d1s06
         iad1=i2e(i2eu)+n1-1+iacto(isx1)*(n3-1+iacto(isx3)              12d19s06
     $       *(n2-1+iacto(isx2)*(n4-1)))                                12d19s06
        else if(icase.eq.2)then                                         9d1s06
         iad1=i2e(i2eu)+n3-1+iacto(isx3)*(n1-1+iacto(isx1)              12d19s06
     $       *(n4-1+iacto(isx4)*(n2-1)))                                12d19s06
        else if(icase.eq.3)then                                         9d1s06
         iad1=i2e(i2eu)+n1-1+iacto(isx1)*(n3-1+iacto(isx3)              12d19s06
     $       *(n4-1+iacto(isx4)*(n2-1)))                                12d19s06
        else if(icase.eq.4)then                                         9d1s06
         iad1=i2e(i2eu)+n3-1+iacto(isx3)*(n1-1+iacto(isx1)              12d19s06
     $       *(n2-1+iacto(isx2)*(n4-1)))                                12d19s06
        end if                                                          9d1s06
       end if                                                           9d1s06
       i2eu=ifind2(isx3,isx2,isx1,isx4,icase)                           9d1s06
       if(isx3.eq.isx2)then                                             9d1s06
        in=min(n3,n2)                                                   9d1s06
        ix=max(n3,n2)                                                   9d1s06
        ii1=((ix*(ix-1))/2)+in-1                                        9d1s06
        na=(iacto(isx3)*(iacto(isx3)+1))/2                              9d1s06
        in=min(n1,n4)                                                   9d1s06
        ix=max(n1,n4)                                                   9d1s06
        ii2=((ix*(ix-1))/2)+in-1                                        9d1s06
        iad2=i2e(i2eu)+ii1+na*ii2                                       9d1s06
       else                                                             9d1s06
        if(icase.eq.1)then                                              9d1s06
         iad2=i2e(i2eu)+n3-1+iacto(isx3)*(n2-1+iacto(isx2)              9d1s06
     $       *(n1-1+iacto(isx1)*(n4-1)))                                9d1s06
        else if(icase.eq.2)then                                         9d1s06
         iad2=i2e(i2eu)+n2-1+iacto(isx2)*(n3-1+iacto(isx3)              9d1s06
     $        *(n4-1+iacto(isx4)*(n1-1)))                               9d1s06
        else if(icase.eq.3)then                                         9d1s06
         iad2=i2e(i2eu)+n3-1+iacto(isx3)*(n2-1+iacto(isx2)              9d1s06
     $        *(n4-1+iacto(isx4)*(n1-1)))                               9d1s06
        else if(icase.eq.4)then                                         9d1s06
         iad2=i2e(i2eu)+n2-1+iacto(isx2)*(n3-1+iacto(isx3)              9d1s06
     $        *(n1-1+iacto(isx1)*(n4-1)))                               9d1s06
        end if                                                          9d1s06
       end if                                                           9d1s06
       tmp=phs(iphdws)*(bc(iad1)-bc(iad2))                              10d6s06
       do jk=1,nhere*nroot                                              3d14s17
        g(jk,i1)=g(jk,i1)+tmp*cvec(jk,i2)                               3d14s17
        g(jk,i2)=g(jk,i2)+tmp*cvec(jk,i1)                               3d14s17
       end do
      end do
      return
      end
