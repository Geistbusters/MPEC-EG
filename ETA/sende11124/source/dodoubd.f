c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dodoubd(idat2,ndoub,cvecb,cveck,numb,nroot,nhere,iacto,5d27s22
     $                  ism,irelo,jdenpt,wgt,bc,ibc)                    11d14s22
      implicit real*8 (a-h,o-z)
      integer*2 idat2(4,ndoub),idwstmp2(2)                              3d21s17
      integer*1 idwstmp1(4)                                             5d30s06
      equivalence (idwstmp2,idwstmp1)                                   5d31s06
      integer*8 ism,irelo                                               12d13s06
      include "common.store"                                            10d6s06
      dimension cvecb(nroot*nhere,numb),cveck(nroot*nhere,numb),        5d27s22
     $     wgt(nroot),jdenpt(*),iacto(8),ism(1),irelo(1),phs(2)         5d27s22
      dimension igoalx(2)
      data phs/1d0,-1d0/
      do i=1,ndoub
       i1signed=idat2(1,i)                                              3d12s17
       i1=iabs(i1signed)                                                3d12s17
       iphdws=(((idat2(1,i)-i1)/(2*i1)))+2                              3d12s17
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
        iad1=jdenpt(i2eu)+ii1+na*ii2                                    3d27s17
        ii1save=ii1
        ii2save=ii2
        nasave=na
       else                                                             9d1s06
        if(icase.eq.1)then                                              9d1s06
         iad1=jdenpt(i2eu)+n1-1+iacto(isx1)*(n3-1+iacto(isx3)           3d27s17
     $       *(n2-1+iacto(isx2)*(n4-1)))                                12d19s06
        else if(icase.eq.2)then                                         9d1s06
         iad1=jdenpt(i2eu)+n3-1+iacto(isx3)*(n1-1+iacto(isx1)           3d27s17
     $       *(n4-1+iacto(isx4)*(n2-1)))                                12d19s06
        else if(icase.eq.3)then                                         9d1s06
         iad1=jdenpt(i2eu)+n1-1+iacto(isx1)*(n3-1+iacto(isx3)           3d27s17
     $       *(n4-1+iacto(isx4)*(n2-1)))                                12d19s06
        else if(icase.eq.4)then                                         9d1s06
         iad1=jdenpt(i2eu)+n3-1+iacto(isx3)*(n1-1+iacto(isx1)           3d27s17
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
        iad2=jdenpt(i2eu)+ii1+na*ii2                                    3d27s17
       else                                                             9d1s06
        if(icase.eq.1)then                                              9d1s06
         iad2=jdenpt(i2eu)+n3-1+iacto(isx3)*(n2-1+iacto(isx2)           3d27s17
     $       *(n1-1+iacto(isx1)*(n4-1)))                                9d1s06
        else if(icase.eq.2)then                                         9d1s06
         iad2=jdenpt(i2eu)+n2-1+iacto(isx2)*(n3-1+iacto(isx3)           3d27s17
     $        *(n4-1+iacto(isx4)*(n1-1)))                               9d1s06
        else if(icase.eq.3)then                                         9d1s06
         iad2=jdenpt(i2eu)+n3-1+iacto(isx3)*(n2-1+iacto(isx2)           3d27s17
     $        *(n4-1+iacto(isx4)*(n1-1)))                               9d1s06
        else if(icase.eq.4)then                                         9d1s06
         iad2=jdenpt(i2eu)+n2-1+iacto(isx2)*(n3-1+iacto(isx3)           3d27s17
     $        *(n1-1+iacto(isx1)*(n4-1)))                               9d1s06
        end if                                                          9d1s06
       end if                                                           9d1s06
       do j=1,nhere                                                     3d27s17
        do k=1,nroot                                                    3d27s17
         jk=k+nroot*(j-1)                                               3d27s17
         ww=phs(iphdws)*wgt(k)*(cvecb(jk,i2)*cveck(jk,i1)               6d16s22
     $        +cveck(jk,i2)*cvecb(jk,i1))                               6d16s22
         bc(iad1)=bc(iad1)+ww                                           6d16s22
         bc(iad2)=bc(iad2)-ww                                           6d16s22
        end do                                                          3d27s17
       end do                                                           3d27s17
      end do
      return
      end
