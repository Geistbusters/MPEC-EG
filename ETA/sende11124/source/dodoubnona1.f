c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dodoubnona1(idat2,ndoub,cvec,g,nroot,nhere,i2e,iacto,  6d13s22
     $                  ism,irelo,ipxder,bc,ibc,igoal)                        11d14s22
      implicit real*8 (a-h,o-z)
      integer*2 idat2(4,ndoub),idwstmp2(2)                              3d21s17
      integer*1 idwstmp1(4)                                             5d30s06
      equivalence (idwstmp2,idwstmp1)                                   5d31s06
      integer*8 ism,irelo                                               12d13s06
      include "common.store"                                            10d6s06
      dimension cvec(nroot*nhere,*),g(nroot*nhere,*),                   6d13s22
     $     i2e(1),iacto(8),ism(1),irelo(1),phs(2),ipxder(4,8,8,8)       6d28s22
      data phs/1d0,-1d0/
      do i=1,ndoub
       i1signed=idat2(1,i)                                              3d12s17
       i1=iabs(i1signed)                                                3d12s17
       iphdws=(((idat2(1,i)-i1)/(2*i1)))+2                              3d12s17
       i2=idat2(2,i)
       idwstmp2(1)=idat2(3,i)                                           5d30s06
       idwstmp2(2)=idat2(4,i)                                           5d30s06
       n1=irelo(idwstmp1(1))-1                                          6d13s22
       n2=irelo(idwstmp1(2))-1                                          6d13s22
       n3=irelo(idwstmp1(3))-1                                          6d13s22
       n4=irelo(idwstmp1(4))-1                                          6d13s22
       isx1=ism(idwstmp1(1))                                             9d1s06
       isx2=ism(idwstmp1(2))                                             9d1s06
       isx3=ism(idwstmp1(3))                                             9d1s06
       isx4=ism(idwstmp1(4))                                             9d1s06
       i2eu=ipxder(1,isx1,isx3,isx2)                                    6d13s22
       ixa=n1+iacto(isx1)*n3                                            6d13s22
       ixb=n3+iacto(isx3)*n1-ixa                                        6d13s22
       irow=ixa+ipxder(2,isx1,isx3,isx2)*ixb                            6d28s22
       nrow=iacto(isx1)*iacto(isx3)                                     6d28s22
       ixc=n2+iacto(isx2)*n4                                            6d13s22
       ixd=n4+iacto(isx4)*n2-ixc                                        6d13s22
       icol=ixc+ipxder(3,isx1,isx3,isx2)*ixd                            6d28s22
       ix12=irow+nrow*icol                                              6d28s22
       nrr=iacto(isx2)*iacto(isx4)                                      6d28s22
       ix21=icol+nrr*irow-ix12                                          6d28s22
       iad1=i2e(i2eu)+ix12+ipxder(4,isx1,isx3,isx2)*ix21                6d28s22
       i2eu=ipxder(1,isx1,isx4,isx2)                                    6d13s22
       ixa=n1+iacto(isx1)*n4                                            6d13s22
       ixb=n4+iacto(isx4)*n1-ixa                                        6d13s22
       irow=ixa+ipxder(2,isx1,isx4,isx2)*ixb                            6d28s22
       nrow=iacto(isx1)*iacto(isx4)                                     6d28s22
       ixc=n2+iacto(isx2)*n3                                            6d13s22
       ixd=n3+iacto(isx3)*n2-ixc                                        6d13s22
       icol=ixc+ipxder(3,isx1,isx4,isx2)*ixd                            6d28s22
       ix12=irow+nrow*icol                                              6d28s22
       nrr=iacto(isx2)*iacto(isx3)                                      6d28s22
       ix21=icol+nrr*irow-ix12                                          6d28s22
       iad2=i2e(i2eu)+ix12+ipxder(4,isx1,isx4,isx2)*ix21                6d28s22
       tmp=phs(iphdws)*(bc(iad1)-bc(iad2))                              10d6s06
       do jk=1,nhere*nroot                                              3d14s17
        g(jk,i1)=g(jk,i1)+tmp*cvec(jk,i2)                               3d14s17
       end do
      end do
      return
      end
