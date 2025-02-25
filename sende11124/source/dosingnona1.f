c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dosingnona1(idat1,n1,cvec,g,nroot,nhere,ih0,           6d13s22
     $     iaorb,nalpha,i2e,iborb,nbeta,kadd,ism,irelo,iacto,ioffdetb,  6d10s22
     $     ipxder,bc,ibc,igoal)                                               11d14s22
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2,ism,irelo                                   12d13s06
      integer*2 idat1(4,n1),idwstmp2                                    5d30s06
      integer*1 iaorb(nalpha,1),iborb(nbeta,1),idwstmp1(2)              5d30s06
      equivalence (idwstmp2,idwstmp1)                                   5d30s06
      include "common.store"                                            10d6s06
      dimension cvec(nroot,nhere,*),g(nroot,nhere,*),                   6d13s22
     $     ih0(1),i2e(1),ism(1),irelo(1),iacto(8),ipxder(4,8,8,8)       6d28s22
      dimension phs(2)
      data phs/1d0,-1d0/
      do i=1,n1
       ib=idat1(1,i)
       ik=idat1(2,i)
       ibo=ib+ioffdetb                                                  6d10s22
       idwstmp2=idat1(3,i)                                              5d30s06
       isbb=ism(idwstmp1(1))                                            8d29s06
       isbk=ism(idwstmp1(2))                                            8d29s06
       inb=irelo(idwstmp1(1))-1                                         6d13s22
       ink=irelo(idwstmp1(2))-1                                         6d13s22
       ii1=inb+iacto(isbb)*ink                                          6d13s22
       ii2=ink+iacto(isbk)*inb-ii1                                      6d13s22
       iad=ih0(isbb)+ii1                                                6d13s22
       iadh=iad
       sum=bc(iad)                                                      8d29s06
       sumh=sum                                                         12d15s06
       na=iacto(isbb)*iacto(isbk)                                       6d13s22
       do j=1,nalpha                                                    9d1s06
        isx=ism(iaorb(j,ibo))                                           6d10s22
        i2eu=ipxder(1,isbb,isbk,isx)                                    6d13s22
        iisx=irelo(iaorb(j,ibo))-1                                      6d13s22
        ii3=iisx*(iacto(isx)+1)                                         6d13s22
        irow=ii1+ipxder(2,isbb,isbk,isx)*ii2                            6d28s22
        ix12=irow+na*ii3                                                6d28s22
        nrr=iacto(isx)**2                                               6d28s22
        ix21=ii3+nrr*irow-ix12                                          6d28s22
        iad1=i2e(i2eu)+ix12+ipxder(4,isbb,isbk,isx)*ix21                6d28s22
        sum=sum+bc(iad1)                                                9d1s06
        i2eu=ipxder(1,isbb,isx,isx)                                     6d13s22
        ixa=inb+iacto(isbb)*iisx                                        6d13s22
        ixb=iisx+iacto(isx)*inb-ixa                                     6d13s22
        irow=ixa+ipxder(2,isbb,isx,isx)*ixb                             6d28s22
        nrow=iacto(isbb)*iacto(isx)                                     6d28s22
        ixc=iisx+iacto(isx)*ink                                         6d13s22
        ixd=ink+iacto(isbk)*iisx-ixc                                    6d13s22
        icol=ixc+ipxder(3,isbb,isx,isx)*ixd                             6d28s22
        ix12=irow+nrow*icol                                             6d28s22
        nrr=iacto(isx)*iacto(isbk)                                      6d28s22
        ix21=icol+nrr*irow-ix12                                         6d28s22
        iad1=i2e(i2eu)+ix12+ipxder(4,isbb,isx,isx)*ix21                 6d28s22
        sum=sum-bc(iad1)                                                9d1s06
       end do                                                           9d1s06
       sumk=sum
       do k=1,nhere                                                     3d14s17
        ka=k+kadd                                                       3d14s17
        sum1=sum                                                        3d14s17
        do l=1,nbeta                                                    3d14s17
         isx=ism(iborb(l,ka))                                           9d1s06
         iisx=irelo(iborb(l,ka))-1                                      6d13s22
         i2eu=ipxder(1,isbb,isbk,isx)                                   6d13s22
         ii3=iisx*(iacto(isx)+1)                                        6d13s22
         irow=ii1+ipxder(2,isbb,isbk,isx)*ii2                           6d28s22
         ix12=irow+na*ii3                                               6d28s22
         nrr=iacto(isx)**2                                              6d28s22
         ix21=ii3+nrr*irow-ix12                                         6d28s22
         iad1=i2e(i2eu)+ix12+ipxder(4,isbb,isbk,isx)*ix21               6d28s22
         sum1=sum1+bc(iad1)                                             9d1s06
        end do                                                          9d1s06
        tmp=sum1*phs(idat1(4,i))                                        5d30s06
        do j=1,nroot                                                    3d14s17
         g(j,k,ib)=g(j,k,ib)+tmp*cvec(j,k,ik)                            3d14s17
        end do                                                          3d14s17
       end do                                                           3d14s17
      end do
      return
      end
