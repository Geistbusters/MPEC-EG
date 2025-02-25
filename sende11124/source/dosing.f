c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dosing(idat1,n1,cvec,g,numb,nroot,nhere,ih0,           3d10s17
     $     iaorb,nalpha,i2e,iborb,nbeta,kadd,ism,irelo,iacto,ioffdet,   11d14s22
     $     bc,ibc)                                                      11d14s22
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2,ism,irelo                                   12d13s06
      integer*2 idat1(4,n1),idwstmp2                                    5d30s06
      integer*1 iaorb(nalpha,1),iborb(nbeta,1),idwstmp1(2)              5d30s06
      equivalence (idwstmp2,idwstmp1)                                   5d30s06
      include "common.store"                                            10d6s06
      dimension cvec(nroot,nhere,numb),g(nroot,nhere,numb),             3d14s17
     $     ih0(1),i2e(1),ism(1),irelo(1),iacto(8)                       9d1s06
      dimension phs(2)
      data phs/1d0,-1d0/
      data icall/0/
      save icall
      do i=1,n1
       i1=idat1(1,i)
       i2=idat1(2,i)
       i1o=i1+ioffdet                                                   12d20s06
       i2o=i2+ioffdet                                                   12d20s06
       idwstmp2=idat1(3,i)                                              5d30s06
       isbx=ism(idwstmp1(1))                                            8d29s06
       isby=ism(idwstmp1(2))                                            8d29s06
       into=irelo(idwstmp1(1))                                          8d29s06
       infrom=irelo(idwstmp1(2))                                        8d29s06
       iad=ih0(isby)+into-1+iacto(isby)*(infrom-1)                      8d29s06
       sum=bc(iad)                                                      8d29s06
       sumh=sum                                                         12d15s06
       ix=max(into,infrom)                                              9d1s06
       in=min(into,infrom)                                              9d1s06
       ii1=((ix*(ix-1))/2)+in                                           9d1s06
       na=(iacto(isbx)*(iacto(isbx)+1))/2                               9d1s06
       do j=1,nalpha                                                    9d1s06
        isx=ism(iaorb(j,i1o))                                           12d20s06
        i2eu=ifind2(isbx,isbx,isx,isx,icase)                            9d1s06
        iisx=irelo(iaorb(j,i1o))                                        12d20s06
        ii2=(iisx*(iisx+1))/2                                           9d1s06
        iad1=i2e(i2eu)+ii1-1+na*(ii2-1)                                 9d1s06
        sum=sum+bc(iad1)                                                9d1s06
        if(isx.eq.isbx)then                                             9d1s06
         i2eu=ifind2(isx,isx,isx,isx,icase)                             9d1s06
         in=min(into,iisx)                                              9d1s06
         ix=max(into,iisx)                                              9d1s06
         ii3=((ix*(ix-1))/2)+in-1                                       9d1s06
         in=min(infrom,iisx)                                            9d1s06
         ix=max(infrom,iisx)                                            9d1s06
         ii2=((ix*(ix-1))/2)+in-1                                       9d1s06
         iad1=i2e(i2eu)+ii3+na*ii2                                      9d1s06
        else                                                            9d1s06
         i2eu=ifind2(isbx,isx,isbx,isx,icase)                           9d1s06
         if(icase.eq.1)then                                             9d1s06
          iad1=i2e(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*         9d1s06
     $        (infrom-1+iacto(isbx)*(iisx-1)))                          9d1s06
         else if(icase.eq.2)then                                        9d1s06
          iad1=i2e(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*         9d1s06
     $        (iisx-1+iacto(isx)*(infrom-1)))                           9d1s06
         else if(icase.eq.3)then                                        9d1s06
          iad1=i2e(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*         9d1s06
     $        (iisx-1+iacto(isx)*(infrom-1)))                           9d1s06
         else if(icase.eq.4)then                                        9d1s06
          iad1=i2e(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*         9d1s06
     $        (infrom-1+iacto(isbx)*(iisx-1)))                          9d1s06
         end if                                                         9d1s06
        end if                                                          9d1s06
        sum=sum-bc(iad1)                                                9d1s06
       end do                                                           9d1s06
       do k=1,nhere                                                     3d14s17
        ka=k+kadd                                                       3d14s17
        sum1=sum                                                        3d14s17
        do l=1,nbeta                                                    3d14s17
         isx=ism(iborb(l,ka))                                           9d1s06
         i2eu=ifind2(isbx,isbx,isx,isx,icase)                           9d1s06
         iisx=irelo(iborb(l,ka))                                        9d1s06
         ii2=(iisx*(iisx+1))/2                                          9d1s06
         iad1=i2e(i2eu)+ii1-1+na*(ii2-1)                                9d1s06
         sum1=sum1+bc(iad1)                                             9d1s06
        end do                                                          9d1s06
        tmp=sum1*phs(idat1(4,i))                                        5d30s06
        do j=1,nroot                                                    3d14s17
         g(j,k,i1)=g(j,k,i1)+tmp*cvec(j,k,i2)                            3d14s17
         g(j,k,i2)=g(j,k,i2)+tmp*cvec(j,k,i1)                            3d14s17
        end do                                                          3d14s17
       end do                                                           3d14s17
      end do
      return
      end
