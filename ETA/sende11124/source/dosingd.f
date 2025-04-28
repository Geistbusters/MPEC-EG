c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dosingd(idat1,n1,cvecb,cveck,numb,nroot,nhere,         5d27s22
     $     iaorb,nalpha,iborb,nbeta,kadd,ism,irelo,iacto,ioffdet,iden1, 3d27s17
     $     jdenpt,wgt,l2e,bc,ibc,mdenoff,igoal)                         3d16s23
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2,ism,irelo                                   12d13s06
      integer*2 idat1(4,n1),idwstmp2                                    5d30s06
      integer*1 iaorb(nalpha,1),iborb(nbeta,1),idwstmp1(2)              5d30s06
      logical l2e                                                       6d2s22
      equivalence (idwstmp2,idwstmp1)                                   5d30s06
      include "common.store"                                            10d6s06
      dimension cvecb(nroot,nhere,numb),cveck(nroot,nhere,numb),        5d27s22
     $     wgt(nroot),igoalx(2),                                                  5d27s22
     $     iden1(8),jdenpt(*),ism(1),irelo(1),iacto(8)                       9d1s06
      dimension phs(2)
      data phs/1d0,-1d0/
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
       if(l2e)then                                                      3d15s23
        iad=iden1(isby)+into-1+iacto(isby)*(infrom-1)                    3d27s17
        iadt=iden1(isbx)+infrom-1+iacto(isbx)*(into-1)                   7d22s22
        ww=phs(idat1(4,i))                                               6d15s22
        do k=1,nhere
         do j=1,nroot                                                     3d27s17
          bc(iad)=bc(iad)+ww*wgt(j)*cveck(j,k,i1)*cvecb(j,k,i2)          7d22s22
          bc(iadt)=bc(iadt)+ww*wgt(j)*cvecb(j,k,i1)*cveck(j,k,i2)        7d22s22
         end do                                                          3d27s17
        end do                                                           3d27s17
       else                                                             3d15s23
        iad=iden1(isby)+into-1+iacto(isby)*(infrom-1                    3d15s23
     $       +iacto(isby)*mdenoff)                                      3d15s23
        iadt=iden1(isbx)+infrom-1+iacto(isbx)*(into-1                   3d15s23
     $       +iacto(isbx)*mdenoff)                                      3d15s23
        ww=phs(idat1(4,i))                                               6d15s22
        nad=iacto(isby)*iacto(isbx)                                     3d15s23
        do k1=1,nroot                                                   3d15s23
         do k2=1,k1                                                     3d15s23
          do k=1,nhere                                                  3d15s23
           bc(iad)=bc(iad)+ww*cveck(k1,k,i1)*cvecb(k2,k,i2)             3d15s23
           bc(iadt)=bc(iadt)+ww*cveck(k1,k,i2)*cvecb(k2,k,i1)           3d15s23
          end do                                                        3d15s23
          iad=iad+nad                                                   3d15s23
          iadt=iadt+nad                                                 3d15s23
         end do                                                         3d15s23
        end do                                                          3d15s23
       end if                                                           3d15s23
       if(l2e)then                                                      6d2s22
        ix=max(into,infrom)                                              9d1s06
        in=min(into,infrom)                                              9d1s06
        ii1=((ix*(ix-1))/2)+in                                           9d1s06
        na=(iacto(isbx)*(iacto(isbx)+1))/2                               9d1s06
        do j=1,nalpha                                                    9d1s06
         isx=ism(iaorb(j,i1o))                                           12d20s06
         i2eu=ifind2(isbx,isbx,isx,isx,icase)                            9d1s06
         iisx=irelo(iaorb(j,i1o))                                        12d20s06
         ii2=(iisx*(iisx+1))/2                                           9d1s06
         iad0=jdenpt(i2eu)+ii1-1+na*(ii2-1)                              3d27s17
         if(isx.eq.isbx)then                                             9d1s06
          i2eu=ifind2(isx,isx,isx,isx,icase)                             9d1s06
          in=min(into,iisx)                                              9d1s06
          ix=max(into,iisx)                                              9d1s06
          ii3=((ix*(ix-1))/2)+in-1                                       9d1s06
          in=min(infrom,iisx)                                            9d1s06
          ix=max(infrom,iisx)                                            9d1s06
          ii2=((ix*(ix-1))/2)+in-1                                       9d1s06
          iad1=jdenpt(i2eu)+ii3+na*ii2                                   3d27s17
         else                                                            9d1s06
          i2eu=ifind2(isbx,isx,isbx,isx,icase)                           9d1s06
          if(icase.eq.1)then                                             9d1s06
           iad1=jdenpt(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      3d27s17
     $        (infrom-1+iacto(isbx)*(iisx-1)))                          9d1s06
          else if(icase.eq.2)then                                        9d1s06
           iad1=jdenpt(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      3d27s17
     $        (iisx-1+iacto(isx)*(infrom-1)))                           9d1s06
          else if(icase.eq.3)then                                        9d1s06
           iad1=jdenpt(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      3d27s17
     $        (iisx-1+iacto(isx)*(infrom-1)))                           9d1s06
          else if(icase.eq.4)then                                        9d1s06
           iad1=jdenpt(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      3d27s17
     $        (infrom-1+iacto(isbx)*(iisx-1)))                          9d1s06
          end if                                                         9d1s06
         end if                                                          9d1s06
         do k=1,nhere                                                    3d27s17
          do l=1,nroot                                                   3d27s17
           term=ww*wgt(l)*(cvecb(l,k,i1)*cveck(l,k,i2)                  6d16s22
     $          +cveck(l,k,i1)*cvecb(l,k,i2))                           6d16s22
           bc(iad0)=bc(iad0)+term                                        3d27s17
           bc(iad1)=bc(iad1)-term                                        3d27s17
          end do
         end do                                                          3d27s17
        end do                                                           9d1s06
        do k=1,nhere                                                     3d14s17
         ka=k+kadd                                                       3d14s17
         do l=1,nbeta                                                    3d14s17
          isx=ism(iborb(l,ka))                                           9d1s06
          i2eu=ifind2(isbx,isbx,isx,isx,icase)                           9d1s06
          iisx=irelo(iborb(l,ka))                                        9d1s06
          ii2=(iisx*(iisx+1))/2                                          9d1s06
          iad1=jdenpt(i2eu)+ii1-1+na*(ii2-1)                             3d27s17
          do m=1,nroot                                                   3d27s17
           term=ww*wgt(m)*(cvecb(m,k,i1)*cveck(m,k,i2)                  6d17s22
     $          +cveck(m,k,i1)*cvecb(m,k,i2))                           6d16s22
           bc(iad1)=bc(iad1)+term                                       6d17s22
          end do                                                         3d27s17
         end do                                                          9d1s06
        end do                                                           3d14s17
       end if                                                           6d2s22
      end do
      return
      end
