c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dosingdbl(idat1,n1,cvec,numb,nhere,iaorb,nalpha,iborb, 7d21s22
     $   nbeta,kadd,ism,irelo,iacto,ioffdet,jdenpt,bc,ibc,nroot,mdenoff,3d15s23
     $     iden1,igoal)                                                 3d15s23
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2,ism,irelo                                   12d13s06
      integer*2 idat1(4,n1),idwstmp2                                    5d30s06
      integer*1 iaorb(nalpha,1),iborb(nbeta,1),idwstmp1(2)              5d30s06
      equivalence (idwstmp2,idwstmp1)                                   5d30s06
      include "common.store"                                            10d6s06
      dimension cvec(nroot,nhere,numb),igoalx(2),jdenpt(*),ism(*),      3d15s23
     $     irelo(*),iacto(8),iden1(*)                                   3d15s23
      dimension phs(2)
      data phs/1d0,-1d0/
      do i=1,n1
       i1=idat1(1,i)
       i2=idat1(2,i)
       i1o=i1+ioffdet                                                   12d20s06
       i2o=i2+ioffdet                                                   12d20s06
       idwstmp2=idat1(3,i)                                              5d30s06
       isbx=ism(idwstmp1(2))                                            8d29s06
       isby=ism(idwstmp1(1))                                            8d29s06
       inx=irelo(idwstmp1(2))-1                                         7d21s22
       iny=irelo(idwstmp1(1))-1                                         7d21s22
       iad=iden1(isbx)+inx+iacto(isbx)*(iny+iacto(isbx)*mdenoff)        3d15s23
       iadt=iden1(isby)+iny+iacto(isby)*(inx+iacto(isbx)*mdenoff)        3d15s23
       nad=iacto(isbx)*iacto(isbx)                                      3d15s23
       do k1=1,nroot                                                    3d15s23
        do k2=1,k1                                                      3d15s23
         do k=1,nhere                                                   3d15s23
          bc(iad)=bc(iad)+cvec(k1,k,i1)*cvec(k2,k,i2)*phs(idat1(4,i))   3d15s23
          bc(iadt)=bc(iadt)+cvec(k1,k,i2)*cvec(k2,k,i1)*phs(idat1(4,i)) 3d15s23
         end do                                                         3d15s23
         iad=iad+nad                                                    3d15s23
         iadt=iadt+nad                                                  3d15s23
        end do                                                          3d15s23
       end do                                                           3d15s23
c
c     -2 sum c cbkc
c
       ww=-2d0*phs(idat1(4,i))                                          3d28s23
       jsbx=isby                                                        3d27s23
       jsby=isbx                                                        3d27s23
       jnx=iny                                                          3d27s23
       jny=inx                                                          3d27s23
       do j=1,nalpha                                                    9d1s06
        if(iaorb(j,i1o).ne.idwstmp1(1).and.                             7d21s22
     $     iaorb(j,i1o).ne.idwstmp1(2))then                             7d21s22
         is14=ism(iaorb(j,i1o))                                         7d21s22
         if14=irelo(iaorb(j,i1o))-1                                     7d21s22
         nad=iacto(is14)*iacto(isbx)*iacto(isby)*iacto(is14)            3d14s23
         i2eu=ifind2(is14,isbx,isby,is14,icase)                         7d21s22
         if(icase.eq.1)then                                             7d21s22
          iad=jdenpt(i2eu)+if14+iacto(is14)*(inx+iacto(isbx)            7d21s22
     $           *(iny+iacto(isby)*if14))                               7d21s22
          wu=ww                                                         7d21s22
         else if(icase.eq.4)then                                        7d21s22
          iad=jdenpt(i2eu)+inx+iacto(isbx)*(if14+iacto(is14)            7d21s22
     $           *(iny+iacto(isby)*if14))                               7d21s22
          wu=-ww                                                        7d21s22
         else if(icase.eq.3)then                                        7d21s22
          iad=jdenpt(i2eu)+if14+iacto(is14)*(inx+iacto(isbx)            7d21s22
     $           *(if14+iacto(is14)*iny))                               7d21s22
          wu=-ww                                                        7d21s22
         else                                                           7d21s22
          iad=jdenpt(i2eu)+inx+iacto(isbx)*(if14+iacto(is14)            7d21s22
     $           *(if14+iacto(is14)*iny))                               7d21s22
          wu=ww                                                         7d21s22
         end if                                                         7d21s22
         iad=iad+nad*mdenoff                                            3d14s23
         do k1=1,nroot                                                  3d14s23
          do k2=1,k1                                                    3d14s23
           do k=1,nhere                                                    3d27s17
            term=wu*cvec(k1,k,i1)*cvec(k2,k,i2)                                 7d21s22
            bc(iad)=bc(iad)+term                                          7d21s22
           end do                                                       3d14s23
           iad=iad+nad                                                  3d14s23
          end do                                                        3d14s23
         end do                                                          3d27s17
         is14=ism(iaorb(j,i1o))                                         7d21s22
         if14=irelo(iaorb(j,i1o))-1                                     7d21s22
         nad=iacto(is14)*iacto(isbx)*iacto(isby)*iacto(is14)            3d14s23
         i2eu=ifind2(is14,jsbx,jsby,is14,icase)                         7d21s22
         if(icase.eq.1)then                                             7d21s22
          iad=jdenpt(i2eu)+if14+iacto(is14)*(jnx+iacto(jsbx)            7d21s22
     $           *(jny+iacto(jsby)*if14))                               7d21s22
          wu=ww                                                         7d21s22
         else if(icase.eq.4)then                                        7d21s22
          iad=jdenpt(i2eu)+jnx+iacto(jsbx)*(if14+iacto(is14)            7d21s22
     $           *(jny+iacto(jsby)*if14))                               7d21s22
          wu=-ww                                                        7d21s22
         else if(icase.eq.3)then                                        7d21s22
          iad=jdenpt(i2eu)+if14+iacto(is14)*(jnx+iacto(jsbx)            7d21s22
     $           *(if14+iacto(is14)*jny))                               7d21s22
          wu=-ww                                                        7d21s22
         else                                                           7d21s22
          iad=jdenpt(i2eu)+jnx+iacto(jsbx)*(if14+iacto(is14)            7d21s22
     $           *(if14+iacto(is14)*jny))                               7d21s22
          wu=ww                                                         7d21s22
         end if                                                         7d21s22
         iad=iad+nad*mdenoff                                            3d14s23
         do k1=1,nroot                                                  3d14s23
          do k2=1,k1                                                    3d14s23
           do k=1,nhere                                                    3d27s17
            term=wu*cvec(k1,k,i2)*cvec(k2,k,i1)                         3d27s23
            bc(iad)=bc(iad)+term                                          7d21s22
           end do                                                       3d14s23
           iad=iad+nad                                                  3d14s23
          end do                                                        3d14s23
         end do                                                          3d27s17
        end if                                                          7d21s22
       end do                                                           3d14s17
      end do
      return
      end
