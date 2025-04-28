c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine doab2dbl(idatac,mca,idatbc,mcb,nadetb,nadetk,nbdetb,      3d16s17
     $     nbdetk,nhereb,ilb,ihb,nherek,ilk,ihk,vecab,vecak,            7d22s22
     $     irelo,ism,iacto,iu1,iu2,jdenpt,bc,ibc,nroot,mdenoff,igoal,ff)3d27s23
      implicit real*8 (a-h,o-z)                                         3d16s17
      integer*2 idatac(4,mca),idatbc(4,mcb),idwstmp2                    3d16s17
      integer*1 idwstmp1(2)                                             7d22s22
      equivalence (idwstmp2,idwstmp1)                                   3d16s17
      integer*8 ism(*),irelo(*)                                         3d16s17
      include "common.store"                                            3d16s17
      dimension iacto(8),phs(2),vecab(nroot,nadetb,nbdetb),             3d15s23
     $     vecak(nroot,nadetk,nbdetk),jdenpt(*),igoalx(2)               3d15s23
      data phs/1d0,-1d0/                                                3d16s17
      common/fnd2cm/inv(2,8,8,8)                                        4d26s18
      do ib=1,mcb
       i1b=idatbc(iu1,ib)
       i2b=idatbc(iu2,ib)
       idwstmp2=idatbc(3,ib)
       nbb=irelo(idwstmp1(iu1))-1                                       7d22s22
       ibb=ism(idwstmp1(iu1))                                           7d22s22
       nbk=irelo(idwstmp1(iu2))-1                                       7d22s22
       ibk=ism(idwstmp1(iu2))                                           7d22s22
       isv1=idwstmp1(iu2)
       isv2=idwstmp1(iu1)
       pb=phs(idatbc(4,ib))                                             3d27s23
       do ia=1,mca
        i1a=idatac(1,ia)
        i2a=idatac(2,ia)
        idwstmp2=idatac(3,ia)
        nab=irelo(idwstmp1(1))-1                                        7d22s22
        iab=ism(idwstmp1(1))                                            7d22s22
        nak=irelo(idwstmp1(2))-1                                        7d22s22
        iak=ism(idwstmp1(2))                                            7d22s22
        nad=iacto(iab)*iacto(iak)*iacto(ibb)*iacto(ibk)                 3d15s23
        if(i1a.ge.ilb.and.i1a.le.ihb)then                               7d22s22
         pa=pb*phs(idatac(4,ia))
         i2eu=ifind2(iab,iak,ibb,ibk,icase)                              3d16s17
         if(icase.eq.1)then                                             7d8s22
          iad=jdenpt(i2eu)+nab+iacto(iab)*(nak+iacto(iak)               7d22s22
     $         *(nbb+iacto(ibb)*nbk))                                   7d22s22
          wu=pa                                                         7d22s22
         else if(icase.eq.4)then                                        7d8s22
          iad=jdenpt(i2eu)+nak+iacto(iak)*(nab+iacto(iab)               7d22s22
     $         *(nbb+iacto(ibb)*nbk))                                   7d22s22
          wu=-pa                                                        7d22s22
         else if(icase.eq.3)then                                        7d8s22
          iad=jdenpt(i2eu)+nab+iacto(iab)*(nak+iacto(iak)               7d22s22
     $         *(nbk+iacto(ibk)*nbb))                                   7d22s22
          wu=-pa                                                        7d22s22
         else                                                           7d8s22
          iad=jdenpt(i2eu)+nak+iacto(iak)*(nab+iacto(iab)               7d22s22
     $         *(nbk+iacto(ibk)*nbb))                                   7d22s22
          wu=pa                                                         7d22s22
         end if                                                         7d8s22
         iad=iad+nad*mdenoff                                            3d15s23
         do k1=1,nroot                                                  3d15s23
          do k2=1,k1-1                                                  3d27s23
           bc(iad)=bc(iad)+wu*(vecak(k1,i2a,i2b)*vecab(k2,i1a,i1b)      3d30s23
     $                        +vecak(k2,i2a,i2b)*vecab(k1,i1a,i1b))
     $            *2d0                                                   3d30s23
           iad=iad+nad                                                  3d15s23
          end do                                                        3d15s23
          bc(iad)=bc(iad)+wu*vecak(k1,i2a,i2b)*vecab(k2,i1a,i1b)*4d0    3d27s23
          iad=iad+nad                                                   3d27s23
         end do                                                         3d15s23
        end if                                                          3d16s17
       end do
      end do
      return
      end
