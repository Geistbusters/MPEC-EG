c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine doab1d(idat1b,nb,idat1a,na,vecab,vecak,numa,numb,nroot,5d27s22
     $     nhere,il,ih,irelo,ism,iacto,jdenpt,wgt,bc,ibc)               11d14s22
      implicit real*8 (a-h,o-z)                                         3d10s17
      integer*2 idat1b(4,nb),idat1a(4,na),idwstmp2                      3d12s17
      integer*1 idwstmp1(4)                                             3d10s17
      equivalence (idwstmp2,idwstmp1)                                   3d10s17
      integer*8 ism,irelo                                               3d10s17
      include "common.store"                                            10d6s06
      dimension vecab(nroot,numa,numb),vecak(nroot,numa,numb),          5d27s22
     $     wgt(nroot),jdenpt(*),iacto(8),ism(1),irelo(1),phs(2)         5d27s22
      dimension igoalx(2)
      data phs/1d0,-1d0/                                                3d10s17
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      do ib=1,nb                                                        3d10s17
       i1b=idat1b(1,ib)                                                 3d10s17
       i2b=idat1b(2,ib)                                                 3d10s17
       idwstmp2=idat1b(3,ib)                                             3d10s17
       ntob=irelo(idwstmp1(1))                                          3d10s17
       istob=ism(idwstmp1(1))                                           3d10s17
       nfromb=irelo(idwstmp1(2))                                        3d10s17
       isfmb=ism(idwstmp1(2))                                           3d10s17
       pb=phs(idat1b(4,ib))                                              3d10s17
       do ia=1,na                                                       3d10s17
        i1a=idat1a(1,ia)                                                3d10s17
        i2a=idat1a(2,ia)                                                3d10s17
        if((i1a.ge.il.and.i1a.le.ih).or.(i2a.ge.il.and.i2a.le.ih))then  3d10s17
         idwstmp2=idat1a(3,ia)                                           3d10s17
         ntoa=irelo(idwstmp1(1))                                        3d10s17
         istoa=ism(idwstmp1(1))                                         3d10s17
         nfroma=irelo(idwstmp1(2))                                      3d10s17
         isfma=ism(idwstmp1(2))                                         3d10s17
         pa=phs(idat1a(4,ia))*pb                                         3d10s17
         i2eu=inv(1,istob,isfmb,istoa)                                  4d26s18
         icase=inv(2,istob,isfmb,istoa)                                 4d26s18
         if(istob.eq.isfmb)then                                         3d10s17
          in=min(nfromb,ntob)                                           3d15s17
          ix=max(nfromb,ntob)                                           3d15s17
          ii1=((ix*(ix-1))/2)+in-1                                      3d10s17
          in=min(nfroma,ntoa)                                           3d15s17
          ix=max(nfroma,ntoa)                                           3d15s17
          ii2=((ix*(ix-1))/2)+in-1                                      3d10s17
          nrow=(iacto(istob)*(iacto(istob)+1))/2                          3d10s17
          iad1=jdenpt(i2eu)+ii1+nrow*ii2                                3d27s17
         else                                                           3d10s17
          if(icase.eq.1)then                                            3d10s17
           iad1=jdenpt(i2eu)+ntob-1+iacto(istob)*(nfromb-1+iacto(isfmb) 3d27s17
     $        *(ntoa-1+iacto(istoa)*(nfroma-1)))                        3d15s17
          else if(icase.eq.2)then                                       3d10s17
           iad1=jdenpt(i2eu)+nfromb-1+iacto(isfmb)*(ntob-1+iacto(istob) 3d27s17
     $        *(nfroma-1+iacto(isfma)*(ntoa-1)))                        3d15s17
          else if(icase.eq.3)then                                       3d10s17
           iad1=jdenpt(i2eu)+ntob-1+iacto(istob)*(nfromb-1+iacto(isfmb) 3d27s17
     $        *(nfroma-1+iacto(isfma)*(ntoa-1)))                        3d15s17
          else if(icase.eq.4)then                                       3d10s17
           iad1=jdenpt(i2eu)+nfromb-1+iacto(isfmb)*(ntob-1+iacto(istob) 3d27s17
     $        *(ntoa-1+iacto(istoa)*(nfroma-1)))                        3d15s17
          end if                                                        3d10s17
         end if                                                         3d10s17
         if(i1a.ge.il.and.i1a.le.ih)then                                3d10s17
          i1ao=i1a+1-il                                                 3d10s17
          do k=1,nroot                                                  3d10s17
           term=pa*wgt(k)*(vecab(k,i2a,i2b)                             6d17s22
     $          *vecak(k,i1a,i1b)                                       5d27s22
     $          +vecab(k,i2a,i1b)*vecak(k,i1a,i2b))                     5d27s22
           bc(iad1)=bc(iad1)+term                                       6d17s22
          end do                                                        3d10s17
         end if                                                         3d10s17
         if(i2a.ge.il.and.i2a.le.ih)then                                3d10s17
          i2ao=i2a+1-il                                                 3d10s17
          do k=1,nroot                                                  3d10s17
           term=pa*wgt(k)*(                                             6d17s22
     $          vecab(k,i1a,i2b)*vecak(k,i2a,i1b)                       5d27s22
     $          +vecab(k,i1a,i1b)*vecak(k,i2a,i2b))                     5d27s22
           bc(iad1)=bc(iad1)+term                                       6d17s22
          end do                                                        3d10s17
         end if                                                         3d10s17
        end if
       end do                                                            3d10s17
      end do
      return                                                            3d10s17
      end                                                               3d10s17
