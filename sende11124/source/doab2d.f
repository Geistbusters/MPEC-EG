c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine doab2d(idatac,mca,idatbc,mcb,nadetb,nadetk,nbdetb,      3d16s17
     $     nbdetk,nroot,nhereb,ilb,ihb,nherek,ilk,ihk,vecab,vecak,      3d27s17
     $     vecbb,vecbk,irelo,ism,iacto,iu1,iu2,jdenpt,wgt,bc,ibc)       11d14s22
      implicit real*8 (a-h,o-z)                                         3d16s17
      integer*2 idatac(4,mca),idatbc(4,mcb),idwstmp2                    3d16s17
      integer*1 idwstmp1(4)                                             3d16s17
      equivalence (idwstmp2,idwstmp1)                                   3d16s17
      integer*8 ism(*),irelo(*)                                         3d16s17
      include "common.store"                                            3d16s17
      dimension iacto(8),phs(2),vecab(nroot,nadetb,nbdetb),                  3d16s17
     $     vecak(nroot,nadetk,nbdetk),jdenpt(*),wgt(nroot),igoalx(2)              3d27s17
      data phs/1d0,-1d0/                                                3d16s17
      common/fnd2cm/inv(2,8,8,8)                                        4d26s18
      do ib=1,mcb
       i1b=idatbc(iu1,ib)
       i2b=idatbc(iu2,ib)
       idwstmp2=idatbc(3,ib)
       nbb=irelo(idwstmp1(iu1))
       ibb=ism(idwstmp1(iu1))
       nbk=irelo(idwstmp1(iu2))
       ibk=ism(idwstmp1(iu2))
       pb=phs(idatbc(4,ib))
       do ia=1,mca
        i1a=idatac(1,ia)
        i2a=idatac(2,ia)
        idwstmp2=idatac(3,ia)
        nab=irelo(idwstmp1(1))
        iab=ism(idwstmp1(1))
        nak=irelo(idwstmp1(2))
        iak=ism(idwstmp1(2))
        if((i1a.ge.ilb.and.i1a.le.ihb).or.(i2a.ge.ilk.and.i2a.le.ihk))
     $       then                                                       3d16s17
         pa=pb*phs(idatac(4,ia))
         i2eu=inv(1,iab,iak,ibb)                                        4d26s18
         icase=inv(2,iab,iak,ibb)                                       4d26s18
         if(iab.eq.iak)then                                             3d16s17
          ix=max(nab,nak)
          in=min(nab,nak)
          ii1=((ix*(ix-1))/2)+in-1                                      3d16s17
          ix=max(nbb,nbk)
          in=min(nbb,nbk)
          ii2=((ix*(ix-1))/2)+in-1                                      3d16s17
          nrow=(iacto(iab)*(iacto(iab)+1))/2
          iad1=jdenpt(i2eu)+ii1+ii2*nrow                                3d27s17
         else                                                           3d16s17
          if(icase.eq.1)then
           iad1=jdenpt(i2eu)+nab-1+iacto(iab)*(nak-1+iacto(iak)*(nbb-1  3d27s17
     $          +iacto(ibb)*(nbk-1)))                                    3d16s17
          else if(icase.eq.2)then
           iad1=jdenpt(i2eu)+nak-1+iacto(iak)*(nab-1+iacto(iab)*(nbk-1  3d27s17
     $          +iacto(ibk)*(nbb-1)))                                    3d16s17
          else if(icase.eq.3)then
           iad1=jdenpt(i2eu)+nab-1+iacto(iab)*(nak-1+iacto(iak)*(nbk-1  3d27s17
     $          +iacto(ibk)*(nbb-1)))                                    3d16s17
          else
           iad1=jdenpt(i2eu)+nak-1+iacto(iak)*(nab-1+iacto(iab)*(nbb-1  3d27s17
     $          +iacto(ibb)*(nbk-1)))                                    3d16s17
          end if
         end if
         if(iad1.lt.1)then
          write(6,*)('iad1 is lt 1!!! '),iad1,ia,ib,i2eu,icase,iab,iak,
     $         ibb,ibk,nab,nak,nbb,nbk,jdenpt(i2eu)
         end if
         if(i1a.ge.ilb.and.i1a.le.ihb)then
          i1ao=i1a-ilb+1
          do i=1,nroot                                                  3d16s17
           term=pa*vecak(i,i2a,i2b)*vecab(i,i1a,i1b)*wgt(i)             6d17s22
           bc(iad1)=bc(iad1)+term                                       6d17s22
          end do                                                        3d16s17
         end if
         if(i2a.ge.ilk.and.i2a.le.ihk)then
          i2ao=i2a-ilk+1
          do i=1,nroot                                                  3d16s17
           term=pa*vecab(i,i1a,i1b)*vecak(i,i2a,i2b)*wgt(i)             6d17s22
           bc(iad1)=bc(iad1)+term                                       6d17s22
          end do                                                        3d16s17
         end if
        end if                                                          3d16s17
       end do
      end do
      return
      end
