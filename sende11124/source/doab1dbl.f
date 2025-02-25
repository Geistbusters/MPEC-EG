c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine doab1dbl(idat1b,nb,idat1a,na,veca,numa,numb,nhere,il,  7d21s22
     $     ih,irelo,ism,iacto,jdenpt,bc,ibc,nroot,mdenoff,igoal)        3d15s23
      implicit real*8 (a-h,o-z)                                         3d10s17
      integer*2 idat1b(4,nb),idat1a(4,na),idwstmp2                      3d12s17
      integer*1 idwstmp1(4)                                             3d10s17
      equivalence (idwstmp2,idwstmp1)                                   3d10s17
      integer*8 ism,irelo                                               3d10s17
      include "common.store"                                            10d6s06
      dimension veca(nroot,numa,numb),jdenpt(*),iacto(*),ism(*),        3d15s23
     $     irelo(*),phs(2)                                              3d15s23
      dimension igoalx(2)
      data phs/1d0,-1d0/                                                3d10s17
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      do ib=1,nb                                                        3d10s17
       i1b=idat1b(1,ib)                                                 3d10s17
       i2b=idat1b(2,ib)                                                 3d10s17
       idwstmp2=idat1b(3,ib)                                             3d10s17
       ntob=irelo(idwstmp1(2))-1                                        7d21s22
       istob=ism(idwstmp1(2))                                           3d10s17
       nfromb=irelo(idwstmp1(1))-1                                      7d21s22
       isfmb=ism(idwstmp1(1))                                           3d10s17
       pb=phs(idat1b(4,ib))                                              3d10s17
       do ia=1,na                                                       3d10s17
        i1a=idat1a(1,ia)                                                3d10s17
        i2a=idat1a(2,ia)                                                3d10s17
        if(i1a.ge.il.and.i1a.le.ih)then                                 7d21s22
         idwstmp2=idat1a(3,ia)                                           3d10s17
         ntoa=irelo(idwstmp1(2))-1                                      7d21s22
         istoa=ism(idwstmp1(2))                                         3d10s17
         nfroma=irelo(idwstmp1(1))-1                                    7d21s22
         isfma=ism(idwstmp1(1))                                         3d10s17
         pa=phs(idat1a(4,ia))*pb                                         3d10s17
         i2eu=ifind2(istob,isfmb,istoa,isfma,icase)                     3d10s17
         nad=iacto(istob)*iacto(isfmb)*iacto(istoa)*iacto(isfma)        3d15s23
         if(icase.eq.1)then                                             7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=pa                                                         7d21s22
         else if(icase.eq.4)then                                        7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=-pa                                                         7d21s22
         else if(icase.eq.3)then                                        7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=-pa                                                         7d21s22
         else                                                           7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=pa                                                         7d21s22
         end if                                                         7d8s22
         iad=iad+nad*mdenoff                                            3d15s23
         do k1=1,nroot                                                  3d15s23
          do k2=1,k1                                                    3d15s23
           term=2d0*wu*veca(k2,i2a,i2b)*veca(k1,i1a,i1b)                3d27s23
           bc(iad)=bc(iad)+term                                           7d21s22
           iad=iad+nad                                                  3d15s23
          end do                                                        3d15s23
         end do                                                         3d15s23
        end if
        if(i2a.ge.il.and.i2a.le.ih)then                                 7d21s22
         idwstmp2=idat1a(3,ia)                                           3d10s17
         ntoa=irelo(idwstmp1(1))-1                                      7d21s22
         istoa=ism(idwstmp1(1))                                         3d10s17
         nfroma=irelo(idwstmp1(2))-1                                    7d21s22
         isfma=ism(idwstmp1(2))                                         3d10s17
         pa=phs(idat1a(4,ia))*pb                                         3d10s17
         i2eu=ifind2(istob,isfmb,istoa,isfma,icase)                     3d10s17
         nad=iacto(istob)*iacto(isfmb)*iacto(istoa)*iacto(isfma)        3d15s23
         if(icase.eq.1)then                                             7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=pa                                                         7d21s22
         else if(icase.eq.4)then                                        7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=-pa                                                         7d21s22
         else if(icase.eq.3)then                                        7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=-pa                                                         7d21s22
         else                                                           7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=pa                                                         7d21s22
         end if                                                         7d8s22
         iad=iad+nad*mdenoff                                            3d15s23
         do k1=1,nroot                                                  3d15s23
          do k2=1,k1                                                    3d15s23
           term=2d0*wu*veca(k2,i1a,i2b)*veca(k1,i2a,i1b)                3d27s23
           bc(iad)=bc(iad)+term                                           7d21s22
           iad=iad+nad                                                  3d15s23
          end do                                                        3d15s23
         end do                                                         3d15s23
        end if
       end do                                                           3d27s23
       i1b=idat1b(2,ib)                                                 3d27s23
       i2b=idat1b(1,ib)                                                 3d10s17
       idwstmp2=idat1b(3,ib)                                            3d27s23
       ntob=irelo(idwstmp1(1))-1                                        3d27s23
       istob=ism(idwstmp1(1))                                           3d10s17
       nfromb=irelo(idwstmp1(2))-1                                      3d27s23
       isfmb=ism(idwstmp1(2))                                           3d27s23
       pb=phs(idat1b(4,ib))                                              3d10s17
       do ia=1,na                                                       3d10s17
        i1a=idat1a(1,ia)                                                3d10s17
        i2a=idat1a(2,ia)                                                3d10s17
        if(i1a.ge.il.and.i1a.le.ih)then                                 7d21s22
         idwstmp2=idat1a(3,ia)                                           3d10s17
         ntoa=irelo(idwstmp1(2))-1                                      7d21s22
         istoa=ism(idwstmp1(2))                                         3d10s17
         nfroma=irelo(idwstmp1(1))-1                                    7d21s22
         isfma=ism(idwstmp1(1))                                         3d10s17
         pa=phs(idat1a(4,ia))*pb                                         3d10s17
         i2eu=ifind2(istob,isfmb,istoa,isfma,icase)                     3d10s17
         nad=iacto(istob)*iacto(isfmb)*iacto(istoa)*iacto(isfma)        3d15s23
         if(icase.eq.1)then                                             7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=pa                                                         7d21s22
         else if(icase.eq.4)then                                        7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=-pa                                                         7d21s22
         else if(icase.eq.3)then                                        7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=-pa                                                         7d21s22
         else                                                           7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=pa                                                         7d21s22
         end if                                                         7d8s22
         iad=iad+nad*mdenoff                                            3d15s23
         do k1=1,nroot                                                  3d15s23
          do k2=1,k1                                                    3d15s23
           term=2d0*wu*veca(k2,i2a,i2b)*veca(k1,i1a,i1b)                3d27s23
           bc(iad)=bc(iad)+term                                           7d21s22
           iad=iad+nad                                                  3d15s23
          end do                                                        3d15s23
         end do                                                         3d15s23
        end if
        if(i2a.ge.il.and.i2a.le.ih)then                                 7d21s22
         idwstmp2=idat1a(3,ia)                                           3d10s17
         ntoa=irelo(idwstmp1(1))-1                                      7d21s22
         istoa=ism(idwstmp1(1))                                         3d10s17
         nfroma=irelo(idwstmp1(2))-1                                    7d21s22
         isfma=ism(idwstmp1(2))                                         3d10s17
         pa=phs(idat1a(4,ia))*pb                                         3d10s17
         i2eu=ifind2(istob,isfmb,istoa,isfma,icase)                     3d10s17
         nad=iacto(istob)*iacto(isfmb)*iacto(istoa)*iacto(isfma)        3d15s23
         if(icase.eq.1)then                                             7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=pa                                                         7d21s22
         else if(icase.eq.4)then                                        7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(ntoa+iacto(istoa)*nfroma))                              7d8s22
          wu=-pa                                                         7d21s22
         else if(icase.eq.3)then                                        7d8s22
          iad=jdenpt(i2eu)+ntob+iacto(istob)*(nfromb+iacto(isfmb)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=-pa                                                         7d21s22
         else                                                           7d8s22
          iad=jdenpt(i2eu)+nfromb+iacto(isfmb)*(ntob+iacto(istob)         7d8s22
     $         *(nfroma+iacto(isfma)*ntoa))                              7d8s22
          wu=pa                                                         7d21s22
         end if                                                         7d8s22
         iad=iad+nad*mdenoff                                            3d15s23
         do k1=1,nroot                                                  3d15s23
          do k2=1,k1                                                    3d15s23
           term=2d0*wu*veca(k2,i1a,i2b)*veca(k1,i2a,i1b)                3d27s23
           bc(iad)=bc(iad)+term                                           7d21s22
           iad=iad+nad                                                  3d15s23
          end do                                                        3d15s23
         end do                                                         3d15s23
        end if
       end do                                                            3d10s17
      end do
      return                                                            3d10s17
      end                                                               3d10s17
