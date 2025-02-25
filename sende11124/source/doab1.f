c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine doab1(idat1b,nb,idat1a,na,veca,g,numa,numb,nroot,nhere,3d15s17
     $     il,ih,i2e,irelo,ism,iacto,msta,mstb,nsymb,bc,ibc)            11d14s22
      implicit real*8 (a-h,o-z)                                         3d10s17
      integer*2 idat1b(4,nb),idat1a(4,na),idwstmp2                      3d12s17
      integer*1 idwstmp1(4)                                             3d10s17
      equivalence (idwstmp2,idwstmp1)                                   3d10s17
      integer*8 ism,irelo                                               3d10s17
      integer msta(8),mstb(8)                                           4d26s18
      logical log1,log2                                                 4d9s18
      include "common.store"                                            10d6s06
      dimension veca(nroot,numa,numb),g(nroot,nhere,numb),              3d14s17
     $     i2e(1),iacto(8),ism(1),irelo(1),phs(2)                       3d10s17
      data phs/1d0,-1d0/                                                3d10s17
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      ioffa=0                                                           4d26s18
      do istoa=1,nsymb                                                    4d26s18
       isfma=istoa                                                      4d26s18
       do ia0=1,msta(istoa)                                             4d26s18
        ia=ia0+ioffa                                                    4d26s18
        i1a=idat1a(1,ia)                                                 3d10s17
        i2a=idat1a(2,ia)                                                 3d10s17
        log1=i1a.ge.il.and.i1a.le.ih                                     4d9s18
        log2=i2a.ge.il.and.i2a.le.ih                                     4d9s18
        fact1=0d0                                                       4d30s18
        if(log1)fact1=1d0                                               4d30s18
        fact2=0d0                                                       4d30s18
        if(log2)fact2=1d0                                               4d30s18
        if(log1.or.log2)then                                             4d9s18
         idwstmp2=idat1a(3,ia)                                           3d10s17
         ntoa=irelo(idwstmp1(1))                                         4d9s18
         nfroma=irelo(idwstmp1(2))                                       4d9s18
         in=min(nfroma,ntoa)                                            4d26s18
         ix=max(nfroma,ntoa)                                            4d26s18
         ii2=((ix*(ix-1))/2)+in-1                                       4d26s18
         pa=phs(idat1a(4,ia))                                            4d9s18
         i1ao=i1a+1-il                                                   4d9s18
         i2ao=i2a+1-il                                                   4d9s18
         nval1=ntoa-1+iacto(istoa)*(nfroma-1)                            4d9s18
         nval2=nfroma-1+iacto(isfma)*(ntoa-1)                            4d9s18
         if(i1ao.ge.1.and.i1ao.le.nhere)then
          ioffb=0
          do istob=1,nsymb                                               4d26s18
           nrow=(iacto(istob)*(iacto(istob)+1))/2                          3d10s17
           isfmb=istob                                                   4d26s18
           i2eu=inv(1,istob,isfmb,istoa)                                  4d9s18
           jad1=i2e(i2eu)+ii2*nrow-1                                     4d26s18
           do ib0=1,mstb(istob)                                          4d26s18
            ib=ib0+ioffb                                                 4d26s18
            i1b=idat1b(1,ib)                                               4d9s18
            i2b=idat1b(2,ib)                                               4d9s18
            idwstmp2=idat1b(3,ib)                                          4d9s18
            ntob=irelo(idwstmp1(1))                                        4d9s18
            nfromb=irelo(idwstmp1(2))                                      4d9s18
            pb=pa*phs(idat1b(4,ib))                                        4d9s18
            in=min(nfromb,ntob)                                           3d15s17
            ix=max(nfromb,ntob)                                           3d15s17
            ii1=((ix*(ix-1))/2)+in                                       4d26s18
            iad1=jad1+ii1                                                4d26s18
            fact=pb*bc(iad1)                                             4d30s18
            facta=fact1*fact
            factb=fact2*fact                                             4d30s18
            do k=1,nroot                                                 4d30s18
             g(k,i1ao,i1b)=g(k,i1ao,i1b)+facta*veca(k,i2a,i2b)           4d30s18
             g(k,i1ao,i2b)=g(k,i1ao,i2b)+facta*veca(k,i2a,i1b)           4d30s18
            end do                                                       4d30s18
           end do                                                        4d26s18
           ioffb=ioffb+mstb(istob)                                       4d26s18
          end do                                                         4d26s18
         end if
         if(i2ao.ge.1.and.i2ao.le.nhere)then
          ioffb=0
          do istob=1,nsymb                                               4d26s18
           nrow=(iacto(istob)*(iacto(istob)+1))/2                          3d10s17
           isfmb=istob                                                   4d26s18
           i2eu=inv(1,istob,isfmb,istoa)                                  4d9s18
           jad1=i2e(i2eu)+ii2*nrow-1                                     4d26s18
           do ib0=1,mstb(istob)                                          4d26s18
            ib=ib0+ioffb                                                 4d26s18
            i1b=idat1b(1,ib)                                               4d9s18
            i2b=idat1b(2,ib)                                               4d9s18
            idwstmp2=idat1b(3,ib)                                          4d9s18
            ntob=irelo(idwstmp1(1))                                        4d9s18
            nfromb=irelo(idwstmp1(2))                                      4d9s18
            pb=pa*phs(idat1b(4,ib))                                        4d9s18
            in=min(nfromb,ntob)                                           3d15s17
            ix=max(nfromb,ntob)                                           3d15s17
            ii1=((ix*(ix-1))/2)+in                                       4d26s18
            iad1=jad1+ii1                                                4d26s18
            fact=pb*bc(iad1)                                             4d30s18
            facta=fact1*fact
            factb=fact2*fact                                             4d30s18
            do k=1,nroot                                                 4d30s18
             g(k,i2ao,i1b)=g(k,i2ao,i1b)+factb*veca(k,i1a,i2b)           4d30s18
             g(k,i2ao,i2b)=g(k,i2ao,i2b)+factb*veca(k,i1a,i1b)           4d30s18
            end do                                                       4d30s18
           end do                                                        4d26s18
           ioffb=ioffb+mstb(istob)                                       4d26s18
          end do                                                         4d26s18
         end if
        end if                                                          4d26s18
       end do                                                           4d26s18
       ioffa=ioffa+msta(istoa)                                          4d26s18
      end do                                                            4d26s18
      return                                                            3d10s17
      end                                                               3d10s17
