c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gen12nona1(ipuse,istng,no,nconf,idat1,ih1,idat2,       6d13s22
     $     ih2,itmp,itab,ndet,nsymb,bc,ibc)                             11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 istng(no,nconf),itmp(no),idwstmp1(4)                    5d30s06
      integer*2 idwstmp2(2),idtmp(4)                                    8d29s06
      integer*8 idtmp8                                                  5d2s18
      integer nfrom(2),nto(2),m1s(8),m2s(8),ioff(8),                    8d8s06
     $     ioff1(8),ioff2(8),itab(8,8),ndet(8),idat1(8),ih1(8),         6d13s22
     $     idat2(8),ih2(8),ioffc(36),                                   6d13s22
     $     nsymb                                                        5d2s18
      equivalence (idwstmp2,idwstmp1)                                   5d30s06
      equivalence (idtmp,idtmp8)                                        8d29s06
      include "common.store"                                            8d29s06
c
c     generate data for singles and doubles differences for
c     this processor.
c
      do i=1,nsymb                                                      8d29s06
       ioff1(i)=0                                                       8d29s06
       ioff2(i)=0                                                       8d29s06
      end do                                                            8d29s06
      nn=(nsymb*(nsymb+1))/2                                            8d29s06
      do i=1,nn                                                         8d29s06
       ioffc(i)=0                                                       8d29s06
      end do                                                            8d29s06
      mxdelt1=0                                                         6d1s06
      mxdelt2=0                                                         6d1s06
      m1=0
      m2=0
      isymb=1                                                           8d8s06
      ioffb=ndet(1)                                                     12d13s06
      ioffb0=0                                                          12d19s06
      do ibra=1,nconf
       ist=isymb                                                        5d7s18
       do ipass=ist,nsymb                                               5d7s18
        if(ibra.gt.ioffb)then                                            12d13s06
         ioffb0=ioffb0+ndet(isymb)                                       12d19s06
         isymb=isymb+1                                                   12d13s06
         ioffb=ioffb+ndet(isymb)                                         12d13s06
        end if                                                           12d13s06
       end do                                                           5d7s18
       isymk=1                                                          8d8s06
       ioffk=ndet(1)                                                    12d13s06
       ioffk0=0                                                         12d19s06
       do iket=1,nconf                                                  6d11s22
        ist=isymk                                                       5d7s18
        do ipass=ist,nsymb                                              5d7s18
         if(iket.gt.ioffk)then                                           12d13s06
          ioffk0=ioffk0+ndet(isymk)                                      12d19s06
          isymk=isymk+1                                                  12d13s06
          ioffk=ioffk+ndet(isymk)                                        12d13s06
         end if                                                          12d13s06
        end do                                                          5d7s18
        i12=itab(isymb,isymk)                                           8d8s06
 3323   format('bra:ket ',2i4,' sym',2i2)
        isum=0
        do i=1,no
         if(istng(i,ibra).ne.istng(i,iket))then
          isum=isum+1
          itmp(isum)=i
         end if
        end do
        if(i12.eq.ipuse)then                                            6d11s22
         if(isum.eq.2)then
          ioff1(isymb)=ioff1(isymb)+1                                   8d29s06
          m1=ioff1(isymb)                                               8d29s06
          if(m1.ge.1.and.m1.le.ih1(isymb))then                          6d13s22
           ii=m1-1                                                      6d13s22
           itest=ibra-ioffb0                                            4d24s07
           idtmp(1)=ibra-ioffb0                                         12d19s06
           if(idtmp(1).ne.itest)then                                    4d24s07
            write(6,*)('in gen12, integer*2 not enough '),itest         4d24s07
            stop                                                        4d24s07
           end if                                                       4d24s07
           itest=iket-ioffk0                                            4d24s07
           idtmp(2)=iket-ioffk0                                         12d19s06
           if(idtmp(2).ne.itest)then                                    4d24s07
            write(6,*)('in gen12, integer*2 not enough '),itest         4d24s07
            stop                                                        4d24s07
           end if                                                       4d24s07
           mxdelt1=max(mxdelt1,iabs(ibra-iket))                          6d1s06
           if(istng(itmp(1),ibra).eq.0)then
            idwstmp1(1)=itmp(1)                                          5d30s06
            idwstmp1(2)=itmp(2)                                          5d30s06
           else
            idwstmp1(1)=itmp(2)                                          5d30s06
            idwstmp1(2)=itmp(1)                                          5d30s06
           end if
           idtmp(3)=idwstmp2(1)                                         8d29s06
           nb=0
           do i=min(idwstmp1(1),idwstmp1(2))+1,                          5d30s06
     $          max(idwstmp1(1),idwstmp1(2))-1                           5d30s06
            if(istng(i,ibra).eq.1)nb=nb+1
           end do
           if(mod(nb,2).eq.0)then
            idtmp(4)=1                                                  8d29s06
           else
            idtmp(4)=2                                                  8d29s06
           end if
           ibc(idat1(isymb)+ii)=idtmp8                                  8d29s06
    1     format('single',i5,5x,5i5)
          end if
         else if(isum.eq.4)then
          ioff2(isymb)=ioff2(isymb)+1                                   8d29s06
          m2=ioff2(isymb)                                               8d29s06
          if(m2.ge.1.and.m2.le.ih2(isymb))then                          6d13s22
           ii=m2-1                                                      6d13s22
           itest=ibra-ioffb0                                            4d24s07
           idtmp(1)=ibra-ioffb0                                         12d21s06
           if(idtmp(1).ne.itest)then                                    4d24s07
            write(6,*)('integer*2 failure in gen12 '),itest             4d24s07
            stop                                                        4d24s07
           end if                                                       4d24s07
           itest=iket-ioffk0                                            4d24s07
           idtmp(2)=iket-ioffk0                                         12d21s06
           if(idtmp(2).ne.itest)then                                    4d24s07
            write(6,*)('integer*2 failure in gen12 '),itest             4d24s07
            stop                                                        4d24s07
           end if                                                       4d24s07
           mxdelt2=max(mxdelt2,iabs(ibra-iket))                          6d1s06
           n1=1                                                          5d30s06
           n2=3                                                          5d30s06
           do i=1,4
            if(istng(itmp(i),ibra).eq.0)then
             idwstmp1(n2)=itmp(i)                                        5d30s06
             n2=n2+1
            else
             idwstmp1(n1)=itmp(i)                                        5d30s06
             n1=n1+1
            end if
           end do
           idtmp(3)=idwstmp2(1)                                         8d29s06
           idtmp(4)=idwstmp2(2)                                         8d29s06
           nfrom(1)=idwstmp1(1)                                          5d30s06
           nfrom(2)=idwstmp1(2)                                          5d30s06
           call phsdet(istng(1,ibra),itmp,nfrom,iphsb,no)                11d14s05
           nto(1)=idwstmp1(3)                                            5d30s06
           nto(2)=idwstmp1(4)                                            5d30s06
           call phsdet(istng(1,iket),itmp,nto,iphsk,no)                  11d14s05
           iphs=iphsb+iphsk                                              11d14s05
           if(mod(iphs,2).eq.0)then
            idtmp(1)=-idtmp(1)                                           8d29s06
           end if
           ibc(idat2(isymb)+ii)=idtmp8                                  8d29s06
          end if
         end if
        end if                                                          8d29s06
       end do
      end do
      return
      end
