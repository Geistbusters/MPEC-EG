c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gen12(istng,no,nconf,idat1,il1,ih1,idat2,il2,ih2,
     $     idatc,ilc,ihc,itmp,m1s,m2s,itab,ndet,nsymb,bc,ibc)           11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 istng(no,nconf),itmp(no),idwstmp1(4)                    5d30s06
      integer*2 idwstmp2(2),idtmp(4)                                    8d29s06
      integer*8 idtmp8                                                  5d2s18
      integer nfrom(2),nto(2),m1s(8),m2s(8),ioff(8),                    8d8s06
     $     ioff1(8),ioff2(8),itab(8,8),ndet(8),idat1(8),il1(8),ih1(8),  8d25s06
     $     idat2(8),il2(8),ih2(8),idatc(36),ilc(36),ihc(36),ioffc(36),  5d2s18
     $     nsymb                                                        5d2s18
      equivalence (idwstmp2,idwstmp1)                                   5d30s06
      equivalence (idtmp,idtmp8)                                        8d29s06
      include "common.store"                                            8d29s06
      mxdtmp=0
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
       do iket=1,ibra-1
        ist=isymk                                                       5d7s18
        do ipass=ist,nsymb                                              5d7s18
         if(iket.gt.ioffk)then                                           12d13s06
          ioffk0=ioffk0+ndet(isymk)                                      12d19s06
          isymk=isymk+1                                                  12d13s06
          ioffk=ioffk+ndet(isymk)                                        12d13s06
         end if                                                          12d13s06
        end do                                                          5d7s18
        i12=itab(isymb,isymk)                                           8d8s06
        isum=0
        ibx=1                                                           7d22s22
        ikx=3                                                           7d22s22
        do i=1,no
         idf=istng(i,ibra)-istng(i,iket)                                7d22s22
         if(idf.gt.0)then                                               7d22s22
          idwstmp1(ibx)=i                                               7d22s22
          isum=isum+1                                                   7d22s22
          ibx=2                                                         7d22s22
         else if(idf.lt.0)then                                          7d22s22
          idwstmp1(ikx)=i                                               7d22s22
          isum=isum+1                                                   7d22s22
          ikx=4                                                         7d22s22
         end if                                                         7d22s22
        end do
        if(isum.eq.2)then                                               7d22s22
         itmp(1)=idwstmp1(1)                                            7d22s22
         itmp(2)=idwstmp1(3)                                            7d22s22
 3323    format('bra:ket ',2i4,' sym',2i2,4x,2i2)
        else if(isum.eq.4)then                                          7d22s22
         itmp(1)=idwstmp1(1)                                            7d22s22
         itmp(2)=idwstmp1(2)                                            7d22s22
         itmp(3)=idwstmp1(3)                                            7d22s22
         itmp(4)=idwstmp1(4)                                            7d22s22
        end if                                                          7d22s22
        if(i12.eq.1)then                                                8d25s06
         if(isum.eq.2)then
          ioff1(isymb)=ioff1(isymb)+1                                   8d29s06
          m1=ioff1(isymb)                                               8d29s06
          if(m1.ge.il1(isymb).and.m1.le.ih1(isymb))then                 8d29s06
           ii=m1-il1(isymb)                                             8d29s06
           itest=ibra-ioffb0                                            4d24s07
           mxdtmp=max(mxdtmp,itest)
           idtmp(1)=ibra-ioffb0                                         12d19s06
           if(idtmp(1).ne.itest)then                                    4d24s07
            write(6,*)('in gen12, integer*2 not enough '),itest         4d24s07
            stop                                                        4d24s07
           end if                                                       4d24s07
           itest=iket-ioffk0                                            4d24s07
           mxdtmp=max(mxdtmp,itest)
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
          if(m2.ge.il2(isymb).and.m2.le.ih2(isymb))then                 8d29s06
           ii=m2-il2(isymb)                                             8d29s06
           itest=ibra-ioffb0                                            4d24s07
           mxdtmp=max(mxdtmp,itest)
           idtmp(1)=ibra-ioffb0                                         12d21s06
           if(idtmp(1).ne.itest)then                                    4d24s07
            write(6,*)('integer*2 failure in gen12 '),itest             4d24s07
            stop                                                        4d24s07
           end if                                                       4d24s07
           itest=iket-ioffk0                                            4d24s07
           mxdtmp=max(mxdtmp,itest)
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
        else                                                            8d29s06
         if(isum.eq.2)then                                              8d29s06
          ii=((isymb*(isymb-1))/2)+isymk                                8d29s06
          ioffc(ii)=ioffc(ii)+1                                         8d29s06
          m1=ioffc(ii)                                                  8d29s06
          if(m1.ge.ilc(ii).and.m1.le.ihc(ii))then                       8d29s06
           itest=ibra-ioffb0                                            4d24s07
           mxdtmp=max(mxdtmp,itest)
           idtmp(1)=ibra-ioffb0                                         12d21s06
           if(idtmp(1).ne.itest)then                                    4d24s07
            write(6,*)('integer*2 failure in gen12 '),itest             4d24s07
            stop                                                        4d24s07
           end if                                                       4d24s07
           itest=iket-ioffk0                                            4d24s07
           mxdtmp=max(mxdtmp,itest)
           idtmp(2)=iket-ioffk0                                         12d21s06
           if(idtmp(2).ne.itest)then                                    4d24s07
            write(6,*)('integer*2 failure in gen12 '),itest             4d24s07
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
           ibc(idatc(ii)+m1-ilc(ii))=idtmp8                             8d29s06
          end if                                                        8d29s06
         end if                                                         8d29s06
        end if                                                          8d29s06
       end do
      end do
      return
      end
