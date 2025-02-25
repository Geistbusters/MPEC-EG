c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dodoubdbl(idat2,ndoub,cvec,numb,nhere,iacto,ism,irelo, 7d21s22
     $     jdenpt,bc,ibc,nroot,mdenoff,igoal)                           3d15s23
      implicit real*8 (a-h,o-z)
      integer*2 idat2(4,ndoub),idwstmp2(2)                              3d21s17
      integer*1 idwstmp1(4)                                             5d30s06
      equivalence (idwstmp2,idwstmp1)                                   5d31s06
      integer*8 ism,irelo                                               12d13s06
      include "common.store"                                            10d6s06
      dimension cvec(nroot,nhere,numb),jdenpt(*),iacto(8),ism(*),       3d14s23
     $     irelo(*),phs(2)                                              3d14s23
      dimension igoalx(2)
      data phs/1d0,-1d0/
      do i=1,ndoub
       i1signed=idat2(1,i)                                              3d12s17
       i1=iabs(i1signed)                                                3d12s17
       iphdws=(((idat2(1,i)-i1)/(2*i1)))+2                              3d12s17
       i2=idat2(2,i)
       idwstmp2(1)=idat2(3,i)                                           5d30s06
       idwstmp2(2)=idat2(4,i)                                           5d30s06
       n1=irelo(idwstmp1(1))-1                                          7d21s22
       n2=irelo(idwstmp1(2))-1                                          7d21s22
       n3=irelo(idwstmp1(3))-1                                          7d21s22
       n4=irelo(idwstmp1(4))-1                                          7d21s22
       isx1=ism(idwstmp1(1))                                             9d1s06
       isx2=ism(idwstmp1(2))                                             9d1s06
       isx3=ism(idwstmp1(3))                                             9d1s06
       isx4=ism(idwstmp1(4))                                             9d1s06
c     2 (b1k1b2k2-b1k2b2k1)
       i2eu=ifind2(isx1,isx3,isx2,isx4,icase)                           12d19s06
       nad=iacto(isx1)*iacto(isx2)*iacto(isx3)*iacto(isx4)              3d14s23
       if(icase.eq.1)then                                               7d21s22
        iadj=jdenpt(i2eu)+n1+iacto(isx1)*(n3+iacto(isx3)                 7d21s22
     $         *(n2+iacto(isx2)*n4))                                    7d21s22
        wj=phs(iphdws)
       else if(icase.eq.4)then                                          7d21s22
        iadj=jdenpt(i2eu)+n3+iacto(isx3)*(n1+iacto(isx1)                7d21s22
     $         *(n2+iacto(isx2)*n4))                                    7d21s22
        wj=-phs(iphdws)
       else if(icase.eq.3)then                                          7d21s22
        iadj=jdenpt(i2eu)+n1+iacto(isx1)*(n3+iacto(isx3)                7d21s22
     $         *(n4+iacto(isx4)*n2))                                    7d21s22
        wj=-phs(iphdws)
       else                                                             7d21s22
        iadj=jdenpt(i2eu)+n3+iacto(isx3)*(n1+iacto(isx1)                7d21s22
     $         *(n4+iacto(isx4)*n2))                                    7d21s22
        wj=phs(iphdws)                                                  7d21s22
       end if                                                           7d21s22
       icasej=icase
       i2eu=ifind2(isx1,isx4,isx2,isx3,icase)                           12d19s06
       if(icase.eq.1)then                                               7d21s22
        iadk=jdenpt(i2eu)+n1+iacto(isx1)*(n4+iacto(isx4)                 7d21s22
     $         *(n2+iacto(isx2)*n3))                                    7d21s22
        wk=phs(iphdws)
       else if(icase.eq.4)then                                          7d21s22
        iadk=jdenpt(i2eu)+n4+iacto(isx4)*(n1+iacto(isx1)                7d21s22
     $         *(n2+iacto(isx2)*n3))                                    7d21s22
        wk=-phs(iphdws)
       else if(icase.eq.3)then                                          7d21s22
        iadk=jdenpt(i2eu)+n1+iacto(isx1)*(n4+iacto(isx4)                7d21s22
     $         *(n3+iacto(isx3)*n2))                                    7d21s22
        wk=-phs(iphdws)
       else                                                             7d21s22
        iadk=jdenpt(i2eu)+n4+iacto(isx4)*(n1+iacto(isx1)                7d21s22
     $         *(n3+iacto(isx3)*n2))                                    7d21s22
        wk=phs(iphdws)                                                  7d21s22
       end if                                                           7d21s22
       wj=wj*2d0                                                        3d27s23
       wk=wk*2d0                                                        3d27s23
       iadj=iadj+mdenoff*nad                                            3d14s23
       iadk=iadk+mdenoff*nad                                            3d14s23
       do k1=1,nroot                                                    3d14s23
        do k2=1,k1                                                      3d14s23
         do j=1,nhere                                                     3d27s17
          vv=cvec(k1,j,i1)*cvec(k2,j,i2)                                3d14s23
          bc(iadj)=bc(iadj)+wj*vv                                         7d21s22
          bc(iadk)=bc(iadk)-wk*vv                                         7d21s22
         end do                                                           3d27s17
         iadj=iadj+nad                                                  3d14s23
         iadk=iadk+nad                                                  3d14s23
        end do                                                          3d14s23
       end do                                                           3d14s23
c     2 (k1b1k2b2-k1b2k2b1)
       i2eu=ifind2(isx3,isx1,isx4,isx2,icase)                           3d27s23
       if(icase.eq.1)then                                               7d21s22
        iadj=jdenpt(i2eu)+n3+iacto(isx3)*(n1+iacto(isx1)                3d27s23
     $         *(n4+iacto(isx4)*n2))                                    3d27s23
        wj=phs(iphdws)
       else if(icase.eq.4)then                                          7d21s22
        iadj=jdenpt(i2eu)+n1+iacto(isx1)*(n3+iacto(isx3)                3d27s23
     $         *(n4+iacto(isx4)*n2))                                    3d27s23
        wj=-phs(iphdws)
       else if(icase.eq.3)then                                          7d21s22
        iadj=jdenpt(i2eu)+n3+iacto(isx3)*(n1+iacto(isx1)                3d27s23
     $         *(n2+iacto(isx2)*n4))                                    3d27s23
        wj=-phs(iphdws)
       else                                                             7d21s22
        iadj=jdenpt(i2eu)+n1+iacto(isx1)*(n3+iacto(isx3)                7d21s22
     $         *(n2+iacto(isx2)*n4))                                    7d21s22
        wj=phs(iphdws)                                                  7d21s22
       end if                                                           7d21s22
       icasej=icase
       i2eu=ifind2(isx3,isx2,isx4,isx1,icase)                           3d27s23
       if(icase.eq.1)then                                               7d21s22
        iadk=jdenpt(i2eu)+n3+iacto(isx3)*(n2+iacto(isx2)                3d27s23
     $         *(n4+iacto(isx4)*n1))                                    3d27s23
        wk=phs(iphdws)
       else if(icase.eq.4)then                                          7d21s22
        iadk=jdenpt(i2eu)+n2+iacto(isx2)*(n3+iacto(isx3)                3d27s23
     $         *(n4+iacto(isx4)*n1))                                    3d27s23
        wk=-phs(iphdws)
       else if(icase.eq.3)then                                          7d21s22
        iadk=jdenpt(i2eu)+n3+iacto(isx3)*(n2+iacto(isx2)                3d27s23
     $         *(n1+iacto(isx1)*n4))                                    3d27s23
        wk=-phs(iphdws)
       else                                                             7d21s22
        iadk=jdenpt(i2eu)+n2+iacto(isx2)*(n3+iacto(isx3)                3d27s23
     $         *(n1+iacto(isx1)*n4))                                    3d27s23
        wk=phs(iphdws)                                                  7d21s22
       end if                                                           7d21s22
       wj=wj*2d0                                                        3d27s23
       wk=wk*2d0                                                        3d27s23
       iadj=iadj+mdenoff*nad                                            3d14s23
       iadk=iadk+mdenoff*nad                                            3d14s23
       do k1=1,nroot                                                    3d14s23
        do k2=1,k1                                                      3d14s23
         do j=1,nhere                                                     3d27s17
          vv=cvec(k1,j,i2)*cvec(k2,j,i1)                                3d27s23
          bc(iadj)=bc(iadj)+wj*vv                                         7d21s22
          bc(iadk)=bc(iadk)-wk*vv                                         7d21s22
         end do                                                           3d27s17
         iadj=iadj+nad                                                  3d14s23
         iadk=iadk+nad                                                  3d14s23
        end do                                                          3d14s23
       end do                                                           3d14s23
      end do
      return
      end
