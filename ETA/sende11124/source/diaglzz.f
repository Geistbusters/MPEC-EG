c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diaglzz(hdig,nconf,iaorb,nalpha,numa,iborb,nbeta,numb,
     $     ixlzze,nsymb,nsbeta,islz,multh,ilc,ihc,nlzz,bc,ibc)          11d14s22
      implicit real*8 (a-h,o-z)
c
c     compute diagonal matrix elements of lz^2 operator
c
      include "common.store"
      include "common.cas"
      integer*1 iaorb(nalpha,numa),iborb(nbeta,numb)
      integer*8 iarg1,iarg2                                             5d7s18
      dimension hdig(nconf),ixlzze(8,*),nsbeta(8),ilc(8),ihc(8),islz(*) 12d31s19
      npass=nlzz/2                                                      12d31s19
      norb=iacto(1)                                                     8d8s06
      do i=2,nsymb                                                      8d8s06
       norb=norb+iacto(i)                                               8d8s06
      end do                                                            8d8s06
      isma=ibcoff                                                        8d8s06
      ireloa=isma+norb                                                  8d8s06
      ibcoff=ireloa+norb                                                8d8s06
      call enough('diaglzz.  1',bc,ibc)
      jsm=isma-1                                                         8d8s06
      jrelo=ireloa-1                                                     8d8s06
      do i=1,nsymb                                                      8d8s06
       do j=1,iacto(i)                                                  8d8s06
        ibc(jsm+j)=i                                                    8d8s06
        ibc(jrelo+j)=j                                                  8d8s06
  153   format(3i5)
       end do
       jsm=jsm+iacto(i)                                                 8d8s06
       jrelo=jrelo+iacto(i)                                             8d8s06
      end do                                                            8d8s06
      ismb=ibcoff                                                        8d8s06
      irelob=ismb+norb                                                  8d8s06
      ibcoff=irelob+norb                                                8d8s06
      call enough('diaglzz.  2',bc,ibc)
      jsm=ismb-1                                                         8d8s06
      jrelo=irelob-1                                                     8d8s06
      do i=1,nsymb                                                      8d8s06
       do j=1,iacto(i)                                                  8d8s06
        ibc(jsm+j)=i                                                    8d8s06
        ibc(jrelo+j)=j                                                  8d8s06
       end do
       jsm=jsm+iacto(i)                                                 8d8s06
       jrelo=jrelo+iacto(i)                                             8d8s06
      end do                                                            8d8s06
c
c     pure alpha contributions
c
      itmpa=ibcoff
      ibcoff=itmpa+numa
      call enough('diaglzz.  3',bc,ibc)
      call diaglzz1(bc(itmpa),numa,nalpha,iaorb,ixlzze,debug,nsymb,     12d22s19
     $     nadet,ibc(isma),ibc(ireloa),islz,multh,npass,bc,ibc)         11d14s22
c
c     pure beta contributions
c
      itmpb=ibcoff
      ibcoff=itmpb+numb
      call enough('diaglzz.  4',bc,ibc)
      call diaglzz1(bc(itmpb),numb,nbeta,iborb,ixlzze,debug,nsymb,      12d22s19
     $     nbdet,ibc(isma),ibc(ireloa),islz,multh,npass,bc,ibc)         11d14s22
c
c     alpha-beta contributions
c
      ii=0
      jreloa=ireloa-1                                                   8d8s06
      jrelob=irelob-1                                                   8d8s06
      jsma=isma-1                                                       8d8s06
      jsmb=ismb-1                                                       8d8s06
      ioffa=0                                                           8d8s06
      do isb=1,nsymb                                                    8d8s06
       ioffb=0                                                          8d8s06
       do is2=1,nsbeta(isb)-1                                           8d8s06
        ioffb=ioffb+nbdet(is2)                                          8d8s06
       end do                                                           8d8s06
       do ia=ilc(isb),ihc(isb)                                          12d22s19
        iap=ia+ioffa                                                    8d8s06
        ia1=itmpa+iap-1                                                 8d8s06
        do ib=1,nbdet(nsbeta(isb))                                      8d8s06
         ibp=ib+ioffb                                                   8d8s06
         ib1=itmpb+ibp-1                                                8d8s06
         ii=ii+1
 3322    format('det ',i4,' beta ',2i4,i2,' alpha ',2i4,i2,1p2e15.7)
         hdig(ii)=bc(ia1)+bc(ib1)                                       12d22s19
 9332    format(i4,2i5,i2,4f18.8)
        end do
       end do
       ioffa=ioffa+nadet(isb)                                           8d8s06
      end do                                                            8d8s06
      ibcoff=isma                                                       8d8s06
      return
      end
