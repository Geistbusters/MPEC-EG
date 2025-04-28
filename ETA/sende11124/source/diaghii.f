c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diaghii(hdig,nconf,iaorb,nalpha,numa,
     $                   iborb,nbeta,numb,ih0e,i2e,shift,ilc,ihc,       8d8s06
     $                   debug,nsymb,nsbeta,ihdig,bc,ibc)               11d14s22
      implicit real*8 (a-h,o-z)
c
c     compute diagonal matrix elements of hamiltonian
c
      include "common.store"
      include "common.cas"
      integer*1 iaorb(nalpha,numa),iborb(nbeta,numb)
      integer*8 iarg1,iarg2                                             5d7s18
      integer*2 ihdig(4,*)                                              6d12s23
      dimension hdig(*),i2e(1),ilc(8),ihc(8),nsbeta(8),ih0e(8)          6d12s23
      data icall/0/
      save icall
      icall=icall+1
      loop=0
      loopx=100000
      norb=iacto(1)                                                     8d8s06
      do i=2,nsymb                                                      8d8s06
       norb=norb+iacto(i)                                               8d8s06
      end do                                                            8d8s06
      isma=ibcoff                                                        8d8s06
      ireloa=isma+norb                                                  8d8s06
      ibcoff=ireloa+norb                                                8d8s06
      call enough('diaghii.  1',bc,ibc)
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
      call enough('diaghii.  2',bc,ibc)
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
      call enough('diaghii.  3',bc,ibc)
      call diaghii1(bc(itmpa),numa,nalpha,iaorb,ih0e,i2e,debug,nsymb,   8d8s06
     $     nadet,ibc(isma),ibc(ireloa),bc,ibc)                          11d14s22
c
c     pure beta contributions
c
      itmpb=ibcoff
      ibcoff=itmpb+numb
      call enough('diaghii.  4',bc,ibc)
      call diaghii1(bc(itmpb),numb,nbeta,iborb,ih0e,i2e,debug,nsymb,    8d8s06
     $     nbdet,ibc(isma),ibc(ireloa),bc,ibc)                          11d14s22
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
       do ia=ilc(isb),ihc(isb)                                          8d8s06
        iap=ia+ioffa                                                    8d8s06
        ia1=itmpa+iap-1                                                 8d8s06
        do ib=1,nbdet(nsbeta(isb))                                      8d8s06
         ibp=ib+ioffb                                                   8d8s06
         ib1=itmpb+ibp-1                                                8d8s06
         ii=ii+1
 3322    format('det ',i4,' beta ',2i4,i2,' alpha ',2i4,i2,1p2e15.7)
         hdig(ii)=bc(ia1)+bc(ib1)+shift*debug                            7d12s06
         ihdig(1,ii)=ia                                                 8d22s06
         ihdig(2,ii)=ib                                                 8d22s06
         ihdig(3,ii)=isb                                                8d22s06
         ihdig(4,ii)=0                                                  8d22s06
         addsum=0d0
         do ja=1,nalpha
          jao=ibc(jreloa+iaorb(ja,iap))                                 8d8s06
          jja=((jao*(jao+1))/2)-1
          na=ibc(jsma+iaorb(ja,iap))                                    8d8s06
          nna=(iacto(na)*(iacto(na)+1))/2                               8d8s06
          do jb=1,nbeta
           jbo=ibc(jrelob+iborb(jb,ibp))                                8d8s06
           jjb=((jbo*(jbo+1))/2)-1
           nb=ibc(jsmb+iborb(jb,ibp))                                   8d8s06
           i2eu=ifind2(nb,nb,na,na,icase)                               8d15s06
           if(icase.eq.1)then                                           4d24s07
            nna=(iacto(nb)*(iacto(nb)+1))/2                              8d14s06
            iadd=i2e(i2eu)+jjb+jja*nna                                  4d24s07
           else if(icase.eq.5)then                                      4d24s07
            nna=(iacto(na)*(iacto(na)+1))/2                             4d24s07
            iadd=i2e(i2eu)+jja+jjb*nna                                  4d24s07
           else                                                         4d24s07
            write(6,*)('icase is not 1 or 5 '),icase,nb,na
            stop
           end if
           hdig(ii)=hdig(ii)+bc(iadd)                                    7d12s06
          end do
         end do
 9332    format(i4,2i5,i2,4f18.8)
        end do
       end do
       ioffa=ioffa+nadet(isb)                                           8d8s06
      end do                                                            8d8s06
      ibcoff=isma                                                       8d8s06
      return
      end
