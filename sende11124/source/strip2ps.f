c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine strip2ps(iptrps,ibasisps,nfcn,nfcnc,ibasisc,iptrcb,    3d24s21
     $     ikeep,nkeep,norbx,iptrqs,ibasisqs,nfcnq,ikeepa,bc,ibc)       11d10s22
      implicit real*8 (a-h,o-z)                                         3d23s21
      dimension iptrps(2,*),ibasisps(3,*),ibasisc(3,*),iptrcb(2,*),     3d24s21
     $     ikeep(*),nab4(2,3),iptrqs(2,*),ibasisqs(3,*),ikeepa(*)       3d31s21
      include "common.store"                                            3d23s21
      nkeep=0                                                           3d24s21
      do ifi=1,nfcnc                                                    3d24s21
       ikeepa(ifi)=ifi                                                  3d25s21
       nclo=ibasisc(1,ifi)
       nclop=nclo+1                                                     3d24s21
       iic=iptrcb(1,nclop)+ibasisc(2,ifi)-1
       iio=iptrcb(2,nclop)+ibasisc(3,ifi)-1
       nopen=popcnt(ibc(iio))
       do if=1,nfcn                                                     3d24s21
        if(iabs(ibasisps(1,if)-ibasisc(1,ifi)).le.1)then                3d24s21
         nclo=ibasisps(1,if)
         nclop=nclo+1
         jjc=iptrps(1,nclop)+ibasisps(2,if)-1                             3d23s21
         jjo=iptrps(2,nclop)+ibasisps(3,if)-1                             3d23s21
         nopenj=popcnt(ibc(jjo))                                        3d24s21
         call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  3d24s21
     $        norbx,nnot,nab4,bc,ibc)                                   11d14s22
         if(nnot.gt.0.and.nnot.le.2)then                                3d24s21
          do ifq=1,nfcnq                                                3d24s21
           if(iabs(ibasisqs(1,ifq)-ibasisc(1,ifi)).le.1)then            3d24s21
            nclo=ibasisqs(1,ifq)                                        3d24s21
            nclop=nclo+1                                                3d24s21
            jjc=iptrqs(1,nclop)+ibasisqs(2,ifq)-1                       3d24s21
            jjo=iptrqs(2,nclop)+ibasisqs(3,ifq)-1                       3d24s21
            nopenj=popcnt(ibc(jjo))                                     3d24s21
            call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,     3d24s21
     $           nopen,norbx,nnot,nab4,bc,ibc)                          11d14s22
            if(nnot.gt.0.and.nnot.le.2)then                                3d24s21
             nkeep=nkeep+1                                                 3d24s21
             ikeep(nkeep)=ifi                                              3d24s21
             go to 1                                                       3d24s21
            end if                                                      3d24s21
           end if                                                       3d24s21
          end do                                                        3d24s21
         end if                                                         3d24s21
        end if                                                          3d24s21
       end do                                                           3d24s21
    1  continue                                                         3d24s21
      end do                                                            3d24s21
      return
      end
