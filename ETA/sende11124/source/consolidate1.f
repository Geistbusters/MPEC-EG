c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine consolidate1(ibasisfi,iptrfi,mdoop,nfcnfi,idorbfi,     3d24s21
     $     isorbfi,nsymb,ibasisc,iptrcb,nfcnc,nec,bc,ibc)               11d14s22
      implicit real*8 (a-h,o-z)                                         3d24s21
      integer*1 idorbfi(*),isorbfi(*)                                   3d24s21
      dimension ibasisfi(*),iptrfi(4,mdoop,*),nfcnfi(*),ibasisc(3,*),   3d24s21
     $     iptrcb(2,*)                                                  3d24s21
      include "common.store"                                            3d24s21
      nfcnc=0                                                            3d24s21
      do isb=1,nsymb                                                    3d24s21
       nfcnc=nfcnc+nfcnfi(isb)                                          3d24s21
      end do                                                            3d24s21
      iorig=ibcoff                                                      3d24s21
      itmp=iorig+2*nfcnc                                                3d24s21
      itmpc=itmp+mdoop                                                  3d24s21
      ibcoff=itmpc+mdoop                                                3d24s21
      jtmp=itmp-1                                                       3d24s21
      jtmpc=itmpc-1                                                     3d24s21
      call enough('consolidate1.  1',bc,ibc)
      do i=1,mdoop                                                      3d24s21
       nhere=0                                                          3d24s21
       do isb=1,nsymb                                                   3d24s21
        call consolidate2(ibc(ibasisfi(isb)),nfcnfi(isb),i,nhere)       3d24s21
       end do                                                           3d24s21
       ibc(jtmp+i)=nhere                                                3d24s21
      end do                                                            3d24s21
      isto=iorig                                                        3d24s21
      do i=1,mdoop                                                      3d24s21
       iptrcb(1,i)=isto                                                 3d24s21
       iptrcb(2,i)=isto+ibc(jtmp+i)                                     3d24s21
       isto=isto+2*ibc(jtmp+i)                                          3d24s21
       ibc(jtmpc+i)=0                                                   3d24s21
      end do                                                            3d24s21
      ioff=1                                                            3d24s21
      do i=1,mdoop                                                      3d24s21
       do isb=1,nsymb                                                   3d24s21
        call consolidate3(ibc(ibasisfi(isb)),nfcnfi(isb),i,             3d24s21
     $       iptrfi(1,1,isb),ibc(itmpc),ibasisc,iptrcb,ioff,idorbfi,    3d24s21
     $       isorbfi,nec,bc,ibc)                                        11d14s22
       end do                                                           3d24s21
      end do                                                            3d24s21
      ibcoff=itmp                                                       3d24s21
      return                                                            3d24s21
      end                                                               3d24s21
