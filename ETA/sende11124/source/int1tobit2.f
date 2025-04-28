c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine int1tobit2(ibasis,ibasisn,iptr,iptrbit,idorb,isorb,
     $     myfcn,nec,mdoop,bc,ibc)                                      11d10s22
      implicit real*8 (a-h,o-z)                                         11d1s22
      integer*1 idorb(*),isorb(*)                                       3d24s21
      dimension ibasis(3,*),iptrbit(2,*),ibasisn(3,*),iptr(4,*)         3d24s21
      include "common.store"                                            3d24s21
c
c     count
c
      iorig=ibcoff                                                      3d24s21
      itmp=iorig+myfcn*2                                                3d24s21
      itmpc=itmp+mdoop                                                  3d24s21
      ibcoff=itmpc+mdoop                                                 3d24s21
      call enough('int1tobit2.  1',bc,ibc)
      jtmp=itmp-1                                                       3d24s21
      jtmpc=itmpc-1                                                     3d24s21
      do i=1,mdoop                                                      3d24s21
       ibc(jtmp+i)=0                                                    3d24s21
       ibc(jtmpc+i)=0                                                   3d24s21
      end do                                                            3d24s21
      do i=1,myfcn                                                      3d24s21
       nclo=ibasis(1,i)                                                 3d24s21
       ibc(itmp+nclo)=ibc(itmp+nclo)+1                                  3d24s21
      end do                                                            3d24s21
      isto=iorig                                                        3d24s21
      do i=1,mdoop                                                      3d24s21
       iptrbit(1,i)=isto                                                3d24s21
       iptrbit(2,i)=isto+ibc(jtmp+i)                                    3d24s21
       isto=isto+ibc(jtmp+i)*2                                          3d24s21
      end do                                                            3d24s21
c
c     store
c
      do i=1,myfcn                                                      3d24s21
       nclo=ibasis(1,i)                                                 3d24s21
       nclop=nclo+1                                                     3d24s21
       nopen=nec-2*nclo                                                 3d24s21
       ic=iptr(2,nclop)+nclo*(ibasis(2,i)-1)                            3d24s21
       io=iptr(4,nclop)+nopen*(ibasis(3,i)-1)                           3d24s21
       ibasisn(1,i)=nclo                                                3d24s21
       iadc=iptrbit(1,nclop)+ibc(itmpc+nclo)
       iado=iptrbit(2,nclop)+ibc(itmpc+nclo)
       ibc(iadc)=0                                                      3d24s21
       do j=0,nclo-1                                                    3d24s21
        ibc(iadc)=ibset(ibc(iadc),idorb(ic+j))                           3d24s21
       end do                                                           3d24s21
       ibc(iado)=0                                                      3d24s21
       do j=0,nopen-1                                                    3d24s21
        ibc(iado)=ibset(ibc(iado),isorb(io+j))                           3d24s21
       end do                                                           3d24s21
       ibc(itmpc+nclo)=ibc(itmpc+nclo)+1                                3d24s21
       ibasisn(2,i)=ibc(itmpc+nclo)                                     3d24s21
       ibasisn(3,i)=ibc(itmpc+nclo)                                     3d24s21
      end do                                                            3d24s21
      ibcoff=itmp                                                       3d24s21
      return                                                            3d24s21
      end                                                               3d24s21
