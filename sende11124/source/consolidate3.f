c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine consolidate3(ibasis,nfcn,igoal,iptrf,itmp,ibasisc,     3d24s21
     $     iptrcb,ioff,idorbfi,isorbfi,nec,bc,ibc)                      11d14s22
      implicit real*8 (a-h,o-z)                                         10d31s22
      integer*1 idorbfi(*),isorbfi(*)                                   3d24s21
      dimension ibasis(3,*),iptrf(4,*),itmp(*),ibasisc(3,*),            3d24s21
     $     iptrcb(2,*)                                                  3d24s21
      include "common.store"                                            3d24s21
      igoalm=igoal-1                                                    3d24s21
      do i=1,nfcn                                                       3d24s21
       if(ibasis(1,i).eq.igoalm)then                                     3d24s21
        nclo=ibasis(1,i)                                                3d24s21
        nopen=nec-nclo*2
        ic=iptrf(2,igoal)+nclo*(ibasis(2,i)-1)                          3d24s21
        io=iptrf(4,igoal)+nopen*(ibasis(3,i)-1)                         3d24s21
        ibasisc(1,ioff)=nclo                                            3d24s21
        iad=iptrcb(1,igoal)+itmp(igoal)                                 3d24s21
        ibc(iad)=0                                                      3d24s21
        do j=0,nclo-1                                                   3d24s21
         ibc(iad)=ibset(ibc(iad),idorbfi(ic+j))                           3d24s21
        end do                                                          3d24s21
        iad=iptrcb(2,igoal)+itmp(igoal)                                 3d24s21
        ibc(iad)=0                                                      3d24s21
        do j=0,nopen-1                                                   3d24s21
         ibc(iad)=ibset(ibc(iad),isorbfi(io+j))                           3d24s21
        end do                                                          3d24s21
        itmp(igoal)=itmp(igoal)+1                                       3d24s21
        ibasisc(2,ioff)=itmp(igoal)                                     3d24s21
        ibasisc(3,ioff)=itmp(igoal)                                     3d24s21
        ioff=ioff+1                                                     3d24s21
       end if                                                           3d24s21
      end do                                                            3d24s21
      return                                                            3d24s21
      end                                                               3d24s21
