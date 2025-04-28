c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine wfetch(nopen,mdon,nclo,is2,ism2,ndet,ncsf,ivec,iaorb,  5d11s21
     $     icsfxv,bc,ibc)                                               11d9s22
      implicit real*8 (a-h,o-z)
      include "common.store"
      integer*8 i18,i28,i38,i48,icsfxv                                  5d11s21
      real*16 fl                                                        7d11s19
      COMMON/FACT16/FL(922),NCALL                                       5d27s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
c
c     nalpha+nbeta=nopen, nalpha-nbeta=ism2, 2*nalpha=nopen+ism2
c
      nalpha=(nopen+ism2)/2                                             3d20s20
      nbeta=nopen-nalpha                                                3d20s20
      itop=nopen+1                                                      3d20s20
      ibot1=nalpha+1                                                    3d20s20
      ibot2=itop-ibot1+1                                                3d20s20
      det=fl(itop)-fl(ibot1)-fl(ibot2)                                  3d20s20
      ndet=nint(exp(det))                                               3d20s20
      iaorb=ibcoff                                                      3d20s20
      ibcoff=iaorb+ndet*nopen                                           5d11s21
      if(ibcoff.lt.0)write(6,*)('enoughb '),ibcoff
      call enough('wfetch.  1',bc,ibc)
      call fetchcsf2(nopen,is2,ncsf,vecadd,idum,ibc(iaorb),ism2,        5d11s21
     $       ndet,bc,ibc)                                               11d15s22
      ivec=nint(vecadd)                                                 5d11s21
      return                                                            3d20s20
      end                                                               3d20s20
