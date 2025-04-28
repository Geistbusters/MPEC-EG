c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumpbas(nfcn,iptr,ibasis,idorb,isorb,nod,nos,mdoop,nec,11d23s19
     $     icode,nfcnx,ncsf,mdon,nwavrec,nrestart,bc,ibc)               11d14s22
      implicit real*8 (a-h,o-z)                                         11d1s22
      integer*1 idorb(nod),isorb(nos)                                   11d23s19
      dimension ibasis(3,nfcn),iptr(4,mdoop),ncsf(*)                    4d28s21
      include "common.store"                                            4d28s21
c
c     nff0(nclop,1) =  no fcns
c                2  =  offset in iff0
c                3  =  offset in vectors (starting from 1)
c     iff0 i0c,i0o pairs
c     these are both integer*4
c
      itmp1=ibcoff                                                      4d28s21
      itmp2=itmp1+mdoop                                                 4d28s21
      itmp3=itmp2+mdoop*3                                               4d28s21
      ibcoff=itmp3+nfcn*2                                               4d28s21
      call dumpbas1(nfcn,iptr,ibasis,idorb,isorb,mdoop,nec,ncsf,        4d28s21
     $     ibc(itmp1),ibc(itmp2),ibc(itmp3),mdon,nwavrec,nfcnx,nrestart)8d11s22
      ibcoff=itmp1                                                      4d28s21
      return                                                            11d23s19
      end                                                               11d23s19
