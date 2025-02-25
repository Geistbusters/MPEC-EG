c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genf0(iptr,ibasis,nfcn,mdon,mdoo,nff,iff,ncsf,bc,ibc)  11d11s22
      implicit real*8 (a-h,o-z)
c
c     take old bit storage and turn it into new storage
c
      dimension iptr(2,*),ibasis(3,*),nff(mdoo+1,3),iff(*),ncsf(*)      11d25s20
      include "common.store"                                            11d25s20
      do i=1,mdoo+1                                                     11d25s20
       nff(i,1)=0                                                       11d25s20
      end do                                                            11d25s20
      ioff=1                                                            11d25s20
      ioffc=1                                                           11d25s20
      last=-1                                                           11d25s20
      do if=1,nfcn
       nclo=ibasis(1,if)
       nclop=nclo+1                                                     11d25s20
       iarg=nclop-mdon                                                  11d25s20
       if(nclo.ne.last)then                                             11d25s20
        nff(nclop,2)=ioff                                               11d25s20
        nff(nclop,3)=ioffc                                              11d25s20
        last=nclo                                                       11d26s20
       end if                                                           11d25s20
       nff(nclop,1)=nff(nclop,1)+1                                      11d25s20
       iic=iptr(1,nclop)+ibasis(2,if)-1                                 11d25s20
       iio=iptr(2,nclop)+ibasis(3,if)-1                                 11d25s20
       iff(ioff)=ibc(iic)                                               11d25s20
       ioff=ioff+1                                                      11d25s20
       iff(ioff)=ibc(iio)                                               11d25s20
       ioff=ioff+1                                                      11d25s20
       ioffc=ioffc+ncsf(iarg)                                           11d25s20
      end do                                                            11d25s20
      return                                                            11d25s20
      end                                                               11d25s20
