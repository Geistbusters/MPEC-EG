c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sortdet(ialpha,iaorb,norb,nconf,itmp,itmp2,nalpha,
     $     ndet,itab,icall,lprint,bc,ibc)                               11d9s22
      implicit real*8 (a-h,o-z)                                         11d1s22
      integer*1 ialpha(norb,nconf),itmp(norb,nconf),
     $     iaorb(nalpha,nconf),itmp2(nalpha,nconf)
c
c     sort dets by symmetry.
c
      include "common.hf"                                               8d8s06
      include "common.cas"
      include "common.store"
      logical lprint                                                    4d20s18
      dimension ndet(8),ioff(8),itab(8,8)                               8d9s14
      isort1=ibcoff
      iso=isort1+nconf
      ibcoff=iso+norb
      call enough('sortdet.  1',bc,ibc)
      jso=iso
      do i=1,nsymb
       ndet(i)=0
       do j=1,iacto(i)
        ibc(jso)=i
        jso=jso+1
 1010   format(2i5)
       end do
      end do
      jso=iso-1
      jsort1=isort1-1
      do i=1,nconf
       isym=1
       do j=1,nalpha
        isym=itab(isym,ibc(jso+iaorb(j,i)))
       end do
       ibc(jsort1+i)=isym
       ndet(isym)=ndet(isym)+1
      end do
      ioff(1)=1                                                         8d8s06
      do i=2,nsymb                                                      8d8s06
       ioff(i)=ioff(i-1)+ndet(i-1)
      end do
      do i=1,nconf
       ioffu=ioff(ibc(jsort1+i))
       do j=1,norb
        itmp(j,ioffu)=ialpha(j,i)
       end do
       do j=1,nalpha
        itmp2(j,ioffu)=iaorb(j,i)
       end do
       ioff(ibc(jsort1+i))=ioff(ibc(jsort1+i))+1
      end do                                                            8d8s06
      do j=1,nconf
       do i=1,norb
        ialpha(i,j)=itmp(i,j)
       end do
      end do
      do i=1,nconf
       do j=1,nalpha
        iaorb(j,i)=itmp2(j,i)
       end do
      end do
      ibcoff=isort1
      return
      end
