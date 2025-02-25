c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine spindig(sdig,iaorb,nalpha,iborb,nbeta,nsymb,nsbeta,    3d22s17
     $     nherec,ilc,bc,ibc)                                           11d9s22
      implicit real*8 (a-h,o-z)
      integer*1 iaorb,iborb
      include "common.cas"
      include "common.store"                                            5d30s06
c
c     explicitly form p-space S**2 diagonal matrix elements that are    3d22s17
c     on this proc.                                                     3d22s17
c
      dimension sdig(1),iaorb(nalpha,1),iborb(nbeta,1),                 3d22s17
     $          nherec(1),ilc(1),nsbeta(8)                              3d22s17
      xms=0.5d0*dfloat(nalpha-nbeta)
      ioff=0
      prevalue=xms*xms+0.5d0*dfloat(nalpha+nbeta)                       3d22s17
      do isb=1,nsymb                                                    3d22s17
       nasum=0                                                          8d15s06
       do i=1,isb-1                                                     3d22s17
        nasum=nasum+nadet(i)                                            3d22s17
       end do                                                           8d15s06
       jsb=nsbeta(isb)                                                  3d22s17
       nbsum=0                                                          8d15s06
       do i=1,jsb-1                                                     3d22s17
        nbsum=nbsum+nbdet(i)                                            3d22s17
       end do                                                           8d15s06
       do ib=1,nbdet(jsb)                                               3d22s17
        ibp=ib+nbsum                                                    3d22s17
        do ia=0,nherec(isb)-1
         iap=ia+ilc(isb)+nasum                                          3d22s17
         isame=0
         do i1=1,nalpha
          do i2=1,nbeta
           if(iborb(i2,ibp).eq.iaorb(i1,iap))then
            isame=isame+1
            go to 10
           end if
          end do
   10     continue
         end do
         ioff=ioff+1
         sdig(ioff)=prevalue-dfloat(isame)                              3d22s17
        end do                                                          3d22s17
       end do                                                           3d22s17
      end do                                                            3d22s17
      return
      end
