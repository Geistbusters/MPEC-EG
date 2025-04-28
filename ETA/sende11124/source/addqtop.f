c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine addqtop(vec,ncsft,nroot,vecq,ncsfq,ipointq,myfcn,ibas, 3d23s22
     $     mdon,ncsf,bc,ibc)                                            11d14s22
      implicit real*8 (a-h,o-z)
      dimension vec(ncsft,*),vecq(ncsfq,*),ipointq(*),ibas(3,*),ncsf(*) 3d23s22
      include "common.store"
      do ir=1,nroot                                                     3d23s22
       iad1=1                                                           3d23s22
       iad2=1                                                           3d23s22
       do if=1,myfcn
        iarg=ibas(1,if)+1-mdon
        if(ipointq(if).ne.0)then
         do i=0,ncsf(iarg)-1
          vec(iad1+i,ir)=vecq(iad2+i,ir)
         end do
         iad2=iad2+ncsf(iarg)
        end if
        iad1=iad1+ncsf(iarg)
       end do
      end do
      return
      end
