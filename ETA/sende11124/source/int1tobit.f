c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine int1tobit(iptr,iptrbit,idorb,isorb,idbit,isbit,
     $     mdoop,nsymb,nec,bc,ibc)                                      11d10s22
      implicit real*8 (a-h,o-z)                                         11d1s22
      integer*1 idorb(*),isorb(*)
      dimension iptr(4,mdoop,*),iptrbit(2,mdoop,*)                      5d7s20
      include "common.store"
      ndq=0
      nsq=0
      do ipass=1,2
       do isb=1,nsymb
        do i=1,mdoop                                                    5d7s20
         im=i-1
         nclo=i-1
         nopen=nec-2*nclo
         ic=iptr(2,i,isb)
         io=iptr(4,i,isb)
         if(nclo.ge.0.and.ipass.eq.2)then
          iptrbit(1,i,isb)=jdbit                                        5d7s20
          do j=1,iptr(1,i,isb)
           ibc(jdbit)=0
           do k=0,nclo-1
            ibit=idorb(ic+k)
            ibc(jdbit)=ibset(ibc(jdbit),ibit)
           end do
           jdbit=jdbit+1
           ic=ic+nclo
          end do
         end if
         if(nopen.ge.0.and.ipass.eq.2)then
          iptrbit(2,i,isb)=jsbit                                        5d7s20
          do j=1,iptr(3,i,isb)
           ibc(jsbit)=0
           do k=0,nopen-1
            ibit=isorb(io+k)
            ibc(jsbit)=ibset(ibc(jsbit),ibit)
           end do
           jsbit=jsbit+1
           io=io+nopen
          end do
         end if
         if(ipass.eq.1)then
          ndq=ndq+iptr(1,i,isb)
          nsq=nsq+iptr(3,i,isb)
         end if
        end do
       end do
       if(ipass.eq.1)then
        idbit=ibcoff
        isbit=idbit+ndq
        ibcoff=isbit+nsq
        call enough('int1tobit.  1',bc,ibc)
        jdbit=idbit
        jsbit=isbit
       end if
      end do                                                            5d7s20
      return                                                            5d7s20
      end                                                               5d7s20
