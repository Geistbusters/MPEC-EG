c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine enough(strng,bc,ibc)                                   11d28s22
      implicit real*8 (a-h,o-z)
      character*(*) strng                                               1d31s21
      include "common.store"                                            8d12s22
      parameter (id=10)
      dimension irecent(id)
      character*20 label(id)
      save                                                              10d8s15
      ienough=ienough+1
      if(ienough.le.id)then
       irecent(ienough)=ibcoff+1-ioffsx                                 10d31s22
       label(ienough)=memneed
      else
       do i=1,id-1
        irecent(i)=irecent(i+1)
        label(i)=label(i+1)
       end do
       irecent(id)=ibcoff+1-ioffsx                                      10d31s22
       label(id)=memneed
      end if
      memneed='                    '
      if(ibcoff-ioffsx.gt.maxbct)then                                   10d31s22
       write(6,*)('in enough, ibcoff =  '),ibcoff-ioffsx,ibcoff,ioffsx                10d31s22
       write(6,*)('which exceeds available memory '),maxbct             1d31s21
       write(6,*)('by factor of '),dfloat(ibcoff-ioffsx)/dfloat(maxbct) 10d31s22
       write(6,*)('ienough = '),ienough
       write(6,*)('len '),len(strng)
        write(6,*)('enough for: '),strng(1:min(20,len(strng)))          4d1s22
       write(6,*)('recent memory marks ')
       if(ienough.lt.id)then
        do i=1,ienough
         write(6,1)i,irecent(i),label(i)
    1    format(i2,i10,1x,a20)
        end do
       else
        do i=1,id
         write(6,1)i,irecent(i),label(i)
        end do
       end if
       stop
      else if(ibcoff-ioffsx.lt.0)then                                   10d31s22
       write(6,*)('in enough, ibcoff =  '),ibcoff-ioffsx,ibcoff,ioffsx                10d31s22
       write(6,*)('len '),len(strng)
        write(6,*)('enough for: '),strng(1:min(20,len(strng)))          4d1s22
       stop                                                             9d26s20
      end if
      if(ibcoff-ioffsx.gt.ihighwtr)ihighwtr=ibcoff-ioffsx               10d31s22
      return
      end
