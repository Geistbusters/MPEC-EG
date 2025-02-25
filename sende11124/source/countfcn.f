c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine countfcn(nff,mdon,mdoop,nfcnf)                         3d18s22
      implicit real*8 (a-h,o-z)                                         3d18s22
      dimension nff(mdoop)
      write(6,*)('in countfcn ')
      nfcnf=0
      do i=mdon+1,mdoop
       write(6,*)('for nclop '),i,nff(i)
       nfcnf=nfcnf+nff(i)
      end do
      write(6,*)('final nfcnf '),nfcnf
      return
      end
