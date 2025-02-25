c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dcbit(istring,norb,type)
      character*(*) type
      dimension iocc(64)
      integer*8 istring
      ii=0
      do i=1,norb
       if(btest(istring,i))then
        ii=ii+1
        iocc(ii)=i
       end if
      end do
      write(6,*)('for '),type,(iocc(i),i=1,ii)
      return
      end
