c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prtsigs(ipack,npack,io,iaorb)
      integer*8 ipack(npack),itry,iaorb(io,npack)
      do i=1,npack
       itry=ipack(i)
       do j=1,io
        iaorb(j,i)=0
        itry=itry/2
        itry4=itry
        if(mod(itry4,2).ne.0)then
         iaorb(j,i)=-1
        else
         iaorb(j,i)=1
        end if
       end do
    1  format(i5,5x,64i3)
      end do
      return
      end
