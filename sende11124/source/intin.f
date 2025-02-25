c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine intin(line,is,ie,ival)
      character*(*) line
c
c     search in line from is to ie for integer.
c     if found, store in ival.
c
      do i=is,ie
       if(line(i:i).eq.' ')then
        i0=i
        go to 1
       end if
      end do
      return
    1 continue
      do i=i0+1,ie
       if(line(i:i).ne.' ')then
        i1=i
        go to 2
       end if
      end do
      return
    2 continue
      do i=i1+1,ie
       if(line(i:i).eq.' ')then
        i2=i-1
        go to 3
       end if
      end do
      i2=ie
    3 continue
c
c     we now have blank delimited string from i1 to i2
c
      read(line(i1:i2),*,err=4)ival
      return
    4 continue
      write(6,*)('error reading in intin '),i1,i2
      write(6,*)line(1:is)
      write(6,*)('"'),line(i1:i2),('"')
      stop
      end
