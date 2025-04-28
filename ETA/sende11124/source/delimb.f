c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine delimb(line,i1,ie)                                      6d1s21
      character*(*) line
c
c     search for blank or comma delimited field in line starting at i1.
c     i1 will be set to start of field and ie will be set to
c     end of field. If rest of line is blank, ie=0.
c
      do i=i1,len(line)
       if(line(i:i).ne.' ')then                                         6d1s21
        is=i
        go to 1
       end if
      end do
      ie=0
      return
    1 continue
      i1=is
      do i=is+1,len(line)
       if(line(i:i).eq.' ')then                                         6d1s21
        ie=i-1
        return
       end if
      end do
      ie=len(line)
      return
      end
