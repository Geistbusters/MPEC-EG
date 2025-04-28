c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine delim(line,i1,ie)                                      6d1s21
      character*(*) line
c
c     search for blank or comma delimited field in line starting at i1.
c     i1 will be set to start of field and ie will be set to
c     end of field. If rest of line is blank, ie=0.
c     # marks the end of the line, with the rest comments
c
      ienew=len(line)                                                   12d2s22
      do i=len(line),i1,-1                                              12d2s22
       if(line(i:i).eq.'#')then                                         12d2s22
        ienew=i-1                                                       12d2s22
        go to 2                                                         12d2s22
       end if                                                           12d2s22
      end do                                                            12d2s22
    2 continue                                                          12d2s22
      do i=i1,ienew                                                     12d2s22
       if(line(i:i).ne.' '.and.line(i:i).ne.','.and.line(i:i).ne.'('    4d10s19
     $      .and.line(i:i).ne.')'.and.line(i:i).ne.'=')then             4d10s19
        is=i
        go to 1
       end if
      end do
      ie=0
      return
    1 continue
      i1=is
      do i=is+1,ienew                                                   12d2s22
       if(line(i:i).eq.' '.or.line(i:i).eq.','.or.line(i:i).eq.'('.or.  4d10s19
     $      line(i:i).eq.')'.or.line(i:i).eq.'=')then                   4d10s19
        ie=i-1
        return
       end if
      end do
      ie=ienew                                                          12d2s22
      return
      end
