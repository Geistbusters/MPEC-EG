c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pdump(x,n,ismult,wname,pname)                          3d30s22
      implicit real*8 (a-h,o-z)                                         3d30s22
      character*(*) wname                                               3d30s22
      character*(*) pname                                               3d30s22
      character*100 line                                                3d30s22
      dimension x(n,*)                                                  3d30s22
      nz=0                                                              3d31s22
      do ik=1,n                                                         3d30s22
       do ib=1,ik                                                       3d30s22
        if(abs(x(ib,ik)).gt.1d-10)then                                  3d30s22
         if(nz.eq.0)write(6,*)(' ')                                     3d31s22
         nz=nz+1                                                        3d31s22
         write(line(1:5),1)ib,ismult                                    3d30s22
    1    format('<',i2,x,i1)                                            3d30s22
         do i=1,len(wname)                                              3d30s22
          line(5+i:5+i)=wname(i:i)                                      3d30s22
         end do                                                         3d30s22
         is=5+len(wname)+1                                              3d30s22
         line(is:is)='|'                                                3d30s22
         in=is                                                          3d30s22
         is=in+len(pname)                                               3d30s22
         ie=is+13+14                                                    3d30s22
         if(ie.gt.100)then                                              3d30s22
          write(6,*)('pdump error!!!! ')                                3d30s22
          stop                                                          3d30s22
         end if                                                         3d30s22
         do i=1,len(pname)                                              3d30s22
          line(in+i:in+i)=pname(i:i)                                    3d30s22
         end do                                                         3d30s22
         in=in+len(pname)+1                                             3d30s22
         ie=in+4                                                        3d30s22
         write(line(in:ie),2)ik,ismult                                  3d30s22
    2    format('|',i2,x,i1)                                            3d30s22
         in=ie                                                          3d30s22
         do i=1,len(wname)                                              3d30s22
          line(in+i:in+i)=wname(i:i)                                    3d30s22
         end do                                                         3d30s22
         is=in+len(wname)+1                                             3d30s22
         ie=is+16                                                       3d30s22
         write(line(is:ie),4)x(ib,ik)                                   3d30s22
    4    format('>=',es15.7)                                            3d30s22
         write(6,*)(line(i:i),i=1,ie-15),x(ib,ik)
    3    format(100a1)                                                  3d30s22
        end if                                                          3d30s22
       end do                                                           3d30s22
      end do                                                            3d30s22
      return                                                            3d30s22
      end                                                               3d30s22
