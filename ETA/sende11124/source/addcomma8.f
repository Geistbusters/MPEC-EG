c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine addcomma8(number,string,njust)                         2d2s21
      character*(*)string                                               2d2s21
      integer*8 number                                                  2d2s21
      character*40 line
      write(string,*)number                                             2d2s21
      is=1                                                              2d2s21
      call delim(string,is,ie)                                          2d2s21
      ndig=ie+1-is                                                      2d2s21
      ncomma=ndig/3                                                     2d2s21
      nlead=ndig-ncomma*3                                               2d2s21
      ngot=0                                                            2d2s21
      ip=is                                                             2d2s21
      do i=1,nlead                                                      2d2s21
       ngot=ngot+1                                                      2d2s21
       line(ngot:ngot)=string(ip:ip)                                    2d2s21
       ip=ip+1                                                          2d2s21
      end do                                                            2d2s21
    1 continue                                                          2d2s21
      if(ip.lt.ie)then                                                  2d2s21
       if(ngot.gt.0)then                                                2d2s21
        ngot=ngot+1                                                      2d2s21
        line(ngot:ngot)=','                                              2d2s21
       end if                                                           2d2s21
       do i=1,3                                                         2d2s21
        ngot=ngot+1                                                     2d2s21
        line(ngot:ngot)=string(ip:ip)                                   2d2s21
        ip=ip+1                                                         2d2s21
       end do                                                           2d2s21
       go to 1                                                          2d2s21
      end if                                                            2d2s21
      if(ngot.gt.len(string))then                                       2d2s21
       write(6,*)('in addcomma8, no. characters = '),ngot               2d2s21
       write(6,*)('exceeds length of input string = '),len(string)      2d2s21
       stop 'addcomma8'                                                 2d2s21
      end if                                                            2d2s21
      if(ngot.lt.njust)then                                             2d2s21
       nlead=njust-ngot                                                 2d2s21
       do i=1,nlead                                                     2d2s21
        string(i:i)=' '                                                 2d2s21
       end do                                                           2d2s21
       is=nlead+1                                                       2d2s21
       ie=nlead+ngot                                                    2d2s21
      else                                                              2d2s21
       is=1                                                             2d2s21
       ie=ngot                                                          2d2s21
      end if                                                            2d2s21
      string(is:ie)=line(1:ngot)                                        2d2s21
      do i=ie+1,len(string)                                             2d2s21
       string(i:i)=' '                                                  2d2s21
      end do                                                            2d2s21
      return                                                            2d2s21
      end                                                               2d2s21
