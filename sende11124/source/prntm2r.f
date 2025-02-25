c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prntm2r(a,ni,nj,id)                                    1d2s23
      implicit real  *8 (a-h,o-z)
c***********************************************************************
c     this subroutine prints the elements of the ni rows and nj columns
c     of array a.  a is dimensioned (id,id) in the calling program
c     dimension statement.
c***********************************************************************
      dimension a(id,*)                                                 1d11s23
      character*132 line
      common/unitcm/iunit                                               11d9s17
      include "common.print"
      if(iunit.gt.0.and.iunit.ne.5.and.iunit.lt.100)then                11d9s17
       iunitu=iunit                                                     11d9s17
      else                                                              11d9s17
       iunitu=6                                                         11d9s17
      end if                                                            11d9s17
      if(iabs(nj).gt.40000)then                                          3d3s21
       write(6,*)('we have a wild column no.: '),nj,loc(nj)
       call dws_synca
       call dws_finalize
       stop
      end if
c     numcol=9                                                          9d13s88
      numcol=5                                                          9d13s88
      jcol=numcol                                                       9d13s88
      jmin=1
      itmp=-nj                                                          2d11s89
      if(iprtr(27).eq.0)then
      tiny=1d-11                                                        3d25s20
      else                                                              5d31s22
       tiny=0d0
      end if                                                            5d31s22
      if(iunitu.eq.6.and.tiny.ne.0d0)then                               6d3s22
       write(6,1066)tiny                                                3d16s20
 1066 format('tiny = ',es7.0)                                           3d16s20
      end if                                                            11d9s17
      if(ni.ge.0)then
      write(iunitu,120)ni,itmp                                          11d9s17
  120 format(1x,2i6,'@')                                                2d11s89
   20 jmax=min(jcol,nj)
      is=6                                                              1d2s23
      line(1:is)='  m-n'                                                1d2s23
      is=7                                                              1d2s23
      do im=jmin,jmax                                                   1d2s23
       if(mod(im,2).eq.0)then                                           1d2s23
        k=2                                                             1d2s23
        ii=im/2
       else
        k=1
        ii=(im+1)/2
       end if
       write(line(is:is+13),30)ii,k                                     1d2s23
   30  format(1x,i5,'k',i1,6x)                                          1d2s23
       is=is+14                                                         1d2s23
      end do                                                            1d2s23
      write(iunitu,64)(line(j:j),j=1,is-1)                              1d2s23
      do 50 i=1,ni
       write(line(1:6),60)i                                             9d1s94
   60  format(i5,1x)                                                    4d10s20
       is=7                                                             9d1s94
       do 61 j=jmin,jmax                                                9d1s94
        if(iprtr(27).eq.0)then                                          5d19s22
         ie=is+13                                                        9d1s94
        else                                                            5d19s22
         ie=is+21                                                        9d1s94
        end if                                                          5d19s22
        if(abs(a(i,j)).lt.tiny)then                                     9d1s94
         write(line(is:ie),63)                                          9d1s94
   63    format(7x,'tiny ')                                             9d1s94
        else if(abs(a(i,j)).ge.tiny)then                                     9d1s94
         if(iprtr(27).eq.0)then                                         5d19s22
          write(line(is:ie),62)a(i,j)                                    9d1s94
         else                                                           5d19s22
          write(line(is:ie),662)a(i,j)                                    9d1s94
         end if                                                         5d19s22
   62    format(1pe14.6)                                                9d1s94
  662    format(1pe22.14)                                               5d19s22
        else                                                            9d1s94
         write(line(is:ie),1063)                                        1d5s99
 1063    format(7x,'NaNQ')                                              1d5s99
        end if                                                          9d1s94
        is=ie+1                                                         9d1s94
   61  continue                                                         9d1s94
       write(iunitu,64)(line(j:j),j=1,is-1)                             11d9s17
   64  format(132a1)                                                    9d1s94
   50 continue
      if(jmax.eq.nj)then                                                2d11s89
       write(iunitu,121)                                                11d9s17
  121  format(1x,'$')                                                   2d11s89
       return                                                           2d11s89
      end if                                                            2d11s89
      jmin=jmax+1
      jcol=jcol+numcol                                                  9d13s88
      go to 20
      else
       write(6,*)('printing transpose of matrix ')
       write(6,120)nj,-ni                                               12d12s94
  720  jmax=min(jcol,ni)                                                12d12s94
       write(6,30)(im,im=jmin,jmax)
       do 750 i=1,nj                                                    12d12s94
        write(line(1:6),60)i                                            12d12s94
        is=7                                                             9d1s94
        do 761 j=jmin,jmax                                                9d1s94
         ie=is+13                                                        9d1s94
         if(abs(a(j,i)).gt.tiny)then                                     9d1s94
          write(line(is:ie),62)a(j,i)                                    9d1s94
         else if(abs(a(j,i)).le.tiny)then                               1d5s99
          write(line(is:ie),63)                                          9d1s94
         else                                                            9d1s94
          write(line(is:ie),1063)                                       1d5s99
         end if                                                          9d1s94
         is=ie+1                                                         9d1s94
  761   continue                                                         9d1s94
        write(6,64)(line(j:j),j=1,is-1)                                  9d1s94
  750  continue
       if(jmax.eq.ni)then                                                2d11s89
        write(6,121)                                                     2d11s89
        return                                                           2d11s89
       end if                                                            2d11s89
       jmin=jmax+1
       jcol=jcol+numcol                                                  9d13s88
      go to 720                                                         12d12s94
      end if
      end
