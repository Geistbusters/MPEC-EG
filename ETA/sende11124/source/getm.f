c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine getm(sym,natom,xmass,inuc,lfn)
      implicit real*8 (a-h,o-z)
      character*(*) lfn
      character*4 sym(natom)
      character*40 line
      parameter (ida=5)
      dimension ns(ida),xmass(natom),ipos(ida),iz(ida),iiu(ida),        7d28s03
     $     prob(ida),igrab(ida)                                         7d28s03
      data ifirst/0/                                                    4d12s19
      save                                                              4d12s19
c
c     given symbol for element followed by mass number, e.g. H1, O16,
c     obtain atomic or nuclear mass from NIST table.
c     also can give mass number first, e.g. 1H, 16O.                    5d17s11
c
      if(ifirst.eq.0)then                                               4d12s19
       write(6,*)('in getm to get atom masses from file atwgts')
       if(inuc.eq.0)then
        write(6,*)('return atomic masses')
       else
        write(6,*)('return nuclear masses ')
       end if
      end if                                                            4d12s19
      if(natom.gt.ida)then
       write(6,*)('in getm, natom = '),i5
       write(6,*)('which is larger than ida ')
       stop 'ida'
      end if
      do 10 i=1,natom
       iiu(i)=1                                                         7d28s03
c
c     isym is where atomic symbol starts in string
c
       if(ichar(sym(i)(1:1)).ge.48.and.ichar(sym(i)(1:1)).le.57)then    5d17s11
c
c     mass number is first
c
        if(ichar(sym(i)(2:2)).ge.48.and.ichar(sym(i)(2:2)).le.57)then   5d17s11
         isym=3                                                         5d18s11
        else                                                            5d17s11
         isym=2                                                         5d18s11
        end if                                                          5d17s11
       else                                                             5d17s11
c
c     mass number is second                                             5d17s11
c
        isym=1                                                          5d17s11
       end if                                                           5d17s11
       isymp=isym+1                                                     5d17s11
       if((ichar(sym(i)(isymp:isymp)).ge.48.and.                        5d17s11
     $      ichar(sym(i)(isymp:isymp)).le.57).or.                       5d17s11
     $      sym(i)(isymp:isymp).eq.' ')then                             5d17s11
        ns(i)=1                                                         5d17s11
       else
        ns(i)=2
       end if
       if(isym.eq.1)then                                                5d17s11
        if(sym(i)(ns(i)+1:ns(i)+1).eq.' ')iiu(i)=0                       7d28s03
       end if                                                           5d17s11
       igrab(i)=0                                                       7d28s03
       ipos(i)=0
       xmass(i)=0d0
       prob(i)=0d0                                                      7d28s03
       iz(i)=0                                                          7d28s03
   10 continue
      open(unit=99,file=lfn,form='formatted')
      ihit=0
      izcount=0
    1 continue
       read(99,2,end=3)line
    2  format(a40)                                                      7d28s03
       if(line(1:3).eq.'___')then
        do 20 i=1,natom
         ipos(i)=0
   20   continue
        izcount=izcount+1                                               7d30s03
       end if
       do 4 i=1,natom
        if(line(1:3).eq.'___')then                                      7d30s03
        if(igrab(i).eq.1.and.iz(i).lt.izcount)then                      7d30s03
         igrab(i)=0                                                     7d28s03
         ihit=ihit+1                                                    7d28s03
        end if                                                          7d28s03
        end if                                                          7d28s03
        if(line(5:4+ns(i)).eq.sym(i)(isym:isym+ns(i)-1).or.             5d18s11
     $       ipos(i).eq.1)then                                          5d17s11
         if(line(1:1).ne.' ')then
          read(line(1:3),*)iznuc
         end if
         if(ns(i).eq.2.or.line(6:6).eq.' ')then
          ipos(i)=1
          if(iiu(i).eq.1)then
           if(isym.eq.1)then                                            5d17s11
            n1=ns(i)+1                                                  5d17s11
            n2=ns(i)+2                                                  5d17s11
           else                                                         5d17s11
            n1=1                                                        5d17s11
            n2=isym-1                                                   5d17s11
           end if                                                       5d17s11
           m2=9+n2-n1                                                   5d17s11
           if(line(9:m2).eq.sym(i)(n1:n2))then                          5d17s11
            do 5 j=20,40
             if(line(j:j).eq.'(')then
              read(line(12:j-1),*)xmass(i)
              iz(i)=iznuc
              ipos(i)=0
              ihit=ihit+1
              go to 6
             end if
    5       continue
           end if
          else                                                          7d28s03
           do 55 j=20,40                                                7d28s03
            if(line(j:j).eq.'(')then                                    7d28s03
             read(line(12:j-1),*)xm                                     7d28s03
             is=j+1                                                     7d28s03
             go to 57                                                   7d28s03
            end if                                                      7d28s03
   55      continue                                                     7d28s03
           go to 56                                                     7d28s03
   57      continue                                                     7d28s03
            do 58 j=is,40                                               7d28s03
             if(line(j:j).eq.')')then                                   7d28s03
              ie=j+1                                                    7d28s03
              go to 59                                                  7d28s03
             end if                                                     7d28s03
   58       continue                                                    7d28s03
            write(6,*)('can''t find ) in getm')                         4d12s19
            stop                                                        7d28s03
   59       continue                                                    7d28s03
            nnz=0                                                       7d28s03
            do 60 j=ie,40                                               7d28s03
             if(line(j:j).ne.' '.and.line(j:j).ne.'#')nnz=nnz+1         7d31s03
             if(line(j:j).ne.' ')nnz=nnz+1                              7d28s03
             if(line(j:j).eq.'(')then                                   7d28s03
              read(line(ie:j-1),*)frac                                  7d28s03
              igrab(i)=1                                                7d28s03
              if(frac.gt.prob(i))then                                   7d28s03
               prob(i)=frac                                             7d28s03
               xmass(i)=xm                                              7d28s03
               iz(i)=iznuc                                              7d28s03
              end if                                                    7d28s03
              go to 56                                                  7d28s03
             end if                                                     7d28s03
   60       continue                                                    7d28s03
            if(nnz.ne.0)then                                            7d28s03
             read(line(ie:40),*)frac                                    7d28s03
             prob(i)=frac                                               7d28s03
             xmass(i)=xm                                                7d28s03
             iz(i)=iznuc                                                7d28s03
             ihit=ihit+1                                                7d28s03
            end if                                                      7d28s03
   56      continue                                                     7d28s03
          end if
         end if
        end if
    6  continue
    4 continue
       if(ihit.lt.natom)then
        go to 1
       else
        go to 30
       end if
    3 continue
      do 31 i=1,natom
       if(xmass(i).eq.0d0)then
        write(6,*)('could not find mass data for "'),sym(i),('"')       4d12s19
        is=1                                                            1d2s24
        call delim(sym(i),is,ie)                                        1d2s24
        if(sym(i)(is:is).eq.'0'.or.sym(i)(ie:ie).eq.'0')then            1d2s24
         write(6,*)('take this to be dummy atom ')                      1d2s24
        end if                                                          1d2s24
       end if
   31 continue
   30 continue
      xme=5.485 799 110 d-4
      if(ifirst.eq.0)write(6,*)('mass of electron in u '),xme           4d12s19
      xmei=1d0/xme
      write(6,*)('symbol          atomic mass in u'),
     $     ('     atomic mass in au    nuclear mass in au')
      do 11 i=1,natom
       au=xmass(i)*xmei
       aun=au-dfloat(iz(i))
       write(6,12)sym(i),xmass(i),au,aun
   12  format(2x,a4,5x,3f22.14)
       if(inuc.eq.0)then
        xmass(i)=au
       else
        xmass(i)=aun                                                    7d26s00
       end if
   11 continue
      close(unit=99)                                                    7d26s00
      ifirst=1                                                          4d12s19
      return
      end
