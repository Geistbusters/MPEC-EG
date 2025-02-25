c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mostpop(anuc,atnum,natom,iretrn,pwd)
      implicit real*8 (a-h,o-z)
      character*80 line
      character*(*) pwd
c
c     return nuclear mass of most likely isotope for
c     atoms.
c
      dimension anuc(natom),atnum(3,natom)
      iretrn=0
      nok=0
      open(unit=1,file=pwd//'/atwgts')
      ns=0
      nf=0
      natu=0                                                            1d3s24
      do i=1,natom                                                      1d3s24
       if(atnum(1,i).ne.0d0)natu=natu+1                                 1d3s24
      end do                                                            1d3s24
    1 continue
      read(1,2,end=3)line
       is=1
    2  format(a80)
       if(line(1:5).eq.'_____')then
        if(nf.eq.1)then
         do i=1,natom
          if(abs(atnum(1,i)-dfloat(nzh)).lt.1d-5)then
           anuc(i)=wgtx
           nok=nok+1
          end if
         end do
        end if
        if(nok.eq.natu)then                                             1d3s24
         close(unit=1)
         return
        end if
        nf=max(nf,1)
        ns=0
        popx=0d0
        go to 1
       end if
       if(nf.eq.1)then
        if(ns.eq.0)then
         is=1
         call delim(line,is,ie)
         read(line(is:ie),*)nzh
         ns=1
         is=ie+1
         call delim(line,is,ie)
         is=ie+1
        else if(nzh.eq.1)then
c
c     atomic symbol if hydrogen
c
         call delim(line,is,ie)
         is=ie+1
        end if
c
c     number of neutrons and protons
c
        call delim(line,is,ie)
        is=ie+1
c
c     atomic mass in amu
c
        call delim(line,is,ie)
c
c     find ( marking start of uncertainty
c
        do i=is,ie
         if(line(i:i).eq.'(')then
          read(line(is:i-1),*)wgt
          go to 22
         end if
        end do
        read(line(is:ie),*)wgt
   22   continue
c
c     nuclear mass in amu
c
        wgt=wgt-dfloat(nzh)*5.485 799 110 d-4
c
c     population or blank if unstable
c
        is=ie+1
        call delim(line,is,ie)
        if(ie.eq.0)go to 1
c
c     find ( marking start of uncertainty
c
       do i=is,ie
        if(line(i:i).eq.'(')then
         read(line(is:i-1),*)pop
         if(pop.gt.popx)then
          popx=pop
          wgtx=wgt
         end if
         go to 1
        end if
       end do
       read(line(is:ie),*)pop
       if(pop.gt.popx)then
        popx=pop
        wgtx=wgt
       end if
      end if
      go to 1
    3 continue
      write(6,*)('could not find nuclear mass data for all atoms')
      iretrn=1
      stop
      end
