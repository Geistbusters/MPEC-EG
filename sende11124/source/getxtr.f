c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine getxtr(iunit,elem,nxtr,xtra,ltra,mxtr,icvra,icv)       12d13s22
      implicit real*8 (a-h,o-z)                                         12d12s19
c
c     scan through file produced by diff to get extra uncontracted      12d12s19
c     functions.                                                        12d12s19
c
      character*(2) elem                                                12d12s19
      integer*8 ltra(mxtr),icvra(*)                                     12d13s22
      dimension xtra(mxtr)                                              12d12s19
      character*20 line                                                 12d12s19
    1 continue
       read(iunit,2,end=3)line
    2  format(a20)
       if(line(1:1).eq.'<')then                                         12d12s19
        is=2                                                            12d12s19
        call delim(line,is,ie)                                          12d12s19
        nh=ie+1-is                                                      12d12s19
        if(nh.gt.2)go to 1                                              12d12s19
        if(line(is:is).ne.elem(1:1))go to 1                             12d12s19
        if(line(is:is+1).ne.elem)go to 1                                12d12s19
        is=ie+1                                                         12d12s19
        call delim(line,is,ie)                                          12d12s19
        if(line(is:is).eq.'S')then
         lgot=0
        else if(line(is:is).eq.'P')then
         lgot=1
        else if(line(is:is).eq.'D')then
         lgot=2
        else if(line(is:is).eq.'F')then
         lgot=3
        else if(line(is:is).eq.'G')then
         lgot=4
        else if(line(is:is).eq.'H')then
         lgot=5
        else if(line(is:is).eq.'I')then
         lgot=6
        else if(line(is:is).eq.'K')then                                 12d12s19
         lgot=7
        else                                                            12d12s19
         write(6,*)('don''t know how to parse "'),line(is:ie),
     $        ('" for l !!!')
         stop
        end if                                                          12d12s19
        read(iunit,2,end=3)line
        read(line(2:20),*)expn
        nxtr=nxtr+1                                                     12d12s19
        if(nxtr.gt.mxtr)then                                            12d12s19
         write(6,*)('toooo many extra fcns: '),nxtr,mxtr                12d12s19
         stop                                                           12d12s19
        end if                                                          12d12s19
        ltra(nxtr)=lgot                                                 12d12s19
        xtra(nxtr)=expn                                                 12d12s19
        icvra(nxtr)=icv                                                 12d13s22
       end if                                                           12d12s19
       go to 1
    3  continue
       return
       end
