c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine oplist(id,name,ipt,npt,data,idata,iosym,iprt,i2e,nop,  12d18s19
     $     opnc,lzzonly)                                                12d20s19
      implicit real*8 (a-h,o-z)
c
c     set up rules for property operators
c     lzzonly = 2 means lz,lzz only.
c     lzzonly = 6 means lx,ly,lz,lxx,lyy,lzz only.
c
      logical lz2only                                                   12d20s19
      character*10 name(id)
      dimension ipt(id),npt(id),data(id),idata(7,id),iosym(id),iprt(id),5d27s21
     $     i2e(6,id),opnc(id),tmp(6)                                    5d20s21
      include "common.basis"
      include "common.input"
      ioff=1
c
c     charge of electron ...
c
      e=-1d0                                                            12d18s19
      nop=0                                                             12d20s19
      if(lzzonly.eq.2)go to 1                                           12d31s19
      if(lzzonly.eq.6)go to 3                                           12d31s19
c
c     mux
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=1
      iosym(nop)=ipropsym(1)                                             12d15s19
      name(nop)='mux       '
      i2e(1,nop)=0                                                      1d13s23
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(1,ioff)=1
      data(ioff)=e                                                      1d3s20
      ioff=ioff+1
      if(ioff.gt.id)stop 'id'
c
c     muy
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=1
      iosym(nop)=ipropsym(2)                                             12d15s19
      name(nop)='muy       '                                            1d3s20
      i2e(1,nop)=0
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(2,ioff)=1
      data(ioff)=e                                                      1d3s20
      ioff=ioff+1
      if(ioff.gt.id)stop 'id'
c
c     muz
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=1
      iosym(nop)=ipropsym(3)                                             12d15s19
      name(nop)='mu0       '
      i2e(1,nop)=0
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(3,ioff)=1
      data(ioff)=1d0*e                                                  12d18s19
      ioff=ioff+1
      if(ioff.gt.id)stop 'id'
c
c     lx/i
c
    3 continue                                                          12d31s19
      nop=nop+1                                                         12d20s19
      lxop=nop                                                          12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=2
      iosym(nop)=ipropsym(6)                                             12d15s19
      name(nop)='lx/i      '                                            1d3s20
      i2e(1,nop)=0
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
       idata(i,ioff+1)=0
      end do
      idata(3,ioff)=1
      idata(5,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(6,ioff)=1
      data(ioff)=-1d0
      ioff=ioff+1
c
c     ly/i
c
      nop=nop+1                                                         12d20s19
      lyop=nop                                                          12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=2
      iosym(nop)=ipropsym(5)                                             12d15s19
      name(nop)='ly/i      '                                            1d3s20
      i2e(1,nop)=0
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
       idata(i,ioff+1)=0
      end do
      idata(3,ioff)=1
      idata(4,ioff)=1
c     per edmonds 2.1.3                                                 5d19s21
      data(ioff)=-1d0                                                   5d19s21
      ioff=ioff+1
      idata(1,ioff)=1
      idata(6,ioff)=1
      data(ioff)=+1d0                                                   5d19s21
      ioff=ioff+1
c
c     lz/i
c
    1 continue                                                          12d20s19
      nop=nop+1                                                         12d20s19
      lzop=nop                                                          12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=2
      iosym(nop)=ipropsym(4)                                             12d15s19
      name(nop)='lz/i      '                                            1d3s20
      i2e(1,nop)=0
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
       idata(i,ioff+1)=0
      end do
      idata(1,ioff)=1
      idata(5,ioff)=1
c     per edmonds 2.1.3                                                 5d19s21
      data(ioff)=-1d0                                                   5d19s21
      ioff=ioff+1
      idata(2,ioff)=1
      idata(4,ioff)=1
      data(ioff)=+1d0                                                   5d19s21
      ioff=ioff+1
c
c     lz^2
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=5
      iosym(nop)=1
      name(nop)='Lz^2     '
      i2e(1,nop)=lzop                                                   12d20s19
      i2e(2,nop)=lzop                                                   12d20s19
      i2e(3,nop)=0                                                      5d20s21
      i2e(4,nop)=0                                                      5d20s21
      i2e(5,nop)=0                                                      5d20s21
      i2e(6,nop)=0                                                      5d20s21
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,7                                                         5d27s21
        idata(i,j)=0
       end do
      end do
      idata(2,ioff)=2
      idata(4,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(5,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(2,ioff)=1
      idata(4,ioff)=1
      idata(5,ioff)=1
      data(ioff)=2d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(4,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(1,ioff)=2
      idata(5,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1                                                       12d16s19
c
      if(lzzonly.eq.2)go to 2                                           12d31s19
      if(lzzonly.eq.6)go to 6                                           12d31s19
c
c     LxLz+LzLx
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=6
      iosym(nop)=ipropsym(5)                                            12d16s19
      name(nop)='LzLx+LzLx '
      i2e(1,nop)=lzop                                                   12d20s19
      i2e(2,nop)=lxop                                                   12d20s19
      i2e(3,nop)=0                                                      5d20s21
      i2e(4,nop)=0                                                      5d20s21
      i2e(5,nop)=0                                                      5d20s21
      i2e(6,nop)=0                                                      5d20s21
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+ipt(nop)-1
       do i=1,7                                                         5d27s21
        idata(i,j)=0
       end do
      end do
      idata(3,ioff)=1
      idata(4,ioff)=1
      data(ioff)=-1d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(3,ioff)=1
      idata(4,ioff)=1
      idata(5,ioff)=1
      data(ioff)=-2d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(3,ioff)=1
      idata(5,ioff)=2                                                   12d15s19
      data(ioff)=2d0
      ioff=ioff+1
      idata(2,ioff)=2
      idata(4,ioff)=1
      idata(6,ioff)=1
      data(ioff)=2d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(2,ioff)=1
      idata(5,ioff)=1
      idata(6,ioff)=1
      data(ioff)=-2d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(6,ioff)=1
      data(ioff)=-1d0
      ioff=ioff+1
c
c     Lx^2
c
    6 continue                                                          12d31s19
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=5
      iosym(nop)=1
      name(nop)='Lx^2      '
      i2e(1,nop)=lxop                                                   12d20s19
      i2e(2,nop)=lxop                                                   12d20s19
      i2e(3,nop)=0                                                      5d20s21
      i2e(4,nop)=0                                                      5d20s21
      i2e(5,nop)=0                                                      5d20s21
      i2e(6,nop)=0                                                      5d20s21
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,7                                                         5d27s21
        idata(i,j)=0
       end do
      end do
      idata(3,ioff)=2
      idata(5,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
      idata(3,ioff)=1
      idata(6,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(3,ioff)=1
      idata(5,ioff)=1
      idata(6,ioff)=1
      data(ioff)=2d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(5,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(2,ioff)=2
      idata(6,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
c
c
c     Ly^2
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=5
      iosym(nop)=1
      name(nop)='Ly^2      '
      i2e(1,nop)=lyop                                                   12d20s19
      i2e(2,nop)=lyop                                                   12d20s19
      i2e(3,nop)=0                                                      5d20s21
      i2e(4,nop)=0                                                      5d20s21
      i2e(5,nop)=0                                                      5d20s21
      i2e(6,nop)=0                                                      5d20s21
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,7                                                         5d27s21
        idata(i,j)=0
       end do
      end do
      idata(1,ioff)=2
      idata(6,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(4,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(3,ioff)=1
      idata(4,ioff)=1
      idata(6,ioff)=1
      data(ioff)=2d0
      ioff=ioff+1
      idata(3,ioff)=1
      idata(6,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(3,ioff)=2
      idata(4,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
      if(lzzonly.eq.6)go to 2                                           12d31s19
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=1
      iosym(nop)=1
      name(nop)='charge   '
      i2e(1,nop)=0
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      data(ioff)=-1d0
      ioff=ioff+1
c
c     quadrapole moments
c
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=2
      iosym(nop)=1
      name(nop)='Re Q2     '
      i2e(1,nop)=0
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,7                                                         5d27s21
        idata(i,j)=0
       end do
      end do
      idata(1,ioff)=2
      data(ioff)=0.5d0*e                                                12d18s19
      ioff=ioff+1
      idata(2,ioff)=2
      data(ioff)=-0.5d0*e                                               12d18s19
      ioff=ioff+1
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=1
      iosym(nop)=ipropsym(4)
      name(nop)='Im Q2    '
      i2e(1,nop)=0
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(1,ioff)=1
      idata(2,ioff)=1
      data(ioff)=1d0*e                                                  12d18s19
      ioff=ioff+1
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=1
      iosym(nop)=ipropsym(5)
      name(nop)='Re Q1    '
      i2e(1,nop)=0
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(1,ioff)=1
      idata(3,ioff)=1
      data(ioff)=-1d0*e                                                 12d18s19
      ioff=ioff+1
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=1
      iosym(nop)=ipropsym(6)
      name(nop)='Im Q1    '
      i2e(1,nop)=0
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(2,ioff)=1
      idata(3,ioff)=1
      data(ioff)=-1d0*e                                                 12d18s19
      ioff=ioff+1
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=3
      iosym(nop)=1
      name(nop)='Q0        '
      i2e(1,nop)=0
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,7                                                         5d27s21
        idata(i,j)=0
       end do
      end do
      idata(1,ioff)=2
      data(ioff)=-1d0*e/sqrt(6d0)                                       12d18s19
      ioff=ioff+1
      idata(2,ioff)=2
      data(ioff)=-1d0*e/sqrt(6d0)                                       12d18s19
      ioff=ioff+1
      idata(3,ioff)=2
      data(ioff)=2d0*e/sqrt(6d0)                                        12d18s19
      ioff=ioff+1
c
c     mass polarization                                                 5d20s21
c
      nop=nop+1                                                         5d20s21
      ipt(nop)=ioff                                                     5d20s21
      npt(nop)=1                                                        5d20s21
      iosym(nop)=ipropsym(1)                                            5d20s21
      name(nop)='d/dx      '
      iddx=nop                                                          5d20s21
      i2e(1,nop)=0                                                      1d13s23
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(4,ioff)=1
      data(ioff)=1d0                                                    5d20s21
      ioff=ioff+1
      if(ioff.gt.id)stop 'id'
      nop=nop+1                                                         5d20s21
      ipt(nop)=ioff                                                     5d20s21
      npt(nop)=1                                                        5d20s21
      iosym(nop)=ipropsym(2)                                            5d20s21
      name(nop)='d/dy      '
      iddy=nop                                                          5d20s21
      i2e(1,nop)=0                                                      1d13s23
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(5,ioff)=1
      data(ioff)=1d0                                                    5d20s21
      ioff=ioff+1
      if(ioff.gt.id)stop 'id'
      nop=nop+1                                                         5d20s21
      ipt(nop)=ioff                                                     5d20s21
      npt(nop)=1                                                        5d20s21
      iosym(nop)=ipropsym(3)                                            5d20s21
      name(nop)='d/dz      '
      iddz=nop                                                          5d20s21
      i2e(1,nop)=0                                                      1d13s23
      iprt(nop)=1
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
      end do
      idata(6,ioff)=1
      data(ioff)=1d0                                                    5d20s21
      ioff=ioff+1
      if(ioff.gt.id)stop 'id'
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=3
      iosym(nop)=1
      name(nop)='kin      '
      i2e(1,nop)=iddx                                                   5d20s21
      i2e(2,nop)=iddx                                                   5d20s21
      i2e(3,nop)=iddy                                                   5d20s21
      i2e(4,nop)=iddy                                                   5d20s21
      i2e(5,nop)=iddz                                                   5d20s21
      i2e(6,nop)=iddz                                                   5d20s21
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,7                                                         5d27s21
        idata(i,j)=0
       end do
      end do
      idata(4,ioff)=2
      data(ioff)=-0.5d0                                                 5d20s21
      ioff=ioff+1
      idata(5,ioff)=2
      data(ioff)=-0.5d0                                                 5d20s21
      ioff=ioff+1
      idata(6,ioff)=2
      data(ioff)=-0.5d0                                                 5d20s21
      ioff=ioff+1
c
c     angular momentum operators wrt cm of target or projectial.
c
      nop=nop+1                                                         12d20s19
      lxop=nop                                                          12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=2
      iosym(nop)=ipropsym(6)                                             12d15s19
      name(nop)='Ax/i      '                                            1d3s20
      i2e(1,nop)=0
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do i=1,7                                                          5d27s21
       idata(i,ioff)=0
       idata(i,ioff+1)=0
      end do
      idata(7,ioff)=1                                                   5d27s21
      idata(7,ioff+1)=1                                                 5d27s21
      idata(3,ioff)=1
      idata(5,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(6,ioff)=1
      data(ioff)=-1d0
      ioff=ioff+1
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=5
      iosym(nop)=1
      name(nop)='Ax^2      '                                            5d27s21
      i2e(1,nop)=lxop                                                   12d20s19
      i2e(2,nop)=lxop                                                   12d20s19
      i2e(3,nop)=0                                                      5d20s21
      i2e(4,nop)=0                                                      5d20s21
      i2e(5,nop)=0                                                      5d20s21
      i2e(6,nop)=0                                                      5d20s21
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,6
        idata(i,j)=0
       end do
       idata(7,j)=1                                                     5d27s21
      end do
      idata(3,ioff)=2
      idata(5,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
      idata(3,ioff)=1
      idata(6,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(3,ioff)=1
      idata(5,ioff)=1
      idata(6,ioff)=1
      data(ioff)=2d0
      ioff=ioff+1
      idata(2,ioff)=1
      idata(5,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(2,ioff)=2
      idata(6,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
      nop=nop+1                                                         12d20s19
      lyop=nop                                                          12d20s19
      ipt(nop)=ioff                                                       12d15s19
      npt(nop)=2
      iosym(nop)=ipropsym(5)                                             12d15s19
      name(nop)='Ay/i      '                                            5d27s21
      i2e(1,nop)=0
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do i=1,6
       idata(i,ioff)=0
       idata(i,ioff+1)=0
      end do
      idata(7,ioff)=1                                                   5d27s21
      idata(7,ioff+1)=1                                                 5d27s21
      idata(3,ioff)=1
      idata(4,ioff)=1
c     per edmonds 2.1.3                                                 5d19s21
      data(ioff)=-1d0                                                   5d19s21
      ioff=ioff+1
      idata(1,ioff)=1
      idata(6,ioff)=1
      data(ioff)=+1d0                                                   5d19s21
      ioff=ioff+1
      nop=nop+1                                                         12d20s19
      ipt(nop)=ioff
      npt(nop)=5
      iosym(nop)=1
      name(nop)='Ay^2      '
      i2e(1,nop)=lyop                                                   12d20s19
      i2e(2,nop)=lyop                                                   12d20s19
      i2e(3,nop)=0                                                      5d20s21
      i2e(4,nop)=0                                                      5d20s21
      i2e(5,nop)=0                                                      5d20s21
      i2e(6,nop)=0                                                      5d20s21
      iprt(nop)=1
      if(ioff+ipt(nop)-1.gt.id)stop 'id'
      do j=ioff,ioff+npt(nop)-1
       do i=1,6
        idata(i,j)=0
       end do
       idata(7,j)=1                                                     5d27s21
      end do
      idata(1,ioff)=2
      idata(6,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(4,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(1,ioff)=1
      idata(3,ioff)=1
      idata(4,ioff)=1
      idata(6,ioff)=1
      data(ioff)=2d0
      ioff=ioff+1
      idata(3,ioff)=1
      idata(6,ioff)=1
      data(ioff)=1d0
      ioff=ioff+1
      idata(3,ioff)=2
      idata(4,ioff)=2
      data(ioff)=-1d0
      ioff=ioff+1
    2 continue                                                          12d20s19
      do iop=1,nop
       sum=0d0                                                          12d19s19
       do ia=1,natom                                                    12d19s19
        zz=e*atnum(1,ia)                                                12d19s19
        do i=1,3
         tmp(i)=xcart(i,ia)                                             12d19s19
        end do                                                          12d19s19
        do i=4,6                                                        12d19s19
         tmp(i)=0d0                                                     12d19s19
        end do                                                          12d19s19
        do i=ipt(iop),ipt(iop)+npt(iop)-1                               12d19s19
         term=zz*data(i)                                                12d19s19
         do j=1,6
          if(idata(j,i).gt.0)then                                       12d19s19
           term=term*(tmp(j)**idata(j,i))                               12d19s19
          end if                                                        12d19s19
         end do                                                         12d19s19
         sum=sum+term                                                   12d19s19
        end do                                                          12d19s19
       end do                                                           12d19s19
       opnc(iop)=sum                                                    12d19s19
      end do
      return
      end
