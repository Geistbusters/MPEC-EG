c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine getcup(irw0,nopen,i2sb,i2smb,i2sk,i2smk,ntype,iout,    10d15s21
     $     imat,mcsf,imap,bc,ibc)                                       11d14s22
      implicit real*8 (a-h,o-z)
      logical ldebug,lcomp                                              3d23s22
      integer*8 ipackc,ipack8                                           10d15s21
      integer*4 ipack48(2)                                              10d15s21
      integer*1 ipackc1(8),imap(*)                                      10d15s21
      equivalence (ipackc,ipackc1),(ipack8,ipack48)                     10d15s21
      include "common.store"                                            10d14s21
      dimension mcsf(2)                                                 10d14s21
      data icall/0/                                                     3d23s22
      save icall                                                        3d23s22
      icall=icall+1
      ldebug=.false.                                                    3d23s22
      sigc=1d0                                                          10d15s21
      isigc=1                                                           10d15s21
      ntype=0                                                           10d15s21
      iout=ibcoff                                                       10d15s21
      if(iabs(i2sb-i2sk).le.2.and.iabs(i2smb-i2smk).le.2)then           10d15s21
       if(ldebug)write(6,*)('Hi, my name is getcup!')
       if(ldebug)write(6,*)('irw0! '),irw0,nopen                                  10d15s21
       if(i2smb.lt.0)then                                               10d15s21
        if(ldebug)
     $      write(6,*)('let us change signs on ms to use positive ms'),
     $              i2sb,i2smb,i2sk,i2smk
        if(i2sb.ne.i2sk)sigc=-sigc                                      10d15s21
        idelm=-(i2smb-i2smk)/2                                          10d15s21
        isigc=-1                                                        10d15s21
       else                                                             10d15s21
        idelm=(i2smb-i2smk)/2                                           10d15s21
       end if                                                           10d15s21
       idelm=idelm+1                                                    10d15s21
       idels=(i2sb-i2sk)/2                                              10d15s21
       idels=idels+1                                                    10d15s21
       ioffb=(i2sb-mod(i2sb,2))*(i2sb-mod(i2sb,2)+2)                    10d15s21
       ioffb=ioffb/8                                                    10d15s21
       if(i2smb.lt.0)then                                               10d15s21
        ioffb=ioffb+((mod(i2sb,2)-i2smb)/2)                             10d15s21
       else                                                             10d15s21
        ioffb=ioffb+((mod(i2sb,2)+i2smb)/2)                             10d15s21
       end if                                                           10d15s21
       ioff=irw0+(nopen-mod(nopen,2))/2                                 10d15s21
       iooa=ibc(ioff)                                                   10d15s21
       if(ldebug)write(6,*)('iooa = '),iooa                                       10d15s21
       ioff=iooa+idels+3*(idelm+3*ioffb)                                10d15s21
       idata=ibc(ioff)                                                  10d15s21
       if(ldebug)write(6,*)('my data starts at '),idata                           10d15s21
       ipack8=ibc(idata)                                                10d15s21
       if(ldebug)write(6,*)('ncsfs '),ipack48                                     10d15s21
       mcsf(1)=ipack48(1)                                               10d15s21
       mcsf(2)=ipack48(2)                                               10d15s21
       nn=mcsf(1)*mcsf(2)                                               10d15s21
       idata=idata+1                                                    10d15s21
       ntype=ibc(idata)                                                 10d15s21
       if(ldebug)write(6,*)('no. types '),ntype                                   10d15s21
       if(ntype.lt.0.or.ntype.gt.1000)then                              10d15s21
        write(6,*)('bad ntype! '),ntype                                 10d15s21
        stop 'getcup'                                                   10d15s21
       end if                                                           10d15s21
       idata=idata+1                                                    10d15s21
       imt=idata                                                        10d15s21
       imat=iout+ntype                                                  10d15s21
       ibcoff=imat+ntype*nn                                             10d15s21
       if(ibcoff.lt.0)write(6,*)('enoughd '),ibcoff
       call enough('getcup.  1',bc,ibc)
       do i=3,8                                                         10d15s21
        ipackc1(i)=0                                                    10d15s21
       end do                                                           10d15s21
       do i=0,ntype-1                                                   10d15s21
        ipack8=ibc(imt)                                                 10d15s21
        im1=imap(iabs(ipack48(1)))                                      10d15s21
        im2=imap(iabs(ipack48(2)))                                      10d15s21
        if(ldebug)write(6,*)('raw labels: '),ipack48,
     $       ('mapped w/o sign: '),im1,im2
        ipackc1(1)=ipack48(1)*isigc                                     10d15s21
        ipackc1(2)=ipack48(2)*isigc                                     10d15s21
        if(ldebug)write(6,*)('after isigc: '),ipackc1(1),ipackc1(2)
        if(ipackc1(1).ge.0)then                                         10d15s21
         ipackc1(1)=im1                                                 10d15s21
        else                                                            10d15s21
         ipackc1(1)=-im1                                                10d15s21
        end if                                                          10d15s21
        if(ipackc1(2).ge.0)then                                         10d15s21
         ipackc1(2)=im2                                                 10d15s21
        else                                                            10d15s21
         ipackc1(2)=-im2                                                10d15s21
        end if                                                          10d15s21
        if(ldebug)write(6,*)('after map: '),ipackc1(1),ipackc1(2)
        ibc(iout+i)=ipackc                                              10d15s21
        imt=imt+1                                                       10d15s21
        lcomp=ibc(imt).ne.0                                             3d23s22
        imt=imt+1                                                       3d23s22
        jmat=imat+i*nn                                                  10d15s21
        if(lcomp)then                                                   3d23s22
         nused=ibc(imt)                                                 3d23s22
         imt=imt+1                                                      3d23s22
         nusedi=nused/2                                                 3d23s22
         if(2*nusedi.ne.nused)nusedi=nusedi+1                           3d23s22
         icmp1=imt                                                      3d23s22
         imt=imt+mcsf(1)                                                3d23s22
         icmp2=imt                                                      3d23s22
         imt=imt+nusedi                                                 3d23s22
         icmp3=imt                                                      3d23s22
         imt=imt+nused                                                  3d23s22
         call uncompxut(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(jmat),       3d23s22
     $        mcsf(1),mcsf(2))                                          3d23s22
         if(sigc.ne.1d0)then                                            3d23s22
          do ii=0,nn-1                                                   3d23s22
           bc(jmat+ii)=bc(jmat+ii)*sigc                                  3d23s22
          end do                                                         3d23s22
         end if                                                         3d23s22
        else                                                            3d23s22
         do ii=0,nn-1                                                    10d15s21
          bc(jmat+ii)=bc(imt+ii)*sigc                                    10d15s21
         end do                                                          10d15s21
         imt=imt+nn                                                      10d15s21
        end if                                                          3d23s22
        if(ldebug)then                                                  10d15s21
         write(6,*)('for type '),i,ipackc1
         call prntm2(bc(jmat),mcsf(1),mcsf(2),mcsf(1))                   10d15s21
        end if
       end do                                                           10d15s21
      end if                                                            10d15s21
      return                                                            10d15s21
      end                                                               10d15s21
