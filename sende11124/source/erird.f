c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine erird(la,zza,xna,xa,ya,za,                              10d4s16
     $               lb,zzb,xnb,xb,yb,zb,                               10d4s16
     $               lc,zzc,xnc,xc,yc,zc,                               10d4s16
     $               ld,zzd,xnd,xd,yd,zd,                               10d4s16
     $               icartt,idoit,fmul,idxa,idya,idza,idxb,idyb,idzb,    3d3s17
     $               idxc,idyc,idzc,idxd,idyd,idzd,bc,ibc)              11d14s22
      implicit real*8 (a-h,o-z)
      character*120 line
      logical log1,log2,log3,log4
c
c     compute electron repulsion integrals via recursion and rys quad   10d4s16
c     for derivatives, non-relativistic case.
c     possible cases:
c     paraerid, one of idxa,idya,idza is 1 or 2.
c     paraeridj, one of idxc,idyc,idzc is 1
c     paraerid2b, one of the pairs idxa,idxb, idya,idyb, idza,idzb is 1,1.
c     paraerid2c, one of the pairs idxb,idxd, idyb,idyd, idzb,idzd is 1,1.
c
      include "common.store"
      include "common.spher"
      include "common.rys"                                              7d5s12
      dimension nhere(8,4),ibest(8),xbest(8,3),ybest(8,2),              1d3s17
     $     mhere(8,4),icoef(8,4),ipow(8,4),ltmp(4),joff(8,4),           12d14s16
     $     isymd2h(3,8),lcent(4),centd(4,4),iperm(4,8),jperm(4,2)       12d23s16
      data isymd2h/0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1,
     $     1,1,1/
      data iperm/1,2,3,4, 1,2,4,3, 2,1,3,4, 2,1,4,3,
     $           3,4,1,2, 3,4,2,1, 4,3,1,2, 4,3,2,1/
      data icall,ifail/0,0/                                                     7d22s14
      save                                                              7d22s14
      if(icall.eq.0)then                                                10d4s16
       pi=acos(-1d0)                                                    10d4s16
      end if                                                            10d4s16
      icall=icall+1                                                     7d22s14
      ida=max(idxa,idya,idza)                                            3d3s17
      idb=max(idxb,idyb,idzb)                                            3d3s17
      idc=max(idxc,idyc,idzc)                                            3d3s17
      idd=max(idxd,idyd,idzd)                                            3d3s17
      idbg=0
      if(idoit.eq.1)idbg=1
 1066 continue
      if(idoit.eq.1.or.idbg.ne.0)then                                   8d15s24
      write(6,1)la,zza,xa,ya,za,lb,zzb,xb,yb,zb,lc,zzc,xc,yc,zc,        1d29s10
     $      ld,zzd,xd,yd,zd,icall,idxa,idya,idza,idxb,idyb,idzb,         3d3s17
     $      idxc,idyc,idzc,idxd,idyd,idzd                                3d3s17
      idbg=1
      end if
    1 format('in erird: ',4(i3,1pe15.7,0p3f8.4),i8,/4(5x,3i2))             3d3s17
      ltmp(1)=la+1                                                      1d13s10
      ltmp(2)=lb+1                                                      1d13s10
      ltmp(3)=lc+1                                                      1d13s10
      ltmp(4)=ld+1                                                      1d13s10
      do i=1,4                                                          1d13s10
       koff=-1                                                          1d14s10
       do isym=1,8                                                      1d13s10
        joff(isym,i)=koff                                               1d14s10
        nhere(isym,i)=0                                                 1d13s10
        mhere(isym,i)=0                                                 1d14s10
        if(ipt(isym,ltmp(i)).gt.0)then                                  1d13s10
         isto=ipt(isym,ltmp(i))                                         1d13s10
         nhere(isym,i)=ibc(isto)                                        1d13s10
         mhere(isym,i)=ibc(isto+1)                                      1d13s10
         icoef(isym,i)=isto+2+mhere(isym,i)+3*nhere(isym,i)             1d13s10
         ipowb=isto+2+ibc(isto+1)                                       1d13s10
         koff=koff+mhere(isym,i)                                        1d14s10
        end if                                                          1d13s10
       end do                                                           1d13s10
      end do                                                            1d13s10
      irnode=ibcoff                                                     1d7s09
      ls=la+lb+lc+ld                                                    11d30s16
      nqr=(ls+2+ida+idb+idc+idd)/2                                       3d3s17
      if(mod(ls+2+ida+idb+idc+idd,2).ne.0)nqr=nqr+1                     8d14s24
      irwgt=irnode+nqr                                                   1d7s09
      ibcoff=irwgt+nqr                                                   1d12s10
      call enough('erird.  1',bc,ibc)
      zzab=zza+zzb                                                         1d13s10
      zzcd=zzc+zzd                                                         1d13s10
      zzabcd=zzab+zzcd
      dxab=xa-xb                                                        12d1s16
      dyab=ya-yb                                                        12d1s16
      dzab=za-zb                                                        12d1s16
      dxcd=xc-xd                                                        12d1s16
      dycd=yc-yd                                                        12d1s16
      dzcd=zc-zd                                                        12d1s16
      pref=(zza*zzb*((dxab)**2+(dyab)**2+(dzab)**2)/zzab)               12d1s16
     $    +(zzc*zzd*((dxcd)**2+(dycd)**2+(dzcd)**2)/zzcd)               12d1s16
      pref=exp(-pref)                                                   1d22s10
      xab=zza*xa+zzb*xb
      yab=zza*ya+zzb*yb
      zab=zza*za+zzb*zb
      xcd=zzc*xc+zzd*xd
      ycd=zzc*yc+zzd*yd
      zcd=zzc*zc+zzd*zd
      bigx=((zzcd*xab-zzab*xcd)**2+(zzcd*yab-zzab*ycd)**2
     $     +(zzcd*zab-zzab*zcd)**2)/(zzab*zzcd*zzabcd)                  1d29s10
   45 format(i5,f10.5)
      if(idbg.ne.0)then
       write(6,*)('pref = '),pref
       write(6,*)('bigx = '),bigx
       write(6,*)('nqr = '),nqr,nqx,ibcrys
      end if
      if(nqr.ge.nqx)then                                                6d22s12
       write(6,*)('asking for more rys quadrature points than are '),   6d22s12
     $     ('available '),nqr,nqx                                       6d22s12
       call dws_finalize                                                6d22s12
       stop                                                             6d22s12
      end if                                                            6d22s12
      iaddq=ibc(ibcrys+2*(nqr-1))                                       6d22s12
      nreg=ibc(ibcrys+1+2*(nqr-1))                                      6d22s12
      if(idbg.ne.0)then
       write(6,*)('iaddq, nreg '),iaddq,nreg,ibcrys+1+2*(nqr-1)
       write(6,*)('whats under iaddq: '),bc(iaddq)
      end if
      ratio=bigx/bc(iaddq)                                              10d29s15
      if(ratio.le.nreg+1)then                                           10d29s15
       ireg=int(ratio)                                                  10d29s15
      else                                                              10d29s15
       ireg=nreg+1                                                      10d29s15
      end if                                                            10d29s15
      ireg=ireg+1                                                       6d22s12
      if(idbg.ne.0)then
       write(6,*)('ireg, vs nreg: '),ireg,nreg
      end if
      rmsx=0d0
      rmsw=0d0
      if(ireg.le.nreg)then                                              6d22s12
       a=dfloat(ireg-1)*bc(iaddq)                                       6d22s12
       b=a+bc(iaddq)                                                    6d22s12
       y=(2d0*bigx-a-b)/bc(iaddq)                                       6d22s12
       y2=2d0*y                                                         6d22s12
       try=bc(iaddq+2*nqr+ireg)                                         6d22s12
       itry=nint(try)                                                   6d22s12
       itry=itry+ibc(ibcrys)-1
       if(idbg.ne.0)then
        write(6,*)('try, itry '),try,itry
       end if
       do i=1,nqr                                                       6d22s12
        ncoef=nint(bc(itry))                                            6d22s12
        itry=itry+1                                                     6d22s12
        dw=0d0                                                          6d22s12
        ddw=0d0                                                         6d22s12
        dn=0d0                                                          6d22s12
        ddn=0d0                                                         6d22s12
        j1=itry                                                         6d25s12
        do j=ncoef,2,-1                                                 6d22s12
         j2=j1+1                                                        6d22s12
         svw=dw                                                         6d22s12
         svn=dn                                                         6d22s12
         dw=y2*dw-ddw+bc(j1)                                            6d22s12
         dn=y2*dn-ddn+bc(j2)                                            6d22s12
         ddw=svw                                                        6d22s12
         ddn=svn                                                        6d22s12
         j1=j2+1                                                        6d25s12
        end do                                                          6d22s12
        j2=j1+1                                                         6d25s12
        wgt=y*dw-ddw+0.5d0*bc(j1)                                       6d25s12
        xnod=y*dn-ddn+0.5d0*bc(j2)                                       6d25s12
        itry=itry+2*ncoef                                               6d22s12
        bc(irwgt+i-1)=wgt                                               7d5s12
        bc(irnode+i-1)=xnod                                              7d5s12
       end do                                                           6d22s12
      else                                                              6d22s12
       j1=iaddq+1                                                       6d25s12
       argi=1d0/sqrt(bigx)
       do i=1,nqr
        ghw=bc(j1)
        wgt=ghw*argi
        j1=j1+1
        ghx=bc(j1)
        j1=j1+1
        r2=ghx**2
        trial=r2/(bigx-r2)
        bc(irwgt+i-1)=wgt                                               7d5s12
        bc(irnode+i-1)=trial                                            7d5s12
       end do
      end if                                                            6d22s12
      if(idbg.ne.0)then
       do i=1,nqr
        write(6,*)i,bc(irnode+i-1),bc(irwgt+i-1)
       end do
      end if
      jrnode=irnode-1                                                   1d11s10
      jrwgt=irwgt-1                                                     1d11s10
      front0=fmul*pref*2d0/sqrt(pi)                                     8d20s15
      lap=la+1                                                          3d6s17
      lbp=lb+1                                                          3d6s17
      lcp=lc+1                                                          3d6s17
      ldp=ld+1                                                          3d6s17
      lapp=lap+max(idxa,idya,idza)                                      3d6s17
      lbpp=lbp+max(idxb,idyb,idzb)                                      3d6s17
      lcpp=lcp+max(idxc,idyc,idzc)                                      3d6s17
      ldpp=ldp+max(idxd,idyd,idzd)                                      3d6s17
      rho=sqrt(zzab*zzcd/zzabcd)                                           1d13s10
      dnorm=xna*xnb*xnc*xnd                                             1d13s10
      front=dnorm*front0*rho/sqrt((zzab*zzcd)**3)                         1d13s10
      iader=idxa+idya+idza                                              3d6s17
      zza2=zza*2d0                                                      3d6s17
      if(iader.gt.0)front=front*(zza2**iader)                           3d6s17
      ibder=idxb+idyb+idzb                                              3d6s17
      zzb2=zzb*2d0                                                      3d6s17
      if(ibder.gt.0)front=front*(zzb2**ibder)                           3d6s17
      icder=idxc+idyc+idzc                                              3d6s17
      zzc2=zzc*2d0                                                      3d6s17
      if(icder.gt.0)front=front*(zzc2**icder)                           3d6s17
      idder=idxd+idyd+idzd                                              3d6s17
      zzd2=zzd*2d0                                                      3d6s17
      if(idder.gt.0)front=front*(zzd2**idder)                           3d6s17
      if(idbg.ne.0)then                                                 4d26s12
       write(6,*)('front: '),front,dnorm,front0,rho
      end if                                                            4d26s12
      nwdsab=lap*lbp                                                    1d14s10
      nwdsabp=lapp*lbpp                                                 10d28s16
      nwdsab3=nwdsab*3                                                  1d29s10
      nwdsabc=nwdsab*lcp                                                1d14s10
      nwdsabcp=nwdsabp*lcpp                                             10d28s16
      nwdscd=lcp*ldp                                                    1d14s10
      nwdsabcd=nwdsab*nwdscd                                            1d13s10
      nwdsabcd3=nwdsabcd*3                                              1d29s10
      ncarti=ibc(ipp(lap))*ibc(ipp(lbp))*ibc(ipp(lcp))*ibc(ipp(ldp))    1d14s10
      if(idbg.ne.0)then
       write(6,*)('computing ncarti from ')
       write(6,*)lap,ibc(ipp(lap))
       write(6,*)lbp,ibc(ipp(lbp))
       write(6,*)lcp,ibc(ipp(lcp))
       write(6,*)ldp,ibc(ipp(ldp))
       write(6,*)ncarti
      end if
      icarti=ibcoff
      ibcoff=icarti+ncarti                                              1d14s10
      icartt=ibcoff                                                     1d14s10
      ibcoff=icartt+ncarti                                              1d14s10
      call enough('erird.  2',bc,ibc)
      do i=0,ncarti-1                                                   1d14s10
       bc(icarti+i)=0d0                                                 1d14s10
      end do                                                            1d14s10
      testax=xa*zzcd-zzc*xc-zzd*xd+zzb*(dxab)                           12d1s16
      testaxl=zzb*dxab*zzcd                                             12d1s16
      testbx=xb*zzcd-zzc*xc-zzd*xd-zza*dxab                             12d1s16
      testbxl=-zza*dxab*zzcd                                            12d1s16
      testcx=xc*zzab-zza*xa-zzb*xb+zzd*dxcd
      testcxl=zzd*dxcd*zzab
      testdx=xd*zzab-zza*xa-zzb*xb-zzc*dxcd
      testdxl=-zzc*dxcd*zzab
      testay=ya*zzcd-zzc*yc-zzd*yd+zzb*dyab
      testayl=zzb*dyab*zzcd
      testby=yb*zzcd-zzc*yc-zzd*yd-zza*dyab                             12d1s16
      testbyl=-zza*dyab*zzcd                                            12d1s16
      testcy=yc*zzab-zza*ya-zzb*yb+zzd*dycd                             12d1s16
      testcyl=zzd*dycd*zzab                                             12d1s16
      testdy=yd*zzab-zza*ya-zzb*yb-zzc*dycd                             12d1s16
      testdyl=-zzc*dycd*zzab                                            12d1s16
      testaz=za*zzcd-zzc*zc-zzd*zd+zzb*dzab                             12d1s16
      testazl=zzb*dzab*zzcd                                             12d1s16
      testbz=zb*zzcd-zzc*zc-zzd*zd-zza*dzab                             12d1s16
      testbzl=-zza*dzab*zzcd                                            12d1s16
      testcz=zc*zzab-zza*za-zzb*zb+zzd*dzcd                             12d1s16
      testczl=zzd*dzcd*zzab                                             12d1s16
      testdz=zd*zzab-zza*za-zzb*zb-zzc*dzcd                             12d1s16
      testdzl=-zzc*dzcd*zzab                                            12d1s16
      fact=zzab+zzcd
      factl=zzcd*zzab
      bgg=-(zza+zzcd)/zzb
      bggl=-zzcd*zza/zzb
      agg=-(zzb+zzcd)/zza
      aggl=-zzcd*zzb/zza
      cgg=-(zzd+zzab)/zzc
      cggl=-zzab*zzd/zzc
      dgg=-(zzc+zzab)/zzd
      dggl=-zzab*zzc/zzd
      iu2i=ibcoff
      irys=iu2i+nqr
      ifactu=irys+nqr
      ifactax=ifactu+nqr
      ifactbx=ifactax+nqr
      ifactcx=ifactbx+nqr
      ifactdx=ifactcx+nqr
      ifactay=ifactdx+nqr
      ifactby=ifactay+nqr
      ifactcy=ifactby+nqr
      ifactdy=ifactcy+nqr
      ifactaz=ifactdy+nqr
      ifactbz=ifactaz+nqr
      ifactcz=ifactbz+nqr
      ifactdz=ifactcz+nqr
      igab=ifactdz+nqr
      igcd=igab+nqr
      iagg=igcd+nqr
      ibgg=iagg+nqr
      icgg=ibgg+nqr
      idgg=icgg+nqr
      idmd=idgg+nqr
      icmd=idmd+nqr
      ibmd=icmd+nqr
      iamd=ibmd+nqr
      ixr=iamd+nqr
      nsz=nqr*lapp*lbpp*lcpp*ldpp                                       10d28s16
      iyr=ixr+nsz
      izr=iyr+nsz
      ibcoff=izr+nsz
      call enough('erird.  3',bc,ibc)
      ixr4=ibcoff                                                       11d4s16
      xnan=-2d0                                                         6d12s24
      do i=0,3*nsz-1
       bc(ixr+i)=xnan                                                   6d12s24
      end do
      iyr4=ixr4+nqr*(ls+1)                                              11d30s16
      izr4=iyr4+nqr*(ls+1)                                              11d30s16
      ibcoff=izr4+nqr*(ls+1)                                            11d30s16
      call enough('erird.  4',bc,ibc)
      do irq=1,nqr
       irqm=irq-1
       u2i=1d0/(rho*rho*bc(jrnode+irq))
       bc(iu2i+irqm)=u2i
       bc(irys+irqm)=bc(jrwgt+irq)
       bc(ixr+irqm)=pi
       bc(iyr+irqm)=pi
       bc(izr+irqm)=pi*bc(irys+irqm)
       bc(ixr4+irqm)=pi
       bc(iyr4+irqm)=pi
       bc(izr4+irqm)=pi*bc(irys+irqm)
       bc(ifactu+irqm)=1d0/(fact+factl*u2i)
       testaxu=testax+testaxl*u2i
       bc(ifactax+irqm)=testaxu*bc(ifactu+irqm)
       testbxu=testbx+testbxl*u2i
       bc(ifactbx+irqm)=testbxu*bc(ifactu+irqm)
       testcxu=testcx+testcxl*u2i
       bc(ifactcx+irqm)=testcxu*bc(ifactu+irqm)
       testdxu=testdx+testdxl*u2i
       bc(ifactdx+irqm)=testdxu*bc(ifactu+irqm)
       testayu=testay+testayl*u2i
       bc(ifactay+irqm)=testayu*bc(ifactu+irqm)
       testbyu=testby+testbyl*u2i
       bc(ifactby+irqm)=testbyu*bc(ifactu+irqm)
       testcyu=testcy+testcyl*u2i
       bc(ifactcy+irqm)=testcyu*bc(ifactu+irqm)
       testdyu=testdy+testdyl*u2i
       bc(ifactdy+irqm)=testdyu*bc(ifactu+irqm)
       testazu=testaz+testazl*u2i
       bc(ifactaz+irqm)=testazu*bc(ifactu+irqm)
       testbzu=testbz+testbzl*u2i
       bc(ifactbz+irqm)=testbzu*bc(ifactu+irqm)
       testczu=testcz+testczl*u2i
       bc(ifactcz+irqm)=testczu*bc(ifactu+irqm)
       testdzu=testdz+testdzl*u2i
       bc(ifactdz+irqm)=testdzu*bc(ifactu+irqm)
       bc(igab+irqm)=(1d0+zzab*u2i)
       bc(igcd+irqm)=(1d0+zzcd*u2i)
       bc(iagg+irqm)=agg+aggl*u2i
       bc(icgg+irqm)=cgg+cggl*u2i
       bc(idgg+irqm)=dgg+dggl*u2i
       bc(ibgg+irqm)=bgg+bggl*u2i
      end do
c
c     compute integrals required w/o differentition
c
      do id=1,ldp
       tmp=0.5d0*dfloat(id-1)
       do iq=0,nqr-1
        bc(idmd+iq)=tmp*bc(ifactu+iq)
       end do
       jdp=id
       jd0=id-1
       jdm=max(0,id-2)
       do ic=1,lcp
        tmp=0.5d0*dfloat(ic-1)
        do iq=0,nqr-1
         bc(icmd+iq)=tmp*bc(ifactu+iq)
        end do
        jcp=ic
        jc0=ic-1
        jcm=max(0,ic-2)
        k00=lbpp*(lcpp*jd0+jc0)                                         10d28s16
        kp0=lbpp*(jcp+lcpp*jd0)                                         10d28s16
        k0p=lbpp*(jc0+lcpp*jdp)                                         10d28s16
        km0=lbpp*(jcm+lcpp*jd0)                                         10d28s16
        k0m=lbpp*(jc0+lcpp*jdm)                                         10d28s16
        do ib=1,lbp
         tmp=0.5d0*dfloat(ib-1)
         do iq=0,nqr-1
          bc(ibmd+iq)=tmp*bc(ifactu+iq)
         end do
         jbp=ib
         jb0=ib-1
         jbm=max(0,ib-2)
         k000=lapp*(jb0+k00)                                            10d28s16
         kp00=lapp*(jbp+k00)                                            10d28s16
         k0p0=lapp*(jb0+kp0)                                            10d28s16
         k00p=lapp*(jb0+k0p)                                            10d28s16
         km00=lapp*(jbm+k00)                                            10d28s16
         k0m0=lapp*(jb0+km0)                                            10d28s16
         k00m=lapp*(jb0+k0m)                                            10d28s16
         do ia=1,lap
          tmp=0.5d0*dfloat(ia-1)
          do iq=0,nqr-1
           bc(iamd+iq)=tmp*bc(ifactu+iq)
          end do
          iad=nqr*(ia-1+k000)
          iadx=ixr+iad
          iady=iyr+iad
          iadz=izr+iad
          iam=nqr*(max(0,ia-2)+k000)
          iamx=ixr+iam
          iamy=iyr+iam
          iamz=izr+iam
          ibm=nqr*(ia-1+km00)
          ibmx=ixr+ibm
          ibmy=iyr+ibm
          ibmz=izr+ibm
          icm=nqr*(ia-1+k0m0)
          icmx=ixr+icm
          icmy=iyr+icm
          icmz=izr+icm
          iem=nqr*(ia-1+k00m)
          idmx=ixr+iem
          idmy=iyr+iem
          idmz=izr+iem
          if(id.lt.ldp)then
           idp=nqr*(ia-1+k00p)
           idpx=ixr+idp
           idpy=iyr+idp
           idpz=izr+idp
           do iq=0,nqr-1
            cmdu=bc(igab+iq)*bc(icmd+iq)
            dmdu=bc(idgg+iq)*bc(idmd+iq)
            bc(idpx+iq)=bc(iamd+iq)*bc(iamx+iq)+bc(ibmd+iq)*bc(ibmx+iq)+
     $      cmdu*bc(icmx+iq)+dmdu*bc(idmx+iq)+bc(ifactdx+iq)*bc(iadx+iq)
            bc(idpy+iq)=bc(iamd+iq)*bc(iamy+iq)+bc(ibmd+iq)*bc(ibmy+iq)+
     $      cmdu*bc(icmy+iq)+dmdu*bc(idmy+iq)+bc(ifactdy+iq)*bc(iady+iq)
            bc(idpz+iq)=bc(iamd+iq)*bc(iamz+iq)+bc(ibmd+iq)*bc(ibmz+iq)+
     $      cmdu*bc(icmz+iq)+dmdu*bc(idmz+iq)+bc(ifactdz+iq)*bc(iadz+iq)
           end do
          end if
          if(ic.lt.lcp)then
           idp=nqr*(ia-1+k0p0)
           idpx=ixr+idp
           idpy=iyr+idp
           idpz=izr+idp
           do iq=0,nqr-1
            dmdu=bc(igab+iq)*bc(idmd+iq)
            cmdu=bc(icgg+iq)*bc(icmd+iq)
            bc(idpx+iq)=bc(iamd+iq)*bc(iamx+iq)+bc(ibmd+iq)*bc(ibmx+iq)
     $     +cmdu*bc(icmx+iq)+dmdu*bc(idmx+iq)+bc(ifactcx+iq)*bc(iadx+iq)
            bc(idpy+iq)=bc(iamd+iq)*bc(iamy+iq)+bc(ibmd+iq)*bc(ibmy+iq)
     $     +cmdu*bc(icmy+iq)+dmdu*bc(idmy+iq)+bc(ifactcy+iq)*bc(iady+iq)
            bc(idpz+iq)=bc(iamd+iq)*bc(iamz+iq)+bc(ibmd+iq)*bc(ibmz+iq)
     $     +cmdu*bc(icmz+iq)+dmdu*bc(idmz+iq)+bc(ifactcz+iq)*bc(iadz+iq)
           end do
          end if
          if(ia.lt.lap)then
           idp=nqr*(ia+k000)
           idpx=ixr+idp
           idpy=iyr+idp
           idpz=izr+idp
           do iq=0,nqr-1
            bmdu=bc(igcd+iq)*bc(ibmd+iq)
            amdu=bc(iagg+iq)*bc(iamd+iq)
            bc(idpx+iq)=amdu*bc(iamx+iq)+bmdu*bc(ibmx+iq)+bc(icmd+iq)
     $   *bc(icmx+iq)+bc(idmd+iq)*bc(idmx+iq)+bc(ifactax+iq)*bc(iadx+iq)
            bc(idpy+iq)=amdu*bc(iamy+iq)+bmdu*bc(ibmy+iq)+bc(icmd+iq)
     $   *bc(icmy+iq)+bc(idmd+iq)*bc(idmy+iq)+bc(ifactay+iq)*bc(iady+iq)
            bc(idpz+iq)=amdu*bc(iamz+iq)+bmdu*bc(ibmz+iq)+bc(icmd+iq)
     $   *bc(icmz+iq)+bc(idmd+iq)*bc(idmz+iq)+bc(ifactaz+iq)*bc(iadz+iq)
           end do
          end if
          if(ib.lt.lbp)then
           idp=nqr*(ia-1+kp00)
           idpx=ixr+idp
           idpy=iyr+idp
           idpz=izr+idp
           do iq=0,nqr-1
            amdu=bc(igcd+iq)*bc(iamd+iq)
            bmdu=bc(ibgg+iq)*bc(ibmd+iq)
            bc(idpx+iq)=amdu*bc(iamx+iq)+bmdu*bc(ibmx+iq)+bc(icmd+iq)
     $   *bc(icmx+iq)+bc(idmd+iq)*bc(idmx+iq)+bc(ifactbx+iq)*bc(iadx+iq)
            bc(idpy+iq)=amdu*bc(iamy+iq)+bmdu*bc(ibmy+iq)+bc(icmd+iq)
     $   *bc(icmy+iq)+bc(idmd+iq)*bc(idmy+iq)+bc(ifactby+iq)*bc(iady+iq)
            bc(idpz+iq)=amdu*bc(iamz+iq)+bmdu*bc(ibmz+iq)+bc(icmd+iq)
     $   *bc(icmz+iq)+bc(idmd+iq)*bc(idmz+iq)+bc(ifactbz+iq)*bc(iadz+iq)
           end do
          end if
         end do
        end do
       end do
      end do
c
c     now compute extra integrals required by differentiation.
c
      if(max(idxa,idya,idza).gt.0)then
c
c     paraerid case or start of paraerid2b case or perhaps onedints case
c     paraerid: der wrt atom a
c     paraerid2b case: wrt atom a and b
c
       if(idxa.gt.0)then
        ixu=ixr
        ifactau=ifactax
        ifactbu=ifactbx
        ifactcu=ifactcx
        ifactdu=ifactdx
       else if(idya.gt.0)then
        ixu=iyr
        ifactau=ifactay
        ifactbu=ifactby
        ifactcu=ifactcy
        ifactdu=ifactdy
       else
        ixu=izr
        ifactau=ifactaz
        ifactbu=ifactbz
        ifactcu=ifactcz
        ifactdu=ifactdz
       end if
       do id=1,ldp
        tmp=0.5d0*dfloat(id-1)
        do iq=0,nqr-1
         bc(idmd+iq)=tmp*bc(ifactu+iq)
        end do
        jdp=id
        jd0=id-1
        jdm=max(0,id-2)
        do ic=1,lcp
         tmp=0.5d0*dfloat(ic-1)
         do iq=0,nqr-1
          bc(icmd+iq)=tmp*bc(ifactu+iq)
         end do
         jcp=ic
         jc0=ic-1
         jcm=max(0,ic-2)
         k00=lbpp*(lcpp*jd0+jc0)                                         10d28s16
         kp0=lbpp*(jcp+lcpp*jd0)                                         10d28s16
         k0p=lbpp*(jc0+lcpp*jdp)                                         10d28s16
         km0=lbpp*(jcm+lcpp*jd0)                                         10d28s16
         k0m=lbpp*(jc0+lcpp*jdm)                                         10d28s16
         do ib=1,lbp
          tmp=0.5d0*dfloat(ib-1)
          do iq=0,nqr-1
           bc(ibmd+iq)=tmp*bc(ifactu+iq)
          end do
          jbp=ib
          jb0=ib-1
          jbm=max(0,ib-2)
          k000=lapp*(jb0+k00)                                            10d28s16
          kp00=lapp*(jbp+k00)                                            10d28s16
          k0p0=lapp*(jb0+kp0)                                            10d28s16
          k00p=lapp*(jb0+k0p)                                            10d28s16
          km00=lapp*(jbm+k00)                                            10d28s16
          k0m0=lapp*(jb0+km0)                                            10d28s16
          k00m=lapp*(jb0+k0m)                                            10d28s16
          do ia=lap,lapp                                                3d6s17
           tmp=0.5d0*dfloat(ia-1)
           do iq=0,nqr-1
            bc(iamd+iq)=tmp*bc(ifactu+iq)
           end do
           iad=nqr*(ia-1+k000)
           iadu=ixu+iad
           iam=nqr*(max(0,ia-2)+k000)
           iamu=ixu+iam
           ibm=nqr*(ia-1+km00)
           ibmu=ixu+ibm
           icm=nqr*(ia-1+k0m0)
           icmu=ixu+icm
           iem=nqr*(ia-1+k00m)
           idmu=ixu+iem
           if(id.lt.ldp.and.ia.ne.lap)then
            idp=nqr*(ia-1+k00p)
            idpu=ixu+idp
            do iq=0,nqr-1
             cmdu=bc(igab+iq)*bc(icmd+iq)
             dmdu=bc(idgg+iq)*bc(idmd+iq)
             bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $            +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $            +dmdu*bc(idmu+iq)+bc(ifactdu+iq)*bc(iadu+iq)
            end do
           end if
           if(ic.lt.lcp.and.ia.ne.lap)then
            idp=nqr*(ia-1+k0p0)
            idpu=ixu+idp
            do iq=0,nqr-1
             dmdu=bc(igab+iq)*bc(idmd+iq)
             cmdu=bc(icgg+iq)*bc(icmd+iq)
             bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $            +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $            +dmdu*bc(idmu+iq)+bc(ifactcu+iq)*bc(iadu+iq)
            end do
           end if
           if(ia.lt.lapp)then
            idp=nqr*(ia+k000)
            idpu=ixu+idp
            do iq=0,nqr-1
             bmdu=bc(igcd+iq)*bc(ibmd+iq)
             amdu=bc(iagg+iq)*bc(iamd+iq)
             bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $            *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $            +bc(ifactau+iq)*bc(iadu+iq)
            end do
           end if
           if(ib.lt.lbp.and.ia.ne.lap)then
            idp=nqr*(ia-1+kp00)
            idpu=ixu+idp
            do iq=0,nqr-1
             amdu=bc(igcd+iq)*bc(iamd+iq)
             bmdu=bc(ibgg+iq)*bc(ibmd+iq)
             bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $            *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $            +bc(ifactbu+iq)*bc(iadu+iq)
            end do
           end if
          end do
         end do
        end do
       end do
       if(max(idxb,idyb,idzb).gt.0)then
c
c     rest of paraerid2b case or part 2 of onedints cas                 8d15s24
c
        do id=1,ldp
         tmp=0.5d0*dfloat(id-1)
         do iq=0,nqr-1
          bc(idmd+iq)=tmp*bc(ifactu+iq)
         end do
         jdp=id
         jd0=id-1
         jdm=max(0,id-2)
         do ic=1,lcp
          tmp=0.5d0*dfloat(ic-1)
          do iq=0,nqr-1
           bc(icmd+iq)=tmp*bc(ifactu+iq)
          end do
          jcp=ic
          jc0=ic-1
          jcm=max(0,ic-2)
          k00=lbpp*(lcpp*jd0+jc0)                                         10d28s16
          kp0=lbpp*(jcp+lcpp*jd0)                                         10d28s16
          k0p=lbpp*(jc0+lcpp*jdp)                                         10d28s16
          km0=lbpp*(jcm+lcpp*jd0)                                         10d28s16
          k0m=lbpp*(jc0+lcpp*jdm)                                         10d28s16
          do ib=lbp,lbpp                                                3d6s17
           tmp=0.5d0*dfloat(ib-1)
           do iq=0,nqr-1
            bc(ibmd+iq)=tmp*bc(ifactu+iq)
           end do
           jbp=ib
           jb0=ib-1
           jbm=max(0,ib-2)
           k000=lapp*(jb0+k00)                                            10d28s16
           kp00=lapp*(jbp+k00)                                            10d28s16
           k0p0=lapp*(jb0+kp0)                                            10d28s16
           k00p=lapp*(jb0+k0p)                                            10d28s16
           km00=lapp*(jbm+k00)                                            10d28s16
           k0m0=lapp*(jb0+km0)                                            10d28s16
           k00m=lapp*(jb0+k0m)                                            10d28s16
           do ia=1,lapp                                                 3d6s17
            tmp=0.5d0*dfloat(ia-1)
            do iq=0,nqr-1
             bc(iamd+iq)=tmp*bc(ifactu+iq)
            end do
            iad=nqr*(ia-1+k000)
            iadu=ixu+iad
            iam=nqr*(max(0,ia-2)+k000)
            iamu=ixu+iam
            ibm=nqr*(ia-1+km00)
            ibmu=ixu+ibm
            icm=nqr*(ia-1+k0m0)
            icmu=ixu+icm
            iem=nqr*(ia-1+k00m)
            idmu=ixu+iem
            if(id.lt.ldp.and.ib.ne.lbp)then
             idp=nqr*(ia-1+k00p)
             idpu=ixu+idp
             do iq=0,nqr-1
              cmdu=bc(igab+iq)*bc(icmd+iq)
              dmdu=bc(idgg+iq)*bc(idmd+iq)
              bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $             +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $             +dmdu*bc(idmu+iq)+bc(ifactdu+iq)*bc(iadu+iq)
             end do
            end if
            if(ic.lt.lcp.and.ib.ne.lbp)then
             idp=nqr*(ia-1+k0p0)
             idpu=ixu+idp
             do iq=0,nqr-1
              dmdu=bc(igab+iq)*bc(idmd+iq)
              cmdu=bc(icgg+iq)*bc(icmd+iq)
              bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $             +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $             +dmdu*bc(idmu+iq)+bc(ifactcu+iq)*bc(iadu+iq)
             end do
            end if
            if(ia.lt.lapp.and.ib.ne.lbp)then
             idp=nqr*(ia+k000)
             idpu=ixu+idp
             do iq=0,nqr-1
              bmdu=bc(igcd+iq)*bc(ibmd+iq)
              amdu=bc(iagg+iq)*bc(iamd+iq)
              bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $             *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $             +bc(ifactau+iq)*bc(iadu+iq)
             end do
            end if
            if(ib.lt.lbpp)then
             idp=nqr*(ia-1+kp00)
             idpu=ixu+idp
             do iq=0,nqr-1
              amdu=bc(igcd+iq)*bc(iamd+iq)
              bmdu=bc(ibgg+iq)*bc(ibmd+iq)
              bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $             *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $             +bc(ifactbu+iq)*bc(iadu+iq)
             end do
            end if
           end do
          end do
         end do
        end do
        if(max(idxc,idyc,idzc).gt.0)then                                8d15s24
c
c     final part of onedints case
c
         if(idxc.gt.0)then                                               8d15s24
          ixu=ixr                                                        8d15s24
          ifactau=ifactax                                                8d15s24
          ifactbu=ifactbx                                                8d15s24
          ifactcu=ifactcx                                                8d15s24
          ifactdu=ifactdx                                                8d15s24
         else if(idyc.gt.0)then                                          8d15s24
          ixu=iyr                                                        8d15s24
          ifactau=ifactay                                                8d15s24
          ifactbu=ifactby                                                8d15s24
          ifactcu=ifactcy                                                8d15s24
          ifactdu=ifactdy                                                8d15s24
         else                                                            8d15s24
          ixu=izr                                                        8d15s24
          ifactau=ifactaz                                                8d15s24
          ifactbu=ifactbz                                                8d15s24
          ifactcu=ifactcz                                                8d15s24
          ifactdu=ifactdz                                                8d15s24
         end if                                                          8d15s24
         do id=1,ldp                                                     8d15s24
          tmp=0.5d0*dfloat(id-1)                                         8d15s24
          do iq=0,nqr-1                                                  8d15s24
           bc(idmd+iq)=tmp*bc(ifactu+iq)                                 8d15s24
          end do                                                         8d15s24
          jdp=id                                                         8d15s24
          jd0=id-1                                                       8d15s24
          jdm=max(0,id-2)                                                8d15s24
          do ic=lcp,lcpp                                                 8d15s24
           tmp=0.5d0*dfloat(ic-1)                                        8d15s24
           do iq=0,nqr-1                                                 8d15s24
            bc(icmd+iq)=tmp*bc(ifactu+iq)                                8d15s24
           end do                                                        8d15s24
           jcp=ic                                                        8d15s24
           jc0=ic-1                                                      8d15s24
           jcm=max(0,ic-2)                                               8d15s24
           k00=lbpp*(lcpp*jd0+jc0)                                       8d15s24
           kp0=lbpp*(jcp+lcpp*jd0)                                       8d15s24
           k0p=lbpp*(jc0+lcpp*jdp)                                       8d15s24
           km0=lbpp*(jcm+lcpp*jd0)                                       8d15s24
           k0m=lbpp*(jc0+lcpp*jdm)                                       8d15s24
           do ib=1,lbpp                                                  8d15s24
            tmp=0.5d0*dfloat(ib-1)                                       8d15s24
            do iq=0,nqr-1                                                8d15s24
             bc(ibmd+iq)=tmp*bc(ifactu+iq)                               8d15s24
            end do                                                       8d15s24
            jbp=ib                                                       8d15s24
            jb0=ib-1                                                     8d15s24
            jbm=max(0,ib-2)                                              8d15s24
            k000=lapp*(jb0+k00)                                          8d15s24
            kp00=lapp*(jbp+k00)                                          8d15s24
            k0p0=lapp*(jb0+kp0)                                          8d15s24
            k00p=lapp*(jb0+k0p)                                          8d15s24
            km00=lapp*(jbm+k00)                                          8d15s24
            k0m0=lapp*(jb0+km0)                                          8d15s24
            k00m=lapp*(jb0+k0m)                                          8d15s24
            do ia=1,lapp                                                 8d15s24
             tmp=0.5d0*dfloat(ia-1)                                      8d15s24
             do iq=0,nqr-1                                               8d15s24
              bc(iamd+iq)=tmp*bc(ifactu+iq)                              8d15s24
             end do                                                      8d15s24
             iad=nqr*(ia-1+k000)                                         8d15s24
             iadu=ixu+iad                                                8d15s24
             iam=nqr*(max(0,ia-2)+k000)                                  8d15s24
             iamu=ixu+iam                                                8d15s24
             ibm=nqr*(ia-1+km00)                                         8d15s24
             ibmu=ixu+ibm                                                8d15s24
             icm=nqr*(ia-1+k0m0)                                         8d15s24
             icmu=ixu+icm                                                8d15s24
             iem=nqr*(ia-1+k00m)                                         8d15s24
             idmu=ixu+iem                                                8d15s24
             if(id.lt.ldp.and.ic.ne.lcp)then                             8d15s24
              idp=nqr*(ia-1+k00p)                                        8d15s24
              idpu=ixu+idp                                               8d15s24
              do iq=0,nqr-1                                              8d15s24
               cmdu=bc(igab+iq)*bc(icmd+iq)                              8d15s24
               dmdu=bc(idgg+iq)*bc(idmd+iq)                              8d15s24
               bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)                       8d15s24
     $              +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)            8d15s24
     $              +dmdu*bc(idmu+iq)+bc(ifactdu+iq)*bc(iadu+iq)         8d15s24
              end do                                                     8d15s24
             end if                                                      8d15s24
             if(ic.lt.lcpp)then                                          8d15s24
              idp=nqr*(ia-1+k0p0)                                        8d15s24
              idpu=ixu+idp                                               8d15s24
              do iq=0,nqr-1                                              8d15s24
               dmdu=bc(igab+iq)*bc(idmd+iq)                              8d15s24
               cmdu=bc(icgg+iq)*bc(icmd+iq)                              8d15s24
               bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)                       8d15s24
     $              +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)            8d15s24
     $              +dmdu*bc(idmu+iq)+bc(ifactcu+iq)*bc(iadu+iq)         8d15s24
              end do                                                     8d15s24
             end if                                                      8d15s24
             if(ia.lt.lapp.and.ic.ne.lcp)then                            8d15s24
              idp=nqr*(ia+k000)                                          8d15s24
              idpu=ixu+idp                                               8d15s24
              do iq=0,nqr-1                                              8d15s24
               bmdu=bc(igcd+iq)*bc(ibmd+iq)                              8d15s24
               amdu=bc(iagg+iq)*bc(iamd+iq)                              8d15s24
               bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq) 8d15s24
     $             *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)                 8d15s24
     $             +bc(ifactau+iq)*bc(iadu+iq)                          8d15s24
              end do                                                     8d15s24
             end if                                                      8d15s24
             if(ib.lt.lbpp.and.ic.ne.lcp)then                            8d15s24
              idp=nqr*(ia-1+kp00)                                        8d15s24
              idpu=ixu+idp                                               8d15s24
              do iq=0,nqr-1                                              8d15s24
               amdu=bc(igcd+iq)*bc(iamd+iq)                              8d15s24
               bmdu=bc(ibgg+iq)*bc(ibmd+iq)                              8d15s24
               bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq) 8d15s24
     $             *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)                 8d15s24
     $             +bc(ifactbu+iq)*bc(iadu+iq)                          8d15s24
              end do                                                     8d15s24
             end if                                                      8d15s24
            end do                                                       8d15s24
           end do                                                        8d15s24
          end do                                                         8d15s24
         end do                                                          8d15s24
        end if                                                           8d15s24
       end if
      else if(max(idxc,idyc,idzc).gt.0)then
c
c     paraeridj case
c     der wrt c
c
       if(idxc.gt.0)then
        ixu=ixr
        ifactau=ifactax
        ifactbu=ifactbx
        ifactcu=ifactcx
        ifactdu=ifactdx
       else if(idyc.gt.0)then
        ixu=iyr
        ifactau=ifactay
        ifactbu=ifactby
        ifactcu=ifactcy
        ifactdu=ifactdy
       else
        ixu=izr
        ifactau=ifactaz
        ifactbu=ifactbz
        ifactcu=ifactcz
        ifactdu=ifactdz
       end if
       do id=1,ldp
        tmp=0.5d0*dfloat(id-1)
        do iq=0,nqr-1
         bc(idmd+iq)=tmp*bc(ifactu+iq)
        end do
        jdp=id
        jd0=id-1
        jdm=max(0,id-2)
        do ic=lcp,lcpp                                                  3d6s17
         tmp=0.5d0*dfloat(ic-1)
         do iq=0,nqr-1
          bc(icmd+iq)=tmp*bc(ifactu+iq)
         end do
         jcp=ic
         jc0=ic-1
         jcm=max(0,ic-2)
         k00=lbpp*(lcpp*jd0+jc0)                                         10d28s16
         kp0=lbpp*(jcp+lcpp*jd0)                                         10d28s16
         k0p=lbpp*(jc0+lcpp*jdp)                                         10d28s16
         km0=lbpp*(jcm+lcpp*jd0)                                         10d28s16
         k0m=lbpp*(jc0+lcpp*jdm)                                         10d28s16
         do ib=1,lbp
          tmp=0.5d0*dfloat(ib-1)
          do iq=0,nqr-1
           bc(ibmd+iq)=tmp*bc(ifactu+iq)
          end do
          jbp=ib
          jb0=ib-1
          jbm=max(0,ib-2)
          k000=lapp*(jb0+k00)                                            10d28s16
          kp00=lapp*(jbp+k00)                                            10d28s16
          k0p0=lapp*(jb0+kp0)                                            10d28s16
          k00p=lapp*(jb0+k0p)                                            10d28s16
          km00=lapp*(jbm+k00)                                            10d28s16
          k0m0=lapp*(jb0+km0)                                            10d28s16
          k00m=lapp*(jb0+k0m)                                            10d28s16
          do ia=1,lap                                                   3d6s17
           tmp=0.5d0*dfloat(ia-1)
           do iq=0,nqr-1
            bc(iamd+iq)=tmp*bc(ifactu+iq)
           end do
           iad=nqr*(ia-1+k000)
           iadu=ixu+iad
           iam=nqr*(max(0,ia-2)+k000)
           iamu=ixu+iam
           ibm=nqr*(ia-1+km00)
           ibmu=ixu+ibm
           icm=nqr*(ia-1+k0m0)
           icmu=ixu+icm
           iem=nqr*(ia-1+k00m)
           idmu=ixu+iem
           if(id.lt.ldp.and.ic.ne.lcp)then
            idp=nqr*(ia-1+k00p)
            idpu=ixu+idp
            do iq=0,nqr-1
             cmdu=bc(igab+iq)*bc(icmd+iq)
             dmdu=bc(idgg+iq)*bc(idmd+iq)
             bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $            +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $            +dmdu*bc(idmu+iq)+bc(ifactdu+iq)*bc(iadu+iq)
            end do
           end if
           if(ic.lt.lcpp)then
            idp=nqr*(ia-1+k0p0)
            idpu=ixu+idp
            do iq=0,nqr-1
             dmdu=bc(igab+iq)*bc(idmd+iq)
             cmdu=bc(icgg+iq)*bc(icmd+iq)
             bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $            +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $            +dmdu*bc(idmu+iq)+bc(ifactcu+iq)*bc(iadu+iq)
            end do
           end if
           if(ia.lt.lap.and.ic.ne.lcp)then
            idp=nqr*(ia+k000)
            idpu=ixu+idp
            do iq=0,nqr-1
             bmdu=bc(igcd+iq)*bc(ibmd+iq)
             amdu=bc(iagg+iq)*bc(iamd+iq)
             bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $            *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $            +bc(ifactau+iq)*bc(iadu+iq)
            end do
           end if
           if(ib.lt.lbp.and.ic.ne.lcp)then
            idp=nqr*(ia-1+kp00)
            idpu=ixu+idp
            do iq=0,nqr-1
             amdu=bc(igcd+iq)*bc(iamd+iq)
             bmdu=bc(ibgg+iq)*bc(ibmd+iq)
             bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $            *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $            +bc(ifactbu+iq)*bc(iadu+iq)
            end do
           end if
          end do
         end do
        end do
       end do
      else if(max(idxb,idyb,idzb).gt.0)then
c
c     start of paraerid2c case
c     der wrt b and d
c
       if(idxb.gt.0)then
        ixu=ixr
        ifactau=ifactax
        ifactbu=ifactbx
        ifactcu=ifactcx
        ifactdu=ifactdx
       else if(idyb.gt.0)then
        ixu=iyr
        ifactau=ifactay
        ifactbu=ifactby
        ifactcu=ifactcy
        ifactdu=ifactdy
       else
        ixu=izr
        ifactau=ifactaz
        ifactbu=ifactbz
        ifactcu=ifactcz
        ifactdu=ifactdz
       end if
       do id=1,ldp
        tmp=0.5d0*dfloat(id-1)
        do iq=0,nqr-1
         bc(idmd+iq)=tmp*bc(ifactu+iq)
        end do
        jdp=id
        jd0=id-1
        jdm=max(0,id-2)
        do ic=1,lcp
         tmp=0.5d0*dfloat(ic-1)
         do iq=0,nqr-1
          bc(icmd+iq)=tmp*bc(ifactu+iq)
         end do
         jcp=ic
         jc0=ic-1
         jcm=max(0,ic-2)
         k00=lbpp*(lcpp*jd0+jc0)                                         10d28s16
         kp0=lbpp*(jcp+lcpp*jd0)                                         10d28s16
         k0p=lbpp*(jc0+lcpp*jdp)                                         10d28s16
         km0=lbpp*(jcm+lcpp*jd0)                                         10d28s16
         k0m=lbpp*(jc0+lcpp*jdm)                                         10d28s16
         do ib=lbp,lbpp                                                 3d6s17
          tmp=0.5d0*dfloat(ib-1)
          do iq=0,nqr-1
           bc(ibmd+iq)=tmp*bc(ifactu+iq)
          end do
          jbp=ib
          jb0=ib-1
          jbm=max(0,ib-2)
          k000=lapp*(jb0+k00)                                            10d28s16
          kp00=lapp*(jbp+k00)                                            10d28s16
          k0p0=lapp*(jb0+kp0)                                            10d28s16
          k00p=lapp*(jb0+k0p)                                            10d28s16
          km00=lapp*(jbm+k00)                                            10d28s16
          k0m0=lapp*(jb0+km0)                                            10d28s16
          k00m=lapp*(jb0+k0m)                                            10d28s16
          do ia=1,lap                                                   3d6s17
           tmp=0.5d0*dfloat(ia-1)
           do iq=0,nqr-1
            bc(iamd+iq)=tmp*bc(ifactu+iq)
           end do
           iad=nqr*(ia-1+k000)
           iadu=ixu+iad
           iam=nqr*(max(0,ia-2)+k000)
           iamu=ixu+iam
           ibm=nqr*(ia-1+km00)
           ibmu=ixu+ibm
           icm=nqr*(ia-1+k0m0)
           icmu=ixu+icm
           iem=nqr*(ia-1+k00m)
           idmu=ixu+iem
           if(id.lt.ldp.and.ib.ne.lbp)then
            idp=nqr*(ia-1+k00p)
            idpu=ixu+idp
            do iq=0,nqr-1
             cmdu=bc(igab+iq)*bc(icmd+iq)
             dmdu=bc(idgg+iq)*bc(idmd+iq)
             bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $            +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $            +dmdu*bc(idmu+iq)+bc(ifactdu+iq)*bc(iadu+iq)
            end do
           end if
           if(ic.lt.lcp.and.ib.ne.lbp)then
            idp=nqr*(ia-1+k0p0)
            idpu=ixu+idp
            do iq=0,nqr-1
             dmdu=bc(igab+iq)*bc(idmd+iq)
             cmdu=bc(icgg+iq)*bc(icmd+iq)
             bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $            +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $            +dmdu*bc(idmu+iq)+bc(ifactcu+iq)*bc(iadu+iq)
            end do
           end if
           if(ia.lt.lap.and.ib.ne.lbp)then
            idp=nqr*(ia+k000)
            idpu=ixu+idp
            do iq=0,nqr-1
             bmdu=bc(igcd+iq)*bc(ibmd+iq)
             amdu=bc(iagg+iq)*bc(iamd+iq)
             bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $            *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $            +bc(ifactau+iq)*bc(iadu+iq)
            end do
           end if
           if(ib.lt.lbpp)then
            idp=nqr*(ia-1+kp00)
            idpu=ixu+idp
            do iq=0,nqr-1
             amdu=bc(igcd+iq)*bc(iamd+iq)
             bmdu=bc(ibgg+iq)*bc(ibmd+iq)
             bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $            *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $            +bc(ifactbu+iq)*bc(iadu+iq)
            end do
           end if
          end do
         end do
        end do
       end do
       if(max(idxd,idyd,idzd).gt.0)then
c
c     rest of paraerid2c case
c
        do id=ldp,ldpp
         tmp=0.5d0*dfloat(id-1)
         do iq=0,nqr-1
          bc(idmd+iq)=tmp*bc(ifactu+iq)
         end do
         jdp=id
         jd0=id-1
         jdm=max(0,id-2)
         do ic=1,lcp                                                    3d6s17
          tmp=0.5d0*dfloat(ic-1)
          do iq=0,nqr-1
           bc(icmd+iq)=tmp*bc(ifactu+iq)
          end do
          jcp=ic
          jc0=ic-1
          jcm=max(0,ic-2)
          k00=lbpp*(lcpp*jd0+jc0)                                         10d28s16
          kp0=lbpp*(jcp+lcpp*jd0)                                         10d28s16
          k0p=lbpp*(jc0+lcpp*jdp)                                         10d28s16
          km0=lbpp*(jcm+lcpp*jd0)                                         10d28s16
          k0m=lbpp*(jc0+lcpp*jdm)                                         10d28s16
          do ib=1,lbpp                                                  3d6s17
           tmp=0.5d0*dfloat(ib-1)
           do iq=0,nqr-1
            bc(ibmd+iq)=tmp*bc(ifactu+iq)
           end do
           jbp=ib
           jb0=ib-1
           jbm=max(0,ib-2)
           k000=lapp*(jb0+k00)                                            10d28s16
           kp00=lapp*(jbp+k00)                                            10d28s16
           k0p0=lapp*(jb0+kp0)                                            10d28s16
           k00p=lapp*(jb0+k0p)                                            10d28s16
           km00=lapp*(jbm+k00)                                            10d28s16
           k0m0=lapp*(jb0+km0)                                            10d28s16
           k00m=lapp*(jb0+k0m)                                            10d28s16
           do ia=1,lap                                                   3d6s17
            tmp=0.5d0*dfloat(ia-1)
            do iq=0,nqr-1
             bc(iamd+iq)=tmp*bc(ifactu+iq)
            end do
            iad=nqr*(ia-1+k000)
            iadu=ixu+iad
            iam=nqr*(max(0,ia-2)+k000)
            iamu=ixu+iam
            ibm=nqr*(ia-1+km00)
            ibmu=ixu+ibm
            icm=nqr*(ia-1+k0m0)
            icmu=ixu+icm
            iem=nqr*(ia-1+k00m)
            idmu=ixu+iem
            if(id.lt.ldpp)then
             idp=nqr*(ia-1+k00p)
             idpu=ixu+idp
             do iq=0,nqr-1
              cmdu=bc(igab+iq)*bc(icmd+iq)
              dmdu=bc(idgg+iq)*bc(idmd+iq)
              bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $             +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $             +dmdu*bc(idmu+iq)+bc(ifactdu+iq)*bc(iadu+iq)
             end do
            end if
            if(ic.lt.lcp.and.id.ne.ldp)then
             idp=nqr*(ia-1+k0p0)
             idpu=ixu+idp
             do iq=0,nqr-1
              dmdu=bc(igab+iq)*bc(idmd+iq)
              cmdu=bc(icgg+iq)*bc(icmd+iq)
              bc(idpu+iq)=bc(iamd+iq)*bc(iamu+iq)
     $             +bc(ibmd+iq)*bc(ibmu+iq)+cmdu*bc(icmu+iq)
     $             +dmdu*bc(idmu+iq)+bc(ifactcu+iq)*bc(iadu+iq)
             end do
            end if
            if(ia.lt.lap.and.id.ne.ldp)then
             idp=nqr*(ia+k000)
             idpu=ixu+idp
             do iq=0,nqr-1
              bmdu=bc(igcd+iq)*bc(ibmd+iq)
              amdu=bc(iagg+iq)*bc(iamd+iq)
              bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $             *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $             +bc(ifactau+iq)*bc(iadu+iq)
             end do
            end if
            if(ib.lt.lbpp.and.id.ne.ldp)then
             idp=nqr*(ia-1+kp00)
             idpu=ixu+idp
             do iq=0,nqr-1
              amdu=bc(igcd+iq)*bc(iamd+iq)
              bmdu=bc(ibgg+iq)*bc(ibmd+iq)
              bc(idpu+iq)=amdu*bc(iamu+iq)+bmdu*bc(ibmu+iq)+bc(icmd+iq)
     $             *bc(icmu+iq)+bc(idmd+iq)*bc(idmu+iq)
     $             +bc(ifactbu+iq)*bc(iadu+iq)
             end do
            end if
           end do
          end do
         end do
        end do
       end if
      end if
      if(idbg.ne.0)then
       itmp=ibcoff
       ibcoff=itmp+lap*lbp*lcp*ldp
       call enough('erird.  5',bc,ibc)
       do iq=1,nqr
        write(6,*)('for rys quadrature point no. '),iq
        write(6,*)('cartx ')
        do id=0,ld
         jd=id+idxd
         do ic=0,lc
          jc=ic+idxc
          do ib=0,lb
           jb=ib+idxb
           do ia=0,la
            ja=ia+idxa
            iad1=itmp+ia+lap*(ib+lbp*(ic+lcp*id))
            iad2=ixr+iq-1+nqr*(ja+lapp*(jb+lbpp*(jc+lcpp*jd)))
            bc(iad1)=bc(iad2)
           end do
          end do
         end do
        end do
        call prntm2(bc(itmp),lap*lbp,lcp*ldp,lap*lbp)
        write(6,*)('carty ')
        do id=0,ld
         jd=id+idyd
         do ic=0,lc
          jc=ic+idyc
          do ib=0,lb
           jb=ib+idyb
           do ia=0,la
            ja=ia+idya
            iad1=itmp+ia+lap*(ib+lbp*(ic+lcp*id))
            iad2=iyr+iq-1+nqr*(ja+lapp*(jb+lbpp*(jc+lcpp*jd)))
            bc(iad1)=bc(iad2)
           end do
          end do
         end do
        end do
        call prntm2(bc(itmp),lap*lbp,lcp*ldp,lap*lbp)
        write(6,*)('cartz ')
        do id=0,ld
         jd=id+idzd
         do ic=0,lc
          jc=ic+idzc
          do ib=0,lb
           jb=ib+idzb
           do ia=0,la
            ja=ia+idza
            iad1=itmp+ia+lap*(ib+lbp*(ic+lcp*id))
            iad2=izr+iq-1+nqr*(ja+lapp*(jb+lbpp*(jc+lcpp*jd)))
            bc(iad1)=bc(iad2)/bc(irys+iq-1)
           end do
          end do
         end do
        end do
        call prntm2(bc(itmp),lap*lbp,lcp*ldp,lap*lbp)
       end do
       ibcoff=itmp
      end if
      jpowd=ipp(ldp)+1                                                  1d14s10
      sz=0d0
      do id=1,ibc(ipp(ldp))                                             1d14s10
       j1=(ibc(jpowd)+idxd)*nwdsabcp                                    3d6s17
       j2=(ibc(jpowd+1)+idyd)*nwdsabcp                                  3d6s17
       j3=(ibc(jpowd+2)+idzd)*nwdsabcp                                  3d6s17
       jpowd=jpowd+3                                                    1d14s10
       jpowc=ipp(lcp)+1                                                 1d14s10
       iadd=(id-1)*ibc(ipp(lcp))                                        1d14s10
       do ic=1,ibc(ipp(lcp))                                            1d14s10
        k1=(ibc(jpowc)+idxc)*nwdsabp+j1                                 3d6s17
        k2=(ibc(jpowc+1)+idyc)*nwdsabp+j2                               3d6s17
        k3=(ibc(jpowc+2)+idzc)*nwdsabp+j3                               3d6s17
        jpowc=jpowc+3                                                   1d14s10
        jpowb=ipp(lbp)+1                                                1d14s10
        iadc=(iadd+ic-1)*ibc(ipp(lbp))                                  1d14s10
        do ib=1,ibc(ipp(lbp))                                           1d14s10
         i1=(ibc(jpowb)+idxb)*lapp+k1+idxa                              3d6s17
         i2=(ibc(jpowb+1)+idyb)*lapp+k2+idya                            3d6s17
         i3=(ibc(jpowb+2)+idzb)*lapp+k3+idza                            3d6s17
         jpowb=jpowb+3                                                  1d14s10
         jpowa=ipp(lap)+1                                               1d14s10
         iadb=(iadc+ib-1)*ibc(ipp(lap))+icarti-1                        1d14s10
         do ia=1,ibc(ipp(lap))                                          1d14s10
          l1=ixr+nqr*(ibc(jpowa)+i1)                                    3d6s17
          l2=iyr+nqr*(ibc(jpowa+1)+i2)                                  10d4s16
          l3=izr+nqr*(ibc(jpowa+2)+i3)                                  10d4s16
          jpowa=jpowa+3                                                 1d21s10
          iada=iadb+ia                                                  1d14s10
          do iq=0,nqr-1                                                 10d4s16
           bc(iada)=bc(iada)+bc(l1+iq)*bc(l2+iq)*bc(l3+iq)
          end do                                                        1d14s10
          sz=sz+bc(iada)**2
         end do                                                         1d14s10
        end do                                                          1d14s10
       end do                                                           1d14s10
      end do                                                            1d14s10
c
c     we have
      if(idbg.ne.0)then
       write(6,*)('after rys quadrature: ')
       call prntm2(bc(icarti),1,ncarti,1)
      end if
      jcarti=icarti                                                     1d14s10
      jcartt=icartt                                                     1d14s10
      nabc=ibc(ipp(lap))*ibc(ipp(lbp))*ibc(ipp(lcp))                    1d14s10
      if(idbg.ne.0)then
       write(6,*)('computing nabc from ')
       write(6,*)lap,ibc(ipp(lap))
       write(6,*)lbp,ibc(ipp(lbp))
       write(6,*)lcp,ibc(ipp(lcp))
       write(6,*)nabc
      end if
      call cart2spher(bc(icarti),bc(icartt),nabc,ld,bc,ibc)             11d14s22
      ldt=2*ld+1                                                        1d14s10
      nab=ibc(ipp(lap))*ibc(ipp(lbp))                                   1d14s10
      do i=1,ldt                                                        1d14s10
       do j=1,ibc(ipp(lcp))                                             1d14s10
        iad1=icartt-1+nab*(j-1+ibc(ipp(lcp))*(i-1))                     1d14s10
        iad2=icarti-1+nab*(i-1+ldt*(j-1))                               1d14s10
        do k=1,nab                                                      1d14s10
         bc(iad2+k)=bc(iad1+k)*front                                    12d1s16
        end do                                                          1d14s10
       end do                                                           1d14s10
      end do                                                            1d14s10
      jcarti=icarti                                                     1d14s10
      jcartt=icartt                                                     1d14s10
      nabd=nab*ldt                                                      1d14s10
      call cart2spher(bc(icarti),bc(icartt),nabd,lc,bc,ibc)             11d14s22
      lct=2*lc+1                                                        1d14s10
      lcdt=lct*ldt
      do i=1,ibc(ipp(lap))                                              1d14s10
       do j=1,ibc(ipp(lbp))                                             1d14s10
        iad1=icarti-1+lcdt*(i-1+ibc(ipp(lap))*(j-1))                    1d14s10
        do k=1,lcdt                                                     1d14s10
         iad2=icartt+i-1+ibc(ipp(lap))*(j-1+ibc(ipp(lbp))*(k-1))        1d14s10
         bc(iad1+k)=bc(iad2)                                            1d14s10
        end do                                                          1d14s10
       end do                                                           1d14s10
      end do                                                            1d14s10
      jcarti=icarti                                                     1d14s10
      jcartt=icartt                                                     1d14s10
      ncda=lcdt*ibc(ipp(lap))                                           1d14s10
      call cart2spher(bc(icarti),bc(icartt),ncda,lb,bc,ibc)             11d14s22
      lbt=2*lb+1                                                        1d14s10
      do i=1,ibc(ipp(lap))                                              1d14s10
       do j=1,lbt                                                       1d14s10
        iad1=icartt-1+lcdt*(i-1+ibc(ipp(lap))*(j-1))                    1d14s10
        iad2=icarti-1+lcdt*(j-1+lbt*(i-1))                              1d14s10
        do k=1,lcdt                                                     1d14s10
         bc(iad2+k)=bc(iad1+k)                                          1d14s10
        end do                                                          1d14s10
       end do                                                           1d14s10
      end do                                                            1d14s10
      jcarti=icarti                                                     1d14s10
      jcartt=icartt                                                     1d14s10
      ncdb=lcdt*lbt                                                     1d29s10
      call cart2spher(bc(icarti),bc(icartt),ncdb,la,bc,ibc)             11d14s22
      lat=2*la+1                                                        1d29s10
      if(idbg.ne.0)then
       write(6,*)('after cart2spher: ')
       call prntm2(bc(icartt),1,lat*lbt*lct*ldt,1)
      end if
   50 format(4i5,1pe26.18)
      jrij=icartt                                                       1d29s10
      rmsz=0d0
      do ia=1,lat                                                       1d29s10
       do ib=1,lbt                                                      1d29s10
        do ic=1,lct                                                     1d29s10
         do id=1,ldt                                                    1d29s10
          rmsz=rmsz+bc(jrij)**2
          if(bc(jrij).ne.bc(jrij))then
           write(6,*)('got nan in eri!! ')
      write(6,*)ia,ib,ic,id
      idbg=1
      go to 1066
          end if
          jrij=jrij+1                                                   1d29s10
         end do                                                         1d29s10
        end do                                                          1d29s10
       end do                                                           1d29s10
      end do                                                            1d29s10
      rmsz=sqrt(rmsz/dfloat(lat*lbt*lct*ldt))
      return
      end
