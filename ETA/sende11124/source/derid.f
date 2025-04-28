c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine derid(la,zza,xna,xa,ya,za,iadda,                       9d8s10
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               icartt,iax,iay,iaz,ibx,iby,ibz,icx,icy,icz,        9d8s10
     $               idx,idy,idz,iflag,ldebp,fmul,bc,ibc)               11d14s22
      implicit real*8 (a-h,o-z)
      external second                                                   3d16s20
c
c     compute electron repulsion integrals.
c     iax etc are order of differentiation of gaussian a with respect   9d8s10
c     to x etc                                                          9d8s10
c     type of operator:
c     iflag=1: 1/rij
c     iflag=2: x12/rij                                                  3d25s19
c     iflag=3: y12/rij                                                  3d25s19
c     iflag=4: z12/rij                                                  3d25s19
c
      include "common.store"
      include "common.spher"
      include "common.rys"                                              7d5s12
      data icall/0/
      save
      logical ldeb,ldebp                                                2d26s20
      dimension iovrm(8),ikinm(8),ioffb(8),ioffk(8),nhere(8,4),         1d13s10
     $     mhere(8,4),icoef(8,4),ipow(8,4),ltmp(4),mult(8,8),joff(8,4)  1d14s10
      dimension breit(10)                                               1d19s23
      common/xcom/iadd1,iadd2,iadd3
      data mult/1,2,3,4,5,6,7,8, 2,1,4,3,6,5,8,7,                       1d13s10
     $          3,4,1,2,7,8,5,6, 4,3,2,1,8,7,6,5,                       1d13s10
     $          5,6,7,8,1,2,3,4, 6,5,8,7,2,1,4,3,                       1d13s10
     $          7,8,5,6,3,4,1,2, 8,7,6,5,4,3,2,1/                       1d13s10
      ldeb=.false.                                                      2d14s23
      icall=icall+1
 1066 continue                                                          2d26s20
      if(ldeb)write(6,1)la,zza,xa,ya,za,lb,zzb,xb,yb,zb,lc,zzc,xc,yc,zc,        1d29s10
     $     ld,zzd,xd,yd,zd                                              1d29s10
      if(ldeb)write(6,335)iax,iay,iaz,ibx,iby,ibz,icx,icy,icz,idx,idy,
     $     idz,iflag
  335 format(4(3i1,3x),5x,i1)
      nex=0
      nexgh=0
      if(iflag.ge.2.and.iflag.le.4)nexgh=nexgh+1
    1 format('in derid: ',4(i3,1pe15.7,0p3f8.4))                          1d29s10
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
      if(iflag.gt.1)then                                                7d13s15
       iadd=3+lc+ld+nexgh                                                     7d13s15
      else                                                              7d13s15
       iadd=2+lc+ld                                                     7d13s15
      end if                                                            7d13s15
      nq1=(la+lb+iadd+max(iax,iay,iaz)+max(ibx,iby,ibz)                 7d13s15
     $     +max(icx,icy,icz)+max(idx,idy,idz))/2                        9d8s10
      nq1=nq1+nex
      nq1=nq1+iadd1                                                     7d14s15
      nq2=(iadd+max(icx,icy,icz)+max(idx,idy,idz))/2                    7d13s15
      nq2=nq2+nex
      nq2=nq2+iadd2                                                     7d14s15
      lap=la+1                                                          1d13s10
      lbp=lb+1                                                          1d13s10
      lcp=lc+1                                                          1d13s10
      ldp=ld+1                                                          1d13s10
      if(lap.gt.40)then
       write(6,*)('lap: '),lap
       call dws_sync
       call dws_finalize
       stop
      end if
      if(lbp.gt.40)then
       write(6,*)('lbp: '),lbp
       call dws_sync
       call dws_finalize
       stop
      end if
      if(lcp.gt.40)then
       write(6,*)('ldp: '),lcp
       call dws_sync
       call dws_finalize
       stop
      end if
      if(ldp.gt.40)then
       write(6,*)('ldp: '),ldp
       call dws_sync
       call dws_finalize
       stop
      end if
      ncarti=ibc(ipp(lap))*ibc(ipp(lbp))*ibc(ipp(lcp))*ibc(ipp(ldp))    1d14s10
      icartt=ibcoff                                                     1d14s10
      ibcoff=icartt+ncarti                                              1d14s10
      inode1=ibcoff
      iwgt1=inode1+nq1
      itmp=iwgt1+nq1
      ibcoff=itmp+nq1
      if(ibcoff-ioffsx.gt.maxbc)write(6,*)('enough no. 1'),inode1,nq1,  11d2s22
     $     nex,iadd,                                                    11d2s22
     $     la,lb,iax,iay,iaz,ibx,iby,ibz,icx,icy,icz,idx,idy,idz
      call enough('derid.  1',bc,ibc)
      ibcoff=itmp                                                       1d27s23
      iad1=(nq1*(nq1-1))+ibcgh                                          1d26s23
      iad2=iad1+nq1                                                     1d26s23
      jnode1=iad2-1                                                     1d26s23
      jwgt1=iad1-1                                                      1d26s23
      ibcoff=itmp
      inode2=ibcoff
      iwgt2=inode2+nq2
      itmp=iwgt2+nq2
      ibcoff=itmp+nq2
      if(ibcoff-ioffsx.gt.maxbc)write(6,*)('enough no. 2')              11d2s22
      call enough('derid.  2',bc,ibc)
      iad1=(nq2*(nq2-1))+ibcgh                                          1d26s23
      iad2=iad1+nq2                                                     1d26s23
      ibcoff=itmp
      jnode2=iad2-1                                                     1d26s23
      jwgt2=iad1-1                                                      1d26s23
      irnode=ibcoff                                                     1d7s09
      iadd=iadd+iax+iay+iaz+ibx+iby+ibz+icx+icy+icz+idx+idy+idz         7d13s15
      if(iflag.gt.1)iadd=iadd+1                                         7d14s15
      nqr=(la+lb+iadd)/2                                                7d13s15
      nqr=nqr+nex
      nqr=nqr+iadd3                                                     7d14s15
      irwgt=irnode+nqr                                                   1d7s09
      ibcoff=irwgt+nqr                                                   1d12s10
      if(ibcoff-ioffsx.gt.maxbc)write(6,*)('enough no. 3')              11d2s22
      call enough('derid.  3',bc,ibc)
      zzab=zza+zzb                                                         1d13s10
      zzcd=zzc+zzd                                                         1d13s10
      zzabcd=zzab+zzcd
      pref=(zza*zzb*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)/zzab)            1d29s10
     $     +(zzc*zzd*((xc-xd)**2+(yc-yd)**2+(zc-zd)**2)/zzcd)            1d29s10
      pref0=pref                                                        4d21s16
      pref=exp(-pref)                                                   1d22s10
      if(pref.ne.pref)then
       write(6,*)('pref is nan '),pref0
       write(6,1)la,zza,xa,ya,za,lb,zzb,xb,yb,zb,lc,zzc,xc,yc,zc,        1d29s10
     $      ld,zzd,xd,yd,zd                                              1d29s10
       write(6,335)iax,iay,iaz,ibx,iby,ibz,icx,icy,icz,idx,idy,
     $     idz
      call dws_sync
            call dws_finalize
       stop
      end if
      xab=zza*xa+zzb*xb
      yab=zza*ya+zzb*yb
      zab=zza*za+zzb*zb
      xcd=zzc*xc+zzd*xd
      ycd=zzc*yc+zzd*yd
      zcd=zzc*zc+zzd*zd
      xabcd=xab+xcd
      yabcd=yab+ycd
      zabcd=zab+zcd
      bigx=((zzcd*xab-zzab*xcd)**2+(zzcd*yab-zzab*ycd)**2
     $     +(zzcd*zab-zzab*zcd)**2)/(zzab*zzcd*zzabcd)                  1d29s10
      if(nqr.ge.nqx)then                                                6d22s12
       write(6,*)('asking for more rys quadrature points than are '),   6d22s12
     $     ('available '),nqr,nqx                                       6d22s12
      end if                                                            6d22s12
      iaddq=ibc(ibcrys+2*(nqr-1))                                       6d22s12
      nreg=ibc(ibcrys+1+2*(nqr-1))                                      6d22s12
      if(ldeb)then
       write(6,*)('iaddq, nreg '),iaddq,nreg,ibcrys+1+2*(nqr-1)
      end if
      ratio=bigx/bc(iaddq)                                              10d29s15
      if(ratio.le.nreg+1)then                                           10d29s15
       ireg=int(ratio)                                                  10d29s15
      else                                                              10d29s15
       ireg=nreg+1                                                      10d29s15
      end if                                                            10d29s15
      ireg=ireg+1                                                       6d22s12
      if(ireg.le.nreg)then                                              6d22s12
       a=dfloat(ireg-1)*bc(iaddq)                                       6d22s12
       b=a+bc(iaddq)                                                    6d22s12
       y=(2d0*bigx-a-b)/bc(iaddq)                                       6d22s12
       y2=2d0*y                                                         6d22s12
       try=bc(iaddq+2*nqr+ireg)                                         6d22s12
       itry=nint(try)                                                   6d22s12
       itry=itry+ibc(ibcrys)-1
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
      jrnode=irnode-1                                                   1d11s10
      jrwgt=irwgt-1                                                     1d11s10
      pi=acos(-1d0)                                                     1d11s10
      front0=fmul*pref*2d0/sqrt(pi)                                     8d20s15
      if(front0.ne.front0)then
       write(6,*)('front0 nan: '),fmul,pref,pi,loc(fmul)
      end if
      bsia=0.5d0/zza                                                     1d13s10
      bsib=0.5d0/zzb                                                     1d13s10
      if(zzc.ne.0d0)then                                                2d17s23
       bsic=0.5d0/zzc                                                     1d13s10
      else                                                              2d17s23
       bsic=0d0                                                         2d17s23
      end if                                                            2d17s23
      if(zzd.ne.0d0)then                                                2d17s23
       bsid=0.5d0/zzd                                                     1d13s10
      else                                                              2d17s23
       bsid=0d0                                                         2d17s23
      end if                                                            2d17s23
      rho=sqrt(zzab*zzcd/zzabcd)                                           1d13s10
      dnorm=xna*xnb*xnc*xnd                                             1d13s10
      front=dnorm*front0*rho/sqrt((zzab*zzcd)**3)                         1d13s10
      if(front.ne.front)then
       write(6,*)('starting front nan!!! '),dnorm,front0,rho,zzab,zzcd
       write(6,*)('xna ... '),xna,xnb,xnc,xnd
       write(6,1)la,zza,xa,ya,za,lb,zzb,xb,yb,zb,lc,zzc,xc,yc,zc,        1d29s10
     $     ld,zzd,xd,yd,zd                                              1d29s10
       write(6,335)iax,iay,iaz,ibx,iby,ibz,icx,icy,icz,idx,idy,
     $     idz,iflag
       stop 'derid'
      end if
      iader=iax+iay+iaz                                                 9d8s10
      zza2=zza*2d0                                                      9d8s10
      if(iader.gt.0)front=front*(zza2**iader)                           9d8s10
      ibder=ibx+iby+ibz                                                 9d8s10
      zzb2=zzb*2d0                                                      9d8s10
      if(ibder.gt.0)front=front*(zzb2**ibder)                           9d8s10
      icder=icx+icy+icz                                                 9d8s10
      zzc2=zzc*2d0                                                      9d8s10
      if(icder.gt.0)front=front*(zzc2**icder)                           9d8s10
      idder=idx+idy+idz                                                 9d8s10
      zzd2=zzd*2d0                                                      9d8s10
      if(idder.gt.0)front=front*(zzd2**idder)                           9d8s10
      if(front.ne.front)then
       write(6,*)('ending front nan!!! '),zza2,zzb2,zzc3,zzd2,iader,
     $      ibder,icder,idder
            call dws_sync
            call dws_finalize
       stop
      end if
      icartab=ibcoff                                                    1d13s10
      kcartab=icartab-1                                                 1d13s10
      nwdsab=lap*lbp                                                    7d6s15
      nwdsab3=nwdsab*3                                                  1d29s10
      nwdsabc=nwdsab*lcp                                                7d6s15
      icartaby=icartab+nwdsab                                           1d29s10
      icartabz=icartaby+nwdsab                                          1d29s10
      kcartaby=icartaby-1                                               1d29s10
      kcartabz=icartabz-1                                               1d29s10
      ibcoff=icartab+nwdsab3                                            1d29s10
      icartabcd=ibcoff                                                    1d13s10
      nwdscd=lcp*ldp                                                    7d6s15
      nwdsabcd=nwdsab*nwdscd                                            1d13s10
      nwdsabcd3=nwdsabcd*3                                              1d29s10
      ibcoff=icartabcd+nwdsabcd3                                        1d29s10
      icartabcdy=icartabcd+nwdsabcd                                     1d29s10
      icartabcdz=icartabcdy+nwdsabcd                                    1d29s10
      ifcna=ibcoff                                                      1d13s10
      ifcnb=ifcna+max(2,lap+max(iax,iay,iaz))                           9d8s10
      ifcnc=ifcnb+max(2,lbp+max(ibx,iby,ibz))                           9d8s10
      ifcnd=ifcnc+max(2,lcp+max(icx,icy,icz))                           9d8s10
      ibcoff=ifcnd+max(2,ldp+max(idx,idy,idz))                          9d8s10
      icarti=ibcoff
      ibcoff=icarti+ncarti                                              1d14s10
      if(ibcoff-ioffsx.gt.maxbc)write(6,*)('enough no. 4')              11d2s22
      call enough('derid.  4',bc,ibc)
      do i=0,ncarti-1                                                   1d14s10
       bc(icarti+i)=0d0                                                 1d14s10
      end do                                                            1d14s10
      do irq=1,nqr                                                      1d14s10
       use2=rho*rho*bc(jrnode+irq)                                      7d5s12
       rwgt=bc(jrwgt+irq)*front                                         1d14s10
       if(rwgt.ne.rwgt)then
        write(6,*)('rwgt is nan for irq = '),irq,('!!!')
        if(idbg.eq.0.or..not.ldeb)then                                  2d26s20
         idbg=1                                                         2d26s20
         ldeb=.true.                                                    2d26s20
         go to 1066                                                      2d26s20
        else                                                            2d26s20
         call dws_sync
         call dws_finalize
         stop
        end if                                                          2d26s20
       end if
       as1=1d0/(use2+zzcd)                                               1d14s10
       as2=use2*as1                                                      1d14s10
       ss2x=xcd*as1                                                     1d29s10
       ss2y=ycd*as1                                                     1d29s10
       ss2z=zcd*as1                                                     1d29s10
       as1=sqrt(as1)                                                    1d14s10
       t3i=1d0/(zzab*zzcd+zzabcd*use2)                                  1d29s10
       as3=sqrt((zzcd+use2)*t3i)                                        1d29s10
       ss3x=(use2*xabcd+zzcd*xab)*t3i                                   1d29s10
       ss3y=(use2*yabcd+zzcd*yab)*t3i                                   1d29s10
       ss3z=(use2*zabcd+zzcd*zab)*t3i                                   1d29s10
   55  format('rysq: ',i5,1p2e15.7)
       do i=0,nwdsabcd3-1                                               1d29s10
        bc(icartabcd+i)=0d0                                             1d14s10
       end do                                                           1d14s10
       if(ldeb)then                                                     3d16s20
        call second(time2)                                              3d16s20
        call second(time1)                                              3d16s20
        tovr=time1-time2                                                3d16s20
        tsto=0d0                                                        3d16s20
       end if                                                           3d16s20
       do iq1=1,nq1                                                     1d14s10
        arg=bc(jnode1+iq1)*as3                                          1d29s10
        arge1=arg                                                       2d12s20
        argx=arg+ss3x                                                   1d29s10
        arga=argx-xa
        argb=argx-xb
        term=bc(jwgt1+iq1)                                              1d14s10
        bc(ifcna)=1d0
        bc(ifcnb)=term
        bc(ifcna+1)=-arga                                               1d29s10
        ff=bsia
        do i=2,la+iax                                                   9d8s10
         bc(ifcna+i)=-arga*bc(ifcna+i-1)-ff*bc(ifcna+i-2)               2d16s10
         ff=ff+bsia
        end do
        bc(ifcnb+1)=-argb*term                                          1d29s10
        ff=bsib
        do i=2,lb+ibx                                                   9d8s10
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
         ff=ff+bsib
        end do
        jcartab=icartab                                                 1d14s10
        do ib=0,lb                                                      7d6s15
         do ia=0,la                                                     7d6s15
          bc(jcartab)=bc(ifcna+ia+iax)*bc(ifcnb+ib+ibx)                 7d6s15
          jcartab=jcartab+1                                             1d14s10
         end do                                                         1d14s10
        end do                                                          1d14s10
        argy=arg+ss3y                                                   1d29s10
        arga=argy-ya
        argb=argy-yb
        term=bc(jwgt1+iq1)                                              1d14s10
        bc(ifcna)=1d0
        bc(ifcnb)=term
        bc(ifcna+1)=-arga                                               1d29s10
        ff=bsia
        do i=2,la+iay                                                   9d8s10
         bc(ifcna+i)=-arga*bc(ifcna+i-1)-ff*bc(ifcna+i-2)               2d16s10
         ff=ff+bsia
        end do
        bc(ifcnb+1)=-argb*term                                          1d29s10
        ff=bsib
        do i=2,lb+iby                                                   9d8s10
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10(
         ff=ff+bsib
        end do
        jcartab=icartaby                                                1d29s10
        do ib=0,lb                                                      7d6s15
         do ia=0,la                                                     7d6s15
          bc(jcartab)=bc(ifcna+ia+iay)*bc(ifcnb+ib+iby)                 7d6s15
          jcartab=jcartab+1                                             1d14s10
         end do                                                         1d14s10
        end do                                                          1d14s10
        argz=arg+ss3z                                                   1d29s10
        arga=argz-za
        argb=argz-zb
        term=bc(jwgt1+iq1)                                              1d14s10
        bc(ifcna)=1d0
        bc(ifcnb)=term
        bc(ifcna+1)=-arga                                               1d29s10
        ff=bsia
        do i=2,la+iaz                                                   9d8s10
         bc(ifcna+i)=-arga*bc(ifcna+i-1)-ff*bc(ifcna+i-2)               2d16s10(
         ff=ff+bsia
        end do
        bc(ifcnb+1)=-argb*term                                          1d29s10
        ff=bsib
        do i=2,lb+ibz                                                   9d8s10
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10(
         ff=ff+bsib
        end do
        jcartab=icartabz                                                1d29s10
        do ib=0,lb                                                      7d6s15
         do ia=0,la                                                     7d6s15
          bc(jcartab)=bc(ifcna+ia+iaz)*bc(ifcnb+ib+ibz)                 7d6s15
          jcartab=jcartab+1                                             1d14s10
         end do                                                         1d14s10
        end do                                                          1d14s10
        if(iflag.ne.0)then                                              3d18s20
         shftx=argx*as2+ss2x                                            3d18s20
         shfty=argy*as2+ss2y                                            3d18s20
         shftz=argz*as2+ss2z                                            3d18s20
        else                                                            3d18s20
c
c     we do this to avoid numberical problems when using finite nucleus.
c     but it seems to screw up breit integrals!
c
         shftx=argx*as2+(ss2x-xc)                                        7d14s15
         deltax=xc-xd                                                    7d14s15
         shfty=argy*as2+(ss2y-yc)                                        7d14s15
         deltay=yc-yd                                                    7d14s15
         shftz=argz*as2+(ss2z-zc)                                        7d14s15
         deltaz=zc-zd                                                    7d14s15
        end if                                                          3d18s20
        do iq2=1,nq2                                                    1d14s10
         arg=bc(jnode2+iq2)*as1
         breit12=arge1-arg                                              2d12s20
         argx2=arg+shftx                                                7d13s15
         if(iflag.ne.0)then                                             3d18s20
          argc=argx2-xc                                                  7d13s15
          argd=argx2-xd                                                  7d13s15
         else                                                           3d18s20
          argc=argx2                                                     7d14s15
          argd=argc+deltax                                               7d14s15
         end if                                                         3d18s20
         breit(1)=1d0
         x12=argx-argx2                                                  7d13s15
         breit(2)=x12                                                   3d25s19
         breit(3)=1d0
         breit(4)=1d0
         term=bc(jwgt2+iq2)*breit(iflag)                                7d14s15
         bc(ifcnc)=1d0                                                  1d14s10
         bc(ifcnd)=term
         bc(ifcnc+1)=-argc                                              1d29s10
         ff=bsic                                                        1d14s10
         do i=2,lc+icx                                                  9d8s10
          bc(ifcnc+i)=-argc*bc(ifcnc+i-1)-ff*bc(ifcnc+i-2)              2d16s10(
          ff=ff+bsic                                                    1d14s10
         end do                                                         1d14s10
         bc(ifcnd+1)=-argd*term                                         7d14s15
         ff=bsid                                                        1d14s10
         do i=2,ld+idx                                                  9d8s10
          bc(ifcnd+i)=-argd*bc(ifcnd+i-1)-ff*bc(ifcnd+i-2)              2d16s10(
          ff=ff+bsid                                                    1d14s10
         end do                                                         1d14s10
         jcartabcd=icartabcd                                            1d14s10
         do id=0,ld                                                     7d6s15
          do ic=0,lc                                                    7d6s15
           ff=bc(ifcnc+ic+icx)*bc(ifcnd+id+idx)                         7d6s15
           do iab=1,nwdsab                                              1d14s10
            bc(jcartabcd)=bc(jcartabcd)+ff*bc(kcartab+iab)              1d14s10
            jcartabcd=jcartabcd+1                                       1d14s10
           end do                                                       1d14s10
          end do                                                        1d14s10
         end do                                                         1d14s10
         argy2=arg+shfty                                                7d13s15
         y12=argy-argy2                                                  7d13s15
         breit(2)=1d0
         breit(3)=y12                                                   3d25s19
         if(iflag.ne.0)then                                             3d18s20
          argc=argy2-yc                                                  7d13s15
          argd=argy2-yd                                                  7d13s15
         else                                                           3d18s20
          argc=argy2                                                     7d14s15
          argd=argc+deltay                                               7d14s15
         end if                                                         3d18s20
         term=bc(jwgt2+iq2)*breit(iflag)                                7d14s15
         bc(ifcnc)=1d0                                                  1d14s10
         bc(ifcnd)=term                                                 1d14s10
         bc(ifcnc+1)=-argc                                              1d29s10
         ff=bsic                                                        1d14s10
         do i=2,lc+icy                                                  9d8s10
          bc(ifcnc+i)=-argc*bc(ifcnc+i-1)-ff*bc(ifcnc+i-2)              2d16s10
          ff=ff+bsic                                                    1d14s10
         end do                                                         1d14s10
         bc(ifcnd+1)=-argd*term                                         1d29s10
         ff=bsid                                                        1d14s10
         do i=2,ld+idy                                                  9d8s10
          bc(ifcnd+i)=-argd*bc(ifcnd+i-1)-ff*bc(ifcnd+i-2)              2d16s10
          ff=ff+bsid                                                    1d14s10
         end do                                                         1d14s10
         jcartabcd=icartabcdy                                           1d29s10
         do id=0,ld                                                     7d6s15
          do ic=0,lc                                                    7d6s15
           ff=bc(ifcnc+ic+icy)*bc(ifcnd+id+idy)                         7d6s15
           do iab=1,nwdsab                                              1d14s10
            bc(jcartabcd)=bc(jcartabcd)+ff*bc(kcartaby+iab)             1d29s10
            jcartabcd=jcartabcd+1                                       1d14s10
           end do                                                       1d14s10
          end do                                                        1d14s10
         end do                                                         1d14s10
         argz2=arg+shftz                                                7d13s15
         z12=argz-argz2                                                  7d13s15
         breit(3)=1d0
         breit(4)=z12                                                   3d25s19
         if(iflag.ne.0)then                                             3d18s20
          argc=argz2-zc                                                  7d13s15
          argd=argz2-zd                                                  7d13s15
         else                                                           3d18s20
          argc=argz2                                                     7d14s15
          argd=argc+deltaz                                               7d14s15
         end if                                                         3d18s20
         term=bc(jwgt2+iq2)*breit(iflag)                                7d14s15
         bc(ifcnc)=1d0                                                  1d14s10
         bc(ifcnd)=term                                                 1d14s10
         bc(ifcnc+1)=-argc                                              1d29s10
         ff=bsic                                                        1d14s10
         do i=2,lc+icz                                                  9d8s10
          bc(ifcnc+i)=-argc*bc(ifcnc+i-1)-ff*bc(ifcnc+i-2)              2d16s10
          ff=ff+bsic                                                    1d14s10
         end do                                                         1d14s10
         bc(ifcnd+1)=-argd*term                                         1d29s10
         ff=bsid                                                        1d14s10
         do i=2,ld+idz                                                  9d8s10
          bc(ifcnd+i)=-argd*bc(ifcnd+i-1)-ff*bc(ifcnd+i-2)              2d16s10
          ff=ff+bsid                                                    1d14s10
         end do                                                         1d14s10
         jcartabcd=icartabcdz                                           1d29s10
         if(ldeb)then
          write(6,*)('for iq1 = '),iq1
          call prntm2(bc(icartab),lap,lbp,lap)
          call prntm2(bc(icartaby),lap,lbp,lap)
          call prntm2(bc(icartabz),lap,lbp,lap)
          write(6,*)('fcnc '),icz,term,bc(jwgt2+iq2),breit(iflag)
          call prntm2(bc(ifcnc+icz),1,lcp,1)
          write(6,*)('fcnd '),idz,argz,argz2
          call prntm2(bc(ifcnd+idz),1,ldp,1)
         end if
         do id=0,ld                                                     7d6s15
          do ic=0,lc                                                    7d6s15
           ff=bc(ifcnc+ic+icz)*bc(ifcnd+id+idz)                         7d6s15
           do iab=1,nwdsab                                              1d14s10
            bc(jcartabcd)=bc(jcartabcd)+ff*bc(kcartabz+iab)             1d29s10
            jcartabcd=jcartabcd+1                                       1d14s10
           end do                                                       1d14s10
          end do                                                        1d14s10
         end do                                                         1d14s10
        end do                                                          1d14s10
 3303   format(7es15.7)
       end do                                                           1d14s10
       rwgtu=rwgt                                                       7d13s15
       if(ldeb)call second(timea)                                       3d16s20
       jpowd=ipp(ldp)+1                                                 1d14s10
       do id=1,ibc(ipp(ldp))                                            1d14s10
        j1=ibc(jpowd)*nwdsabc+icartabcd                                 7d6s15
        j2=ibc(jpowd+1)*nwdsabc+icartabcdy                              7d6s15
        j3=ibc(jpowd+2)*nwdsabc+icartabcdz                              7d6s15
        jpowd=jpowd+3                                                   1d14s10
        jpowc=ipp(lcp)+1                                                1d14s10
        iadd=(id-1)*ibc(ipp(lcp))                                       1d14s10
        do ic=1,ibc(ipp(lcp))                                           1d14s10
         k1=ibc(jpowc)*nwdsab+j1                                        7d6s15
         k2=ibc(jpowc+1)*nwdsab+j2                                      7d6s15
         k3=ibc(jpowc+2)*nwdsab+j3                                      7d6s15
         jpowc=jpowc+3                                                  1d14s10
         jpowb=ipp(lbp)+1                                               1d14s10
         iadc=(iadd+ic-1)*ibc(ipp(lbp))                                 1d14s10
         do ib=1,ibc(ipp(lbp))                                          1d14s10
          i1=ibc(jpowb)*lap+k1                                          7d6s15
          i2=ibc(jpowb+1)*lap+k2                                        7d6s15
          i3=ibc(jpowb+2)*lap+k3                                        7d6s15
          jpowb=jpowb+3                                                 1d14s10
          jpowa=ipp(lap)+1                                              1d14s10
          iadb=(iadc+ib-1)*ibc(ipp(lap))+icarti-1                       1d14s10
          do ia=1,ibc(ipp(lap))                                         1d14s10
           l1=ibc(jpowa)+i1                                             7d6s15
           l2=ibc(jpowa+1)+i2                                           7d6s15
           l3=ibc(jpowa+2)+i3                                           7d6s15
           jpowa=jpowa+3                                                1d21s10
           iada=iadb+ia                                                 1d14s10
           bc(iada)=bc(iada)+bc(l1)*bc(l2)*bc(l3)*rwgtu                 7d13s15
          end do                                                        1d14s10
         end do                                                         1d14s10
        end do                                                          1d14s10
       end do                                                           1d14s10
       if(ldeb)then                                                     3d16s20
        call second(timeb)                                              3d16s20
        telap=timeb-timea-tovr                                          3d16s20
        tsto=tsto+telap                                                 3d16s20
       end if                                                           3d16s20
      end do                                                            1d14s10
      if(ldeb)then                                                      3d16s20
       call second(time2)                                               3d16s20
       telap=time2-time1-tovr                                           3d16s20
       write(6,*)('total time for nqr: '),telap                         3d16s20
       write(6,*)('time for store: '),tsto                              3d16s20
      end if                                                            3d16s20
      sz0=0d0                                                           1d21s10
      do i=1,ncarti                                                     1d21s10
       sz0=sz0+bc(icarti+i-1)**2                                        1d21s10
      end do                                                            1d21s10
      sz0=sqrt(sz0/dfloat(ncarti))                                      1d21s10
      if(sz0.ne.sz0)then
       write(6,*)('not a number!!! '),sz0,ncarti
            call dws_sync
            call dws_finalize
       stop
      end if
      jcarti=icarti                                                     1d14s10
      jcartt=icartt                                                     1d14s10
      nabc=ibc(ipp(lap))*ibc(ipp(lbp))*ibc(ipp(lcp))                    1d14s10
      if(ldeb)call second(time1)                                        3d16s20
      if(iusecart.eq.0)then                                             2d21s20
       do isymd=1,8                                                      1d14s10
        if(nhere(isymd,4).gt.0)then                                      1d14s10
         call dgemm('n','n',nabc,mhere(isymd,4),nhere(isymd,4),1d0,      1d14s10
     $      bc(jcarti),nabc,bc(icoef(isymd,4)),nhere(isymd,4),0d0,      1d14s10
     $      bc(jcartt),nabc,                                            1d14s10
     d' derid.  1')
         jcarti=jcarti+nabc*nhere(isymd,4)                               1d14s10
         jcartt=jcartt+nabc*mhere(isymd,4)                               1d14s10
        end if                                                           1d14s10
       end do                                                            1d14s10
       nab=ibc(ipp(lap))*ibc(ipp(lbp))                                   1d14s10
       ldt=2*ld+1                                                        1d14s10
       do i=1,ldt                                                        1d14s10
        do j=1,ibc(ipp(lcp))                                             1d14s10
         iad1=icartt-1+nab*(j-1+ibc(ipp(lcp))*(i-1))                     1d14s10
         iad2=icarti-1+nab*(i-1+ldt*(j-1))                               1d14s10
         do k=1,nab                                                      1d14s10
          bc(iad2+k)=bc(iad1+k)                                          1d14s10
         end do                                                          1d14s10
        end do                                                           1d14s10
       end do                                                            1d14s10
       jcarti=icarti                                                     1d14s10
       jcartt=icartt                                                     1d14s10
       nabd=nab*ldt                                                      1d14s10
       do isymc=1,8                                                      1d14s10
        if(nhere(isymc,3).gt.0)then                                      1d14s10
         call dgemm('n','n',nabd,mhere(isymc,3),nhere(isymc,3),1d0,      1d14s10
     $       bc(jcarti),nabd,bc(icoef(isymc,3)),nhere(isymc,3),0d0,      1d14s10
     $       bc(jcartt),nabd,                                            1d14s10
     d' derid.  2')
         jcarti=jcarti+nabd*nhere(isymc,3)                               1d14s10
         jcartt=jcartt+nabd*mhere(isymc,3)                               1d14s10
        end if                                                           1d14s10
       end do                                                            1d14s10
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
       do isymb=1,8                                                      1d14s10
        if(nhere(isymb,2).gt.0)then                                      1d14s10
         call dgemm('n','n',ncda,mhere(isymb,2),nhere(isymb,2),1d0,      1d14s10
     $       bc(jcarti),ncda,bc(icoef(isymb,2)),nhere(isymb,2),0d0,      1d14s10
     $       bc(jcartt),ncda,                                            1d14s10
     d' derid.  3')
         jcarti=jcarti+ncda*nhere(isymb,2)                               1d14s10
         jcartt=jcartt+ncda*mhere(isymb,2)                               1d14s10
        end if                                                           1d14s10
       end do                                                            1d14s10
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
       do isyma=1,8                                                      1d14s10
        if(nhere(isyma,1).gt.0)then                                      1d14s10
         call dgemm('n','n',ncdb,mhere(isyma,1),nhere(isyma,1),1d0,      1d29s10
     $       bc(jcarti),ncdb,bc(icoef(isyma,1)),nhere(isyma,1),0d0,      1d29s10
     $      bc(jcartt),ncdb,                                            1d29s10
     d' derid.  4')
         jcarti=jcarti+ncdb*nhere(isyma,1)                               1d29s10
         jcartt=jcartt+ncdb*mhere(isyma,1)                               1d29s10
        end if                                                           1d14s10
       end do                                                            1d14s10
       if(ldeb)then
        call second(time2)                                              3d16s20
        telap=time2-time1-tovr                                          3d16s20
        write(6,*)('time for cart to sphere '),telap                    3d16s20
        if(nbasall.eq.10)then
         call dws_sync
         call dws_finalize
         stop
        end if
       end if
      else                                                              2d21s20
       icartt=icarti                                                    2d21s20
      end if                                                            2d21s20
      ibcoff=inode1
      return
      end
