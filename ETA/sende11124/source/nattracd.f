c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine nattracd(lb,zbra,xnbra,xb,yb,zb,lk,zket,xnket,xk,yk,zk,2d4s16
     $                   atnum,xi,yi,zi,icarti2,                        2d4s16
     $                  ibx,iby,ibz,ikx,iky,ikz,ifcn,ipx,ipb,ipk,xmult, 6d30s16
     $                  ldeb3,xib,bc,ibc)                               11d10s22
c
c     like nattrac except compute derivative of bra (ibx,iby,ibz),
c     and/or ket(ikx,iky,ikz), and/or 1/r
c     ifcn=0 means no der of 1/r
c     with ifcn=1 means der wrt xi, 2 means der wrt yi, and 3 means    6d30s16
c     der wrt zi.                                                      6d30s16
c     ipx=0 means no der of bigx part, ipx=1 means first der of bigx     6d30s16
c     part, and ipx=2 means 2nd der of bigx part.                       6d30s16
c     ipb=0 means no der of bra s part, ipb=1 means first der of bra   6d30s16
c     s part, and ipb=2 means second der of bra s part                 6d30s16
c     ipk=0 means no der of bra s part, ipk=1 means first der of ket   6d30s16
c     s part, and ipk=2 means second der of ket s part                 6d30s16
c
c
      implicit real*8 (a-h,o-z)
      logical ldeb,ldeb2,ldeb3                                          4d25s16
      include "common.store"
      include "common.spher"
      include "common.rys"                                              6d22s12
      ldeb=ldeb3                                                        4d25s16
      ibder=ibx+iby+ibz
      ikder=ikx+iky+ikz
      lbp=lb+1
      lkp=lk+1
      zz=-atnum                                                         1d19s10
      ldeb2=.false.
      ldeb2=ldeb3
      ldeb=ldeb3                                                        6d24s16
      if(ldeb)write(6,1)atnum,xi,yi,zi,lb,zbra,xb,yb,zb,lk,zket,xk,yk,zk    1d26s10
    1 format('in nattracd: ',f5.1,3f8.4,i3,1pe15.7,0p3f8.4,i3,1pe15.7,   1d26s10
     $         0p3f8.4)                                                 1d26s10
      if(ldeb)write(6,*)('xmult = '),xmult,ibx,iby,ibz,ikx,iky,ikz,ifcn,
     $     ipx,ipb,ipk
      if(ldeb)write(6,*)('atnum '),atnum,xnbra,xnket
      iadd=ibder+ikder
      nq=(lb+lk+2+iadd)/2
      nq=nq+4                                                           11d29s22
      nq0=nq
      npb=ibc(ipp(lbp))                                                 1d26s10
      npk=ibc(ipp(lkp))                                                 1d26s10
      icarti2=ibcoff                                                    1d22s10
      ibcoff=icarti2+npb*npk                                            1d26s10
      call enough('nattracd.  1',bc,ibc)
      inode=ibcoff
      iwgt=inode+nq
      itmp=iwgt+nq
      ibcoff=itmp+nq
      call enough('nattracd.  2',bc,ibc)
      ibcoff=itmp
      iad1=(nq*(nq-1))+ibcgh                                            1d19s23
      iad2=iad1+nq
      jnode=iad2-1                                                      1d19s23
      jwgt=iad1-1                                                       1d19s23
      irnode=ibcoff                                                     1d7s09
      nqr=nq0                                                            1d12s10
      if(ifcn.ne.0)nqr=nqr+1                                            6d30s16
      if(ldeb3)write(6,*)('nq,nqr: '),nq,nqr
      irwgt=irnode+nqr                                                   1d7s09
      ibcoff=irwgt+nqr                                                   1d12s10
      call enough('nattracd.  3',bc,ibc)
      zbki=1d0/(zbra+zket)                                              1d22s10
      sftx=(zbra*xb+zket*xk)*zbki                                       1d22s10
      sfty=(zbra*yb+zket*yk)*zbki                                       1d22s10
      sftz=(zbra*zb+zket*zk)*zbki                                       1d22s10
      pref=zbra*zket*zbki*((xk-xb)**2+(yk-yb)**2+(zk-zb)**2)            1d22s10
      pref=exp(-pref)                                                   1d22s10
      pref=pref*xmult                                                   6d29s16
      if(ibder.gt.0)then
       zbra2=zbra*2d0
       pref=pref*(zbra2**ibder)
      end if
      if(ikder.gt.0)then
       zket2=zket*2d0
       pref=pref*(zket2**ikder)
      end if
      bigx=(zbra+zket)*((sftx-xi)**2+(sfty-yi)**2+(sftz-zi)**2)         1d26s10
      if(ipx.eq.0)then                                                  6d30s16
       dbigx=0d0                                                        6d30s16
       dbigx2=0d0                                                       6d30s16
      else if(ipx.eq.1)then                                             6d30s16
       if(ifcn.eq.1)then                                                 6d29s16
        dbigx=(zbra+zket)*2d0*(sftx-xi)                                 6d29s16
       else if(ifcn.eq.2)then                                            6d29s16
        dbigx=(zbra+zket)*2d0*(sfty-yi)                                 6d29s16
       else if(ifcn.eq.3)then                                            6d29s16
        dbigx=(zbra+zket)*2d0*(sftz-zi)                                 6d29s16
       end if                                                            6d29s16
      else if(ipx.eq.2)then                                             6d30s16
       dbigx=-2d0*(zbra+zket)                                           6d30s16
       if(ifcn.eq.1)then                                                 6d29s16
        dbigx2=((zbra+zket)*2d0*(sftx-xi))**2                           6d30s16
       else if(ifcn.eq.2)then                                            6d29s16
        dbigx2=((zbra+zket)*2d0*(sfty-yi))**2                           6d30s16
       else if(ifcn.eq.3)then                                            6d29s16
        dbigx2=((zbra+zket)*2d0*(sftz-zi))**2                           6d30s16
       end if                                                            6d29s16
       if(ldeb2)write(6,*)('dbigx = '),dbigx
       if(ldeb2)write(6,*)('dbigx2 = '),dbigx2
      end if                                                            6d30s16
      if(ldeb3)write(6,*)('bigx: '),bigx
      if(nqr.ge.nqx)then                                                6d22s12
       write(6,*)('asking for more rys quadrature points than are '),   6d22s12
     $     ('available '),nqr,nqx                                       6d22s12
       call dws_finalize                                                6d22s12
       stop                                                             6d22s12
      end if                                                            6d22s12
      iaddq=ibc(ibcrys+2*(nqr-1))                                       6d22s12
      nreg=ibc(ibcrys+1+2*(nqr-1))                                      6d22s12
      reg=min(dfloat(nreg*2+1),bigx/bc(iaddq))                          1d19s23
      ireg=reg                                                          4d4s23
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
        xnd=y*dn-ddn+0.5d0*bc(j2)                                       6d25s12
        itry=itry+2*ncoef                                               6d22s12
        bc(irwgt+i-1)=wgt                                               7d5s12
        bc(irnode+i-1)=xnd                                              7d5s12
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
      sftx=zbra*xb+zket*xk                                              1d27s10
      sfty=zbra*yb+zket*yk                                              1d27s10
      sftz=zbra*zb+zket*zk                                              1d27s10
      pi=acos(-1d0)                                                     1d11s10
      front0=zz*pref*2d0/sqrt(pi)                                       1d26s10
      icartx=ibcoff
      nwds=lbp*lkp
      icarty=icartx+nwds                                                1d26s10
      icartz=icarty+nwds                                                1d26s10
      nwds3=nwds*3                                                      1d26s10
      ibcoff=icartx+nwds3                                               1d26s10
      ifcnb=ibcoff
      ifcnk=ifcnb+max(3,lbp+ibder)                                            1d14s10
      ibcoff=ifcnk+max(3,lkp+ikder)                                           1d14s10
      if(ipb.gt.0)then                                                  6d30s16
       idfcnb=ibcoff                                                    6d30s16
       ibcoff=idfcnb+ipb*max(3,lbp+ibder)                               6d30s16
       idfcnb2=idfcnb+max(3,lbp+ibder)                                  6d30s16
      end if                                                            6d30s16
      if(ipk.gt.0)then                                                  6d30s16
       idfcnk=ibcoff                                                    6d30s16
       ibcoff=idfcnk+ipk*max(3,lkp+ikder)                               6d30s16
       idfcnk2=idfcnk+max(3,lkp+ikder)                                  6d30s16
      end if                                                            6d30s16
      icarti=ibcoff                                                     1d14s10
      nwdi=npb*npk                                                      1d26s10
      ibcoff=icarti+nwdi                                                1d14s10
      call enough('nattracd.  4',bc,ibc)
      bsik=0.5d0/zket
      bsib=0.5d0/zbra
      rho=sqrt(zket+zbra)                                               1d14s10
      dnorm=xnbra*xnket
      front=dnorm*front0/(rho*rho)                                      1d14s10
  100 continue
      do i=0,nwdi-1                                                     1d14s10
       bc(icarti+i)=0d0                                                 1d14s10
      end do                                                            1d14s10
      idcarti=npb                                                       1d26s10
      do irq=1,nqr                                                      1d14s10
       use2=rho*rho*bc(jrnode+irq)                                      7d5s12
       t2=bc(jrnode+irq)/(1d0+bc(jrnode+irq))                           6d29s16
       asi2=1d0/(use2+zbra+zket)                                        1d26s10
       asi=sqrt(asi2)                                                   1d26s10
       trial=use2*asi2
       rwgt=bc(jrwgt+irq)*front                                         1d14s10
       if(ldeb2)write(6,*)('u = '),use2,asi,rho,bc(jrnode+irq),rwgt,
     $      bc(jrwgt+irq),t2,trial,t2-trial,jrnode+irq,irq,ireg,nreg,
     $      iaddq,ibcrys,nqr,ibcrys+2*(nqr-1)
       sx=(sftx+use2*xib)*asi2                                           1d26s10
       sy=(sfty+use2*yi)*asi2                                           1d26s10
       sz=(sftz+use2*zi)*asi2                                           1d26s10
       do i=0,nwds3-1                                                   1d26s10
        bc(icartx+i)=0d0
       end do
       do iq=1,nq
        arg=bc(jnode+iq)*asi
        argx=arg+sx                                                     1d26s10
        argb=argx-xb                                                    1d26s10
        argk=argx-xk                                                    1d26s10
        term=bc(jwgt+iq)                                                1d14s10
        bc(ifcnb)=1d0
        bc(ifcnk)=term
        bc(ifcnb+1)=-argb
        ff=bsib
        if(ipb.eq.1.and.ifcn.eq.1)then
         bc(idfcnb)=0d0                                                 6d30s16
         bc(idfcnb+1)=-trial                                            6d30s16
         do i=2,lb+ibx
          bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
          bc(idfcnb+i)=-trial*bc(ifcnb+i-1)-argb*bc(idfcnb+i-1)         6d30s16
     $         -ff*bc(idfcnb+i-2)                                       6d30s16
          ff=ff+bsib
         end do
         ifcnbu=idfcnb                                                  6d30s16
         if(ldeb2)then
          write(6,*)('fcn ')
          call prntm2(bc(ifcnb),1,max(2,lb+ibx),1)
          write(6,*)('dfcn ')
          call prntm2(bc(idfcnb),1,max(2,lb+ibx),1)
          write(6,*)('ufcn ')
          call prntm2(bc(ifcnbu),1,max(2,lb+ibx),1)
         end if
        else if(ipb.eq.2.and.ifcn.eq.1)then                             6d30s16
         bc(idfcnb)=0d0                                                 6d30s16
         bc(idfcnb+1)=-trial                                            6d30s16
         bc(idfcnb2)=0d0                                                6d30s16
         bc(idfcnb2+1)=0d0                                              6d30s16
         do i=2,lb+ibx
          bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
          bc(idfcnb+i)=-trial*bc(ifcnb+i-1)-argb*bc(idfcnb+i-1)         6d30s16
     $         -ff*bc(idfcnb+i-2)                                       6d30s16
          bc(idfcnb2+i)=-trial*bc(idfcnb+i-1)*2d0-argb*bc(idfcnb2+i-1)  6d30s16
     $         -ff*bc(idfcnb2+i-2)                                      6d30s16
          ff=ff+bsib
         end do
         ifcnbu=idfcnb2                                                 6d30s16
         if(ldeb2)then
          call prntm2(bc(ifcnb),1,max(2,lb+ibx),1)
          call prntm2(bc(idfcnb),1,max(2,lb+ibx),1)
          call prntm2(bc(idfcnb2),1,max(2,lb+ibx),1)
          call prntm2(bc(ifcnbu),1,max(2,lb+ibx),1)
         end if
        else
         do i=2,lb+ibx
          bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
          ff=ff+bsib
         end do
         ifcnbu=ifcnb                                                   6d30s16
        end if                                                          6d30s16
        bc(ifcnk+1)=-argk*term
        ff=bsik
        if(ipk.eq.1.and.ifcn.eq.1)then                                  6d30s16
         bc(idfcnk)=0d0                                                 6d30s16
         bc(idfcnk+1)=-term*trial                                       6d30s16
         do i=2,lk+ikx
          bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
          bc(idfcnk+i)=-trial*bc(ifcnk+i-1)-argk*bc(idfcnk+i-1)         6d30s16
     $         -ff*bc(idfcnk+i-2)                                       6d30s16
          ff=ff+bsik
         end do
         ifcnku=idfcnk                                                  6d30s16
         if(ldeb2)then
          write(6,*)('fcn ')
          call prntm2(bc(ifcnk),1,max(2,lk+ikx),1)
          write(6,*)('dfcn ')
          call prntm2(bc(idfcnk),1,max(2,lk+ikx),1)
          write(6,*)('ufcn ')
          call prntm2(bc(ifcnku),1,max(2,lk+ikx),1)
         end if
        else if(ipk.eq.2.and.ifcn.eq.1)then                             6d30s16
         bc(idfcnk)=0d0                                                 6d30s16
         bc(idfcnk+1)=-term*trial                                       6d30s16
         bc(idfcnk2)=0d0                                                6d30s16
         bc(idfcnk2+1)=0d0                                              6d30s16
         do i=2,lk+ikx
          bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
          bc(idfcnk+i)=-trial*bc(ifcnk+i-1)-argk*bc(idfcnk+i-1)         6d30s16
     $         -ff*bc(idfcnk+i-2)                                       6d30s16
          bc(idfcnk2+i)=-2d0*trial*bc(idfcnk+i-1)-argk*bc(idfcnk2+i-1)  6d30s16
     $         -ff*bc(idfcnk2+i-2)                                      6d30s16
          ff=ff+bsik
         end do
         ifcnku=idfcnk2                                                 6d30s16
        else                                                            6d30s16
         do i=2,lk+ikx
          bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
          ff=ff+bsik
         end do
         ifcnku=ifcnk                                                   6d30s16
        end if                                                          6d30s16
        jcartx=icartx
        do ik=0,lk
         do ib=0,lb
          bc(jcartx)=bc(jcartx)+bc(ifcnbu+ib+ibx)*bc(ifcnku+ik+ikx)
          jcartx=jcartx+1
         end do
        end do
        if(ldeb2)then
         write(6,*)('cartx ')
         call prntm2(bc(icartx),lb+1,lk+1,lb+1)
        end if
        argy=arg+sy                                                     1d26s10
        argb=argy-yb                                                    1d26s10
        argk=argy-yk                                                    1d26s10
        term=bc(jwgt+iq)                                                1d14s10
        bc(ifcnb)=1d0
        bc(ifcnk)=term
        bc(ifcnb+1)=-argb
        ff=bsib
        if(ipb.eq.1.and.ifcn.eq.2)then                                  6d30s16
         bc(idfcnb)=0d0                                                 6d30s16
         bc(idfcnb+1)=-trial                                            6d30s16
         do i=2,lb+iby
          bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
          bc(idfcnb+i)=-trial*bc(ifcnb+i-1)-argb*bc(idfcnb+i-1)         6d30s16
     $         -ff*bc(idfcnb+i-2)                                       6d30s16
          ff=ff+bsib
         end do
         ifcnbu=idfcnb                                                  6d30s16
        else if(ipb.eq.2.and.ifcn.eq.2)then                             6d30s16
         bc(idfcnb)=0d0                                                 6d30s16
         bc(idfcnb+1)=-trial                                            6d30s16
         bc(idfcnb2)=0d0                                                6d30s16
         bc(idfcnb2+1)=0d0                                              6d30s16
         do i=2,lb+iby
          bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
          bc(idfcnb+i)=-trial*bc(ifcnb+i-1)-argb*bc(idfcnb+i-1)         6d30s16
     $         -ff*bc(idfcnb+i-2)                                       6d30s16
          bc(idfcnb2+i)=-trial*bc(idfcnb+i-1)*2d0-argb*bc(idfcnb2+i-1)  6d30s16
     $         -ff*bc(idfcnb2+i-2)                                      6d30s16
          ff=ff+bsib
         end do
         ifcnbu=idfcnb2                                                 6d30s16
        else
        do i=2,lb+iby
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
         ff=ff+bsib
        end do
        ifcnbu=ifcnb
        end if                                                          6d30s16
        bc(ifcnk+1)=-argk*term
        ff=bsik
        if(ipk.eq.1.and.ifcn.eq.2)then                                  6d30s16
         bc(idfcnk)=0d0                                                 6d30s16
         bc(idfcnk+1)=-term*trial                                       6d30s16
         do i=2,lk+iky
          bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
          bc(idfcnk+i)=-trial*bc(ifcnk+i-1)-argk*bc(idfcnk+i-1)         6d30s16
     $         -ff*bc(idfcnk+i-2)                                       6d30s16
          ff=ff+bsik
         end do
         ifcnku=idfcnk                                                  6d30s16
        else if(ipk.eq.2.and.ifcn.eq.2)then                             6d30s16
         bc(idfcnk)=0d0                                                 6d30s16
         bc(idfcnk+1)=-term*trial                                       6d30s16
         bc(idfcnk2)=0d0                                                6d30s16
         bc(idfcnk2+1)=0d0                                              6d30s16
         do i=2,lk+iky
          bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
          bc(idfcnk+i)=-trial*bc(ifcnk+i-1)-argk*bc(idfcnk+i-1)         6d30s16
     $         -ff*bc(idfcnk+i-2)                                       6d30s16
          bc(idfcnk2+i)=-2d0*trial*bc(idfcnk+i-1)-argk*bc(idfcnk2+i-1)  6d30s16
     $         -ff*bc(idfcnk2+i-2)                                      6d30s16
          ff=ff+bsik
         end do
         ifcnku=idfcnk2                                                 6d30s16
        else                                                            6d30s16
        do i=2,lk+iky
         bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
         ff=ff+bsik
        end do
        ifcnku=ifcnk
        end if                                                          6d30s16
        jcarty=icarty
        do ik=0,lk
         do ib=0,lb
          bc(jcarty)=bc(jcarty)+bc(ifcnbu+ib+iby)*bc(ifcnku+ik+iky)     6d30s16
          jcarty=jcarty+1
         end do
        end do
        if(ldeb2)then
         write(6,*)('carty ')
         call prntm2(bc(icarty),lb+1,lk+1,lb+1)
        end if
        argz=arg+sz                                                     1d26s10
        argb=argz-zb                                                    1d26s10
        argk=argz-zk                                                    1d26s10
        term=bc(jwgt+iq)                                                1d14s10
        bc(ifcnb)=1d0
        bc(ifcnk)=term
        bc(ifcnb+1)=-argb
        ff=bsib
        if(ipb.eq.1.and.ifcn.eq.3)then                                  6d30s16
         bc(idfcnb)=0d0                                                 6d30s16
         bc(idfcnb+1)=-trial                                            6d30s16
         do i=2,lb+ibz
          bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
          bc(idfcnb+i)=-trial*bc(ifcnb+i-1)-argb*bc(idfcnb+i-1)         6d30s16
     $         -ff*bc(idfcnb+i-2)                                       6d30s16
          ff=ff+bsib
         end do
         ifcnbu=idfcnb                                                  6d30s16
        else if(ipb.eq.2.and.ifcn.eq.3)then                             6d30s16
         bc(idfcnb)=0d0                                                 6d30s16
         bc(idfcnb+1)=-trial                                            6d30s16
         bc(idfcnb2)=0d0                                                6d30s16
         bc(idfcnb2+1)=0d0                                              6d30s16
         do i=2,lb+ibz
          bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
          bc(idfcnb+i)=-trial*bc(ifcnb+i-1)-argb*bc(idfcnb+i-1)         6d30s16
     $         -ff*bc(idfcnb+i-2)                                       6d30s16
          bc(idfcnb2+i)=-trial*bc(idfcnb+i-1)*2d0-argb*bc(idfcnb2+i-1)  6d30s16
     $         -ff*bc(idfcnb2+i-2)                                      6d30s16
          ff=ff+bsib
         end do
         ifcnbu=idfcnb2                                                 6d30s16
        else
        do i=2,lb+ibz
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
         ff=ff+bsib
        end do
        ifcnbu=ifcnb
        end if                                                          6d30s16
        bc(ifcnk+1)=-argk*term
        ff=bsik
        if(ipk.eq.1.and.ifcn.eq.3)then                                  6d30s16
         bc(idfcnk)=0d0                                                 6d30s16
         bc(idfcnk+1)=-term*trial                                       6d30s16
         do i=2,lk+ikz
          bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
          bc(idfcnk+i)=-trial*bc(ifcnk+i-1)-argk*bc(idfcnk+i-1)         6d30s16
     $         -ff*bc(idfcnk+i-2)                                       6d30s16
          ff=ff+bsik
         end do
         ifcnku=idfcnk                                                  6d30s16
        else if(ipk.eq.2.and.ifcn.eq.3)then                             6d30s16
         bc(idfcnk)=0d0                                                 6d30s16
         bc(idfcnk+1)=-term*trial                                       6d30s16
         bc(idfcnk2)=0d0                                                6d30s16
         bc(idfcnk2+1)=0d0                                              6d30s16
         do i=2,lk+ikz
          bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
          bc(idfcnk+i)=-trial*bc(ifcnk+i-1)-argk*bc(idfcnk+i-1)         6d30s16
     $         -ff*bc(idfcnk+i-2)                                       6d30s16
          bc(idfcnk2+i)=-2d0*trial*bc(idfcnk+i-1)-argk*bc(idfcnk2+i-1)  6d30s16
     $         -ff*bc(idfcnk2+i-2)                                      6d30s16
          ff=ff+bsik
         end do
         ifcnku=idfcnk2                                                 6d30s16
        else                                                            6d30s16
        do i=2,lk+ikz
         bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
         ff=ff+bsik
        end do
        ifcnku=ifcnk
        end if
        jcartz=icartz
        do ik=0,lk
         do ib=0,lb
          if(ldeb2)write(6,*)bc(jcartz),bc(ifcnb+ib+ibz),
     $         bc(ifcnk+ik+ikz)
          bc(jcartz)=bc(jcartz)+bc(ifcnbu+ib+ibz)*bc(ifcnku+ik+ikz)
          jcartz=jcartz+1
         end do
        end do
       end do
        if(ldeb2)then
         write(6,*)('cartz '),ikz
         call prntm2(bc(icartz),lb+1,lk+1,lb+1)
        end if
       jpowk=ipp(lkp)+1                                                 1d26s10
       rwgtu=rwgt
       if(ipx.eq.1)rwgtu=rwgtu*t2*dbigx                                 6d30s16
       if(ipx.eq.2)rwgtu=rwgtu*t2*(dbigx+t2*dbigx2)                     6d30s16
       if(ldeb2)write(6,*)('rwgtu: '),rwgtu
       do ik=1,npk                                                      1d26s10
        j1=ibc(jpowk)*lbp+icartx                                        1d26s10
        j2=ibc(jpowk+1)*lbp+icarty                                      1d26s10
        j3=ibc(jpowk+2)*lbp+icartz                                      1d26s10
        jpowk=jpowk+3
        jpowb=ipp(lbp)+1                                                1d14s10
        do ib=1,npb                                                     1d26s10
         i1=ibc(jpowb)+j1                                               1d26s10
         i2=ibc(jpowb+1)+j2                                             1d26s10
         i3=ibc(jpowb+2)+j3                                             1d26s10
         jpowb=jpowb+3
         iad4=icarti+ib-1+idcarti*(ik-1)                                1d14s10
         if(ldeb2)write(6,*)bc(iad4),bc(i1),bc(i2),bc(i3),
     $        bc(i1)*bc(i2)*bc(i3),rwgtu,rwgtu/rwgt
         bc(iad4)=bc(iad4)+bc(i1)*bc(i2)*bc(i3)*rwgtu                    1d26s10
        end do
       end do
      end do                                                            1d14s10
      if(ldeb)write(6,*)('before sphericalization '),bc(icarti)
      jcarti=icarti                                                     1d22s10
      jcarti2=icarti2                                                   1d22s10
      do isym=1,8                                                       1d22s10
       if(ipt(isym,lkp).gt.0)then
        istok=ipt(isym,lkp)
        nherek=ibc(istok)
        mherek=ibc(istok+1)
        icoefk=istok+2+mherek+3*nherek
        call dgemm('n','n',npb,mherek,nherek,1d0,bc(jcarti),npb,
     $       bc(icoefk),nherek,0d0,bc(jcarti2),npb,                     1d26s10
     d' nattracd.  1')
        jcarti=jcarti+npb*nherek
        jcarti2=jcarti2+npb*mherek
       end if
      end do                                                            1d22s10
      lkt=2*lk+1                                                        1d22s10
      do i=1,npb                                                        1d22s10
       do j=1,lkt                                                       1d22s10
        iad1=j-1+lkt*(i-1)                                              1d22s10
        iad2=i-1+npb*(j-1)                                              1d22s10
        bc(icarti+iad1)=bc(icarti2+iad2)                                1d22s10
       end do                                                           1d22s10
      end do                                                            1d22s10
      jcarti=icarti                                                     1d22s10
      jcarti2=icarti2                                                   1d22s10
      do isym=1,8                                                       1d22s10
       if(ipt(isym,lbp).gt.0)then
        istob=ipt(isym,lbp)
        nhereb=ibc(istob)
        mhereb=ibc(istob+1)
        icoefb=istob+2+mhereb+3*nhereb
        call dgemm('n','n',lkt,mhereb,nhereb,1d0,bc(jcarti),lkt,
     $       bc(icoefb),nhereb,0d0,bc(jcarti2),lkt,                     1d26s10
     d' nattracd.  2')
        jcarti=jcarti+lkt*nhereb
        jcarti2=jcarti2+lkt*mhereb
       end if
      end do                                                            1d22s10
      lbt=2*lb+1
      if(ldeb)then
       write(6,*)('front, output'),front
       call prntm2(bc(icarti2),lkt,lbt,lkt)
      end if
      ibcoff=inode
      return
      end
