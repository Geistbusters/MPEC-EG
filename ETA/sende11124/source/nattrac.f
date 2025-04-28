c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine nattrac(lb,zbra,xnbra,xb,yb,zb,iaddb,                  1d26s10
     $                   lk,zket,xnket,xk,yk,zk,iaddk,                  1d26s10
     $                   atnum,xi,yi,zi,buff,nsymb,icarti2,bc,ibc)      11d9s22
      implicit real*8 (a-h,o-z)
      include "common.store"
      include "common.spher"
      include "common.rys"                                              6d22s12
      dimension buff(1)                                                 1d26s10
      data icall/0/
      save
      lbp=lb+1
      lkp=lk+1
      zz=-atnum                                                         1d19s10
      ihead=0                                                           1d21s10
    1 format('in nattrac: ',f5.1,3f8.4,i3,1pe15.7,0p3f8.4,i3,1pe15.7,   1d26s10
     $         0p3f8.4)                                                 1d26s10
      nq=(lb+lk+2)/2
      inode=ibcoff
      iwgt=inode+nq
      itmp=iwgt+nq
      ibcoff=itmp+nq
      call enough('nattrac.  1',bc,ibc)
      iad1=(nq*(nq-1))+ibcgh                                            1d19s23
      iad2=iad1+nq
      ibcoff=itmp
      jnode=iad2-1                                                      1d19s23
      jwgt=iad1-1                                                       1d19s23
      irnode=ibcoff                                                     1d7s09
      nqr=nq                                                            1d12s10
      irwgt=irnode+nqr                                                   1d7s09
      ibcoff=irwgt+nqr                                                   1d12s10
      call enough('nattrac.  2',bc,ibc)
      zbki=1d0/(zbra+zket)                                              1d22s10
      sftx=(zbra*xb+zket*xk)*zbki                                       1d22s10
      sfty=(zbra*yb+zket*yk)*zbki                                       1d22s10
      sftz=(zbra*zb+zket*zk)*zbki                                       1d22s10
      pref=zbra*zket*zbki*((xk-xb)**2+(yk-yb)**2+(zk-zb)**2)            1d22s10
      pref=exp(-pref)                                                   1d22s10
      bigx=(zbra+zket)*((sftx-xi)**2+(sfty-yi)**2+(sftz-zi)**2)         1d26s10
      if(nqr.ge.nqx)then                                                6d22s12
       write(6,*)('asking for more rys quadrature points than are '),   6d22s12
     $     ('available '),nqr,nqx                                       6d22s12
       call dws_finalize                                                6d22s12
       stop                                                             6d22s12
      end if                                                            6d22s12
      iaddq=ibc(ibcrys+2*(nqr-1))                                       6d22s12
      nreg=ibc(ibcrys+1+2*(nqr-1))                                      6d22s12
      reg=min(dfloat(nreg*2+1),bigx/bc(iaddq))                          1d19s23
      ireg=reg                                                          4d9s18
      ireg=ireg+1                                                       4d9s18
      if(ireg.le.nreg.and.reg.lt.dfloat(nreg*2))then                    4d9s18
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
      ifcnk=ifcnb+max(3,lbp)                                            1d14s10
      ibcoff=ifcnk+max(3,lkp)                                           1d14s10
      icarti=ibcoff                                                     1d14s10
      npb=ibc(ipp(lbp))                                                 1d26s10
      npk=ibc(ipp(lkp))                                                 1d26s10
      nwdi=npb*npk                                                      1d26s10
      ibcoff=icarti+nwdi                                                1d14s10
      call enough('nattrac.  3',bc,ibc)
      bsik=0.5d0/zket
      bsib=0.5d0/zbra
      rho=sqrt(zket+zbra)                                               1d14s10
      dnorm=xnbra*xnket
      front=dnorm*front0/(rho*rho)                                      1d14s10
      do i=0,nwdi-1                                                     1d14s10
       bc(icarti+i)=0d0                                                 1d14s10
      end do                                                            1d14s10
      idcarti=npb                                                       1d26s10
      do irq=1,nqr                                                      1d14s10
       use2=rho*rho*bc(jrnode+irq)                                      7d5s12
       asi2=1d0/(use2+zbra+zket)                                        1d26s10
       asi=sqrt(asi2)                                                   1d26s10
       rwgt=bc(jrwgt+irq)*front                                         1d14s10
       sx=(sftx+use2*xi)*asi2                                           1d26s10
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
        do i=2,lb
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
         ff=ff+bsib
        end do
        bc(ifcnk+1)=-argk*term
        ff=bsik
        do i=2,lk
         bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
         ff=ff+bsik
        end do
        jcartx=icartx
        do ik=0,lk
         do ib=0,lb
          bc(jcartx)=bc(jcartx)+bc(ifcnb+ib)*bc(ifcnk+ik)
          jcartx=jcartx+1
         end do
        end do
        argy=arg+sy                                                     1d26s10
        argb=argy-yb                                                    1d26s10
        argk=argy-yk                                                    1d26s10
        term=bc(jwgt+iq)                                                1d14s10
        bc(ifcnb)=1d0
        bc(ifcnk)=term
        bc(ifcnb+1)=-argb
        ff=bsib
        do i=2,lb
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
         ff=ff+bsib
        end do
        bc(ifcnk+1)=-argk*term
        ff=bsik
        do i=2,lk
         bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
         ff=ff+bsik
        end do
        jcarty=icarty
        do ik=0,lk
         do ib=0,lb
          bc(jcarty)=bc(jcarty)+bc(ifcnb+ib)*bc(ifcnk+ik)
          jcarty=jcarty+1
         end do
        end do
        argz=arg+sz                                                     1d26s10
        argb=argz-zb                                                    1d26s10
        argk=argz-zk                                                    1d26s10
        term=bc(jwgt+iq)                                                1d14s10
        bc(ifcnb)=1d0
        bc(ifcnk)=term
        bc(ifcnb+1)=-argb
        ff=bsib
        do i=2,lb
         bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)               2d16s10
         ff=ff+bsib
        end do
        bc(ifcnk+1)=-argk*term
        ff=bsik
        do i=2,lk
         bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)               2d16s10
         ff=ff+bsik
        end do
        jcartz=icartz
        do ik=0,lk
         do ib=0,lb
          bc(jcartz)=bc(jcartz)+bc(ifcnb+ib)*bc(ifcnk+ik)
          jcartz=jcartz+1
         end do
        end do
       end do
       jpowk=ipp(lkp)+1                                                 1d26s10
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
         bc(iad4)=bc(iad4)+bc(i1)*bc(i2)*bc(i3)*rwgt                    1d26s10
        end do
       end do
      end do                                                            1d14s10
      icarti2=ibcoff                                                    1d22s10
      ibcoff=icarti2+npb*npk                                            1d26s10
      call enough('nattrac.  4',bc,ibc)
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
     d' nattrac.  1')
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
     d' nattrac.  2')
        jcarti=jcarti+lkt*nhereb
        jcarti2=jcarti2+lkt*mhereb
       end if
      end do                                                            1d22s10
      lbt=2*lb+1
      ibcoff=inode
      return
      end
