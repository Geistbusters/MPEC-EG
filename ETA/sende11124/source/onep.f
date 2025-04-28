c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine onep(lb,zbra,xnbra,xb,yb,zb,                           12d15s19
     $                lk,zket,xnket,xk,yk,zk,                           12d15s19
     $             nsymb,icarti2,idxt,idyt,idzt,ipx,ipy,ipz,idx,idy,idz,11d9s22
     $     bc,ibc)                                                      11d9s22
      implicit real*8 (a-h,o-z)
      include "common.store"
      include "common.spher"
      include "common.rys"                                              1d27s23
      logical ldebug                                                    12d19s19
      dimension store1(1),store2(1)                                     2d22s10
c
c     compute matrix element of
c     (d/dx**idxt)*(d/dy**idyt)*(d/dz**idzt) on bra and
c     x**ipx y**ipy z**ipz d/dx**idx d/dy**idy d/dz**idz on ket.
c
      ldebug=.false.                                                    12d19s19
      lbp=lb+1
      lkp=lk+1
      lbpp=lbp+1
      lkpp=lkp+1
      lbppp=lbpp+1                                                      1d22s10
      lkppp=lkpp+1                                                      1d22s10
      if(ldebug)
     $     write(6,1)lb,zbra,xnbra,xb,yb,zb,lk,zket,xnket,xk,yk,zk
    1 format('in onep for: ',i3,2es15.7,3f12.4,i3,2es15.7,3f12.4)       1d22s10
      if(ldebug)
     $     write(6,*)('args: '),idxt,idyt,idzt,ipx,ipy,ipz,idx,idy,
     $     idz
      nq=(lb+lk+2+idxt+idyt+idzt+ipx+ipy+ipz+idx+idy+idz)/2             12d14s19
      nq=nq+1
      if(ldebug)write(6,*)('nq = '),nq
      inode=ibcoff
      iwgt=inode+nq
      itmp=iwgt+nq
      ibcoff=itmp+nq
      call enough('onep.  1',bc,ibc)
      ibcoff=itmp
      iad1=(nq*(nq-1))+ibcgh                                            1d27s23
      iad2=iad1+nq                                                      1d27s23
      jnode=iad2-1                                                      1d27s23
      jwgt=iad1-1                                                       1d27s23
      icartx=ibcoff
      nwds=(lbp+1)*(lkp+1)
      ibcoff=icartx+nwds*3                                              1d22s10
      ifcnb=ibcoff
      ifcnk=ifcnb+lbpp+max(idxt,idyt,idzt)                              12d19s19
      ibcoff=ifcnk+lkpp+max(idx,idy,idz)                                12d19s19
      call enough('onep.  2',bc,ibc)
      zbki=1d0/(zbra+zket)                                              1d22s10
      sftx=(zbra*xb+zket*xk)*zbki                                       1d22s10
      sfty=(zbra*yb+zket*yk)*zbki                                       1d22s10
      sftz=(zbra*zb+zket*zk)*zbki                                       1d22s10
      pref=zbra*zket*zbki*((xk-xb)**2+(yk-yb)**2+(zk-zb)**2)            1d22s10
      pref=exp(-pref)                                                   1d22s10
      bsik=0.5d0/zket
      bsib=0.5d0/zbra
      asi=1d0/sqrt(zbra+zket)
      zket2=zket*2d0                                                    12d14s19
      dnorm=xnbra*xnket*pref                                            1d25s10
      idxyz=idx+idy+idz                                                 12d14s19
      if(idxyz.gt.0)dnorm=dnorm*(zket2**idxyz)                          12d14s19
      zbra2=zbra*2d0                                                    12d14s19
      idtxyz=idxt+idyt+idzt                                             12d14s19
      if(idtxyz.gt.0)dnorm=dnorm*(zbra2**idtxyz)                        12d14s19
      if(ldebug)write(6,*)('dnorm: '),dnorm
      jcartx=icartx-1
      icarty=icartx+nwds                                                1d22s10
      icartz=icarty+nwds                                                1d22s10
      do i=1,nwds*3                                                     1d22s10
       bc(jcartx+i)=0d0
      end do
      do iq=1,nq
       arg=bc(jnode+iq)*asi
       argx=arg+sftx
       argb=argx-xb
       argk=argx-xk
       term=bc(jwgt+iq)*asi
       bc(ifcnb)=1d0
       if(ipx.ne.0)then
        powr=argx**ipx                                                   12d14s19
        bc(ifcnk)=term*powr                                              12d14s19
       else                                                             12d19s19
        bc(ifcnk)=term                                                  12d19s19
       end if                                                           12d19s19
       bc(ifcnb+1)=-argb
       ff=bsib
       do i=2,lbp+idxt                                                   12d14s19
        bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)                2d16s10
        ff=ff+bsib
       end do
       bc(ifcnk+1)=-argk*bc(ifcnk)                                      12d16s19
       ff=bsik
       do i=2,lkp+idx                                                   12d14s19
        bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)                2d16s10
        ff=ff+bsik
       end do
       jcartx=icartx
       do ik=0,lkp
        do ib=0,lbp
         bc(jcartx)=bc(jcartx)+bc(ifcnb+ib+idxt)*bc(ifcnk+ik+idx)       12d14s19
         jcartx=jcartx+1
        end do
       end do
       argy=arg+sfty
       argb=argy-yb
       argk=argy-yk
       term=bc(jwgt+iq)*asi
       bc(ifcnb)=1d0
       bc(ifcnb+1)=-argb
       if(ipy.gt.0)then                                                 12d19s19
        powr=argy**ipy                                                   12d14s19
        bc(ifcnk)=term*powr                                              12d14s19
       else                                                             12d19s19
        bc(ifcnk)=term                                                  12d19s19
       end if                                                           12d19s19
       ff=bsib
       do i=2,lbp+idyt                                                  12d14s19
        bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)                2d16s10
        ff=ff+bsib
       end do
       bc(ifcnk+1)=-argk*bc(ifcnk)                                      12d16s19
       ff=bsik
       do i=2,lkp+idy                                                   12d14s19
        bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)                2d16s10
        ff=ff+bsik
       end do
       jcarty=icarty
       do ik=0,lkp
        do ib=0,lbp
         bc(jcarty)=bc(jcarty)+bc(ifcnb+ib+idyt)*bc(ifcnk+ik+idy)       12d15s19
         jcarty=jcarty+1
        end do
       end do
       argz=arg+sftz
       argb=argz-zb
       argk=argz-zk
       term=bc(jwgt+iq)*asi
       bc(ifcnb)=1d0
       bc(ifcnb+1)=-argb
       if(ipz.gt.0)then                                                 12d19s19
        powr=argz**ipz                                                   12d16s19
        bc(ifcnk)=term*powr                                              12d14s19
       else                                                             12d19s19
        bc(ifcnk)=term                                                  12d19s19
       end if                                                           12d19s19
       ff=bsib
       do i=2,lbp+idzt                                                  12d15s19
        bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)                2d16s10
        ff=ff+bsib
       end do
       bc(ifcnk+1)=-argk*bc(ifcnk)                                      12d16s19
       ff=bsik
       do i=2,lkp+idz                                                   12d15s19
        bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)                2d16s10
        ff=ff+bsik
       end do
       jcartz=icartz
       do ik=0,lkp
        do ib=0,lbp
         bc(jcartz)=bc(jcartz)+bc(ifcnb+ib+idzt)*bc(ifcnk+ik+idz)       12d15s19
         jcartz=jcartz+1
        end do
       end do
      end do
      ibcoff=ifcnb
      npb=ibc(ipp(lbp))
      npk=ibc(ipp(lkp))
      icarti=ibcoff                                                     1d22s10
      ibcoff=icarti+npb*npk                                             1d22s10
      call enough('onep.  3',bc,ibc)
      jpowk=ipp(lkp)+1                                                  1d22s10
      do ik=1,npk                                                       1d22s10
       jx=ibc(jpowk)*lbpp+icartx
       jy=ibc(jpowk+1)*lbpp+icarty
       jz=ibc(jpowk+2)*lbpp+icartz
       jpowk=jpowk+3
       jpowb=ipp(lbp)+1
       do ib=1,npb                                                      1d22s10
        ix=ibc(jpowb)+jx                                                1d25s10
        iy=ibc(jpowb+1)+jy                                              1d25s10
        iz=ibc(jpowb+2)+jz                                              1d25s10
        jpowb=jpowb+3
        iad=ib-1+npb*(ik-1)                                             1d22s10
        bc(iad+icarti)=bc(ix)*bc(iy)*bc(iz)*dnorm                       1d25s10
       end do                                                           1d22s10
      end do                                                            1d22s10
      icarti2=ibcoff                                                    1d22s10
      ibcoff=icarti2+npb*npk                                            12d15s19
      call enough('onep.  4',bc,ibc)
      if(iusecart.eq.0)then                                             2d21s20
       jcarti=icarti                                                     1d22s10
       jcarti2=icarti2                                                   1d22s10
       do isym=1,8                                                       1d22s10
        if(ipt(isym,lkp).gt.0)then
         istok=ipt(isym,lkp)
         nherek=ibc(istok)
         nherek2=nherek*2
         mherek=ibc(istok+1)
         icoefk=istok+2+mherek+3*nherek
         call dgemm('n','n',npb,mherek,nherek,1d0,bc(jcarti),npb,
     $       bc(icoefk),nherek,0d0,bc(jcarti2),npb,                     1d26s10
     d' onep.  1')
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
         nhereb2=nhereb*2
         mhereb=ibc(istob+1)
         icoefb=istob+2+mhereb+3*nhereb
         call dgemm('n','n',lkt,mhereb,nhereb,1d0,bc(jcarti),lkt,
     $       bc(icoefb),nhereb,0d0,bc(jcarti2),lkt,                     1d26s10
     d' onep.  2')
         jcarti=jcarti+lkt*nhereb
         jcarti2=jcarti2+lkt*mhereb
        end if
       end do                                                            1d22s10
      else                                                              2d21s20
       icarti2=icarti                                                   2d21s20
      end if                                                            2d21s20
      ibcoff=inode
      return
      end
