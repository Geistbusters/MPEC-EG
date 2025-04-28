c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine onei(lb,zbra,xnbra,xb,yb,zb,iaddb,                     1d22s10
     $                lk,zket,xnket,xk,yk,zk,iaddk,store1,store2,       5d3s10
     $                nsymb,icarti2,icartte2,bc,ibc)                    11d9s22
      implicit real*8 (a-h,o-z)
      include "common.store"
      include "common.spher"
      include "common.rys"                                              1d27s23
      dimension store1(1),store2(1)                                     2d22s10
c
c     compute overlap and kinetic energy
c
      lbp=lb+1
      lkp=lk+1
      lbpp=lbp+1
      lkpp=lkp+1
      lbppp=lbpp+1                                                      1d22s10
      lkppp=lkpp+1                                                      1d22s10
    1 format('in onei for: ',i3,1pe15.7,0p3f12.4,i3,1pe15.7,0p3f12.4)   1d22s10
      nq=(lb+lk+2)/2
      nq=nq+1
      inode=ibcoff
      iwgt=inode+nq
      itmp=iwgt+nq
      ibcoff=itmp+nq
      call enough('onei.  1',bc,ibc)
      ibcoff=itmp
      iad1=(nq*(nq-1))+ibcgh                                            1d27s23
      iad2=iad1+nq                                                      1d27s23
      jnode=iad2-1                                                      1d27s23
      jwgt=iad1-1                                                       1d27s23
      icartx=ibcoff
      nwds=(lbp+1)*(lkp+1)
      ibcoff=icartx+nwds*3                                              1d22s10
      ifcnb=ibcoff
      ifcnk=ifcnb+lbpp
      ibcoff=ifcnk+lkpp
      call enough('onei.  2',bc,ibc)
      zbki=1d0/(zbra+zket)                                              1d22s10
      sftx=(zbra*xb+zket*xk)*zbki                                       1d22s10
      sfty=(zbra*yb+zket*yk)*zbki                                       1d22s10
      sftz=(zbra*zb+zket*zk)*zbki                                       1d22s10
      pref=zbra*zket*zbki*((xk-xb)**2+(yk-yb)**2+(zk-zb)**2)            1d22s10
      pref=exp(-pref)                                                   1d22s10
      bsik=0.5d0/zket
      bsib=0.5d0/zbra
      asi=1d0/sqrt(zbra+zket)
      dnorm=xnbra*xnket*pref                                            1d25s10
      dnorm2=dnorm*2d0*zket*zbra
      fderb=1d0/zbra                                                    2s12s10
      fderk=1d0/zket                                                    2s12s10
      fderbk=0.5d0*fderb*fderk                                          2s12s10
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
       bc(ifcnk)=term
       bc(ifcnb+1)=-argb
       ff=bsib
       do i=2,lbp
        bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)                2d16s10
        ff=ff+bsib
       end do
       bc(ifcnk+1)=-argk*term
       ff=bsik
       do i=2,lkp
        bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)                2d16s10
        ff=ff+bsik
       end do
       jcartx=icartx
       do ik=0,lkp
        do ib=0,lbp
         bc(jcartx)=bc(jcartx)+bc(ifcnb+ib)*bc(ifcnk+ik)
         jcartx=jcartx+1
        end do
       end do
       argy=arg+sfty
       argb=argy-yb
       argk=argy-yk
       term=bc(jwgt+iq)*asi
       bc(ifcnb)=1d0
       bc(ifcnk)=term
       bc(ifcnb+1)=-argb
       ff=bsib
       do i=2,lbp
        bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)                2d16s10
        ff=ff+bsib
       end do
       bc(ifcnk+1)=-argk*term
       ff=bsik
       do i=2,lkp
        bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)                2d16s10
        ff=ff+bsik
       end do
       jcarty=icarty
       do ik=0,lkp
        do ib=0,lbp
         bc(jcarty)=bc(jcarty)+bc(ifcnb+ib)*bc(ifcnk+ik)
         jcarty=jcarty+1
        end do
       end do
       argz=arg+sftz
       argb=argz-zb
       argk=argz-zk
       term=bc(jwgt+iq)*asi
       bc(ifcnb)=1d0
       bc(ifcnk)=term
       bc(ifcnb+1)=-argb
       ff=bsib
       do i=2,lbp
        bc(ifcnb+i)=-argb*bc(ifcnb+i-1)-ff*bc(ifcnb+i-2)                2d16s10
        ff=ff+bsib
       end do
       bc(ifcnk+1)=-argk*term
       ff=bsik
       do i=2,lkp
        bc(ifcnk+i)=-argk*bc(ifcnk+i-1)-ff*bc(ifcnk+i-2)                2d16s10
        ff=ff+bsik
       end do
       jcartz=icartz
       do ik=0,lkp
        do ib=0,lbp
         bc(jcartz)=bc(jcartz)+bc(ifcnb+ib)*bc(ifcnk+ik)
         jcartz=jcartz+1
        end do
       end do
      end do
      ibcoff=ifcnb
      npb=ibc(ipp(lbp))
      npk=ibc(ipp(lkp))
      icarti=ibcoff                                                     1d22s10
      ibcoff=icarti+npb*npk                                             1d22s10
      icartte=ibcoff                                                    1d22s10
      ibcoff=icartte+npb*npk                                            1d22s10
      call enough('onei.  3',bc,ibc)
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
        ixp=ix+lbppp
        iyp=iy+lbppp
        izp=iz+lbppp
        bc(iad+icartte)=(bc(ixp)*bc(iy)*bc(iz)                           1d22s10
     $      +bc(ix)*bc(iyp)*bc(iz)+bc(ix)*bc(iy)*bc(izp))*dnorm2        1d25s10
       end do                                                           1d22s10
      end do                                                            1d22s10
      icarti2=ibcoff                                                    1d22s10
      icartte2=icarti2+npb*npk                                          1d22s10
      ibcoff=icartte2+npb*npk                                           1d22s10
      call enough('onei.  4',bc,ibc)
      jcarti=icarti                                                     1d22s10
      jcartte=icartte                                                   1d22s10
      jcarti2=icarti2                                                   1d22s10
      jcartte2=icartte2                                                 1d22s10
      do isym=1,8                                                       1d22s10
       if(ipt(isym,lkp).gt.0)then
        istok=ipt(isym,lkp)
        nherek=ibc(istok)
        nherek2=nherek*2
        mherek=ibc(istok+1)
        icoefk=istok+2+mherek+3*nherek
        call dgemm('n','n',npb,mherek,nherek,1d0,bc(jcarti),npb,
     $       bc(icoefk),nherek,0d0,bc(jcarti2),npb,                     1d26s10
     d' onei.  1')
        call dgemm('n','n',npb,mherek,nherek,1d0,bc(jcartte),npb,
     $       bc(icoefk),nherek,0d0,bc(jcartte2),npb,                    1d26s10
     d' onei.  2')
        jcarti=jcarti+npb*nherek
        jcartte=jcartte+npb*nherek
        jcarti2=jcarti2+npb*mherek
        jcartte2=jcartte2+npb*mherek
       end if
      end do                                                            1d22s10
      lkt=2*lk+1                                                        1d22s10
      do i=1,npb                                                        1d22s10
       do j=1,lkt                                                       1d22s10
        iad1=j-1+lkt*(i-1)                                              1d22s10
        iad2=i-1+npb*(j-1)                                              1d22s10
        bc(icarti+iad1)=bc(icarti2+iad2)                                1d22s10
        bc(icartte+iad1)=bc(icartte2+iad2)                              1d22s10
       end do                                                           1d22s10
      end do                                                            1d22s10
      jcarti=icarti                                                     1d22s10
      jcartte=icartte                                                   1d22s10
      jcarti2=icarti2                                                   1d22s10
      jcartte2=icartte2                                                 1d22s10
      do isym=1,8                                                       1d22s10
       if(ipt(isym,lbp).gt.0)then
        istob=ipt(isym,lbp)
        nhereb=ibc(istob)
        nhereb2=nhereb*2
        mhereb=ibc(istob+1)
        icoefb=istob+2+mhereb+3*nhereb
        call dgemm('n','n',lkt,mhereb,nhereb,1d0,bc(jcarti),lkt,
     $       bc(icoefb),nhereb,0d0,bc(jcarti2),lkt,                     1d26s10
     d' onei.  3')
        call dgemm('n','n',lkt,mhereb,nhereb,1d0,bc(jcartte),lkt,
     $       bc(icoefb),nhereb,0d0,bc(jcartte2),lkt,                    1d26s10
     d' onei.  4')
        jcarti=jcarti+lkt*nhereb
        jcartte=jcartte+lkt*nhereb
        jcarti2=jcarti2+lkt*mhereb
        jcartte2=jcartte2+lkt*mhereb
       end if
      end do                                                            1d22s10
      ibcoff=inode
      return
      end
