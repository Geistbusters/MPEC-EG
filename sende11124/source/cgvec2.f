c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cgvec2(hdig,ih0e,i4oa,irhs,ie,iv,idva,idv,idvt,nroot,
     $     nadet,nbdet,nsbeta,nherec,ilc,ihc,nherect,ilct,ihct,ncona,
     $     iaorb,iborb,nalpha,nbeta,m12sym,norb,nsymb,idata1,idata2,
     $     idatb1,idatb2,m1c,idatac,idatbc,mst,mnz,nz,tul,icode,lprt,bc,11d9s22
     $     ibc,tolb)                                                    4d27s23
      implicit real*8 (a-h,o-z)
c
c     solve (H-e+e*v*vt)dv=rhs for dv
c
      logical ldebug,lprt,lnew                                               7d28s22
      dimension hdig(*),irhs(*),iv(*),idv(*),nsbeta(*),nadet(*),        5d27s22
     $     nbdet(*),idva(*),idvt(*),iga(8),ig(8),igt(8),idp(8),         5d27s22
     $     nherec(*),nherect(*),ilc(*),ilct(*)                          5d27s22
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data icall/0/
      data loop,loopx/0,1000/
      save icall
      common/singcm/iuse,nff
      icall=icall+1
      tol=tul                                                           4d9s23
      tolb=tul                                                          4d27s23
      lnew=.true.                                                       4d13s23
      do isb=1,nsymb
       jsb=nsbeta(isb)                                                  5d27s22
       if(min(nadet(isb),nbdet(jsb)).gt.0)then                          5d27s22
        ncol=nadet(isb)*nbdet(jsb)                                       5d27s22
        if(ncona.eq.1)then                                              5d27s22
         bc(idva(isb))=0d0                                              5d27s22
        end if                                                          5d27s22
       end if                                                           5d27s22
      end do
      if(ncona.eq.1)then                                                  5d20s22
       return                                                           5d20s22
      end if                                                            5d20s22
      nconam=ncona-1
      nwds=nroot*ncona                                                  5d27s22
      ix=ibcoff                                                         4d4s23
      ir=ix+nwds                                                        4d4s23
      iz=ir+nwds*2                                                        4d4s23
      ip=iz+nwds                                                        4d4s23
      igdiag=ip+nwds                                                    4d4s23
      iq=igdiag+nwds                                                    4d4s23
      iqq=iq+nwds                                                       4d4s23
      ibeta=iqq+nwds                                                    4d4s23
      iscale=ibeta+nroot                                                  5d27s22
      idotrz=iscale+nroot                                                 5d27s22
      idotrznew=idotrz+nroot                                               5d27s22
      isz=idotrznew+nroot                                                  5d27s22
      ialpha=isz+nroot
      ibcoff=ialpha+nroot                                               4d4s23
      iqqq=ibcoff                                                       4d7s23
      ippp=iqqq+nwds                                                    4d7s23
      ibcoff=ippp+nwds                                                  4d7s23
      do isb=1,nsymb                                                    5d27s22
       jsb=nsbeta(isb)                                                  5d27s22
       ig(isb)=ibcoff                                                   5d27s22
       igt(isb)=ig(isb)+nroot*nbdet(jsb)*nherec(isb)                    5d27s22
       ibcoff=igt(isb)+nroot*nadet(jsb)*nherect(isb)                    5d27s22
      end do                                                            5d27s22
      call enough('cgvec2.  1',bc,ibc)
      nrootm=nroot-1                                                    5d27s22
      nconarm=nwds-1                                                    4d4s23
      do i=0,nconarm                                                     5d27s22
       bc(igdiag+i)=0d0                                                 5d27s22
      end do                                                            5d27s22
      jgdiag=igdiag                                                     5d27s22
      ioff=1                                                            5d27s22
      do isb=1,nsymb                                                    5d27s22
       jsb=nsbeta(isb)                                                  5d27s22
       do ia=0,nherec(isb)-1                                            7d25s22
        iap=ia+ilc(isb)-1                                               7d25s22
        do ib=0,nbdet(jsb)-1                                             5d27s22
         iad1=jgdiag+nroot*(iap+nadet(isb)*ib)                           7d25s22
         do irt=0,nrootm                                                 5d27s22
          bc(iad1+irt)=hdig(ioff+ib)                                    4d4s23
         end do
        end do                                                          5d27s22
        ioff=ioff+nbdet(jsb)                                            7d25s22
       end do                                                           5d27s22
       jgdiag=jgdiag+nroot*nbdet(jsb)*nadet(isb)                        5d27s22
      end do                                                            5d27s22
      call dws_gsumf(bc(igdiag),nwds)                                   4d4s23
      call mtimesh2(bc(idva(1)),bc(iqq),nsymb,nroot,nadet,nbdet,nsbeta, 4d10s23
     $     nherec,nherect,icode,ilc,ihc,ilct,ihct,hdig,ih0e,ibc(iaorb), 4d7s23
     $     ibc(iborb),                                                  4d7s23
     $     nalpha,nbeta,i4oa,m12sym,norb,idata1,idata2,idatb1,idatb2,   4d7s23
     $     m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc,iv,bc(ie))         4d7s23
      jqq=iqq                                                           4d4s23
      jr=ir
      do isb=1,nsymb                                                    4d7s23
       jsb=nsbeta(isb)                                                  4d7s23
       ncol=nroot*nadet(isb)*nbdet(jsb)                                 4d7s23
       ncolm=ncol-1                                                     4d7s23
       iad1=irhs(isb)                                                   4d7s23
       do i=0,ncolm
        bc(jr+i)=bc(iad1+i)-bc(jqq+i)                                   4d7s23
       end do                                                           4d7s23
       jr=jr+ncol                                                       4d7s23
       jqq=jqq+ncol                                                     4d7s23
      end do                                                            4d7s23
      jgdiag=igdiag                                                     5d27s22
      jr=ir                                                             4d4s23
      jz=iz                                                             4d4s23
      jqq=iqq                                                           4d4s23
      xtest=0d0                                                         4d4s23
      do irt=0,nrootm                                                   4d4s23
       bc(iscale+irt)=0d0                                               4d4s23
       bc(idotrz+irt)=0d0                                               4d4s23
       bc(isz+irt)=0d0                                                  4d4s23
      end do                                                            4d4s23
      do isb=1,nsymb                                                    5d27s22
       jsb=nsbeta(isb)                                                  5d27s22
       if(min(nadet(isb),nbdet(jsb)).gt.0)then                          5d27s22
        ncol=nadet(isb)*nbdet(jsb)                                       5d27s22
        ncolm=ncol-1                                                     5d27s22
        iad1=irhs(isb)                                                   5d27s22
        if(icode.eq.1)then                                              6d15s22
         iad3=iv(isb)                                                     5d27s22
         do i=0,ncolm                                                     5d27s22
          do irt=0,nrootm                                                  5d27s22
           bc(iscale+irt)=bc(iscale+irt)+bc(iad1+irt)**2                       5d27s22
           bc(jgdiag+irt)=bc(jgdiag+irt)-bc(ie+irt)                     4d4s23
     $          *(1d0-bc(iad3+irt)**2)                                  4d4s23
           bc(jr+irt)=bc(iad1+irt)-bc(jqq+irt)                          4d4s23
           bc(jgdiag+irt)=1d0/bc(jgdiag+irt)                            4d4s23
           bc(jz+irt)=bc(jr+irt)*bc(jgdiag+irt)                         4d4s23
           bc(idotrz+irt)=bc(idotrz+irt)+bc(jr+irt)*bc(jz+irt)          4d4s23
           bc(isz+irt)=bc(isz+irt)+bc(jr+irt)**2                        4d4s23
          end do                                                          5d27s22
          iad1=iad1+nroot                                                 5d27s22
          jqq=jqq+nroot                                                   4d4s23
          iad3=iad3+nroot                                                 5d27s22
          jgdiag=jgdiag+nroot                                             5d27s22
          jr=jr+nroot                                                   4d4s23
          jz=jz+nroot                                                   4d4s23
         end do                                                           5d27s22
        else                                                            6d15s22
         iad3=iv(isb)                                                     5d27s22
         do i=0,ncolm                                                     5d27s22
          do irt=0,nrootm                                                  5d27s22
           bc(iscale+irt)=bc(iscale+irt)+bc(iad1+irt)**2                       5d27s22
           bc(jgdiag+irt)=bc(jgdiag+irt)-bc(ie+irt)                     4d4s23
           bc(jr+irt)=bc(iad1+irt)-bc(jqq+irt)                          4d4s23
           bc(jgdiag+irt)=1d0/bc(jgdiag+irt)                            4d4s23
           bc(jz+irt)=bc(jr+irt)*bc(jgdiag+irt)                         4d4s23
           bc(idotrz+irt)=bc(idotrz+irt)+bc(jr+irt)*bc(jz+irt)          4d4s23
           bc(isz+irt)=bc(isz+irt)+bc(jr+irt)**2                        4d4s23
          end do                                                          5d27s22
          iad1=iad1+nroot                                                 5d27s22
          jqq=jqq+nroot                                                 4d4s23
          iad3=iad3+nroot                                                 5d27s22
          jgdiag=jgdiag+nroot                                             5d27s22
          jr=jr+nroot                                                   4d4s23
          jz=jz+nroot                                                   4d4s23
         end do                                                           5d27s22
        end if                                                          6d15s22
       end if                                                           5d27s22
      end do                                                            5d27s22
      xtest=0d0                                                         4d4s23
      do irt=0,nrootm                                                   4d4s23
       bc(iscale+irt)=1d0/max(1d-12,bc(iscale+irt))                     4d7s23
       xtest=xtest+bc(isz+irt)*bc(iscale+irt)                           4d4s23
      end do                                                            4d4s23
      xtest=sqrt(xtest/dfloat(nroot))                                   4d4s23
      if(lnew)then
      maxdiis=40                                                        4d27s23
      idiis1=ibcoff                                                     4d12s23
      idiis3=idiis1+maxdiis*nwds                                        4d12s23
      ibcoff=idiis3+maxdiis*nwds                                        4d12s23
      call enough('cgvec2.diis',bc,ibc)
      do irt=0,nrootm                                                    4d12s23
       bc(isz+irt)=0d0                                                  4d12s23
      end do                                                            4d12s23
      jz=iz                                                             4d12s23
      do i=0,nconam                                                     4d12s23
       do irt=0,nrootm                                                  4d12s23
        bc(isz+irt)=bc(isz+irt)+bc(jz+irt)**2                           4d12s23
       end do                                                           4d12s23
       jz=jz+nroot                                                      4d12s23
      end do                                                            4d12s23
      do irt=0,nrootm                                                    4d12s23
       bc(isz+irt)=1d0/sqrt(bc(isz+irt))                                4d12s23
      end do                                                            4d12s23
      do i=0,nconarm                                                    4d12s23
       bc(idiis1+i)=bc(idva(1)+i)                                       4d12s23
       bc(idiis3+i)=bc(iqq+i)                                           4d12s23
      end do                                                            4d12s23
      jdii=idiis1+nwds                                                  4d12s23
      jz=iz                                                             4d12s23
      do i=0,nconam                                                     4d12s23
       do irt=0,nrootm                                                  4d12s23
        bc(jdii+irt)=bc(jz+irt)*bc(isz+irt)                             4d12s23
       end do                                                           4d12s23
       jdii=jdii+nroot                                                  4d12s23
       jz=jz+nroot                                                      4d12s23
      end do                                                            4d12s23
      iterm=1                                                           4d12s23
 1000 continue                                                          4d12s23
      if(lprt)write(6,*)('starting iterm '),iterm
      iterm=iterm+1                                                     4d12s23
      if(iterm.gt.maxdiis)then
       call dws_synca
       call dws_finalize
       stop 'maxdiis:cgvec2'                                            4d13s23
      end if
      iin=idiis1+nwds*(iterm-1)                                         4d12s23
      iout=idiis3+nwds*(iterm-1)                                        4d12s23
      call mtimesh2(bc(iin),bc(iout),nsymb,nroot,nadet,nbdet,           4d12s23
     $       nsbeta,                                                    4d12s23
     $     nherec,nherect,icode,ilc,ihc,ilct,ihct,hdig,ih0e,ibc(iaorb), 4d7s23
     $     ibc(iborb),                                                  4d7s23
     $     nalpha,nbeta,i4oa,m12sym,norb,idata1,idata2,idatb1,idatb2,   4d7s23
     $     m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc,iv,bc(ie))         4d7s23
      nfit=iterm                                                        4d12s23
      idpt=idiis1+nwds*nfit                                             4d12s23
      ifmat=ibcoff                                                      4d12s23
      icoef=ifmat+ncona*nfit                                             4d12s23
      iwgt=icoef+nfit                                                   4d12s23
      iscr=iwgt+ncona                                                    4d12s23
      idat=iscr+2*nfit+nfit*nfit+2*nfit*ncona+ncona                       4d12s23
      ibcoff=idat+ncona                                                  4d12s23
      call enough('cgvec2.ls',bc,ibc)                                   4d12s23
      nok=0                                                             4d13s23
      do irt=0,nrootm                                                   4d12s23
       do i=0,nfit-1                                                     4d12s23
        iad=ifmat+i*ncona                                               4d12s23
        iad3=idiis3+nwds*i+irt                                          4d12s23
        do j=0,nconam
         bc(iad+j)=bc(iad3)                                             4d12s23
         iad3=iad3+nroot                                                4d12s23
        end do                                                          4d12s23
       end do                                                           4d12s23
       iad4=irhs(1)+irt                                                 4d12s23
       do i=0,nconam
        bc(idat+i)=bc(iad4)                                             4d12s23
        iad4=iad4+nroot                                                 4d12s23
        bc(iwgt+i)=1d0                                                   4d12s23
       end do                                                            4d12s23
       iuse=0
       call lsqfit2(bc(ifmat),ncona,nfit,bc(idat),ncona,1,ncona,        4d12s23
     $      bc(icoef),nfit,bc(iscr),bc(iwgt),0,rms,bc,ibc)              4d12s23
       call dgemm('n','n',ncona,1,nfit,1d0,bc(ifmat),ncona,bc(icoef),   4d12s23
     $      nfit,0d0,bc(iscr),ncona,'cgvec2.dgfit')                     4d12s23
       xx=0d0                                                           4d12s23
       jrhs=irhs(1)+irt                                                 4d13s23
       jgdiag=igdiag+irt                                                4d13s23
       do i=0,nconam                                                    4d12s23
        diff=bc(jrhs)-bc(iscr+i)                                        4d12s23
        jrhs=jrhs+nroot                                                 4d12s23
        bc(iscr+i)=bc(jgdiag)*diff                                      4d13s23
        jgdiag=jgdiag+nroot                                             4d13s23
        xx=xx+diff**2                                                   4d12s23
       end do                                                           4d12s23
       xtest=sqrt(xx*bc(iscale+irt))                                    4d13s23
       call dws_bcast(xtest,1)                                          6d22s23
       if(lprt)write(6,*)('for root '),irt+1,xtest,tol,xx,bc(iscale+irt)
       if(xtest.lt.tol)then
        nok=nok+1
       else if(iterm.eq.maxdiis-1)then                                  4d27s23
        nok=nok+1                                                       4d27s23
        tolb=max(tolb,xtest)                                            4d27s23
       end if
       do i=0,nfit-1                                                     4d12s23
        iad=ifmat+i*ncona                                               4d12s23
        iad1=idiis1+nwds*i+irt                                          4d12s23
        do j=0,nconam
         bc(iad+j)=bc(iad1)                                             4d12s23
         iad1=iad1+nroot                                                4d12s23
        end do                                                          4d12s23
       end do                                                           4d12s23
       if(xtest.lt.tol.or.iterm.eq.maxdiis-1)then                       4d27s23
        jppp=ippp+ncona*irt                                             4d13s23
        call dgemm('n','n',ncona,1,nfit,1d0,bc(ifmat),ncona,bc(icoef),  4d13s23
     $       nfit,0d0,bc(jppp),ncona,'cgvec2.jppp')                     4d13s23
       end if                                                           4d13s23
       do i=0,nfit-1                                                    4d12s23
        dot=0d0                                                         4d12s23
        iad=ifmat+ncona*i                                               4d12s23
        do j=0,nconam                                                   4d12s23
         dot=dot+bc(iad+j)*bc(iscr+j)                                   4d12s23
        end do                                                          4d12s23
        do j=0,nconam
         bc(iscr+j)=bc(iscr+j)-dot*bc(iad+j)                            4d12s23
        end do                                                          4d12s23
       end do                                                           4d12s23
       sz=0d0                                                           4d12s23
       do j=0,nconam                                                    4d12s23
        sz=sz+bc(iscr+j)**2                                             4d12s23
       end do
       sz=1d0/sqrt(sz)                                                  4d12s23
       do j=0,nconam                                                    4d12s23
        bc(iscr+j)=bc(iscr+j)*sz                                        4d12s23
       end do                                                           4d12s23
       jdpt=idpt+irt                                                    4d12s23
       do i=0,nconam                                                    4d12s23
        bc(jdpt)=bc(iscr+i)                                             4d12s23
        jdpt=jdpt+nroot                                                 4d12s23
       end do                                                           4d12s23
      end do                                                            4d12s23
      ibcoff=ifmat                                                      4d12s23
      if(lprt)write(6,*)('nok = '),nok
      if(nok.eq.nroot)then                                              4d13s23
       if(lprt)then
        write(6,*)('ippp ')
        call prntm2(bc(ippp),ncona,nroot,ncona)                          4d13s23
       end if
       ixx=idva(1)                                                      4d13s23
       do i=0,nconam                                                    4d13s23
        jppp=ippp+i                                                     4d13s23
        do irt=0,nrootm                                                  4d13s23
         bc(ixx+irt)=bc(jppp)                                           4d13s23
         jppp=jppp+ncona                                                4d13s23
        end do                                                          4d13s23
        ixx=ixx+nroot                                                   4d13s23
       end do                                                           4d13s23
       if(lprt)then
        write(6,*)('cgvec2 diis iterations converged after '),
     $      iterm,('iterations to tolerence '),tolb                     4d27s23
       end if                                                           4d27s23
       go to 1001                                                       4d13s23
      end if                                                            4d13s23
      go to 1000
 1001 continue
      else                                                              4d13s23
      ifirst=0                                                          5d31s22
      macit=0
      itmax=5
      ierror=ibcoff                                                     5d19s22
      ibcoff=ierror+itmax                                               5d19s22
      call enough('cgvec2.  2',bc,ibc)
      ratio=0d0                                                         4d10s23
  100 continue
      macit=macit+1
      if(macit.gt.100)then
       write(6,*)('too many restarts!!!')
       call dws_synca
       call dws_finalize
       stop 'cgvec2'
      end if
      do i=0,nconarm                                                     5d27s22
       bc(ix+i)=0d0
       bc(ip+i)=bc(iz+i)                                                4d4s23
      end do
      iter=0
    1 continue
       iter=iter+1
       if(iter.gt.itmax)then                                            4d4s23
         if(mynowprog.eq.0.and..not.lprt)then                                4d12s23
          write(6,52)(bc(ierror+i),i=0,iter-2)                           5d19s22
         end if
        ratio=bc(ierror+iter-2)/bc(ierror)                              4d10s23
        jx=ix                                                           4d4s23
        do isb=1,nsymb                                                  5d27s22
         jsb=nsbeta(isb)                                                5d27s22
         if(min(nadet(isb),nbdet(jsb)).gt.0)then                        5d27s22
          ncol=nadet(isb)*nbdet(jsb)                                     5d27s22
          iad1=idva(isb)                                                5d27s22
          do i=0,ncol*nroot-1                                           4d4s23
           bc(iad1+i)=bc(iad1+i)+bc(jx+i)                               4d4s23
          end do                                                        5d27s22
          jx=jx+ncol*nroot                                              4d4s23
         end if                                                         5d27s22
        end do                                                          5d27s22
        call mtimesh2(bc(idva(1)),bc(iqq),nsymb,nroot,nadet,nbdet,      4d12s23
     $       nsbeta,                                                    4d12s23
     $     nherec,nherect,icode,ilc,ihc,ilct,ihct,hdig,ih0e,ibc(iaorb), 4d7s23
     $     ibc(iborb),                                                  4d7s23
     $     nalpha,nbeta,i4oa,m12sym,norb,idata1,idata2,idatb1,idatb2,   4d7s23
     $     m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc,iv,bc(ie))         4d7s23
        jqq=iqq                                                         4d4s23
        jr=ir
        jz=iz                                                           4d12s23
        jgdiag=igdiag                                                   4d12s23
        do i=0,nrootm                                                   4d12s23
         bc(idotrz+i)=0d0                                               4d12s23
         bc(isz+i)=0d0                                                  4d12s23
        end do                                                          4d12s23
        do isb=1,nsymb                                                  5d27s22
         jsb=nsbeta(isb)                                                5d27s22
         ncol=nadet(isb)*nbdet(jsb)                                     5d27s22
         ncolm=ncol-1                                                   5d27s22
         iad1=iga(isb)                                                  5d27s22
         iad3=idva(isb)                                                 5d27s22
         iad4=irhs(isb)                                                 5d27s22
         do i=0,ncolm                                                   4d12s23
          do irt=0,nrootm                                                4d12s23
           bc(jr+irt)=bc(iad4+irt)-bc(jqq+irt)                                 4d12s23
           bc(isz+irt)=bc(isz+irt)+bc(jr+irt)**2                        4d12s23
           bc(jz+irt)=bc(jz+irt)*bc(jgdiag+irt)                                4d12s23
           bc(idotrz+irt)=bc(idotrz+irt)+bc(jr+irt)*bc(jz+irt)          4d12s23
          end do                                                        4d12s23
          jr=jr+nroot                                                   4d12s23
          jqq=jqq+nroot                                                 4d12s23
          jz=jz+nroot                                                   4d12s23
          jgdiag=jgdiag+nroot                                           4d12s23
          iad4=iad4+nroot
         end do                                                         4d12s23
        end do                                                          4d12s23
        xx=0d0
        do irt=0,nrootm
         xx=xx+bc(isz+irt)*bc(iscale+irt)
        end do
        write(6,*)('new sz'),xx,sqrt(xx/dfloat(nroot))
        call prntm2(bc(isz),1,nroot,1)
        write(6,*)('new dotrz ')
        call prntm2(bc(idotrz),1,nroot,1)
        xtest=0d0                                                       4d4s23
        do irt=0,nrootm                                                 4d4s23
         xtest=xtest+bc(isz+irt)*bc(iscale+irt)                         4d4s23
        end do                                                          4d4s23
        xtest=sqrt(xtest/dfloat(nroot))                                 4d4s23
        write(6,*)('new xtest '),xtest
        go to 100                                                       5d9s22
       end if
       if(iter.ne.1)then                                                4d4s23
        do irt=0,nrootm                                                 4d4s23
         bc(ibeta+irt)=bc(idotrznew+irt)/bc(idotrz+irt)                 4d4s23
         bc(idotrz+irt)=bc(idotrznew+irt)                               4d4s23
        end do                                                          4d4s23
        write(6,*)('beta ')
        call prntm2(bc(ibeta),1,nroot,1)
        jp=ip                                                           4d4s23
        jz=iz                                                           4d4s23
        do i=0,nconam                                                   4d4s23
         do irt=0,nrootm                                                4d4s23
          bc(jp+irt)=bc(jz+irt)+bc(ibeta+irt)*bc(jp+irt)                4d4s23
         end do                                                         4d4s23
         jp=jp+nroot                                                    4d4s23
         jz=jz+nroot                                                    4d4s23
        end do                                                          4d4s23
       end if                                                           4d4s23
       call mtimesh2(bc(ip),bc(iqq),nsymb,nroot,nadet,nbdet,nsbeta,     4d10s23
     $     nherec,nherect,icode,ilc,ihc,ilct,ihct,hdig,ih0e,ibc(iaorb), 4d7s23
     $     ibc(iborb),                                                  4d7s23
     $     nalpha,nbeta,i4oa,m12sym,norb,idata1,idata2,idatb1,idatb2,   4d7s23
     $     m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc,iv,bc(ie))         4d7s23
       do irt=0,nrootm                                                  4d4s23
        bc(ialpha+irt)=0d0                                                 4d4s23
       end do                                                           4d4s23
       jqq=iqq                                                            5d27s22
       jp=ip                                                            5d27s22
       do i=0,nconam
        do irt=0,nrootm
         bc(ialpha+irt)=bc(ialpha+irt)+bc(jp+irt)*bc(jqq+irt)
        end do
        jp=jp+nroot
        jqq=jqq+nroot
       end do
       write(6,*)('alpha '),(bc(ialpha+irt),irt=0,nrootm)
       do irt=0,nrootm
        bc(ialpha+irt)=bc(idotrz+irt)/bc(ialpha+irt)
        bc(isz+irt)=0d0                                                 4d4s23
        bc(idotrznew+irt)=0d0                                           4d4s23
       end do
       write(6,*)('alpha')
       call prntm2(bc(ialpha),1,nroot,1)
       jx=ix                                                            4d4s23
       jp=ip                                                            4d4s23
       jr=ir                                                            4d4s23
       jz=iz                                                            4d4s23
       jqq=iqq
       jgdiag=igdiag
       do i=0,nconam
        do irt=0,nrootm                                                 4d4s23
         bc(jx+irt)=bc(jx+irt)+bc(ialpha+irt)*bc(jp+irt)                4d4s23
         bc(jr+irt)=bc(jr+irt)-bc(ialpha+irt)*bc(jqq+irt)
         bc(jz+irt)=bc(jr+irt)*bc(jgdiag+irt)                           4d4s23
         bc(isz+irt)=bc(isz+irt)+bc(jr+irt)**2
         bc(idotrznew+irt)=bc(idotrznew+irt)+bc(jr+irt)*bc(jz+irt)
        end do
        jx=jx+nroot
        jp=jp+nroot
        jr=jr+nroot
        jqq=jqq+nroot
        jz=jz+nroot
        jgdiag=jgdiag+nroot
       end do                                                           4d4s23
       write(6,*)('sz now ')
       call prntm2(bc(isz),1,nroot,1)
       write(6,*)('dotrznew ')
       call prntm2(bc(idotrznew),1,nroot,1)
       xdotnew=0d0                                                      4d4s23
       do irt=0,nrootm                                                  4d4s23
        xdotnew=xdotnew+bc(isz+irt)*bc(iscale+irt)                      4d4s23
       end do                                                           4d4s23
       call dws_bcast(xdotnew,1)                                            5d27s22
       xtest=sqrt(xdotnew/dfloat(nroot))                                4d4s23
       bc(ibcoff)=xtest                                                 4d10s23
       bc(ibcoff+1)=ratio                                               4d10s23
       call dws_bcast(bc(ibcoff),2)                                     4d10s23
       xtest=bc(ibcoff)                                                 4d10s23
       ratio=bc(ibcoff+1)                                               4d10s23
       bc(ierror+iter-1)=xtest                                          4d4s23
    2  format('after iteration ',i2,' error is ',es10.3)
       write(6,*)('ratio '),ratio,xtest,tol
       if(xtest.gt.tol.and.ratio.lt.0.9d0)go to 1                         4d10s23
       if(mynowprog.eq.0.and.lprt)then                                  4d12s23
        if(xtest.le.tol)then                                            4d12s23
         write(6,52)(bc(ierror+i),i=0,iter-1)                            5d19s22
         write(6,*)('cg iterations converged after '),iter,
     $      (' iterations and '),macit-1,(' restarts')
   52    format(10es8.1)
        else                                                            4d12s23
         write(6,*)('cg iterations stagnated after '),iter,
     $      (' iterations and '),macit-1,(' restarts')
        end if                                                          4d12s23
      end if                                                             3d3s17
      end if                                                            4d13s23
      do isb=1,nsymb                                                    5d27s22
       jsb=nsbeta(isb)                                                  5d27s22
       if(min(nadet(isb),nbdet(jsb)).gt.0)then                          5d27s22
        ncol=nadet(isb)*nbdet(jsb)                                      5d27s22
        if(.not.lnew)then                                               4d13s23
         ncolm=ncol-1                                                    5d27s22
         iad1=idva(isb)                                                  5d27s22
         do i=0,ncolm                                                    5d27s22
          do irt=0,nrootm                                                 5d27s22
           bc(iad1+irt)=bc(iad1+irt)+bc(jx+irt)                          4d4s23
          end do                                                         5d27s22
          iad1=iad1+nroot                                                5d27s22
          jx=jx+nroot                                                    4d4s23
         end do                                                          5d27s22
        end if                                                          4d13s23
        if(lprt)then                                                    7d28s22
         write(6,*)('solution for symmetry '),isb,idva(isb)
         call prntm2(bc(idva(isb)),nroot*ncol,1,nroot*ncol)
        end if                                                          7d28s22
       end if                                                           5d27s22
      end do                                                            5d27s22
      ibcoff=ix                                                         4d4s23
      return
      end
