c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cgvec(ham,rhs,istart,nhere,nps,e,nroot,v,dv,tolvx,lsym, 7d20s22
     $     lprnt,bc,ibc)                                                11d9s22
      implicit real*8 (a-h,o-z)
c
c     solve (H-e+e*v*vt)dv=rhs for dv
c
      logical lsym,lprnt,ldebug,lnew                                    4d13s23
      dimension ham(nhere,*),rhs(nps,*),v(nps,*),dv(nps,*),e(*)         7d20s22
      include "common.store"
      common/singcm/iuse,nff
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data icall/0/
      save icall
      icall=icall+1
      tolv=tolvx                                                        4d9s23
      lnew=.true.
      if(lprnt)then                                                     7d14s22
       write(6,*)('Hi, my name is cgvec ...'),icall,ibcoff
       if(lnew)then
        write(6,*)('preconditioned diis gmres ')
       else
        write(6,*)('conjugate gradient method')
       end if
       write(6,*)('initial guess for dv: ')
       write(6,*)('lsym = '),lsym
       call prntm2(dv,nps,nroot,nps)
      end if                                                            7d14s22
      if(nps.eq.1)then                                                  5d20s22
       dv(1,1)=0d0                                                      7d20s22
       return                                                           5d20s22
      end if                                                            5d20s22
      npsm=nps-1
      nwds=nps*nroot                                                    4d10s23
      npsrm=nwds-1                                                      4d10s23
      ix=ibcoff                                                         4d4s23
      ir=ix+nwds                                                        4d10s23
      iz=ir+nwds                                                        4d10s23
      iq=iz+nwds                                                        4d10s23
      iqq=iq+nwds                                                       4d10s23
      ipxr=iqq+nwds                                                     4d10s23
      igdiag=ipxr+nwds                                                  4d10s23
      ibeta=igdiag+nwds                                                 4d10s23
      iscale=ibeta+nroot                                                4d3s23
      idotrz=iscale+nroot                                               4d4s23
      idotrznew=idotrz+nroot                                            4d4s23
      ibest1=idotrznew+nroot                                            4d10s23
      ibest2=ibest1+nroot                                               4d10s23
      ibcoff=ibest2+nwds                                                4d10s23
      call enough('cgvec.  1',bc,ibc)
      do i=0,npsm
       bc(igdiag+i)=0d0
      end do
      jbest1=ibest1-1                                                   4d10s23
      do irt=1,nroot                                                     7d20s22
       scale=0d0                                                         4d3s23
       bc(jbest1+irt)=1d10                                              4d10s23
       do i=1,nps                                                       4d4s23
        scale=scale+rhs(i,irt)**2                                       4d4s23
       end do
       bc(iscale+irt-1)=1d0/max(1d-4,scale)                             4d10s23
       if(scale.ne.scale)then                                           10d30s24
        write(6,*)('ooops got NaN in rhs! ')                            10d30s24
        call prntm2(rhs,nps,nroot,nps)                                  10d30s24
        call dws_synca                                                  10d30s24
        call dws_finalize                                               10d30s24
        stop 'cgvec'                                                    10d30s24
       end if                                                           10d30s24
       if(scale.lt.1d-10)bc(iscale+irt-1)=0d0                           10d30s24
      end do                                                            7d20s22
      jgdiag=igdiag+istart-2
      do i=1,nhere
       ip=i+istart-1                                                    7d14s22
       bc(jgdiag+i)=ham(i,ip)                                           7d14s22
      end do
      call dws_gsumf(bc(igdiag),nps)
      if(lprnt)then
       write(6,*)('diagonal elements of ham: ')
       call prntm2(bc(igdiag),1,nps,1)
       write(6,*)('diagonal elements of H-e+e*v*vt ')
      end if                                                            7d14s22
      jgdiag=igdiag+nps                                                 4d4s23
      do irt=1,nroot-1                                                   7d20s22
       do i=0,npsm                                                      4d4s23
        bc(jgdiag+i)=bc(igdiag+i)                                       7d20s22
       end do                                                           7d20s22
       jgdiag=jgdiag+nps                                                4d4s23
      end do                                                            7d20s22
      if(lsym)then                                                      6d23s22
       jgdiag=igdiag                                                    7d20s22
       do irt=1,nroot                                                    7d20s22
        do i=0,npsm
         ip=i+1
         bc(jgdiag+i)=bc(jgdiag+i)-e(irt)*(1d0-v(ip,irt)*v(ip,irt))        7d20s22
        end do
        jgdiag=jgdiag+nps                                               4d4s23
       end do                                                           7d20s22
      else                                                              6d23s22
       jgdiag=igdiag                                                    7d20s22
       do irt=1,nroot                                                    7d20s22
        do i=0,npsm
         ip=i+1
         bc(jgdiag+i)=bc(jgdiag+i)-e(irt)                                7d20s22
        end do
        jgdiag=jgdiag+nps                                               4d4s23
       end do                                                           7d20s22
      end if                                                            6d23s22
      if(lprnt)then
       write(6,*)('final diagonals ')
       call prntm2(bc(igdiag),nps,nroot,nps)                            7d20s22
      end if                                                            7d14s22
      do i=0,npsrm                                                      4d4s23
       bc(igdiag+i)=1d0/bc(igdiag+i)                                    5d19s22
      end do                                                            5d19s22
      if(lnew)then                                                      4d13s23
       maxdiis=120                                                        4d12s23
       idiis1=ibcoff                                                     4d12s23
       idiis3=idiis1+maxdiis*nwds                                        4d12s23
       ibcoff=idiis3+maxdiis*nwds                                        4d12s23
       call enough('cgvec.diis',bc,ibc)
       irest=0                                                          4d14s23
 1010  continue
       jdiis1=idiis1-1                                                  4d13s23
       irest=irest+1
       if(irest.gt.100)then
        write(6,*)('too many restarts in cgvec')
        call dws_synca
        call dws_finalize
        stop
       end if
       do irt=1,nroot                                                   4d13s23
        sz=0d0
        do i=1,nps                                                       4d13s23
         sz=sz+dv(i,irt)**2                                             4d13s23
        end do                                                          4d13s23
        if(sz.ne.0d0)then                                               4d13s23
         sz=1d0/sqrt(sz)                                                4d13s23
         do i=1,nps                                                      4d13s23
          bc(jdiis1+i)=sz*dv(i,irt)                                      4d13s23
         end do                                                          4d13s23
        else                                                            4d13s23
         sz=0d0                                                         4d13s23
         jgdiag=igdiag-1+nps*(irt-1)                                    4d13s23
         do i=1,nps                                                     4d13s23
          bc(jdiis1+i)=bc(jgdiag+i)*rhs(i,irt)                          4d13s23
          sz=sz+bc(jdiis1+i)**2                                         4d13s23
         end do                                                         4d13s23
         sz=1d0/sqrt(sz)                                                4d13s23
         do i=1,nps                                                     4d13s23
          bc(jdiis1+i)=sz*bc(jdiis1+i)                                  4d13s23
         end do                                                         4d13s23
        end if                                                          4d13s23
        jdiis1=jdiis1+nps                                               4d13s23
       end do                                                           4d13s23
       iterm=0                                                          4d13s23
       nrootm=nroot-1                                                   4d13s23
 1000  continue                                                         4d13s23
        iterm=iterm+1                                                   4d13s23
        if(iterm.gt.maxdiis)then                                        4d13s23
         if(lprnt)then
          write(6,*)('ipxr ')
          call prntm2(bc(ipxr),nps,nroot,nps)                               4d13s23
         end if
         jp=ipxr-1                                                        4d13s23
         do irt=1,nroot                                                   4d13s23
          do i=1,nps                                                      4d13s23
           dv(i,irt)=bc(jp+i)                                             4d13s23
          end do                                                          4d13s23
          jp=jp+nps                                                       4d13s23
         end do                                                           4d13s23
         go to 1010                                                     4d14s23
        end if                                                          4d13s23
        iin=idiis1+nwds*(iterm-1)                                        4d13s23
        iout=idiis3+nwds*(iterm-1)                                      4d13s23
        call mtimesh(bc(iin),bc(iout),istart,nhere,nroot,nps,ham,lsym,v,4d13s23
     $       e,nwds)                                                    4d13s23
        nfit=iterm                                                      4d13s23
        ifmat=ibcoff                                                    4d13s23
        icoef=ifmat+nps*nfit                                            4d13s23
        iwgt=icoef+nfit                                                   4d12s23
        iscr=iwgt+nps                                                    4d12s23
        ibcoff=iscr+2*nfit+nfit*nfit+2*nfit*nps+nps                       4d13s23
        call enough('cgvec2.ls',bc,ibc)                                   4d12s23
        nok=0                                                             4d13s23
        do irt=0,nrootm                                                   4d12s23
         irp=irt+1                                                      4d13s23
         do i=0,nfit-1                                                     4d12s23
          iad=ifmat+i*nps                                               4d13s23
          iad3=idiis3+nwds*i+nps*irt                                    4d13s23
          do j=0,npsm
           bc(iad+j)=bc(iad3+j)                                         4d13s23
          end do                                                          4d12s23
         end do                                                           4d12s23
         do i=0,npsm                                                    4d13s23
          bc(iwgt+i)=1d0                                                   4d12s23
         end do                                                            4d12s23
         iuse=0
         call lsqfit2(bc(ifmat),nps,nfit,rhs(1,irp),nps,1,nps,          4d13s23
     $      bc(icoef),nfit,bc(iscr),bc(iwgt),0,rms,bc,ibc)              4d12s23
         call dgemm('n','n',nps,1,nfit,1d0,bc(ifmat),nps,bc(icoef),     4d13s23
     $      nfit,0d0,bc(iscr),nps,'cgvec.dgfit')                        4d13s23
         xx=0d0                                                           4d12s23
         jgdiag=igdiag+nps*irt-1                                        4d13s23
         jscr=iscr-1                                                    4d13s23
         do i=1,nps                                                     4d13s23
          diff=rhs(i,irp)-bc(jscr+i)                                        4d12s23
          bc(jscr+i)=bc(jgdiag+i)*diff                                      4d13s23
          xx=xx+diff**2                                                   4d12s23
         end do                                                           4d12s23
         xtest=sqrt(xx*bc(iscale+irt))                                    4d13s23
         call dws_bcast(xtest,1)                                        6d22s23
         if(xtest.lt.tolv)then                                          4d13s23
          nok=nok+1
         end if
         do i=0,nfit-1                                                     4d12s23
          iad=ifmat+i*nps                                               4d13s23
          iad1=idiis1+nwds*i+irt*nps                                    4d13s23
          do j=0,npsm                                                   4d13s23
           bc(iad+j)=bc(iad1+j)                                         4d13s23
          end do                                                          4d12s23
         end do                                                           4d12s23
       if(xtest.lt.tolv.or.iterm.eq.maxdiis)then                        4d14s23
        jp=ipxr+nps*irt                                                 4d13s23
        call dgemm('n','n',nps,1,nfit,1d0,bc(ifmat),nps,bc(icoef),      4d13s23
     $       nfit,0d0,bc(jp),nps,'cgvec.jp')                            4d13s23
       end if                                                           4d13s23
       do i=0,nfit-1                                                    4d12s23
        dot=0d0                                                         4d12s23
        iad=ifmat+nps*i                                                 4d12s23
        do j=0,npsm                                                     4d12s23
         dot=dot+bc(iad+j)*bc(iscr+j)                                   4d12s23
        end do                                                          4d12s23
        do j=0,npsm
         bc(iscr+j)=bc(iscr+j)-dot*bc(iad+j)                            4d12s23
        end do                                                          4d12s23
       end do                                                           4d12s23
       sz=0d0                                                           4d12s23
       do j=0,npsm                                                      4d12s23
        sz=sz+bc(iscr+j)**2                                             4d12s23
       end do
       sz=1d0/sqrt(sz)                                                  4d12s23
       do j=0,npsm                                                      4d12s23
        bc(iscr+j)=bc(iscr+j)*sz                                        4d12s23
       end do                                                           4d12s23
       jdpt=idiis1+nwds*nfit+nps*irt                                    4d13s23
       do i=0,npsm                                                      4d12s23
        bc(jdpt+i)=bc(iscr+i)                                             4d12s23
       end do                                                           4d12s23
      end do                                                            4d12s23
      ibcoff=ifmat                                                      4d12s23
      if(lprnt)write(6,*)('nok = '),nok
      if(nok.eq.nroot)then                                              4d13s23
       if(lprnt)then
        write(6,*)('ipxr ')
        call prntm2(bc(ipxr),nps,nroot,nps)                               4d13s23
       end if
       jp=ipxr-1                                                        4d13s23
       do irt=1,nroot                                                   4d13s23
        do i=1,nps                                                      4d13s23
         dv(i,irt)=bc(jp+i)                                             4d13s23
        end do                                                          4d13s23
        jp=jp+nps                                                       4d13s23
       end do                                                           4d13s23
       if(lprnt)write(6,*)('cgvec diis iterations converged after '),
     $      iterm,('iterations and '),irest-1,('restarts')
        else                                                            4d13s23
         go to 1000                                                     4d13s23
        end if                                                          4d13s23
      else                                                              4d13s23
      call mtimesh(dv,bc(iq),istart,nhere,nroot,nps,ham,lsym,v,e,nwds)  4d10s23
      jr=ir-1                                                           4d4s23
      jq=iq-1                                                           4d4s23
      jz=iz-1                                                           4d4s23
      jgdiag=igdiag-1                                                   4d4s23
      xtest=0d0                                                         4d4s23
      do irt=0,nroot-1                                                  4d4s23
       irp=irt+1                                                        4d4s23
       dotrz=0d0                                                        4d4s23
       sz=0d0                                                           4d4s23
       jbest2=ibest2-1+nps*irt                                          4d10s23
       do i=1,nps                                                       4d4s23
        bc(jbest2+i)=dv(i,irp)                                          4d10s23
        bc(jr+i)=rhs(i,irp)-bc(jq+i)                                    4d10s23
        bc(jz+i)=bc(jr+i)*bc(jgdiag+i)                                  4d4s23
        dotrz=dotrz+bc(jr+i)*bc(jz+i)                                   4d4s23
        sz=sz+bc(jr+i)**2                                               4d4s23
       end do                                                           4d4s23
       xtest=xtest+sz*bc(iscale+irt)                                    4d4s23
       bc(idotrz+irt)=dotrz                                             4d4s23
       bc(ibest1+irt)=sz                                                4d10s23
       jr=jr+nps                                                        4d4s23
       jq=jq+nps                                                        4d4s23
       jz=jz+nps                                                        4d4s23
       jgdiag=jgdiag+nps                                                4d4s23
      end do                                                            4d4s23
      xtest=sqrt(xtest/dfloat(nroot))                                   4d4s23
      macit=0
      itmax=10
      ifirst=0                                                          7d18s22
      ierror=ibcoff                                                     5d19s22
      ibcoff=ierror+itmax                                               5d19s22
      call enough('cgvec.  2',bc,ibc)
      ratio=0d0                                                         4d10s23
  100 continue
      macit=macit+1
      if(macit.gt.100)then
       write(6,*)('too many restarts!!!'),tolv
       call dws_synca
       call dws_finalize
       stop 'cgvec'
      end if
      do i=0,npsrm                                                      4d4s23
       bc(ix+i)=0d0
       bc(ipxr+i)=bc(iz+i)                                              4d4s23
      end do
      iter=0
      tol=tolv                                                          6d23s22
    1 continue
       iter=iter+1
       if(iter.gt.itmax)then                                            4d10s23
        if(mynowprog.eq.0.and.lprnt)then                                4d10s23
         write(6,52)(bc(ierror+i),i=0,iter-2)                           5d19s22
        end if
        ratio=bc(ierror+iter-2)/bc(ierror)                              4d10s23
        jx=ix-1                                                         4d4s23
        do irt=1,nroot                                                   3d22s23
         do i=1,nps                                                     4d4s23
          dv(i,irt)=dv(i,irt)+bc(jx+i)                                  4d4s23
         end do                                                         3d22s23
         jx=jx+nps                                                      4d4s23
        end do                                                          3d22s23
        call mtimesh(dv,bc(iq),istart,nhere,nroot,nps,ham,lsym,v,e,     4d10s23
     $       nwds)                                                      4d10s23
        jr=ir-1
        jz=iz-1
        jq=iq-1
        xtest=0d0                                                       4d4s23
        jgdiag=igdiag-1                                                 4d4s23
        do irt=1,nroot                                                   3d22s23
         irm=irt-1                                                      4d4s23
         dotrz=0d0
         sz=0d0
         do i=1,nps                                                     4d4s23
          bc(jr+i)=rhs(i,irt)-bc(jq+i)                                  4d4s23
          bc(jz+i)=bc(jr+i)*bc(jgdiag+i)
          dotrz=dotrz+bc(jr+i)*bc(jz+i)
          sz=sz+bc(jr+i)**2
         end do                                                         3d22s23
         if(sz.lt.bc(jbest1+irt))then                                   4d10s23
          bc(jbest1+irt)=sz                                             4d10s23
          jbest2=ibest2-1+nps*(irt-1)                                   4d10s23
          do i=1,nps                                                    4d10s23
           bc(jbest2+i)=dv(i,irt)                                       4d10s23
          end do                                                        4d10s23
         end if                                                         4d10s23
         bc(idotrz+irm)=dotrz
         xtest=xtest+sz*bc(iscale+irm)                                  4d4s23
         jr=jr+nps                                                      4d4s23
         jz=jz+nps                                                      4d4s23
         jq=jq+nps                                                      4d4s23
         jgdiag=jgdiag+nps                                              4d4s23
        end do                                                          3d22s23
        xtest=sqrt(xtest/dfloat(nroot))                                 4d4s23
        go to 100                                                       5d9s22
       end if
       if(iter.ne.1)then
        jpxr=ipxr                                                       4d4s23
        jz=iz                                                           4d4s23
        do irt=0,nroot-1                                                4d4s23
         bc(ibeta+irt)=bc(idotrznew+irt)/bc(idotrz+irt)                 4d4s23
         bc(idotrz+irt)=bc(idotrznew+irt)                               4d4s23
         do i=0,npsm                                                    4d4s23
          bc(jpxr+i)=bc(jz+i)+bc(ibeta+irt)*bc(jpxr+i)                  4d4s23
         end do                                                         4d4s23
         jpxr=jpxr+nps                                                  4d4s23
         jz=jz+nps                                                      4d4s23
        end do                                                          4d4s23
       end if                                                           4d4s23
       call mtimesh(bc(ipxr),bc(iq),istart,nhere,nroot,nps,ham,lsym,v,e,4d10s23
     $       nwds)                                                      4d10s23
       jp=ipxr                                                          4d4s23
       jq=iq                                                            4d4s23
       jx=ix
       jr=ir
       jz=iz                                                            4d4s23
       jgdiag=igdiag                                                    4d4s23
       xdotnew=0d0                                                      4d4s23
       iupd=0
       do irt=0,nroot-1                                                 4d4s23
        irtp=irt+1                                                      4d10s23
        alpha=0d0                                                       4d4s23
        do i=0,npsm                                                     4d4s23
         alpha=alpha+bc(jp+i)*bc(jq+i)                                  4d4s23
        end do                                                          4d4s23
        alpha=bc(idotrz+irt)/alpha                                      4d4s23
        dotrznew=0d0                                                    4d4s23
        sz=0d0                                                          4d4s23
        do i=0,npsm                                                     4d4s23
         bc(jx+i)=bc(jx+i)+alpha*bc(jp+i)                               4d4s23
         bc(jr+i)=bc(jr+i)-alpha*bc(jq+i)                               4d4s23
         bc(jz+i)=bc(jr+i)*bc(jgdiag+i)                                 4d4s23
         sz=sz+bc(jr+i)**2                                              4d4s23
         dotrznew=dotrznew+bc(jr+i)*bc(jz+i)                            4d4s23
        end do                                                          4d4s23
        if(sz.lt.bc(ibest1+irt))then                                    4d10s23
         iupd=iupd+1
         bc(ibest1+irt)=sz                                              4d10s23
         jbest2=ibest2-1+nps*irt                                        4d10s23
         jxm=jx-1                                                       4d10s23
         do i=1,nps                                                     4d10s23
          bc(jbest2+i)=dv(i,irtp)+bc(jxm+i)                             4d10s23
         end do                                                         4d10s23
        end if                                                          4d10s23
        jp=jp+nps                                                       4d4s23
        jq=jq+nps                                                       4d4s23
        jx=jx+nps                                                       4d4s23
        jr=jr+nps                                                       4d4s23
        jz=jz+nps                                                       4d4s23
        jgdiag=jgdiag+nps                                               4d4s23
        xdotnew=xdotnew+sz*bc(iscale+irt)                               4d4s23
        bc(idotrznew+irt)=dotrznew                                      4d4s23
       end do                                                           4d4s23
       xtest=sqrt(xdotnew/dfloat(nroot))                                4d4s23
       bc(ierror+iter-1)=xtest                                          4d4s23
    2  format('after iteration ',i2,' error is ',es10.3)
       bc(ibcoff)=xtest                                                 4d10s23
       bc(ibcoff+1)=ratio                                               4d10s23
       call dws_bcast(bc(ibcoff),2)                                     4d10s23
       xtest=bc(ibcoff)                                                 4d10s23
       ratio=bc(ibcoff+1)                                               4d10s23
       if(xtest.gt.tol.and.ratio.lt.0.9d0)go to 1                         4d10s23
       if(mynowprog.eq.0.and.lprnt)then                                 7d14s22
        if(xtest.le.tol)then                                            4d10s23
         write(6,52)(bc(ierror+i),i=0,iter-1)                            5d19s22
   52    format(10es8.1)                                                 5d27s22
         write(6,*)('cg iterations converged after '),iter,
     $      (' iterations and '),macit-1,(' restarts')                  4d10s23
        else                                                            4d10s23
         write(6,*)('cg iterations stagnating after '),iter,            4d10s23
     $      (' iterations and '),macit-1,(' restarts')                  4d10s23
        end if                                                          4d10s23
      end if                                                             3d3s17
      jx=ix-1                                                           4d4s23
      do irt=1,nroot                                                     3d22s23
       jbest2=ibest2-1+nps*(irt-1)
       do i=1,nps                                                       4d4s23
        dv(i,irt)=bc(jbest2+i)                                          4d10s23
       end do                                                           3d22s23
       jx=jx+nps                                                        4d4s23
      end do                                                            3d22s23
      end if                                                            4d13s23
      if(lprnt)then                                                     7d14s22
       write(6,*)('solution: ')
       call prntm2(dv,nps,nroot,nps)                                    7d20s22
      end if                                                            7d14s22
      ibcoff=ix                                                         4d4s23
      return
      end
