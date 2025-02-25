c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine orbdercas(idamat,ihessa,ihessb,ihessc,ihessd,idoub,    5d9s22
     $     iacto,noc,ipuse,multh,idwsdeb,soln,morb,idarot,tol,bc,ibc)   11d10s22
      implicit real*8 (a-h,o-z)
c
c     solve for orbital ders given damat and hess.
c
      dimension idamat(*),ihessa(8,*),ihessb(8,*),ihessc(8,*),
     $     ihessd(*),idoub(*),iacto(*),noc(*),multh(8,8),ipt(8,2),      5d9s22
     $     soln(*),morb(*),idarot(*)                                    11d28s22
      include "common.store"
      include "common.hf"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      if(idwsdeb.gt.10)write(6,*)('Hi, my name is orbdercas! '),tol
      ibcoffo=ibcoff
c
c     solve linear equations via conjugate gradient
c
      nun=0
      do isb=1,nsymb
       isk=multh(isb,ipuse)
       if(idwsdeb.gt.10)write(6,*)('for symmetry blocks '),isb,isk
       ipt(isb,1)=nun                                                   4d3s23
       if(min(iacto(isb),idoub(isk)).gt.0)then                          5d9s22
        if(idwsdeb.gt.10)then
         write(6,*)('active-doub part '),idamat(isb),
     $      idamat(isb)+iacto(isb)*idoub(isk)
         call prntm2(bc(idamat(isb)),iacto(isb),idoub(isk),iacto(isb))   5d9s22
         if(nsymb.eq.1)call printa(bc(idamat(isb)),iacto,idoub,1,0,
     $        idoub,0,1,0,bc(ibcoff))
        end if
        nun=nun+iacto(isb)*idoub(isk)
       end if
       ipt(isb,2)=nun                                                   4d3s23
       jamat=idamat(isb)+iacto(isb)*idoub(isk)                          5d9s22
       if(min(noc(isk),nvirt(isb)).gt.0)then                            5d9s22
        if(idwsdeb.gt.10)then
         write(6,*)('occupied-virt part '),jamat,
     $      jamat+noc(isk)*nvirt(isb)
         call prntm2(bc(jamat),nvirt(isb),noc(isk),nvirt(isb))            5d9s22
         if(nsymb.eq.1)call printa(bc(jamat),nvirt,noc,1,0,noc,0,1,0,
     $        bc(ibcoff))
        end if
        nun=nun+nvirt(isb)*noc(isk)                                     4d3s23
       end if                                                           5d9s22
      end do
      nunm=nun-1
      ix=ibcoff
      ir=ix+nun
      iz=ir+nun                                                         4d3s23
      iq=iz+nun                                                         4d3s23
      ip=iq+nun
      igdiag=ip+nun
      ibcoff=igdiag+nun
      ipp=ibcoff
      ibcoff=ipp+nun
      call enough('orbdercas.  1',bc,ibc)
      jgdiag=igdiag                                                     5d9s22
      do isb=1,nsymb                                                    5d9s22
       if(ihessd(isb).gt.0)then                                         5d9s22
        isbo=multh(isb,ipuse)                                           5d9s22
        nhessd=idoub(isbo)*iacto(isb)+noc(isbo)*nvirt(isb)              5d9s22
        do i=0,nhessd-1                                                 5d9s22
         bc(jgdiag+i)=bc(ihessd(isb)+i)                                 5d9s22
        end do                                                          5d9s22
        jgdiag=jgdiag+nhessd                                            5d9s22
       end if                                                           5d9s22
      end do                                                            5d9s22
      if(idwsdeb.gt.10)then                                             12d17s16
       write(6,*)('diagonal elements of hessian: ')
       call prntm2(bc(igdiag),1,nun,1)
      end if                                                            12d17s16
      do i=0,nunm                                                       4d3s23
       if(abs(bc(igdiag+i)).gt.1d-10)then                               11d29s22
        bc(igdiag+i)=1d0/bc(igdiag+i)                                    5d19s22
       else                                                             5d20s22
        bc(igdiag+i)=0d0                                                5d20s22
       end if                                                           5d20s22
      end do                                                            5d19s22
      xdot=0d0                                                          4d3s23
      if(idwsdeb.gt.10)then                                             5d19s22
       write(6,*)('initial guess for solution: ')
       call prntm2(soln,nun,1,nun)                                      4d3s23
      end if                                                            5d19s22
      call mtimesx(bc(iq),soln,ihessa,ihessb,ihessc,nsymb,              4d3s23
     $      multh,ipuse,ipt,idoub,iacto,noc,nvirt,nun,bc,ibc)           4d3s23
      if(idwsdeb.gt.10)then                                             5d19s22
       write(6,*)('product: ')                                           5d17s22
       call prntm2(bc(iq),nun,1,nun)                                    4d3s23
      end if                                                            5d19s22
      scale=0d0                                                         4d3s23
      dotrz=0d0
      do i=0,nunm                                                       4d3s23
       bc(ir+i)=bc(idamat(1)+i)-bc(iq+i)                                4d3s23
       bc(iz+i)=bc(ir+i)*bc(igdiag+i)                                   4d3s23
       xdot=xdot+bc(ir+i)**2                                            4d3s23
       dotrz=dotrz+bc(iz+i)*bc(ir+i)                                    4d3s23
       scale=scale+bc(idamat(1)+i)**2                                   4d3s23
      end do                                                            5d9s22
      scale=1d0/scale                                                   4d3s23
      xtest=sqrt(xdot*scale)                                            4d4s23
      if(idwsdeb.gt.10)                                                 5d19s22
     $     write(6,*)('rms size of residual '),xtest,xdot,1d0/scale                    4d3s23
      macit=0                                                           5d9s22
      itmax=10
      ierror=ibcoff                                                     5d19s22
      ibcoff=ierror+itmax                                               5d19s22
  100 continue                                                          5d9s22
      macit=macit+1                                                     5d9s22
      if(macit.gt.100)then                                               5d9s22
       write(6,*)('too many restarts!!!')
       call dws_synca
       call dws_finalize
       stop 'orbdercas'
      end if
      do i=0,nunm
       bc(ix+i)=0d0
      end do
      do i=0,nunm
       bc(ip+i)=bc(iz+i)                                                4d3s23
      end do
      iter=0
      nrepreive=0                                                       4d19s16
    1 continue
       iter=iter+1
       if(iter.gt.itmax)then
        if(mynowprog.eq.0.and.idwsdeb.ne.0)then                         7d14s22
         write(6,52)(bc(ierror+izz),izz=0,iter-2)                          5d19s22
        end if
        do i=0,nunm                                                     4d3s23
         soln(i+1)=soln(i+1)+bc(ix+i)                                   4d3s23
        end do
        call mtimesx(bc(iq),soln,ihessa,ihessb,ihessc,nsymb,            4d3s23
     $      multh,ipuse,ipt,idoub,iacto,noc,nvirt,nun,bc,ibc)           4d3s23
        xdot=0d0
        dotrz=0d0
        do i=0,nunm                                                     4d3s23
         bc(ir+i)=bc(idamat(1)+i)-bc(iq+i)                              4d3s23
         bc(iz+i)=bc(ir+i)*bc(igdiag+i)                                 4d3s23
         xdot=xdot+bc(ir+i)**2
         dotrz=dotrz+bc(ir+i)*bc(iz+i)                                  4d3s23
        end do
        xtest=sqrt(xdot*scale)                                          4d4s23
        go to 100                                                       5d9s22
       end if
       if(iter.ne.1)then                                                4d3s23
        beta=dotrznew/dotrz                                             4d3s23
        dotrz=dotrznew                                                  4d3s23
        do i=0,nunm                                                     4d3s23
         bc(ip+i)=bc(iz+i)+beta*bc(ip+i)                                4d3s23
        end do                                                          4d3s23
       end if                                                           4d3s23
       call mtimesx(bc(iq),bc(ip),ihessa,ihessb,ihessc,nsymb,multh,
     $      ipuse,ipt,idoub,iacto,noc,nvirt,nun,bc,ibc)                 4d3s23
       alpha=0d0                                                        4d3s23
       do i=0,nunm
        alpha=alpha+bc(ip+i)*bc(iq+i)                                   4d3s23
       end do                                                           4d3s23
       alpha=dotrz/alpha                                                4d3s23
       xdotnew=0d0                                                      4d3s23
       dotrznew=0d0                                                     4d3s23
       do i=0,nunm                                                      4d3s23
        bc(ix+i)=bc(ix+i)+alpha*bc(ip+i)                                4d3s23
        bc(ir+i)=bc(ir+i)-alpha*bc(iq+i)                                4d3s23
        bc(iz+i)=bc(ir+i)*bc(igdiag+i)                                  4d3s23
        xdotnew=xdotnew+bc(ir+i)**2                                     4d3s23
        dotrznew=dotrznew+bc(ir+i)*bc(iz+i)                             4d3s23
       end do                                                           4d3s23
       xtest=sqrt(xdotnew*scale)                                        4d4s23
       bc(ierror+iter-1)=xtest                                          4d3s23
    2  format('after iteration ',i2,' error is ',es10.3)
       if(iter.gt.1)then                                                11d28s22
        rratio=bc(ierror+iter-2)/xtest                                  4d3s23
       else                                                             11d28s22
        rratio=10d0                                                      11d28s22
       end if                                                           11d28s22
       if(xtest.gt.tol.and.abs(rratio-1d0).gt.1d-7)go to 1                1d4s23
       if(mynowprog.eq.0.and.idwsdeb.ne.0)then                          7d14s22
        write(6,52)(bc(ierror+iz),iz=0,iter-1)                          5d19s22
        write(6,*)('cg iterations converged after '),iter,
     $      (' iterations and '),macit-1,(' restarts')
   52   format(10es8.1)                                                 5d19s22
      end if                                                             3d3s17
      do i=0,nunm                                                       4d3s23
       soln(i+1)=soln(i+1)+bc(ix+i)                                     4d3s23
      end do                                                            5d16s22
      do isb=1,nsymb                                                    6d16s22
       isk=multh(isb,ipuse)
       nn=noc(isb)+nvirt(isb)
       nnk=noc(isk)+nvirt(isk)                                          6d16s22
       do iz=0,nn*nnk-1                                                  5d16s22
        bc(idarot(isb)+iz)=0d0                                          5d16s22
        bc(idarot(isk)+iz)=0d0                                          6d16s22
       end do                                                           5d16s22
      end do                                                            6d16s22
      do isb=1,nsymb
       isk=multh(isb,ipuse)
       nn=noc(isb)+nvirt(isb)
       nnk=noc(isk)+nvirt(isk)                                          6d16s22
       ix1=ibcoff
       ix2=ix1+nn*nn
       ix3=ix2+nn*nnk                                                   6d16s22
       ibcoff=ix3+nn*nn
       call enough('orbdercas.  3',bc,ibc)
       do iz=ix1,ibcoff-1
        bc(iz)=0d0
       end do
       if(idwsdeb.gt.10)then
        write(6,*)('orbitals in ob basis ')
        call prntm2(bc(morb(isb)),nn,nn,nn)
        write(6,*)('solution for symmetry '),isb,isk
       end if
       if(min(iacto(isb),idoub(isk)).gt.0)then                          5d9s22
        if(idwsdeb.gt.10)then
         write(6,*)('active-doub part: '),ipt(isb,1)
         call prntm2(soln(1+ipt(isb,1)),iacto(isb),idoub(isk),           5d9s22
     $       iacto(isb))
        end if
        ii=1+ipt(isb,1)
        do id=0,idoub(isk)-1
         do ia=0,iacto(isb)-1
          iap=ia+idoub(isb)
          iao=idarot(isb)+iap+nn*id                                     5d16s22
          ioa=idarot(isk)+id+nnk*iap                                    6d16s22
          bc(ioa)=soln(ii)
          bc(iao)=-soln(ii)
          ii=ii+1
         end do
        end do
       end if                                                           5d9s22
       if(min(noc(isk),nvirt(isb)).gt.0)then                            5d9s22
        if(idwsdeb.gt.10)then
         write(6,*)('occupied-virt part: ')
         call prntm2(soln(1+ipt(isb,2)),nvirt(isb),noc(isk),             5d9s22
     $       nvirt(isb))                                                5d9s22
        end if
        ii=1+ipt(isb,2)
        do io=0,noc(isk)-1                                              5d9s22
         do iv=0,nvirt(isb)-1
          ivp=iv+noc(isb)
          ivo=idarot(isb)+ivp+nn*io                                     5d16s22
          iov=idarot(isk)+io+nnk*ivp                                    6d16s22
          bc(iov)=soln(ii)
          bc(ivo)=-soln(ii)
          ii=ii+1
         end do
        end do
       end if                                                           5d9s22
       if(idwsdeb.gt.10)then
        write(6,*)('ix1=idarot matrix ')
        call prntm2(bc(idarot(isb)),nn,nnk,nn)
        if(isk.ne.isb)then
         write(6,*)('for isk = '),isk,('rather than isb ')
         call prntm2(bc(idarot(isk)),nnk,nn,nnk)                        6d16s22
        end if                                                          6d16s22
       end if                                                           5d19s22
       ibcoff=ix1
      end do                                                            5d9s22
      ibcoff=ibcoffo                                                    5d16s22
      return
      end
