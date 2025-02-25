c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diagx(nbas,x,eig,vec,isym,bc,ibc)                      11d14s22
      implicit real*8 (a-h,o-z)
      integer*8 isym                                                    3d3s10
c
c     diagonalization of matrix x.
c     locate decoupled blocks, identical blocks,
c     and enforce symmetry and no mixing.
c
      include "common.store"
      include "common.hf"
      dimension x(nbas,nbas),eig(nbas),vec(nbas,nbas),isym(nbas)        9d10s07
      data icall/0/
      save
      idwsdeb=0
      thrs=1d-12
      thrsdd=thrs*0.01d0                                                9d1s23
      if(idwsdeb.gt.10)then
       write(6,*)('in diagx '),nbas
       write(6,*)('threshold = '),thrs
      end if
      iarg1=nbas
      if(idwsdeb.gt.100)then                                             3d16s12
       write(6,*)('matrix to be diagonalized ')
       call prntm2(x,iarg1,iarg1,iarg1)
      end if
      ieigsav=ibcoff                                                    8d8s7
      ibcoff=ieigsav+nbas                                               8d8s7
      call enough('diagx.  1',bc,ibc)
      ihit=ibcoff
      ibcoff=ihit+nbas
      ihit2=ibcoff                                                      8d18s14
      ibcoff=ihit2+nbas                                                 8d18s14
      call enough('diagx.  2',bc,ibc)
      jhit=ihit-1
      jhit2=ihit2-1
      do 1 i=1,nbas
       ibc(jhit+i)=0
    1 continue
      ngrp=1                                                            8d18s14
      is=1                                                              8d18s14
      ibc(ihit)=1                                                       8d18s14
      ibc(ihit2)=1                                                      8d18s14
      nsofar=1                                                          8d18s14
 1908 continue                                                          8d18s14
      nadd=0
      do j=1,nsofar                                                     8d18s14
       do i=is+1,nbas                                                    8d18s14
        if(ibc(jhit+i).eq.0)then                                         8d18s14
         if(abs(x(i,ibc(jhit2+j))).gt.thrs.and.i.ne.ibc(jhit2+j))then   8d18s14
          ibc(jhit+i)=ngrp                                              8d18s14
          ibc(ihit2+nsofar)=i                                           8d18s14
          nsofar=nsofar+1                                               8d18s14
          nadd=nadd+1                                                   8d18s14
         end if                                                         8d18s14
        end if                                                           8d18s14
       end do                                                           8d18s14
      end do                                                            8d18s14
      if(nadd.ne.0)go to 1908                                           8d18s14
      do i=is,nbas                                                      8d18s14
       if(ibc(jhit+i).eq.0)then
        ngrp=ngrp+1                                                     8d18s14
        is=i                                                            8d18s14
        ibc(jhit+is)=ngrp                                               8d18s14
        ibc(ihit2)=is                                                   8d18s14
        nsofar=1                                                        8d18s14
        go to 1908                                                      8d18s14
       end if                                                           8d18s14
      end do                                                            8d18s14
      igroup=ibcoff
      ibcoff=igroup+ngrp
      igroup2=ibcoff
      ibcoff=igroup2+ngrp
      isymg=ibcoff                                                      9d10s07
      isyme=isymg+ngrp                                                  9d10s07
      ibcoff=isyme+ngrp                                                 9d10s07
      call enough('diagx.  3',bc,ibc)
      do 6 i=1,ngrp
       if(idwsdeb.gt.10)write(6,*)('block '),i
       nz=0
       do 7 j=1,nbas
        if(ibc(jhit+j).eq.i)nz=nz+1
    7  continue
       if(idwsdeb.gt.10)write(6,*)('is size '),nz
       ibc(igroup2+i-1)=nz
       ibck=ibcoff
       ibcoff=ibck+nz*(nz+1)
       ibc(igroup+i-1)=ibck
       call enough('diagx.  4',bc,ibc)
       jj=0
       do 8 j=1,nbas
        if(ibc(jhit+j).eq.i)then
         jj=jj+1
         if(idwsdeb.gt.0)write(6,33)j,i
   33    format(2i5)
         kk=0
         do 9 k=1,nbas
          if(ibc(jhit+k).eq.i)then
           bc(ibck+kk+(jj-1)*nz)=x(k,j)
           kk=kk+1
          end if
    9    continue
        end if
    8  continue
       if(idwsdeb.gt.10)then
       write(6,*)('matrix ')
       iarg1=nz
       call prntm2(bc(ibck),iarg1,iarg1,iarg1)
       end if
       do 10 j=1,i-1
        if(ibc(igroup2+j-1).eq.nz)then
         rms=0d0
         jbck=ibc(igroup+j-1)
         do 11 k=1,nz*nz
          rms=rms+(bc(jbck+k-1)-bc(ibck+k-1))**2
   11    continue
         rms=sqrt(rms/dfloat(nz*nz))
         if(rms.lt.thrs)then
          ibc(igroup2+i-1)=-j
         end if
        end if
   10  continue
    6 continue
      do 12 i=1,ngrp
       nz=ibc(igroup2+i-1)
       if(nz.gt.0)then
        ieig=ibcoff
        ivec=ieig+nz
        lwork=nz*8
        iwork=ivec+nz*nz
        jwork=iwork+lwork
        ifail=jwork+nz*5
        ibcoff=ifail+nz
        call enough('diagx.  5',bc,ibc)
        ibck=ibc(igroup+i-1)
        tol=0d0
        if(idwsdeb.gt.10)then
         write(6,*)('diagonalizing ')
         call prntm2(bc(ibck),nz,nz,nz)
        end if
        call dsyevx('V','A','L',nz,bc(ibck),nz,dum,dum,idum,idum,
     $       tol,neig,bc(ieig),bc(ivec),nz,bc(iwork),lwork,
     $       bc(jwork),bc(ifail),info)
        bc(isyme+i-1)=bc(ieig)                                          9d10s07
        call phasv(bc(ivec),nz,nz,nz)
        if(info.ne.0)then
         write(6,*)('on return from dsyevx, info = '),info
         stop
        end if
        if(idwsdeb.gt.10)then
        write(6,*)('eigenvalues ')
        iarg1=1
        iarg2=nz
        call prntm2(bc(ieig),iarg1,iarg2,iarg1)
        write(6,*)('eigen vectors ')
        call prntm2(bc(ivec),iarg2,iarg2,iarg2)
        end if
        do 13 j=1,nz*(nz+1)
         bc(ibck+j-1)=bc(ieig+j-1)
   13   continue
        ibcoff=ieig
       end if
   12 continue
      do i=1,ngrp                                                       1d19s23
       nz=ibc(igroup2+i-1)                                              1d19s23
       if(nz.lt.0)then                                                  1d19s23
        bc(isyme+i-1)=bc(isyme-nz-1)                                    1d19s23
       end if                                                           1d19s23
      end do                                                            1d19s23
      iarg1=ngrp
      call dsortdws(bc(isyme),ibc(isymg),iarg1)                         1d18s23
      do 14 i=1,nbas
       do 15 j=1,nbas
        vec(j,i)=0d0
   15  continue
   14 continue
      do 16 i=1,ngrp
       nz=ibc(igroup2+i-1)
       ipt=i
       if(nz.lt.0)then
        ipt=-nz
        nz=ibc(igroup2-nz-1)
       end if
       istuff=ibcoff
       ibcoff=istuff+nz
       call enough('diagx.  6',bc,ibc)
       jstuff=istuff
       do 17 j=1,nbas
        if(ibc(jhit+j).eq.i)then
         ibc(jstuff)=j
         jstuff=jstuff+1
        end if
   17  continue
       ieig=ibc(igroup+ipt-1)
       ivec=ieig+nz
       do 18 j=1,nz
        isym(ibc(istuff+j-1))=ibc(isymg+i-1)                            9d10s07
        eig(ibc(istuff+j-1))=bc(ieig+j-1)+dfloat(i)*thrsdd              9d1s23
        iadd=ieigsav+ibc(istuff+j-1)-1                                  8d8s7
        bc(iadd)=bc(ieig+j-1)                                           8d8s7
        do 19 k=1,nz
         vec(ibc(istuff+k-1),ibc(istuff+j-1))=bc(ivec+k-1+(j-1)*nz)
   19   continue
   18  continue
   16 continue
      ibcoff=ihit
      isort=ibcoff
      ibcoff=isort+nbas
      call enough('diagx.  7',bc,ibc)
      iarg2=nbas
      call dsortdws(eig,ibc(isort),iarg2)                               1d18s23
      itmp=ibcoff
      ibcoff=itmp+nbas*nbas
      do 20 i=1,nbas
       ii=ibc(isort+i-1)
       iad=ieigsav+ii-1                                                 8d8s7
       eig(i)=bc(iad)                                                   8d8s7
       do 32 j=1,nbas
        bc(itmp+j-1+(i-1)*nbas)=vec(j,ii)
   32  continue
   20 continue
      jtmp=itmp-1
      do i=1,nbas
       do j=1,nbas
        vec(j,i)=bc(jtmp+j)
       end do
       jtmp=jtmp+nbas
      end do
      do i=1,nbas                                                       9d10s07
       ibc(itmp+i-1)=isym(ibc(isort+i-1))                               9d10s07
      end do                                                            9d10s07
      do i=1,nbas                                                       9d10s07
       isym(i)=ibc(itmp+i-1)                                            9d10s07
      end do                                                            9d10s07
      ibcoff=ieigsav                                                    9d10s07
      return
      end
