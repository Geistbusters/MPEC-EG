c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diagy(nbas,x,eig,vec,isym,bc,ibc,lambdaf,nlambda,      9d1s23
     $     drangcut)                                                    9d1s23
      implicit real*8 (a-h,o-z)
      integer*8 isym                                                    4d19s21
c
c     diagonalization of matrix x.                                      4d19s21
c     enforce symmetry and no mixing based on isym.                     4d19s21
c     if lambdaf ne 0, only use the first nlambda(lambda+1) functions   7d3s23
c     in the diagonalizaton, where lambda is the linear molecule        7d3s23
c     orbital symmetry. if nlambda(lambda+1) lt 0, then keep all fcns.  7d3s23
c
      include "common.store"
      include "common.hf"
      dimension x(nbas,nbas),eig(nbas),vec(nbas,nbas),isym(nbas),       7d3s23
     $     nlambda(*)                                                   7d3s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data icall/0/
      save
      idwsdeb=0
      if(idwsdeb.gt.10)then
       write(6,*)('in diagy '),nbas
       write(6,*)('what we have for isym: ')
       do i=1,nbas
        write(6,*)i,isym(i)
       end do
      end if
      iarg1=nbas
      if(idwsdeb.gt.100)then                                             3d16s12
       write(6,*)('matrix to be diagonalized ')
       call prntm2(x,iarg1,iarg1,iarg1)
      end if
      ieigsav=ibcoff                                                    8d8s7
      ibcoff=ieigsav+nbas                                               8d8s7
      ipt=ibcoff                                                        4d23s21
      ibcoff=ipt+nbas                                                   4d23s21
      iunq=ibcoff                                                       4d19s21
      ibcoff=iunq+nbas                                                  4d19s21
      junq=iunq-1                                                       4d19s21
      call enough('diagy.  1',bc,ibc)
      ngrp=0                                                            4d19s21
      do i=1,nbas                                                       4d19s21
       do j=1,ngrp                                                      4d19s21
        if(isym(i).eq.ibc(junq+j))go to 1908                            4d19s21
       end do                                                           4d19s21
       ngrp=ngrp+1                                                      4d19s21
       ibc(junq+ngrp)=isym(i)                                           4d19s21
 1908  continue                                                         4d19s21
      end do                                                            4d19s21
      igroup=ibcoff
      ibcoff=igroup+ngrp
      igroup2=ibcoff
      ibcoff=igroup2+ngrp
      isymg=ibcoff                                                      9d10s07
      isyme=isymg+ngrp                                                  9d10s07
      ibcoff=isyme+ngrp                                                 9d10s07
      call enough('diagy.  2',bc,ibc)
      rms0=0d0                                                          4d19s21
      szx=0d0                                                           4d19s21
      nrms0=0                                                           4d19s21
      rms1=x(nbas,nbas)**2                                              4d23s21
      do j=1,nbas-1                                                     4d19s21
       rms1=rms1+x(j,j)**2                                              4d23s21
       do i=j+1,nbas                                                    4d19s21
        if(isym(i).ne.isym(j))then                                      4d19s21
         if(abs(x(i,j)).gt.szx)then                                     4d19s21
          szx=abs(x(i,j))                                               4d19s21
          iszx=i                                                        4d19s21
          jszx=j                                                        4d19s21
         end if                                                         4d19s21
         rms0=rms0+x(i,j)**2                                            4d19s21
         nrms0=nrms0+1                                                  4d19s21
        end if                                                          4d19s21
       end do                                                           4d19s21
      end do                                                            4d19s21
      if(nrms0.gt.0)then                                                4d19s21
       rms1=sqrt(rms1/dfloat(nbas))                                     4d23s21
       rms0=sqrt(rms0/dfloat(nrms0))                                    4d19s21
       rmsr=rms0/rms1                                                   4d23s21
       if(rmsr.gt.1d-10)then                                            4d20s21
        write(6,*)('rms offdiagonal for diagy '),rmsr,rms1,rms0
        write(6,*)('largest offdiagonal '),szx,('at '),iszx,jszx
        write(6,*)('but go on anyway ')
       end if                                                           4d20s21
      end if                                                            4d19s21
      ioff=0                                                            4d23s21
      do 6 i=1,ngrp
       if(idwsdeb.gt.10)write(6,*)('block '),i
       nz=0
       do 7 j=1,nbas
        if(isym(j).eq.ibc(junq+i))nz=nz+1                               4d19s21
    7  continue
       if(idwsdeb.gt.10)write(6,*)('is size '),nz
       ibc(igroup2+i-1)=nz
       ibck=ibcoff
       ibcoff=ibck+nz*(nz+1)
       ibc(igroup+i-1)=ibck
       call enough('diagy.  3',bc,ibc)
       jj=0
       do 8 j=1,nbas
        if(isym(j).eq.ibc(junq+i))then                                  4d19s21
         ibc(ipt+ioff+jj)=j                                             4d23s21
         jj=jj+1
         if(idwsdeb.gt.0)write(6,33)j,i
   33    format(2i5)
         kk=0
         do 9 k=1,nbas
          if(isym(k).eq.ibc(junq+i))then                                4d19s21
           bc(ibck+kk+(jj-1)*nz)=x(k,j)
           kk=kk+1
          end if
    9    continue
        end if
    8  continue
       ioff=ioff+nz                                                     4d23s21
       iarg1=nz
       if(idwsdeb.gt.10)then
        write(6,*)('matrix ')
        call prntm2(bc(ibck),iarg1,iarg1,iarg1)
        write(6,*)('nlzz = '),nlzz
        write(6,*)('isym = '),ibc(junq+i)                                7d3s23
       end if                                                           7d3s23
       if(lambdaf.eq.2)then                                                7d3s23
        if(mynowprog.eq.0)write(6,*)('linear molecule')                 7d3s23
        lambda=ibc(junq+i)                                              7d3s23
        if(ibc(junq+i).ge.50)then                                       7d3s23
         lambda=lambda-50                                               7d3s23
        end if                                                          7d3s23
        if(mynowprog.eq.0)                                              7d3s23
     $       write(6,*)('this symmetry block has lambda = '),lambda     7d3s23
        lambdap=lambda+1                                                7d3s23
       end if                                                           7d3s23
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
        call enough('diagy.  4',bc,ibc)
        ibck=ibc(igroup+i-1)
        tol=0d0
        if(idwsdeb.gt.10)then
         write(6,*)('diagonalizing ')
         call prntm2(bc(ibck),nz,nz,nz)
        end if
        if(lambdaf.eq.2)then                                            8d21s23
         ibcopy=ibcoff                                                  8d21s23
         ibcoff=ibcopy+nz*nz                                            8d21s23
         call enough('diagy.bcopy',bc,ibc)                              8d21s23
         do ic=0,nz*nz-1                                                8d21s23
          bc(ibcopy+ic)=bc(ibck+ic)                                     8d21s23
         end do                                                         8d21s23
        end if                                                          8d21s23
        knockout=0                                                      8d21s23
 2022   continue                                                        8d21s23
        knockout=knockout+1                                             8d21s23
        call dsyevx('V','A','L',nz,bc(ibck),nz,dum,dum,idum,idum,
     $       tol,neig,bc(ieig),bc(ivec),nz,bc(iwork),lwork,
     $       bc(jwork),bc(ifail),info)
        bc(isyme+i-1)=bc(ieig)                                          9d10s07
        call phasv(bc(ivec),nz,nz,nz)
        if(lambdaf.eq.2)then
         call dws_bcast(bc(ieig),nz)                                    8d21s23
         if(knockout.gt.15)then                                          8d21s23
          call dws_synca                                                8d21s23
          call dws_finalize                                             8d21s23
          stop "knockout"                                               8d21s23
         end if                                                         8d21s23
         do iz=0,nz-1                                                   8d21s23
          if(abs(bc(ieig+iz)).gt.1d-14)then                             8d21s23
           dyrange=bc(ieig+iz)/bc(ieig+nz-1)                            8d21s23
           izuse=iz                                                     8d21s23
           go to 2023                                                   8d21s23
          else                                                          8d21s23
           bc(ieig+iz)=0d0                                              8d21s23
          end if                                                        8d21s23
         end do                                                         8d21s23
 2023    continue                                                       8d21s23
         lambda=ibc(junq+i)                                              7d3s23
         if(ibc(junq+i).ge.50)then                                       7d3s23
          lambda=lambda-50                                               7d3s23
         end if                                                          7d3s23
         lambdap=lambda+1                                               9d1s23
         if(dyrange.lt.drangcut.or.knockout.le.nlambda(lambdap))then    9d5s23
          if(mynowprog.eq.0)then                                        8d21s23
           write(6,*)('for lambda = '),lambda                           8d21s23
           write(6,*)('dynamic range = '),1d0/dyrange,                  8d21s23
     $          (' is toooo large: delete a fcn')                       8d21s23
          end if
          lvec=ivec+nz*izuse                                            8d21s23
          vx=abs(bc(lvec))                                              8d21s23
          ix=0
          do iq=1,nz-1                                                  8d21s23
           if(abs(bc(lvec+iq)).gt.vx)then                               8d21s23
            vx=abs(bc(lvec+iq))                                         8d21s23
            ix=iq                                                       8d21s23
           end if                                                       8d21s23
          end do                                                        8d21s23
          if(mynowprog.eq.0)write(6,*)                                  8d21s23
     $         ('largest component of vector for lowest eigenvalue is'),
     $         vx,('at '),ix,(', so delete this one ')                  8d21s23
          do iq=0,nz-1                                                  8d21s23
           iad=ibcopy+iq+nz*ix                                          8d21s23
           bc(iad)=0d0                                                  8d21s23
           iad=ibcopy+ix+nz*iq                                          8d21s23
           bc(iad)=0d0                                                  8d21s23
          end do                                                        8d21s23
          do iq=0,nz*nz-1                                               8d21s23
           bc(ibck+iq)=bc(ibcopy+iq)                                    8d21s23
          end do                                                        8d21s23
          go to 2022                                                    8d21s23
         end if                                                         8d21s23
         ibcoff=ibcopy                                                  8d21s23
        end if
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
      iarg1=ngrp
      call dsortdws(bc(isyme),ibc(isymg),iarg1)
      do 14 i=1,nbas
       do 15 j=1,nbas
        vec(j,i)=0d0
   15  continue
   14 continue
      ioff=1                                                            4d23s21
      do i=1,ngrp                                                       4d23s21
       nz=ibc(igroup2+i-1)                                              4d23s21
       ieig=ibc(igroup+i-1)                                             4d23s21
       ivec=ieig+nz                                                     4d23s21
       do j=0,nz-1                                                      4d23s21
        jp=j+ioff                                                       4d23s21
        isym(jp)=ibc(junq+i)                                            4d23s21
        eig(jp)=bc(ieig+j)                                              4d23s21
        do k=0,nz-1                                                     4d23s21
         kp=k+ioff                                                      4d23s21
         kk=ibc(ipt+kp-1)                                               4d23s21
         vec(kk,jp)=bc(ivec+k)                                          4d23s21
        end do                                                          4d23s21
        ivec=ivec+nz                                                    4d23s21
       end do                                                           4d23s21
       ioff=ioff+nz                                                     4d23s21
      end do                                                            4d23s21
      ibcoff=iunq                                                       4d19s21
   88 format(10(10i1,x))
      isort=ibcoff
      ibcoff=isort+nbas
      call enough('diagy.  5',bc,ibc)
      iarg2=nbas
      call dsortdws(eig,ibc(isort),iarg2)
      itmp=ibcoff
      ibcoff=itmp+nbas*nbas
      do 20 i=1,nbas
       ii=ibc(isort+i-1)
       do 32 j=1,nbas
        bc(itmp+j-1+(i-1)*nbas)=vec(j,ii)
   32  continue
   20 continue
      jtmp=itmp-1                                                       1d11s23
      do i=1,nbas                                                       1d11s23
       do j=1,nbas                                                      1d11s23
        vec(j,i)=bc(jtmp+j)                                             1d11s23
       end do                                                           1d11s23
       jtmp=jtmp+nbas                                                   1d11s23
      end do                                                            1d11s23
      do i=1,nbas                                                       9d10s07
       ibc(itmp+i-1)=isym(ibc(isort+i-1))                               9d10s07
      end do                                                            9d10s07
      do i=1,nbas                                                       9d10s07
       isym(i)=ibc(itmp+i-1)                                            4d23s21
      end do                                                            9d10s07
      ibcoff=ieigsav                                                    9d10s07
      return
      end
