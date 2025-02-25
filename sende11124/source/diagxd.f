c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diagxd(nbas,x,eig,vec,isym,bc,ibc,nder,thrse,nok,xskip)9d1s23
      implicit real*8 (a-h,o-z)
      integer*8 isym                                                    3d3s10
c
c     diagonalization of matrix x to form orthogonalization matrix      8d31s23
c     plus derivatives.                                                 8d31s23
c     locate decoupled blocks, identical blocks,
c     and enforce symmetry and no mixing.
c
      include "common.store"
      include "common.hf"
      dimension x(nbas,*),eig(nbas),vec(nbas,*),isym(nbas)              8d31s23
      data icall/0/
      save
      xskip=0d0                                                         9d1s23
      idwsdeb=-1000
      nderp=nder+1                                                      8d31s23
      thrs=1d-12
      thrsdd=thrs*0.01d0                                                9d1s23
      if(idwsdeb.gt.10)then
       write(6,*)('in diagxd '),nbas
       write(6,*)('threshold = '),thrs,thrse
      end if
      iarg1=nbas
      if(idwsdeb.gt.100)then                                             3d16s12
       write(6,*)('matrix to be diagonalized ')
       call prntm2(x,iarg1,iarg1*nderp,iarg1)
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
       ibcoff=ibck+nz*(nz+1)*nderp                                      8d31s23
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
           jbck=ibck+kk+(jj-1)*nz                                       8d31s23
           jl=j                                                         8d31s23
           do idr=0,nder                                                8d31s23
            bc(jbck)=x(k,jl)
            jbck=jbck+nz*nz                                             8d31s23
            jl=jl+nbas                                                  8d31s23
           end do                                                       8d31s23
           kk=kk+1
          end if
    9    continue
        end if
    8  continue
       if(idwsdeb.gt.10)then
       write(6,*)('matrix ')
       iarg1=nz
       call prntm2(bc(ibck),iarg1,iarg1*nderp,iarg1)                    8d31s23
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
        ibcpy=ibcoff                                                    8d31s23
        ibcoff=ibcpy+nz*nz                                              8d31s23
        call enough('diagxd.bcpy',bc,ibc)                               8d31s23
        do j=0,nz*nz-1                                                  8d31s23
         bc(ibcpy+j)=bc(ibck+j)                                         8d31s23
        end do                                                          8d31s23
        call dsyevx('V','A','L',nz,bc(ibck),nz,dum,dum,idum,idum,
     $       tol,neig,bc(ieig),bc(ivec),nz,bc(iwork),lwork,
     $       bc(jwork),bc(ifail),info)
        bc(isyme+i-1)=bc(ieig)                                          9d10s07
        call phasv(bc(ivec),nz,nz,nz)                                   6d27s24
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
        iskip=0                                                         8d31s23
        do j=0,nz-1                                                     8d31s23
         if(bc(ieig+j).lt.thrse)then                                    9d1s23
          iskip=iskip+1                                                 9d1s23
          xskip=max(xskip,bc(ieig+j))                                   9d1s23
         end if                                                         9d1s23
        end do                                                          8d31s23
        itmph=ibcoff                                                    8d31s23
        itmpf=itmph+nz*nz                                               8d31s23
        ideig=itmpf+nz*nz                                               8d31s23
        ibcoff=ideig+nz                                                 8d31s23
        ipvt=ibcoff                                                     8d31s23
        ibcoff=ipvt+nz                                                  8d31s23
        call enough('diagxd.tmpd',bc,ibc)                               8d31s23
        idmt=ibck+nz*nz                                                 8d31s23
        do idr=1,nder                                                   8d31s23
         call dgemm('n','n',nz,nz,nz,1d0,bc(idmt),nz,bc(ivec),nz,0d0,   8d31s23
     $        bc(itmph),nz,'diagxd.tmph')                               8d31s23
         do j=0,nz-1                                                    8d31s23
          jvec=ivec+j*nz                                                8d31s23
          jtmph=itmph+j*nz                                              8d31s23
          sum=0d0                                                       8d31s23
          do k=0,nz-1                                                   8d31s23
           sum=sum+bc(jtmph+k)*bc(jvec+k)                               8d31s23
          end do                                                        8d31s23
          bc(ideig+j)=sum                                               8d31s23
         end do                                                         8d31s23
         if(nz.gt.1)then                                                8d31s23
          if(iskip.eq.0)then                                            8d31s23
           do j=0,nz-1                                                  8d31s23
            do k=0,nz-1                                                 8d31s23
             jk=itmph+j+nz*k                                            8d31s23
             kj=itmpf+k+nz*j                                            8d31s23
             bc(kj)=bc(jk)                                              8d31s23
            end do                                                      8d31s23
           end do                                                       8d31s23
           call dgemm('n','n',nz,nz,nz,1d0,bc(itmpf),nz,bc(ivec),nz,0d0,8d31s23
     $          bc(itmph),nz,'diagxd.tmpf')                             8d31s23
           do j=0,nz-1                                                  8d31s23
            jtmph=itmph+nz*j                                            8d31s23
            do k=0,nz-1                                                 8d31s23
             bot=bc(ieig+k)-bc(ieig+j)                                  8d31s23
             if(abs(bot).lt.thrs)then                                   8d31s23
              if(j.ne.k)
     $            write(6,*)('small denominator! '),bot,k,j,bc(jtmph+k)
              bc(jtmph+k)=0d0                                           8d31s23
             else                                                       8d31s23
              bc(jtmph+k)=-bc(jtmph+k)/bot                              8d31s23
             end if                                                     8d31s23
            end do                                                      8d31s23
           end do                                                       8d31s23
           call dgemm('n','n',nz,nz,nz,1d0,bc(ivec),nz,bc(itmph),nz,0d0,8d31s23
     $          bc(itmpf),nz,'diagxd.dv')                               8d31s23
           icpy=itmpf                                                   8d31s23
           itmpf=itmph                                                  8d31s23
           itmph=icpy                                                   8d31s23
          else                                                          8d31s23
           nkeep=nz-iskip                                               8d31s23
           do j=0,nkeep-1                                               8d31s23
            jp=j+iskip                                                  8d31s23
            ff=bc(ideig+jp)                                             8d31s23
            jtmph=itmph+nz*jp                                           8d31s23
            jvec=ivec+nz*jp                                             8d31s23
            do k=0,nz-1                                                 8d31s23
             bc(jtmph+k)=bc(jvec+k)*ff-bc(jtmph+k)                      8d31s23
            end do                                                      8d31s23
           end do                                                       8d31s23
           jtmph=itmph+nz*iskip
           do lhs=0,nkeep-1                                             8d31s23
            lhsp=lhs+iskip                                              8d31s23
            xx=bc(ieig+lhsp)                                            8d31s23
            jtmph=itmph+nz*lhsp                                         8d31s23
            jvec=ivec+nz*lhsp                                           8d31s23
            do j=0,nz-1                                                 8d31s23
             jtmpf=itmpf+nz*j                                           8d31s23
             jbcpy=ibcpy+nz*j
             do k=0,nz-1                                                8d31s23
              bc(jtmpf+k)=bc(jbcpy+k)+bc(jvec+k)*bc(jvec+j)*xx          8d31s23
             end do                                                     8d31s23
             bc(jtmpf+j)=bc(jtmpf+j)-xx                                 8d31s23
            end do                                                      8d31s23
            call lusolv(bc(itmpf),nz,nz,bc(jtmph),nz,1,ibc(ipvt),ierr,3)8d31s23
           end do                                                       8d31s23
          end if                                                        8d31s23
         else                                                           8d31s23
          bc(itmph)=0d0                                                 8d31s23
         end if                                                         8d31s23
         jtmph=itmph+nz*iskip
         nkeep=nz-iskip                                                 8d31s23
         do j=0,nkeep-1                                                 8d31s23
          jp=j+iskip                                                    8d31s23
          xx=1d0/sqrt(bc(ieig+jp))                                      8d31s23
          dxx=-0.5d0*bc(ideig+jp)*(xx**3)                               8d31s23
          jtmph=itmph+nz*jp                                             8d31s23
          jvec=ivec+nz*jp                                               8d31s23
          jtmpf=itmpf+nz*jp                                             8d31s23
          do k=0,nz-1                                                   9d1s23
           bc(jtmpf+k)=xx*bc(jvec+k)                                    8d31s23
           bc(jtmph+k)=xx*bc(jtmph+k)+dxx*bc(jvec+k)                    8d31s23
          end do                                                        8d31s23
         end do                                                         8d31s23
         do j=0,iskip-1                                                 9d1s23
          jtmpf=itmpf+nz*j                                              9d1s23
          jtmph=itmph+nz*j                                              9d1s23
          do k=0,nz-1                                                   9d1s23
           bc(jtmpf+k)=0d0                                              9d1s23
           bc(jtmph+k)=0d0                                              9d1s23
          end do                                                        9d1s23
         end do                                                         9d1s23
         jtmpf=itmpf+nz*iskip
         jtmph=itmph+nz*iskip
         do j=0,nz*nz-1                                                 8d31s23
          bc(idmt+j)=bc(itmph+j)                                        8d31s23
         end do                                                         8d31s23
         idmt=idmt+nz*nz                                                8d31s23
        end do                                                          8d31s23
        do j=0,nz*nz-1                                                  8d31s23
         bc(ibck+j)=bc(itmpf+j)                                         8d31s23
        end do                                                          8d31s23
        jbck=ibck+nz                                                    9d1s23
        nmvt=nz*nz*nderp                                                9d1s23
        do j=nmvt-1,0,-1                                                9d1s23
         bc(jbck+j)=bc(ibck+j)                                          9d1s23
        end do                                                          9d1s23
        do j=0,nz-1                                                     8d31s23
         bc(ibck+j)=bc(ieig+j)                                          8d31s23
        end do                                                          8d31s23
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
      do 14 i=1,nbas*nderp                                              8d31s23
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
        ibcol=ibc(istuff+j-1)                                           8d31s23
        ilcol=j-1                                                       8d31s23
        do idr=0,nder                                                   9d1s23
         do 19 k=1,nz                                                   8d31s23
          vec(ibc(istuff+k-1),ibcol)=bc(ivec+k-1+ilcol*nz)              8d31s23
   19    continue                                                       8d31s23
         ibcol=ibcol+nbas                                               8d31s23
         ilcol=ilcol+nz                                                 8d31s23
        end do                                                          8d31s23
   18  continue
   16 continue
      ibcoff=ihit
      isort=ibcoff
      ibcoff=isort+nbas
      call enough('diagx.  7',bc,ibc)
      iarg2=nbas
      call dsortdws(eig,ibc(isort),iarg2)                               1d18s23
      itmp=ibcoff
      ibcoff=itmp+nbas*nbas*nderp                                       8d31s23
      do 20 i=1,nbas
       ii=ibc(isort+i-1)
       iad=ieigsav+ii-1                                                 8d8s7
       eig(i)=bc(iad)                                                   8d8s7
       icol=i-1                                                         8d31s23
       do idr=1,nderp                                                   8d31s23
        do 32 j=1,nbas
         bc(itmp+j-1+icol*nbas)=vec(j,ii)                               8d31s23
   32   continue
        icol=icol+nbas                                                  8d31s23
        ii=ii+nbas                                                      8d31s23
       end do                                                           8d31s23
   20 continue
      jtmp=itmp-1
      do i=1,nbas*nderp                                                 8d31s23
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
      iskip=0                                                           8d31s23
      do i=1,nbas                                                       8d31s23
       if(eig(i).lt.thrse)iskip=iskip+1                                 8d31s23
      end do                                                            8d31s23
      if(iskip.gt.0)then
       nok=nbas-iskip                                                   8d31s23
       do i=1,nok                                                       8d31s23
        ip=i+iskip                                                      8d31s23
        eig(i)=eig(ip)                                                  8d31s23
        isym(i)=isym(ip)                                                8d31s23
       end do                                                           8d31s23
       i1off=iskip+1
       i2off=1                                                          8d31s23
       do idr=0,nder                                                    8d31s23
        do i=1,nok                                                      8d31s23
         do j=1,nbas                                                    8d31s23
          vec(j,i2off)=vec(j,i1off)                                     8d31s23
         end do                                                         8d31s23
         i2off=i2off+1                                                  8d31s23
         i1off=i1off+1                                                  8d31s23
        end do                                                          8d31s23
        i1off=i1off+iskip                                               8d31s23
       end do                                                           8d31s23
      end if
      return
      end
