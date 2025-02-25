c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine phouse(iham,ndim,elow,ehi,nroot,ieig,ivec,ipassback,   10d17s14
     $     istart,bc,ibc)                                               11d9s22
      implicit real*8 (a-h,o-z)                                         2d27s12
      external second                                                   10d17s14
c
c     parallel householder etc diagonalization of real symmetric matrix 2d27s12
c     stored decimated upper triangle by columns                        2d27s12
c     ipassback=0 means only return eigenvalues and vectors computed    11d8s12
c     on this proc                                                      11d8s12
c     if ipassback=0, istart is the index of the first eigenvalue       11d9s12
c      returned                                                         11d9s12
c     if we hit an error, return nroot as -1.                           1d8s19
c
      integer*8 iarg8,iarg8b,nnu,nhereu,istartu,iendu,i0,ipy,ipz,ip,ix, 8d30s12
     $     jpre,jx,kpre,is,i0s,ise,ie                                   8d30s12
      include "common.store"
      logical tryrac                                                    3d9s19
      parameter (id8=2000)                                              3d6s12
      integer dws_me,dws_np                                             4d18s12
      dimension ibc8(id8*2)                                             4d18s12
      common/ddidwscm/dws_me,dws_np                                     10d23s14
      common/phpcm/time(3)
      data iseed/935315885/                                             5d3s12
      save                                                              5d3s12
      data icall/0/
      ibcoffo=ibcoff
      icall=icall+1
      nthres=20
      if(ndim.lt.nthres)then                                            2d20s20
       itri=ibcoff
       ntri=(ndim*(ndim+1))/2
       ibcoff=itri+ntri
       call enough('phouse.  1',bc,ibc)
       do i=itri,ibcoff-1
        bc(i)=0d0
       end do
       jham=iham
       do i=0,ntri-1
        if(mod(i,dws_np).eq.dws_me)then
         bc(itri+i)=bc(jham)
         jham=jham+1
        end if
       end do
       call dws_gsumf(bc(itri),ntri)
       ibcoff=itri+ndim*ndim                                            2d20s20
       call square(bc(itri),ndim)
       ieig=ibcoff
       ivec=ieig+ndim
       isym=ivec+ndim*ndim
       ibcoff=isym+ndim
       call enough('phouse.  2',bc,ibc)
       call diagx(ndim,bc(itri),bc(ieig),bc(ivec),ibc(isym),bc,ibc)     11d14s22
       if(elow.le.ehi)then                                              2d20s20
        ilow=0                                                          2d20s20
        nroot=0                                                         2d20s20
        do i=0,ndim-1                                                   2d20s20
         if(bc(ieig+i).lt.elow)then                                     2d20s20
          ilow=i+1                                                      2d20s20
         else if(bc(ieig+i).le.ehi)then                                 2d20s20
          nroot=nroot+1                                                 2d20s20
         end if                                                         2d20s20
        end do                                                          2d20s20
        ieig=ieig+ilow                                                  2d20s20
        ivec=ivec+ilow*ndim                                             2d20s20
       else                                                             2d20s20
       end if                                                           2d20s20
       do i=0,nroot-1
        bc(itri+i)=bc(ieig+i)                                           2d20s20
       end do
       ieig=itri
       itri=itri+nroot                                                  2d20s20
       do i=0,nroot*ndim-1                                              2d20s20
        bc(itri+i)=bc(ivec+i)                                           2d20s20
       end do                                                           2d20s20
       ivec=itri                                                        2d20s20
       ibcoff=ivec+nroot*ndim                                           2d20s20
       if(ipassback.eq.0)then
        if(dws_me.ne.0)then
         ibcoff=ibcoffo
         nroot=0
        end if
        istart=1
        return
       end if
       call dws_bcast(bc(ieig),nroot*(1+ndim))                          2d20s20
       return
      end if                                                            2d20s20
      nproc=dws_np                                                      2d28s12
      nowpro=dws_me                                                     2d28s12
      ibc8(1)=1
      ibc8(2)=0
      ibcoffo=ibcoff                                                    2d27s12
      ndimm=ndim-1
      iarg8=ndim                                                        8d30s12
      iarg8b=ndim-1                                                     8d30s12
      nnu=(iarg8*iarg8b)/2                                              8d30s12
      call distit8(nnu,nhereu,istartu)                                  8d6s12
      istoreu=ibcoff
      ibcoff=istoreu+nhereu
      if(nhereu.gt.0)then                                               5d2s18
       do i=0,nhereu-1                                                  3d15s19
        bc(istoreu+i)=0d0                                               3d15s19
       end do                                                           3d15s19
      end if                                                            5d2s18
      idig=ibcoff
      iodig=idig+ndim
      ibcoff=iodig+ndim
      idig2=ibcoff
      ibcoff=idig2+ndim
      iuvec=ibcoff
      ipvec=iuvec+ndim
      ibcoff=ipvec+ndim
      iendu=istartu+nhereu-1
      call enough('phouse.  3',bc,ibc)
c
c     tridiagonalize via householder reflections
c     from Numerical recipies,2nd ed pg 467
c     plus my || and modified for decimated packed storage
c     I also only store the householder vectors and not form
c     Q directly.
c
      xnpi=1d0/dfloat(dws_np)
      xadd=2d0-0.1d0*xnpi
      xadd2=xadd-1d0
      jham=iham-1
      juvec=iuvec-1
      jpvec=ipvec-1
      jodig=iodig-1
      jdig=idig-1
      jdig2=idig2-1
      call second(time1)
      do i=ndim,2,-1
       do j=1,i
        bc(iuvec+j-1)=0d0
       end do
       im=i-1
       if(im.gt.1)then
        iarg8=i                                                         8d30s12
        iarg8b=im                                                       8d30s12
        ix=(iarg8*iarg8b)/2                                             8d30s12
        xn=dfloat(ix-nowpro)*xnpi+xadd
        ln=max(1,int(xn))
        xx=dfloat(ix+im-nowpro)*xnpi+xadd2
        lx=int(xx)
        jpre=nowpro+1-nproc-ix
        scale=0d0
        do l=ln,lx
         j=jpre+l*nproc
         bc(juvec+j)=bc(jham+l)
         scale=scale+abs(bc(jham+l))
        end do
        bc(juvec+i)=scale
        iarg8=i
        call dws_gsumf(bc(iuvec),i)
        scale=bc(juvec+i)
        if(scale.eq.0d0)then
         bc(jodig+i)=bc(juvec+im)
        else
         scalei=1d0/bc(juvec+i)
         h=0d0
         do k=1,im
          bc(juvec+k)=bc(juvec+k)*scalei
          h=h+bc(juvec+k)**2
         end do
         f=bc(juvec+im)
         g=-sign(sqrt(h),f)
         bc(jodig+i)=scale*g
         h=h-f*g
         bc(juvec+im)=f-g
         srhi=1d0/sqrt(h)
         do j=1,im
          bc(jpvec+j)=bc(juvec+j)*srhi
         end do
         iarg8=im                                                       8d30s12
         iarg8b=i-2                                                     8d30s12
         i0=(iarg8*iarg8b)/2                                            8d30s12
         ipy=max(i0+1,istartu)
         ipz=min(i0+im,iendu)
         do ip=ipy,ipz
          iad1=istoreu+ip-istartu
          iad2=jpvec+ip-i0
          bc(iad1)=bc(iad2)
         end do
         do j=1,im
          bc(jpvec+j)=0d0
         end do
         do j=1,im
          iarg8=j                                                       8d30s12
          iarg8b=j-1                                                    8d30s12
          jx=(iarg8*iarg8b)/2                                           8d30s12
          xn=dfloat(jx-nowpro)*xnpi+xadd
          ln=max(1,int(xn))
          xx=dfloat(jx+j-1-nowpro)*xnpi+xadd2
          lx=int(xx)
          kpre=nowpro+1-nproc-jx
          do l=ln,lx
           k=kpre+l*nproc
           bc(jpvec+j)=bc(jpvec+j)+bc(jham+l)*bc(juvec+k)
           bc(jpvec+k)=bc(jpvec+k)+bc(jham+l)*bc(juvec+j)
          end do
          xx=dfloat(jx+j-nowpro)*xnpi+xadd2
          lxn=int(xx)
          if(lxn.ne.lx)then
           k=kpre+lxn*nproc
           bc(jpvec+j)=bc(jpvec+j)+bc(jham+l)*bc(juvec+k)
          end if
         end do
         iarg8=i-1
         call dws_gsumf(bc(ipvec),i-1)
         hi=1d0/h
         xk=0d0
         do j=1,im
          bc(jpvec+j)=bc(jpvec+j)*hi
          xk=xk+bc(jpvec+j)*bc(juvec+j)
         end do
         xk=xk*0.5d0*hi
         do j=1,im
          bc(jpvec+j)=bc(jpvec+j)-xk*bc(juvec+j)
         end do
         do j=1,im
          iarg8=j                                                       8d30s12
          iarg8b=j-1                                                    8d30s12
          jx=(iarg8*iarg8b)/2                                           8d30s12
          xn=dfloat(jx-nowpro)*xnpi+xadd
          ln=max(1,int(xn))
          xx=dfloat(jx+j-nowpro)*xnpi+xadd2
          lx=int(xx)
          kpre=nowpro+1-nproc-jx
          do l=ln,lx
           k=kpre+l*nproc
           bc(jham+l)=bc(jham+l)-bc(juvec+k)*bc(jpvec+j)
     $          -bc(juvec+j)*bc(jpvec+k)
          end do
         end do
        end if
       else
        if(nowpro.eq.1)then                                             5d15s12
         bc(jodig+i)=bc(iham)                                           5d15s12
        else if(nproc.eq.1)then                                         5d15s12
         bc(jodig+i)=bc(iham+1)                                         5d15s12
        else                                                            5d15s12
         bc(jodig+i)=0d0
        end if                                                          5d15s12
        iarg8=1
        call dws_gsumf(bc(jodig+i),1)
       end if
       bc(jdig+i)=h
      end do
      do i=1,ndim
       bc(jdig2+i)=0d0
       im=i-1
       iarg8=i                                                          8d30s12
       iarg8b=i+1                                                       8d30s12
       iarg8=(iarg8*iarg8b)/2                                           8d30s12
       xx=dfloat(iarg8-nowpro)*xnpi+xadd2                               8d30s12
       lx=int(xx)
       rm=abs(xx-dfloat(lx))
       if(rm.lt.xnpi)then
        bc(jdig2+i)=bc(jham+lx)
       end if
      end do
      iarg8=ndim
      call dws_gsumf(bc(idig2),ndim)
      call second(time2)
      time(1)=time2-time1
      iodig2=ibcoff
      ibcoff=iodig2+ndim
      jodig2=iodig2-1
      ind=ibcoff
      irv1=ind+ndim                                                     5d9s12
      ibcoff=irv1+ndim                                                  5d9s12
      call enough('phouse.  4',bc,ibc)
      indi=ind*intmul                                                   5d9s12
      call enough('phouse.  5',bc,ibc)
      nan=0                                                             5d2s18
      do i=1,ndim
       bc(jdig+i)=bc(jdig2+i)
       bc(jodig2+i)=bc(jodig+i)
       bc(jodig+i)=bc(jodig+i)**2
       if(bc(jdig2+i).ne.bc(jdig2+i))nan=1                              5d2s18
       if(i.gt.1)then                                                   12d14s18
        if(bc(jodig+i).ne.bc(jodig+i))nan=1                              5d2s18
       end if                                                           12d14s18
 3352  format(i5,2es15.7)
      end do
      if(nan.ne.0)then                                                  5d2s18
       write(6,*)('in phouse, got at least one nan for imtqlv input')    5d2s18
       do i=1,ndim
        write(6,3352)i,bc(jdig+i),bc(jodig2+i)
       end do
       call dws_sync                                                    5d2s18
       call dws_finalize                                                5d2s18
       stop                                                             5d2s18
      end if                                                            5d2s18
      idothr=ibcoff                                                      3d27s23
      ieothr=idothr+ndim                                                3d27s23
      ibcoff=ieothr+ndim                                                3d27s23
      call enough('phouse.dothr',bc,ibc)                                3d27s23
      jdothr=idothr-1                                                   3d27s23
      jeothr=ieothr-1                                                   3d27s23
      do i=1,ndim                                                       3d27s23
       bc(jdothr+i)=bc(jdig2+i)                                         3d27s23
       bc(jeothr+i)=bc(iodig2+i)                                         3d27s23
      end do                                                            3d27s23
      lwork=ndim*5                                                      3d27s23
      liwork=ndim*5                                                     3d27s23
      ieig=ibcoff                                                       3d27s23
      iwork=ieig+ndim                                                   3d27s23
      iiwork=iwork+lwork                                                3d27s23
      ifail=iiwork+liwork                                               3d27s23
      ibcoff=ifail+ndim                                                 3d27s23
      call enough('phouse.eig',bc,ibc)                                  3d27s23
      call dstevx('N','A',ndim,bc(idothr),bc(ieothr),dum,dum,idum,idum, 3d27s23
     $     1d-14,ngot,bc(ieig),dum,ndim,bc(iwork),ibc(iiwork),          3d27s23
     $     ibc(ifail),info)                                             3d27s23
      do i=0,ndim-1                                                     3d27s23
       bc(idig+i)=bc(ieig+i)                                            3d27s23
      end do                                                            3d27s23
      ibcoff=ieig                                                       3d27s23
      err=0d0                                                           1d9s18
      if(ierr.ne.0)then
       if(dws_me.eq.0)                                                  12d2s22
     $      write(6,*)('on return from dstevx in phouse, ierr = '),ierr 12d2s22
       err=1d0                                                          1d9s18
      end if
      call dws_gsumf(err,1)                                             1d9s18
      if(err.ne.0d0)then                                                1d9s18
       if(dws_me.eq.0)                                                  12d2s22
     $     write(6,*)('at least one proc had an error in phoused2p')    12d2s22
       nroot=-1                                                         1d9s18
       ibcoff=ibcoffo                                                   1d9s18
       return                                                           1d9s18
      end if                                                            1d9s18
      call dws_bcast(bc(idig),ndim)                                     1d24s20
c
c     option 1: generate roots up to eeig
c     option 2: generate roots up to nroot
c     option 3: generate roots between elow and ehi
c
      if(elow.gt.ehi)then                                               10d23s14
      else
       nzero=-1                                                         4d3s12
       nroot=0                                                          10d17s14
       ilow=0                                                           10d17s14
       do i=1,ndim                                                      2d27s12
        if(bc(jdig+i).lt.elow)ilow=i                                    10d17s14
        if(bc(jdig+i).gt.elow.and.bc(jdig+i).le.ehi)nroot=nroot+1       10d17s14
       end do                                                           2d27s12
       if(nroot.eq.0)then                                               7d3s18
        ibcoff=ibcoffo                                                  7d3s18
        return                                                          7d3s18
       end if                                                           7d3s18
       jdig=jdig+ilow                                                   10d17s14
       idig=idig+ilow                                                   10d17s14
      end if                                                            10d23s14
      dgtst=1d-5
      ndgen=0                                                           5d3s12
      isomap=ibcoff                                                     5d3s12
      ibcoff=isomap+nroot                                               5d3s12
      isomap2=ibcoff                                                    5d1s18
      ibcoff=isomap2+nroot                                              5d1s18
      call enough('phouse.  6',bc,ibc)
      isomapi=isomap*intmul                                             5d3s12
      isomap2=isomap2*intmul                                            5d1s18
      jsomapi=isomapi                                                   5d3s12
      is=1                                                              5d3s12
      nunq=0                                                            5d1s18
      maxdg=0                                                           5d3s12
 3322 continue                                                          5d3s12
       nunq=nunq+1                                                      5d1s18
       in=nroot+1                                                       5d3s12
       ibc(isomap2+is-1)=nunq                                           5d1s18
       do i=is+1,nroot                                                  5d3s12
        diff=bc(idig+i-1)-bc(idig+is-1)                                 5d3s12
        if(diff.gt.dgtst)then                                           5d3s12
         in=i                                                           5d3s12
         go to 3323                                                     5d3s12
        else                                                            5d1s18
         ibc(isomap2+i-1)=nunq                                          5d1s18
        end if                                                          5d3s12
        ndgen=ndgen+1                                                   5d3s12
       end do                                                           5d3s12
 3323  continue                                                         5d3s12
       maxdg=max(maxdg,in-is)                                           5d3s12
       is=in                                                            5d3s12
      if(is.le.nroot)go to 3322                                         5d3s12
      nrootx=nunq                                                       5d2s18
      call distit(nrootx,nheret,istart)                                 5d3s12
      if(ndgen.ne.0)then                                                5d3s12
       isw=nroot+1                                                      5d1s18
       nh=0                                                             5d1s18
       do iw=1,nheret                                                   5d1s18
        igoal=istart+iw-1                                               5d1s18
        do i=0,nroot-1                                                  5d1s18
         if(ibc(isomap2+i).eq.igoal)then                                5d1s18
          isw=min(isw,i+1)                                              5d1s18
          nh=nh+1                                                       5d1s18
         end if                                                         5d1s18
        end do                                                          5d1s18
       end do                                                           5d1s18
       istart=isw                                                       5d1s18
       nheret=nh                                                        5d1s18
      end if                                                            5d3s12
      ivect=ibcoff                                                      3d17s17
      ivect=ibcoff                                                      5d2s18
      ibcoff=ivect+ndim*nheret                                          5d2s18
      err=0d0                                                           1d8s19
      if(nheret.gt.0)then
      irv1=ibcoff                                                       5d3s12
      irv2=irv1+ndim                                                    5d3s12
      irv3=irv2+ndim                                                    5d3s12
      irv4=irv3+ndim                                                    5d3s12
      irv6=irv4+ndim                                                    5d3s12
      ibcoff=irv6+ndim                                                  5d3s12
      call enough('phouse.  7',bc,ibc)
      ngot=nheret                                                       5d14s12
c
c     need a wrapper here because ibc is integer*8 but imtql2 and
c     tinvit address ibc(indi) as integer*4
c
      lwork=ndim*5                                                      3d27s23
      liwork=ndim*5                                                     3d27s23
      ieig=ibcoff                                                       3d27s23
      iwork=ieig+ndim                                                   3d27s23
      iiwork=iwork+lwork                                                3d27s23
      ifail=iiwork+liwork                                               3d27s23
      ivecnew=ifail+ndim                                                3d27s23
      ibcoff=ivecnew+ndim*ngot                                          3d27s23
      call enough('phouse.eigb',bc,ibc)                                  3d27s23
      istarte=istart+ngot-1                                             3d27s23
      call dstevx('V','I',ndim,bc(idothr),bc(ieothr),dum,dum,istart,    3d27s23
     $     istarte,1d-14,ngotx,bc(ieig),bc(ivecnew),ndim,bc(iwork),     3d27s23
     $     ibc(iiwork),ibc(ifail),info)                                 3d27s23
      if(ngot.ne.ngotx)err=1d0                                          6d26s23
      do i=0,ndim*ngot-1                                                3d27s23
       bc(ivect+i)=bc(ivecnew+i)                                        3d27s23
      end do                                                            3d27s23
      if(info.ne.0)then                                                 5d3s12
       write(6,*)('on return from dstevx, info = '),info                3d27s23
       err=1d0                                                          1d8s19
      end if                                                            5d3s12
      ibcoff=irv1                                                       5d2s18
      end if                                                            10d17s14
      call dws_gsumf(err,1)                                             1d8s19
      if(err.ne.0d0)then                                                1d8s19
       nroot=-1                                                         1d8s19
       ibcoff=ibcoffo                                                   1d8s19
       return                                                           1d8s19
      end if                                                            1d8s19
      call second(time3)
      time(2)=time3-time2
c
c     back transform to original basis
c
      nhereux=nnu/nproc
      nhu=nhereux
      nleftu=nnu-nhereux*nproc
      if(nleftu.ne.0)nhereux=nhereux+1
      istoreu2=ibcoff
      ibcoff=istoreu2+nhereux
      znorm=0d0                                                         3d15s19
      call enough('phouse.  8',bc,ibc)
      ipass=0
      is=1
      ist=3
  101 continue
      call second(time1)
       if(nowpro.eq.ipass)then
        do i=0,nhereu-1
         bc(istoreu2+i)=bc(istoreu+i)
        end do
       end if
       iarg8=nhereux
       ipass8=ipass
       call dws_sync
       call dws_12all(bc(istoreu2),nhereux,ipass)
       nhereg=nhu
       if(ipass.lt.nleftu)nhereg=nhereg+1
       ise=is+nhereg-1
       do i=ist,ndim
        im=i-1
        iarg8=im                                                        8d30s12
        iarg8b=i-2                                                      8d30s12
        i0=(iarg8*iarg8b)/2                                             8d30s12
        i0s=max(i0+1,is)
        ie=min(i0+im,ise)
        ngh=ie+1-i0s
        if(i0s.eq.i0+1)then
         itmp=ibcoff
         ibcoff=itmp+nheret
         itmp2=ibcoff
         ibcoff=itmp2+im
         call enough('phouse.  9',bc,ibc)
         do j=0,nheret-1
          bc(itmp+j)=0d0
         end do
        end if
        jlow=i0s-i0
        jhigh=ie-i0
        do j=jlow,jhigh
         bc(itmp2+j-1)=bc(istoreu2+i0s-is+j-jlow)
        end do
        iarg2=istoreu2+i0s-is-jlow
        do k=1,nheret
         iarg1=itmp+k-1
         iarg3=ivect-1+ndim*(k-1)
         do j=jlow,jhigh
          bc(iarg1)=bc(iarg1)+bc(iarg2+j)*bc(iarg3+j)
         end do
        end do
        if(ie.ne.i0+i-1)then
         ibg=ie
         ist=i
         go to 102
        else if(jhigh.ge.jlow)then                                      5d2s18
         do k=1,nheret
          iarg1=ivect-1+ndim*(k-1)
          iarg2=itmp2-1
          val=bc(itmp+k-1)
          do j=1,im
           bc(iarg1+j)=bc(iarg1+j)-bc(iarg2+j)*val
          end do
         end do
         ibcoff=itmp
        end if
       end do
  102  continue
       ipass=ipass+1
       is=is+nhereg
      if(ipass.lt.nproc)go to 101
      call phasv(bc(ivect),ndim,ndim,nheret)                            2d27s12
      call dws_sync
      if(ipassback.eq.0)then                                            11d8s12
       ieig=idig+istart-1                                               11d9s12
       ivec=ivect                                                       11d8s12
       nroot=nheret                                                     11d8s12
       ibcoff=ivec+ndim*nroot                                           11d8s12
c
c     move full list of eigenvalues to be right after eigenvectors
c
       do i=1,ndim                                                      6d20s14
        bc(ibcoff+i-1)=bc(idig+i-1)                                     6d20s14
       end do                                                           6d20s14
       ivnew=ieig+nheret                                                12d5s12
       nmove=ndim*(nroot+1)                                             6d20s14
       do i=0,nmove-1                                                   12d5s12
        bc(ivnew+i)=bc(ivec+i)                                          12d5s12
       end do                                                           12d5s12
       ivec=ivnew                                                       12d5s12
       ibcoff=ivec+ndim*nheret                                          12d5s12
       return                                                           11d8s12
      end if                                                            11d8s12
c
c     we have all vectors - pass them about to everyone                 2d27s12
c                                                                       2d27s12
      ieig=ibcoff                                                       2d27s12
      ibcoff=ieig+nroot                                                 6d7s12
      ivall=ibcoff                                                      2d27s12
      ibcoff=ivall+ndim*nroot                                           2d27s12
      nblock=1                                                          3d6s12
      if(nproc.gt.id8)then                                              3d6s12
       write(6,*)('in phouse, nproc = '),nproc                          3d6s12
       write(6,*)('exceeds id8 = '),id8                                 3d6s12
       stop                                                             3d6s12
      end if                                                            3d6s12
      ioff=nblock+nproc                                                 3d6s12
      call enough('phouse. 10',bc,ibc)
      nhere=nrootx/nproc                                                5d3s12
      nleft=nrootx-nhere*nproc                                          5d3s12
      ibc8(ioff)=0                                                      3d6s12
      ist=1                                                             5d3s12
      do ip=0,nproc-1                                                   2d27s12
       nherex=nhere                                                     5d3s12
       if(ip.lt.nleft)nherex=nherex+1                                   5d3s12
       nh=0                                                             5d2s18
       isw=nroot+1                                                      5d2s18
       do iw=1,nherex                                                   5d2s18
        igoal=ist+iw-1                                                  5d2s18
        do i=0,nroot-1                                                  5d2s18
         if(ibc(isomap2+i).eq.igoal)then                                5d2s18
          isw=min(isw,i+1)                                              5d2s18
          nh=nh+1                                                       5d2s18
         end if                                                         5d2s18
        end do                                                          5d2s18
       end do                                                           5d2s18
       ibc8(nblock+ip)=nh                                               5d2s18
       ibc8(nblock+ip)=ibc8(nblock+ip)*ndim                             3d6s12
       ist=ist+nherex                                                   5d3s12
       if(ip.gt.0)then                                                  2d27s12
        ipm=ip-1                                                        2d27s12
        ibc8(ioff+ip)=ibc8(ioff+ipm)+ibc8(nblock+ipm)                   3d6s12
       end if                                                           2d27s12
      end do                                                            2d27s12
      call dws_allgv2(bc(ivall),ibc8(nblock),ibc8(ioff),bc(ivect))      3d6s12
      call second(time4)
      time(3)=time4-time3
      do i=0,nroot-1                                                    2d27s12
       bc(ieig+i)=bc(idig+i)                                            2d27s12
      end do                                                            2d27s12
      nmove=nroot*(ndim+1)                                              2d27s12
      do i=0,nmove-1                                                    2d27s12
       bc(ibcoffo+i)=bc(ieig+i)                                         2d27s12
      end do                                                            2d27s12
      ieig=ibcoffo                                                      2d27s12
      ivec=ieig+nroot                                                   2d27s12
      ibcoff=ivec+nroot*ndim                                            2d27s12
      return                                                            2d27s12
      end
