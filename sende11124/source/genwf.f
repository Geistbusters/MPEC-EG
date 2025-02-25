c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genwf(iwavedat,ixw1,ixw2,ncsf,nec,nspc,ixmt,opdata,    5d12s21
     $     iopdata,iosym,nop,opname,nbasdws,nsymb,mdon,mdoop,           5d12s21
     $     ism,irel,irefo,norb,multh,idoubo,isadd,nvirt,maxbx,maxbxd,   7d27s21
     $     srh,sr2,npadddi,lprint,genwfeps,bc,ibc)                      11d10s22
      implicit real*8 (a-h,o-z)                                         5d12s21
      external second
      integer*1 ipack1(4)                                               5d12s21
      integer*4 ipack4(2)                                               5d12s21
      integer*8 ipack8                                                  5d12s21
      character*10 opname(*)                                            5d12s21
      character*6 wname                                                 5d13s21
      logical lri,lprint                                                3d2s22
      equivalence (ipack1,npack4)                                       5d12s21
      equivalence (ipack4,ipack8)                                       5d12s21
      dimension iwavedat(nspc,*),ncsf(*),ixmt(8,*),opdata(*),           5d12s21
     $     iopdata(7,*),iosym(*),nbasdws(*),ism(*),irel(*),irefo(*),    5d27s21
     $     multh(8,8),idoubo(*),nvirt(*),ixmtf(8),idum(2)               12d21s22
      data idum/2*1/                                                    12d21s22
      data dum/1d0/                                                     12d21s22
      include "common.store"                                            5d12s21
      npack4=iwavedat(6,1)                                              5d12s21
      ibcoffo=ibcoff                                                    5d26s22
      do i=1,6                                                          5d13s21
       if(iwavedat(13+i,1).ne.0)then                                    5d13s21
        wname(i:i)=char(iwavedat(13+i,1))                                5d13s21
        nnn=i                                                           5d13s21
       end if                                                           5d13s21
      end do                                                            5d13s21
      if(lprint)                                                        3d2s22
     $     write(6,446)iwavedat(1,1),(wname(i:i),i=1,nnn)               12d19s22
  446 format(' genwf for ',i1,20a1)                                     12d19s22
      nlzz=ipack1(2)                                                    5d12s21
      if(ipack1(2).eq.2)then                                            5d12s21
       npass=1                                                          5d12s21
      else                                                              5d12s21
       npass=ipack1(3)*2                                                5d12s21
      end if                                                            5d12s21
      if(npass.eq.1)then                                                5s26s21
       if(lprint)write(6,*)('linear system, lambda = '),ipack1(3)       3d2s22
       ll=ipack1(3)                                                     5d13s21
      else                                                              5d13s21
       if(lprint)write(6,*)('spherical system, L = '),ipack1(3)         3d2s22
       ll=ipack1(3)                                                     5d13s21
       ipack1(4)=0                                                      5d24s21
       iwavedat(6,1)=npack4                                             5d24s21
       isadd=2*ll                                                       5d24s21
       if(mod(ll,2).eq.0)then                                           5d13s21
        naturs=0                                                        5d13s21
       else
        naturs=1                                                        5d13s21
       end if                                                           5d13s21
      end if                                                            5d13s21
      do i=1,nop                                                        5d12s21
       if(opname(i)(1:3).eq.'lx/')llxx=i                                5d12s21
       if(opname(i)(1:3).eq.'ly/')llyy=i                                5d12s21
       if(opname(i)(1:3).eq.'lz/')llzz=i                                5d12s21
      end do                                                            5d12s21
      nroot=iwavedat(3,1)                                               5d12s21
      ilook=iwavedat(4,1)+iwavedat(13,1)                                5d12s21
      ipack8=ibc(ilook)                                                 5d12s21
      if(npass.eq.1)then                                                5d14s21
c
c     linear system: generate other lz component.
c
       itmp=ibcoff                                                      7d27s21
       ibcoff=itmp+nroot*nroot                                          7d27s21
       call enough('genwf.  1',bc,ibc)
       call psioppsi(iwavedat(1,2),nroot,llzz,1d0,ixmt,iosym,0,idum,dum,7d27s21
     $      iwavedat,nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,ism,irel,irefo,7d27s21
     $      norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb, 7d27s21
     $      iloob,ncsf,1,npadddi,bc,ibc)                                11d10s22
c
c     we have generated uwave(1,2)=Lz*wave(1,1), where uwave is
c     unnormalized. normalized <wave(1,2)|Lz|wave(1,1)>=Lambda.
c     thus wave(1,2)=uwave(1,2)/Lambda.
c
c
c     lz*psi=i*lambda*psi=i*lambda*(re+i*im)=-lambda*im+i*lambda*re.
c     <psi|lz|psi>=i*lambda
c     =<re|lz|re>+i<re|lz|im>-i<im|lz|re>+<im|lz|im>.
c     thus <re|lz|re>+<im|lz|im>=0 and
c     <re|lz|im>-<im|lz|re>=lambda
c     <re|lz|im>=+i*lambda, so if k is im, use + sign.
c     re is b1,a1 = 2,1
c     we compute matrix elements of lzu=lz/i, or lz=i*lzu,
c     Thus <re|lz|im>=<re|i*lzu|im>, and <re|lzu|im>=lambda.            8d21s21
c     Now consider <im|lz|re>:
       if(iwavedat(2,1).eq.1.or.iwavedat(2,1).eq.2.or.                  3d10s22
     $      iwavedat(2,1).eq.5.or.iwavedat(2,1).eq.6)then               3d10s22
        xnorm=-1d0/dfloat(ll)                                           3d10s22
       else                                                             8d21s21
        xnorm=+1d0/dfloat(ll)                                           3d10s22
       end if                                                           8d21s21
       call wnorm(iwavedat(1,2),mdon,mdoop,nvirt,multh,nsymb,xnorm,bc,  11d10s22
     $      ibc)                                                        11d10s22
       if(iwavedat(2,1).eq.1.or.iwavedat(2,1).eq.2.or.
     $      iwavedat(2,1).eq.5.or.iwavedat(2,1).eq.6)then               3d10s22
        ipack1(4)=-ipack1(3)                                            5d28s21
        iwavedat(6,2)=npack4                                            5d28s21
       else                                                             8d21s21
        do iswap=1,nspc                                                 8d21s21
         icpy=iwavedat(iswap,1)                                         8d21s21
         iwavedat(iswap,1)=iwavedat(iswap,2)                            8d21s21
         iwavedat(iswap,2)=icpy                                         8d21s21
        end do                                                          8d21s21
        write(6,*)('swapping, so ipack1(4) of first fcn is '),ipack1(4)
        iwavedat(6,1)=npack4                                            8d21s21
        ipack1(4)=-ipack1(3)                                            5d28s21
        write(6,*)('and ipack1(4) of second fcn is '),ipack1(4)
        iwavedat(6,2)=npack4                                            5d28s21
       end if                                                           8d21s21
       do ic=1,2
        npack4=iwavedat(6,ic)
       end do
       call dws_synca                                                   7d19s21
c
c     let us test expectation value of lz...
c
       xnan=-2d0
       nll=min(2,2*ll)                                                   5d28s21
       idot=ibcoff                                                       5d17s21
       nroot2=nroot*nroot                                               3d9s22
       ibcoff=idot+nroot2                                               3d9s22
       call enough('genwf.  2',bc,ibc)
       ip=1                                                              5d28s21
       do mlq=1,nll                                                      5d28s21
        psr=1d0                                                         3d9s22
        psi=1d0                                                         3d9s22
        if(mod(ll,2).eq.0)then                                          3d9s22
         if(mlq.eq.1)then                                               3d9s22
          ml=-ll                                                        3d9s22
         else                                                           3d9s22
          ml=ll                                                         3d9s22
          psi=-1d0                                                      3d9s22
         end if                                                         3d9s22
        else                                                            3d9s22
         psi=-1d0                                                       3d10s22
         if(mlq.eq.1)then                                               3d9s22
          ml=-ll                                                        3d9s22
          psr=-psr                                                      3d9s22
         else                                                           3d9s22
          ml=ll                                                         3d9s22
         end if                                                         3d9s22
        end if                                                          3d9s22
        if(lprint)write(6,*)('for ml = '),ml
        ipp=ip+1                                                         5d17s21
        npack4=iwavedat(6,ip)                                            5d17s21
        mlp=ipack1(4)                                                    5d24s21
        npack4=iwavedat(6,ipp)                                           5d17s21
        mlpp=ipack1(4)                                                   5d24s21
c
c     (real-i*imag)*(i*lzu)*(real+i*imag)=                              5d17s21
c     (real-i*imag)*(-lzu*imag+i*lzu*real)=                             5d17s21
c     -real*lzu*imag+imag*lzu*real+i*(imag*lzu*imag+real*lzu*real)      5d17s21
c
        if(mlp.gt.0)then                                                 5d17s21
         nzr=ip                                                          5d17s21
         nzi=ipp                                                         5d17s21
        else                                                             5d17s21
         nzr=ipp                                                         5d17s21
         nzi=ip                                                          5d17s21
        end if                                                           5d17s21
        itestr=multh(iwavedat(2,nzr),iosym(llzz))                        5d17s21
        if(itestr.ne.iwavedat(2,nzi))then
         write(6,*)('symmetry test: '),itestr,iwavedat(2,nzi)
         write(6,*)('failed!!')
         call dws_synca
         call dws_finalize
         stop 'genwf'                                                   11d9s22
        end if                                                           5d17s21
        itmp=ibcoff                                                     8d21s21
        ibcoff=itmp+nroot*nroot                                         8d21s21
        call enough('genwf.  3',bc,ibc)
        call psioppsi(iwavedat(1,nzr),nroot,llzz,1d0,ixmt,iosym,0,idum, 8d21s21
     $       isum,iwavedat(1,nzi),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,  8d21s21
     $       ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh, 8d21s21
     $       sr2,multh,nsymb,idum2,ncsf,0,npadddi,bc,ibc)               12d21s22
        pp=psr*psi                                                      3d9s22
        do ir=0,nroot2-1                                                3d9s22
         bc(idot+ir)=-pp*0.5d0*bc(itmp+ir)                              3d9s22
        end do                                                          3d9s22
        call psioppsi(iwavedat(1,nzi),nroot,llzz,1d0,ixmt,iosym,0,idum, 8d21s21
     $       isum,iwavedat(1,nzr),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,  8d21s21
     $       ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh, 8d21s21
     $       sr2,multh,nsymb,idum2,ncsf,0,npadddi,bc,ibc)               12d21s22
        do ir=0,nroot2-1                                                3d9s22
         bc(idot+ir)=bc(idot+ir)+0.5d0*pp*bc(itmp+ir)                   3d9s22
        end do                                                          3d9s22
        ibcoff=itmp                                                     8d21s21
        if(lprint)then                                                  3d2s22
         write(6,*)('lz=-real*lzu*imag+imag*lzu*real: ')                  5d17s21
         call prntm2(bc(idot),nroot,nroot,nroot)                        3d9s22
        end if                                                          3d2s22
        xl=dfloat(ml)                                                   5d28s21
        nok=0                                                           5d28s21
        do i=0,nroot-1                                                  5d28s21
         ii=i*(nroot+1)                                                 3d9s22
         if(abs(bc(idot+ii)-xl).gt.1d-5)then                            3d9s22
          if(lprint)write(6,*)('bad lz!!')                              3d2s22
         else                                                           5d28s21
          nok=nok+1                                                     5d28s21
         end if                                                         5d28s21
         bc(idot+ii)=0d0                                                3d9s22
        end do                                                          5d28s21
        sz=0d0                                                          3d9s22
        do i=0,nroot2-1                                                 3d9s22
         sz=sz+bc(idot+i)**2                                            3d9s22
        end do                                                          3d9s22
        sz=sqrt(sz/dfloat(nroot2))                                      3d9s22
        xnok=dfloat(nok)                                                5d28s21
        call dws_bcast(xnok,1)                                          5d28s21
        call dws_bcast(sz,1)                                            3d9s22
        nok=nint(xnok)                                                  5d28s21
        if(nok.ne.nroot.or.sz.gt.genwfeps)then                          11d9s22
         write(6,*)('sz = '),sz,(' is too large in genwf! ')
         call dws_synca                                                 5d28s21
         call dws_finalize                                              5d28s21
         stop 'genwf'                                                   11d9s22
        end if                                                          5d28s21
        ibcoff=itmp                                                      5d17s21
       end do                                                            5d17s21
       ibcoff=idot                                                      3d9s22
       ibcoff=ibcoffo                                                   5d26s22
       return                                                           5d14s21
      end if                                                            5d14s21
c
c     -1 raising operator
c     +1 lowering operator
c
      psy=-1d0                                                          5d17s21
c                                                                       5d12s21
c     we have lz=0, now use ladder operators.                           5d14s21
c     I phase so lz is real.
c     l+=lx+i*ly. Since I return lux=lx/i and luy=ly/i,                 5d14s21
c     l+=lux*i-luy
c                                                                       5d12s21
      nzu=1                                                             5d14s21
      iloo0=iwavedat(4,nzu)+iwavedat(13,nzu)                            5d14s21
      ipack8=ibc(iloo0)                                                 5d12s21
      if(iwavedat(2,nzu).eq.1.or.iwavedat(2,nzu).eq.4)then              5d17s21
       nacts=0                                                          5d17s21
      else                                                              5d17s21
       nacts=1                                                          5d17s21
      end if                                                            5d17s21
      if(mod(nacts+naturs,2).eq.0)then                                  5d17s21
       lri=.false.                                                      5d17s21
       itestr=multh(iwavedat(2,nzu),iosym(llyy))                         5d14s21
       llru=llyy
       lliu=llxx                                                        5d17s21
       ssru=psy                                                         5d17s21
       ssiu=1d0                                                         5d17s21
      else                                                              5d17s21
       lri=.true.                                                       5d17s21
       itestr=multh(iwavedat(2,nzu),iosym(llxx))                        5d17s21
c     l+=(lux*i-luy)*i*f0=(-lux-i*luy)*f0                               5d17s21
c
       llru=llxx                                                        5d17s21
       lliu=llyy                                                        5d17s21
       ssru=-1d0                                                        5d17s21
       ssiu=psy                                                         5d17s21
      end if                                                            5d17s21
      if(itestr.eq.iwavedat(2,nzu+1))then                               5d14s21
       nzr=nzu+1                                                        5d14s21
       nzi=nzu+2                                                        5d14s21
      else                                                              5d14s21
       nzr=nzu+2                                                        5d14s21
       nzi=nzu+1                                                        5d14s21
      end if                                                            5d14s21
      ipack1(4)=0                                                       5d24s21
      iwavedat(6,1)=npack4                                              5d14s21
      ipack1(4)=1                                                       5d24s21
      iwavedat(6,nzr)=npack4                                            5d14s21
      ipack1(4)=-1                                                      5d24s21
      iwavedat(6,nzi)=npack4                                            5d14s21
      itmp=ibcoff                                                       8d23s21
      ibcoff=itmp+nroot*nroot                                           8d23s21
      call enough('genwf.  4',bc,ibc)
      call psioppsi(iwavedat(1,nzr),nroot,llru,ssru,ixmt,iosym,0,idum,   8d23s21
     $     dum,iwavedat,nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,ism,irel,   8d23s21
     $     irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,multh,  8d23s21
     $     nsymb,idum2,ncsf,1,npadddi,bc,ibc)                            11d10s22
      xnorm=1d0/sqrt(dfloat(ll*(ll+1)))                                 8d23s21
      call wnorm(iwavedat(1,nzr),mdon,mdoop,nvirt,multh,nsymb,xnorm,bc, 11d10s22
     $     ibc)                                                         11d10s22
      call psioppsi(iwavedat(1,nzi),nroot,lliu,ssiu,ixmt,iosym,0,idum,  8d23s21
     $     dum,iwavedat,nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,ism,irel,   8d23s21
     $     irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,multh,  8d23s21
     $     nsymb,idum2,ncsf,1,npadddi,bc,ibc)                            11d10s22
      call wnorm(iwavedat(1,nzi),mdon,mdoop,nvirt,multh,nsymb,xnorm,bc, 11d10s22
     $     ibc)                                                         11d10s22
      ibcoff=itmp                                                       8d23s21
      do ml=2,ll                                                        5d14s21
c
c     now we have complex functions to deal with:                       5d14s21
c     psim+=l+psim=(lux*i-luy)*[Re(psim)+i*Im(psim)]                   5d14s21
c          =-luy*Re(psim)-lux*Im(psim)+i*[-luy*Im(psim)+lux*Re(psim)]
c
       nzq=max(nzr,nzi)+1                                               5d14s21
       nzqq=max(nzr,nzi)+2                                              5d14s21
       itest1=multh(iwavedat(2,nzr),iosym(llyy))                        5d14s21
       itest2=multh(iwavedat(2,nzi),iosym(llxx))                        5d14s21
       if(itest1.ne.itest2)then                                         5d14s21
        write(6,*)('oh no!!! itest1 = '),itest1,(' ne itest2 = '),      5d14s21
     $      itest2                                                      5d14s21
        call dws_synca                                                  5d14s21
        call dws_finalize                                               5d14s21
        stop 'genwf'                                                    11d9s22
       end if                                                           5d14s21
       if(itest1.eq.iwavedat(2,nzq))then                                5d14s21
        nzrp=nzq                                                        5d14s21
        nzip=nzqq                                                       5d14s21
       else                                                             5d14s21
        nzrp=nzqq                                                        5d14s21
        nzip=nzq                                                        5d14s21
       end if                                                           5d14s21
       ipack1(4)=ml                                                     5d24s21
       iwavedat(6,nzrp)=npack4                                          5d14s21
       ipack1(4)=-ml                                                    5d24s21
       iwavedat(6,nzip)=npack4                                          5d14s21
       itmp=ibcoff                                                      8d24s21
       ibcoff=itmp+nroot*nroot                                          8d24s21
       call enough('genwf.  5',bc,ibc)
       call psioppsi(iwavedat(1,nzrp),nroot,llyy,-1d0,ixmt,iosym,0,idum,8d23s21
     $     dum,iwavedat(1,nzr),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,ism, 8d23s21
     $     irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,   8d23s21
     $     multh,nsymb,idum2,ncsf,1,npadddi,bc,ibc)                      11d10s22
       call psioppsi(iwavedat(1,nzrp),nroot,llxx,-1d0,ixmt,iosym,0,idum,8d23s21
     $     dum,iwavedat(1,nzi),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,ism, 8d23s21
     $     irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,   8d23s21
     $     multh,nsymb,idum2,ncsf,1,npadddi,bc,ibc)                      11d10s22
       xnorm=1d0/sqrt(dfloat(ll*(ll+1)-ml*(ml-1)))                      8d23s21
       call wnorm(iwavedat(1,nzrp),mdon,mdoop,nvirt,multh,nsymb,xnorm,  11d10s22
     $      bc,ibc)                                                     11d10s22
       call psioppsi(iwavedat(1,nzip),nroot,llyy,-1d0,ixmt,iosym,0,idum,8d23s21
     $     dum,iwavedat(1,nzi),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,ism, 8d23s21
     $     irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,   8d23s21
     $     multh,nsymb,idum2,ncsf,1,npadddi,bc,ibc)                      11d10s22
       call psioppsi(iwavedat(1,nzip),nroot,llxx,+1d0,ixmt,iosym,0,idum,8d23s21
     $     dum,iwavedat(1,nzr),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,ism, 8d23s21
     $     irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,   8d23s21
     $     multh,nsymb,idum2,ncsf,1,npadddi,bc,ibc)                      11d10s22
       call wnorm(iwavedat(1,nzip),mdon,mdoop,nvirt,multh,nsymb,xnorm,  11d10s22
     $      bc,ibc)                                                     11d10s22
       ibcoff=itmp                                                      8d24s21
       nzr=nzrp                                                         5d14s21
       nzi=nzip                                                         5d14s21
      end do                                                            5d14s21
c
c     let us test expectation value of lz...
c
      nll=2*ll                                                          5d14s21
      ip=2                                                              5d17s21
      idot=ibcoff                                                       5d17s21
      ibcoff=idot+nroot*nroot                                           2d26s22
      call enough('genwf.  6',bc,ibc)
      do ml=1,ll                                                        5d17s21
       ipp=ip+1                                                         5d17s21
       npack4=iwavedat(6,ip)                                            5d17s21
       mlp=ipack1(4)                                                    5d24s21
       npack4=iwavedat(6,ipp)                                           5d17s21
       mlpp=ipack1(4)                                                   5d24s21
c
c     (real-i*imag)*(i*lzu)*(real+i*imag)=                              5d17s21
c     (real-i*imag)*(-lzu*imag+i*lzu*real)=                             5d17s21
c     -real*lzu*imag+imag*lzu*real+i*(imag*lzu*imag+real*lzu*real)      5d17s21
c
       if(mlp.gt.0)then                                                 5d17s21
        nzr=ip                                                          5d17s21
        nzi=ipp                                                         5d17s21
       else                                                             5d17s21
        nzr=ipp                                                         5d17s21
        nzi=ip                                                          5d17s21
       end if                                                           5d17s21
       itestr=multh(iwavedat(2,nzr),iosym(llzz))                        5d17s21
       if(itestr.ne.iwavedat(2,nzi))then
        write(6,*)('symmetry test: '),itestr,iwavedat(2,nzi)
        write(6,*)('failed!!')
        call dws_synca
        call dws_finalize
        stop 'genwf'                                                    11d9s22
       end if                                                           5d17s21
       itmp=ibcoff                                                      5d17s21
       ibcoff=itmp+nroot*nroot                                          8d23s21
       call enough('genwf.  7',bc,ibc)
       call psioppsi(iwavedat(1,nzr),nroot,llzz,1d0,ixmt,iosym,0,idum,  8d21s21
     $       isum,iwavedat(1,nzi),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,  8d21s21
     $       ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh, 8d21s21
     $       sr2,multh,nsymb,idum2,ncsf,0,npadddi,bc,ibc)               12d21s22
        do ir=0,nroot-1                                                 7d19s21
         do jr=0,nroot-1                                                2d26s22
          iad=itmp+jr+ir*nroot                                          2d26s22
          jad=idot+jr+ir*nroot                                          2d26s22
          bc(jad)=-bc(iad)                                              2d26s22
         end do                                                         2d26s22
        end do                                                          7d19s21
        call psioppsi(iwavedat(1,nzi),nroot,llzz,1d0,ixmt,iosym,0,idum, 8d21s21
     $       isum,iwavedat(1,nzr),nroot,bc(itmp),mdon,mdoop,ixw1,ixw2,  8d21s21
     $       ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh, 8d21s21
     $       sr2,multh,nsymb,idum2,ncsf,0,npadddi,bc,ibc)               12d21s22
        do ir=0,nroot-1                                                 7d19s21
         do jr=0,nroot-1                                                2d26s22
          iad=itmp+jr+ir*nroot                                          2d26s22
          jad=idot+jr+ir*nroot                                          2d26s22
          bc(jad)=bc(jad)+bc(iad)                                       2d26s22
         end do                                                         2d26s22
        end do                                                          7d19s21
        ibcoff=itmp                                                     8d21s21
        if(lprint)then                                                  3d2s22
         write(6,*)('lz=-real*lzu*imag+imag*lzu*real: ')                  5d17s21
         call prntm2(bc(idot),nroot,nroot,nroot)                                  5d17s21
        end if                                                          3d2s22
       ip=ip+2                                                          5d17s21
      end do                                                            5d17s21
c
c     note: we have normalized the complex function to unity rather than5d14s21
c     the individual real and imaginary parts                           5d14s21
c
      if(lprint)                                                        3d2s22
     $     write(6,*)('individual parts normalized to unity '),ll,nll   3d2s22
      do i=1,nll
       ip=i+1                                                           5d14s21
       call wnorm(iwavedat(1,ip),mdon,mdoop,nvirt,multh,nsymb,sr2,bc,   11d10s22
     $      ibc)                                                        11d10s22
      end do
      ibcoff=ibcoffo                                                    5d26s22
      return                                                            5d12s21
      end                                                               5d12s21
