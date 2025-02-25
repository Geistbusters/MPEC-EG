c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine putwf(fname,iwavedat,idorel,                           2d13s22
     $     nsymb,nbasp,idoubo,irefo,nbasdws,nvirt,irel,ism,multh,maxbx, 2d14s22
     $     maxbxd,iapair,hdd,mdond,mdoopd,bc,ibc)                       11d9s22
      implicit real*8 (a-h,o-z)
      character*(*) fname
      parameter (id=100)                                                11d25s19
      integer*8 ipack8                                                  11d23s19
      integer*4 ipack4(2),jpack4                                        10d7s24
      integer*2 jpack2(2)                                               10d7s24
      integer*1 ipack1(4)                                               5d12s21
      equivalence (ipack8,ipack4)                                       11d23s19
      equivalence (npack4,ipack1)                                       5d12s21
      equivalence (jpack4,jpack2)                                       10d7s24
      include "common.store"
      include "common.basis"
      dimension multhx(64),nbasdws(8),nbasc(8),nbasp(8),isymx(3,8),
     $     idata(id),iunc(2),iwavedat(*),iorb(8),nameci(6),             2d14s22
     $     idumx18(18),nvirt(*),multh(8,8),idoubo(*),irefo(*),          2d14s22
     $     irel(*),ism(*),iapair(3,*),hdd(*),idoubx(8),iacto(8),        5d3s23
     $     icanog(8)                                                    5d3s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      inorm=1                                                           8d19s22
      mxopn=0                                                           1d20s23
      mdoox=0                                                           1d20s23
      if(idorel.eq.0)then                                               1d11s20
       ncomp=1                                                          1d11s20
      else                                                              1d11s20
       ncomp=2                                                          1d11s20
      end if                                                            1d11s20
      if(mynowprog.eq.0)then                                            2d3s22
       write(6,*)('in putwf for file '),fname
      end if                                                            2d3s22
      ibcoffo=ibcoff                                                    5d7s21
      if(mynowprog.eq.0)then                                            7d8s21
       open(unit=1,file=fname(2:len(fname)),form='unformatted')         2d13s22
       open(unit=2,file=fname,form='unformatted')                       2d13s22
       nwavrec=1                                                         5d4s21
       read(1)ismult,isymmrci,nroot,norb,nsb,nameci                      12d27s19
       write(2)ismult,isymmrci,nroot,norb,nsb,nameci                      12d27s19
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)idum
       write(2)(idoubo(i),irefo(i),nbasdws(i),i=1,nsymb)                2d14s22
       if(norb.gt.0)then                                                 5d4s21
        nwavrec=nwavrec+1                                                 5d4s21
        read(1)idum                                                      5d4s21
        write(2)(irel(i),ism(i),i=1,norb)                               2d14s22
       end if                                                            5d4s21
       if(mynowprog.eq.0)then                                           2d3s22
        write(6,*)('spin multiplicity '),ismult
        write(6,*)('symmetry of wavefunction '),isymmrci
        write(6,*)('number of roots = '),nroot
       end if                                                           2d3s22
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)nsymb,idorelx,ngaus,natom,nwcont,numeminus,lmax,nbasallp, 8d30s21
     $     multhx,ascale,potdws,idumx18                                 5d25s21
       write(2)nsymb,idorelx,ngaus,natom,nwcont,numeminus,lmax,nbasallp,2d14s22
     $     multhx,ascale,potdws,idumx18                                 5d25s21
       nextradatatag=idumx18(18)+1                                      2d14s22
       nvec=natom-1                                                      5d26s21
       nexpdata=0                                                        5d26s21
       if(nvec.eq.1)then                                                 5d26s21
        nexpdata=1                                                       5d26s21
       else if(nvec.gt.1)then                                            5d26s21
        nexpdata=3*nvec-3                                                5d26s21
       end if                                                            5d26s21
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(nbasdws(isb),isb=1,nsymb)                                 11d20s19
       write(2)(nbasdws(isb),isb=1,nsymb)                               2d14s22
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(nbasc(isb),isb=1,nsymb)                                   11d20s19
       write(2)(nbasc(isb),isb=1,nsymb)                                   11d20s19
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(nbasp(isb),isb=1,nsymb)                                   11d20s19
       write(2)(nbasp(isb),isb=1,nsymb)                                   11d20s19
       nwavrec=nwavrec+1                                                 5d4s21
       nbasallp=nbasp(1)                                                 12d24s19
       do isb=2,nsymb                                                    12d24s19
        nbasallp=nbasallp+nbasp(isb)                                     12d24s19
       end do                                                            12d24s19
       itmp1=ibcoff
       itmp2=itmp1+natom*3
       iextrad=itmp2+natom*3
       ibcoff=iextrad+nextradatatag                                      5d27s21
       call enough('putwf.  1',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)isymx,(bc(itmp1+i),bc(itmp2+i),i=0,natom*3-1),             5d27s21
     $     (bc(iextrad+i),i=0,nextradatatag-1),                         2d14s22
     $      (idoubx(isb),iacto(isb),isb=1,nsymb),nlzzq                   2d14s22
       nexpt=1+2*natom*3+nextradatatag+2*nsymb+1
       write(2)isymx,(bc(itmp1+i),bc(itmp2+i),i=0,natom*3-1),             5d27s21
     $     (bc(iextrad+i),i=0,nextradatatag-1),                         2d14s22
     $      (idoubx(isb),iacto(isb),isb=1,nsymb),nlzzq                  2d14s22
       ibcoff=iextrad                                                    5d27s21
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)dum
       write(2)((iapair(j,i),j=1,3),i=1,natom)                          2d14s22
       ibcoff=itmp1
       nbdat=ngaus*9+nwcont
       ibcoff=itmp1+nbdat
       call enough('putwf.  2',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(bc(itmp1+i),i=0,nbdat-1)
       write(2)(bc(itmp1+i),i=0,nbdat-1)                                2d14s22
       ibcoff=itmp1+nbasallp*6
       call enough('putwf.  3',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(ibc(itmp1+i),i=0,nbasallp*6-1)
       write(2)(ibc(itmp1+i),i=0,nbasallp*6-1)                          2d14s22
       iout=ibcoff                                                      2d14s22
       jout=iout                                                         11d22s19
       do isb=1,nsymb
        need=nbasp(isb)*nbasdws(isb)*ncomp                               1d11s20
        ibcoff=jout+need                                                 11d22s19
        call enough('putwf.  4',bc,ibc)
        if(need.gt.0)then
         nwavrec=nwavrec+1                                                 5d4s21
         read(1)(bc(jout+i),i=0,need-1)                                  5d3s21
         write(2)(bc(jout+i),i=0,need-1)                                2d14s22
        end if                                                           5d3s21
        jout=ibcoff                                                      11d22s19
        need=nbasdws(isb)*nbasdws(isb)                                   11d25s19
        ibcoff=jout+need                                                 11d25s19
        if(need.gt.0)then
         nwavrec=nwavrec+1                                                 5d4s21
         read(1)(bc(jout+i),i=0,need-1)                                  5d3s21
         write(2)(bc(jout+i),i=0,need-1)                                2d14s22
        end if                                                           5d3s21
        jout=ibcoff                                                      11d25s19
       end do
       jout=iout                                                        2d14s22
       ibcoff=iout                                                      2d14s22
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)nec,jff0,mdon,mdoop,nct                                    5d12s21
       write(2)nec,jff0,mdon,mdoop,nct                                  2d14s22
       nopenx=nec-2*mdon                                                 5d7s21
       mxopn=max(mxopn,nopenx)                                           5d7s21
       mdoo=mdoop-1                                                      5d3s21
       ipack4(1)=mdoo                                                    11d23s19
       ipack4(2)=mdon                                                    11d23s19
       mdoox=max(mdoo,mdoox)                                             11d25s19
       mdoop=mdoo+1                                                      5d3s21
       nhere=0                                                           5d12s21
       if(nhere.gt.id)stop 'id,getwf'
       nff0=jout                                                        2d14s22
       iff0=nff0+3*mdoop                                                2d14s22
       call getbas0(ibc(nff0),mdoop,ibc(iff0),jff0)                       4d30s21
       call dumpbas0(ibc(nff0),mdoop,ibc(iff0),jff0)                    2d14s22
       read(1)iunc,ethred                                                11d23s19
       write(2)iunc,ethred                                                11d23s19
       read(1)nctf,ntot                                                  11d24s19
       write(2)nctf,ntot                                                2d14s22
       read(1)(bc(jout+i),i=0,nctf*nroot-1)                              11d24s19
       nrootw=iwavedat(3)                                               2d14s22
       iveci=iwavedat(4)+iwavedat(13)                                   2d14s22
       ipack8=ibc(iveci)                                                2d14s22
       iveci=iveci+1                                                    2d14s22
       ncsft=ipack4(1)
       write(2)(bc(iveci+i),i=0,ncsft*nrootw-1)                         2d14s22
       imff1=iveci+nrootw*(ncsft+1)                                     2d14s22
       mff1=ibc(imff1)                                                  2d14s22
       inext=imff1+1                                                    2d14s22
       if(mff1.gt.0)then                                                 8d18s21
        ihsdiag=imff1+1                                                  7d9s21
        nff1=ihsdiag+nsymb*mdoop*2                                       7d9s21
        iff1=nff1+nsymb*mdoop                                            7d9s21
        icsf=iff1+ibc(imff1)                                             7d9s21
        inext=icsf+mdoop-mdon                                            7d21s21
       end if                                                           2d14s22
       mff2=ibc(inext)                                                   8d19s21
       if(mff2.gt.0)then                                                 8d19s21
        ihddiagd=inext+1                                                  7d21s21
        nff2d=ihddiagd+nsymb*mdoop*2                                       7d21s21
        iff2d=nff2d+nsymb*mdoop                                            7d21s21
        icsfd=iff2d+mff2                                                   7d21s21
        icsfd2=icsfd+mdoop-mdon                                         2d15s22
       else if(mff2.lt.0)then                                            8d19s21
        mdoubstore=inext+1
        ipack8=ibc(mdoubstore)                                           8d3s21
        ndoub=ipack4(1)                                                  8d12s21
        mdoub=ipack4(2)                                                  8d3s21
        mff2a=iabs(mff2)                                                 8d3s21
        iff2d=mdoubstore+1                                                8d3s21
        nff2d=iff2d+mff2a                                                  8d3s21
        nfdatd=nff2d+mdoop*nsymb                                           8d3s21
        ivdk=nfdatd+10*nsymb                                               8d3s21
        ivdknon=ivdk+nroot*ndoub                                         8d19s21
        icsfd=ivdknon+nroot*mdoub                                       2d14s22
        icsfd2=icsfd+mdoop-mdon                                         2d14s22
        ivdd=ivdk                                                       2d14s22
       end if                                                           2d14s22
       ibcoff=jout                                                       11d24s19
       ieig=ibcoff                                                       3d20s20
       ibcoff=ieig+nroot                                                 3d20s20
       call enough('putwf.  5',bc,ibc)
       read(1)(bc(ieig+i),i=0,nroot-1)                                   3d20s20
       if(mynowprog.eq.0)write(6,*)('diabatic diagonals: '),            2d15s22
     $      (hdd(i),i=1,nroot)                                          2d15s22
       write(2)(hdd(i),i=1,nroot)                                       2d14s22
       read(1)nlzz,lambdaci                                              5d12s21
       jpack4=nlzz                                                      10d7s24
       nder=jpack2(2)                                                   10d7s24
       write(2)nlzz,lambdaci                                              5d12s21
       read(1)mff1                                                      7d8s21
       write(6,*)('mff1 = '),mff1
       if(mff1.gt.0)then                                                7d8s21
        ihsdiag=jout                                                    7d8s21
        nff1=ihsdiag+nsymb*mdoop*2                                      7d9s21
        iff1=nff1+nsymb*mdoop                                           7d8s21
        icsf=iff1+mff1                                                  7d8s21
        jout=icsf+mdoop-mdon                                            7d8s21
        call loadsinga(ibc(nff1),ibc(iff1),mff1,nsymb,mdon,mdoo,        7d8s21
     $       ibc(icsf),1)                                               8d18s22
        call dumpsing(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(icsf),       2d14s22
     $       mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,dum,0,bc,ibc)   11d14s22
        jout=ihsdiag                                                    2d14s22
       else                                                             1d10s24
        write(2)mff1                                                    1d10s24
       end if                                                           7d8s21
       read(1)mff2                                                      7d21s21
       write(6,*)('mff2 '),mff2
       if(mff2.gt.0)then                                                7d8s21
        ihddiag=jout                                                    7d8s21
        nff2=ihddiag+nsymb*mdoop*2                                      7d9s21
        iff2=nff2+nsymb*mdoop                                           7d8s21
        icsf=iff2+mff2                                                  7d8s21
        icsf2=icsf+mdoop-mdon                                            7d8s21
        jout=icsf2+mdoop-mdon                                           7d21s21
        call loaddouba(ibc(nff2),ibc(iff2),mff2,nsymb,mdon,mdoo,        7d8s21
     $       ibc(icsf),ibc(icsf2),inorm)                                8d19s22
        call dumpdoubuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ibc(icsfd),    2d14s22
     $       ibc(icsfd2),mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,0,1,11d14s22
     $       bc,ibc)                                                    11d14s22
        jout=ihddiag                                                    2d14s22
       else if(mff2.lt.0)then                                           8d3s21
        read(1)ndoub,mdoub                                              8d3s21
        mff2a=-mff2                                                     8d3s21
        ipack4(1)=ndoub                                                 8d3s21
        ipack4(2)=mdoub                                                 8d3s21
        mdoubstore=jout                                                 8d3s21
        ibc(jout)=ipack8                                                8d3s21
        iff22=jout+1                                                    8d3s21
        nff22=iff22+mff2a                                               8d3s21
        nfdat=nff22+mdoop*nsymb                                         8d3s21
        ivd=nfdat+10*nsymb                                              8d3s21
        icsf=ivd+nroot*(mdoub+ndoub)                                    8d12s21
        icsf2=icsf+mdoop-mdon                                           8d3s21
        jout=icsf2+2*(mdoop-mdon)                                       8d3s21
        nall=ndoub*nroot                                                8d3s21
        iff2dm=iff2d-1                                                  2d14s22
        call dumpdoubc(ibc(nff2d),ibc(nfdatd),bc(ivdd),nsymb,ndoub,     2d14s22
     $       mdoub,nroot,mdon,mdoo,ibc(icsfd),ibc(icsfd2),iff2dm,bc,ibc,10d7s24
     $       nder)                                                      10d7s24
        call loaddoubc(ibc(iff22),ibc(nff22),ibc(nfdat),                8d3s21
     $       ibc(icsf),ibc(icsf2),mdon,mdoo,nsymb,mff2a,bc(ivd),nall,   8d3s21
     $       isymmrci,nroot,multh,nvirt,inorm)                          8d18s22
        jout=mdoubstore                                                 2d14s22
       else if(mff2.eq.0)then                                           11d9s22
        write(2)mff2                                                     11d9s22
       end if                                                           7d8s21
       ibcoff=jout                                                      7d8s21
      end if                                                            7d8s21
      call dws_synca                                                    2d14s22
      isymmrci=iwavedat(2)
      nroot=iwavedat(3)
      nff0=iwavedat(4)+iwavedat(9)                                            8d18s21
      jff0=iwavedat(4)+iwavedat(10)                                           8d18s21
      nec=iwavedat(7)                                                      7d27s21
      iveci=iwavedat(4)+iwavedat(13)                                          7d27s21
      ipack8=ibc(iveci)                                                 5d12s21
      ncsft=ipack4(1)                                                   8d18s21
      iveci=iveci+1                                                     8d18s21
      imff1=iveci+nroot*(ncsft+1)                                       8d18s21
      mff1=ibc(imff1)                                                   7d9s21
      inext=imff1+1                                                     7d21s21
      if(mff1.gt.0)then                                                 8d18s21
       mdood=mdoopd-1                                                   2d14s22
       ihsdiag=imff1+1                                                  7d9s21
       nff1=ihsdiag+nsymb*mdoopd*2                                      2d14s22
       iff1=nff1+nsymb*mdoopd                                           2d14s22
       icsf=iff1+ibc(imff1)                                             7d9s21
       if(mynowprog.eq.0)then                                           2d14s22
        call dumpsing(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(icsf),mdon,   2d14s22
     $      mdoo,nsymb,multh,nvirt,isymmrci,nroot,dum,1,bc,ibc)         11d14s22
       end if                                                           2d14s22
       mddilow=0                                                        1d20s23
       mddihig=0                                                        1d20s23
       call loadsingb(ibc(ihsdiag),ibc(nff1),nsymb,mdond,mdood,           7d8s21
     $     nvirt,multh,isymmrci,ibc(icsf),nroot,mynowprog,maxbx,inorm,  11d10s22
     $      bc,ibc,mddilow,mddihig)                                     1d19s23
       inext=icsf+mdoopd-mdond                                          2d14s22
      end if                                                            8d18s21
      mff2=ibc(inext)                                                   8d19s21
      if(mff2.gt.0)then                                                 8d19s21
       ihddiag=inext+1                                                  7d21s21
       nff2=ihddiag+nsymb*mdoopd*2                                       7d21s21
       iff2=nff2+nsymb*mdoopd                                            7d21s21
       icsf=iff2+mff2                                                   7d21s21
       icsf2=icsf+mdoopd-mdond                                            7d21s21
       mdood=mdoopd-1                                                   2d14s22
       if(mynowprog.eq.0)then                                           2d14s22
        call dumpdoubuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ibc(icsf),      2d14s22
     $      ibc(icsf2),mdond,mdood,nsymb,multh,nvirt,isymmrci,nroot,1,1,11d14s22
     $      bc,ibc)                                                     11d14s22
       end if                                                           2d14s22
       call loaddoubb(ibc(ihddiag),ibc(nff2),nsymb,mdond,mdood,           7d8s21
     $     nvirt,multh,isymmrci,ibc(icsf),ibc(icsf2),nroot,             2d14s22
     $      mynowprog,maxbxd,inorm,bc,ibc)                              11d10s22
      else if(mff2.lt.0)then                                            8d19s21
       mdoubstore=inext+1
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoub=ipack4(1)                                                  8d12s21
       mdoub=ipack4(2)                                                  8d3s21
       mff2a=iabs(mff2)                                                 8d3s21
       iff2=mdoubstore+1                                                8d3s21
       nff2=iff2+mff2a                                                  8d3s21
       nfdat=nff2+mdoopd*nsymb                                          1d19s23
       ivdk=nfdat+10*nsymb                                               8d3s21
       ivdknon=ivdk+nroot*ndoub                                         8d19s21
      end if                                                            8d19s21
      npass=0                                                           5d12s21
      npack4=iwavedat(6)                                                2d14s22
      nlzz=ipack1(2)                                                    2d14s22
      lambdaci=ipack1(3)                                                2d14s22
      if(nlzz.eq.2)then                                                 5d12s21
       npass=1                                                          5d12s21
       if(lambdaci.eq.0)npass=0                                         8d30s21
      else if(nlzz.eq.6)then                                            5d12s21
       npass=2*lambdaci                                                 5d12s21
      end if                                                            5d12s21
      if(npass.gt.0)then                                                5d12s21
       do ipass=1,npass                                                 5d12s21
        iout=ibcoff                                                      5d12s21
        jout=iout                                                        5d12s21
        ipp=ipass+1                                                     5d12s21
        if(mynowprog.eq.0)then                                          7d8s21
         read(1)mysym                                                    5d12s21
         write(2)mysym                                                  2d14s22
         read(1)necz,jff0z,mdonz,mdoopz,nctz                             5d12s21
         write(2)necz,jff0z,mdonz,mdoopz,nctz                           2d14s22
         nff0z=jout                                                      5d12s21
         iff0z=nff0z+3*mdoopd                                           1d19s23
         jout=iff0z+jff0z                                                5d12s21
         call getbas0(ibc(nff0z),mdoopd,ibc(iff0z),jff0z)               1d19s23
         call dumpbas0(ibc(nff0z),mdoopd,ibc(iff0z),jff0z)              1d19s23
        end if                                                          7d8s21
       end do                                                           5d12s21
      end if                                                            5d12s21
      if(mynowprog.eq.0)then                                            7d8s21
       close(unit=1)
       close(unit=2)
       write(6,*)('we are done writing this file!!')
      end if                                                            7d8s21
      return
      end
