c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine psisopsi(iwbra,nrootb,i2smb,isymc,iwket,nrootk,i2smk,  8d30s21
     $     ih0,iall,xout,mdon,ism,irel,irefo,norb,nvirt,maxbx,           8d30s21
     $     maxbxd,srh,sr2,multh,nsymb,irori,irw0,irw1,irw2,ih0n,nh0,    10d13s21
     $     isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,iifmx,    11d15s21
     $     ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,  12d20s20
     $         isstor,isym,ascale,idorel,iorb,l2e,bc,ibc,n4vso)         3d29s24
c
      implicit real*8 (a-h,o-z)                                         7d27s21
      external second                                                   5d4s22
      integer*4 ipack4(2)                                               7d27s21
      integer*1 ipack1(8)                                               3d14s24
      integer*8 ipack8,i18,i28,i38,i48                                                  7d27s21
      logical lpr                                                       2d17s22
      equivalence (ipack8,ipack4)                                       7d27s21
      equivalence (ipack8,ipack1)                                       3d14s24
      dimension iwbra(*),ih0(*),iwket(*),xout(nrootb,*),ism(*),irel(*), 10d8s21
     $     irefo(*),nvirt(*),multh(8,8),iall(8,8,8,*),isopt(4,4),       11d15s21
     $     ih0n(8,*),i4or(8,8,8,*),ionexr(8,8,8,*),jmatsr(8,8,8,*),     11d15s21
     $     nh0(8),kmatsr(8,8,8,*),kmatsrb(8,8,8,*),i3xr(8,8,8,*),       11d15s21
     $     iifmx(32,4),ntype(4),idv4(2,4,8),nbasp(*),nbaspc(*),         12d20s20
     $     iaddr(36,16),nfcn(36,3),isstor(*)                            12d23s21
      include "common.store"                                            7d27s21
      include "common.print"                                            3d2s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/cpucom/tovr,top(10),tso(11)                                5d4s22
      data loop,loopx/0,13/
      data icall/0/                                                     8d23s21
      save icall                                                        8d23s21
      icall=icall+1                                                     8d23s21
      lpr=icall.eq.-2
      if(l2e.ne.0)then                                                  2d17s22
       itu=l2e                                                          2d17s22
       itt=1                                                            2d17s22
      else                                                              2d17s22
       itu=-1
       idels=i2smb-i2smk                                                3d9s22
       irif=0                                                           3d29s24
       if(iabs(idels).eq.2)then                                         3d9s22
        idels=1                                                         3d9s22
       else if(iabs(idels).eq.4)then                                    3d27s24
        idels=0                                                         3d27s24
        irif=1                                                          3d29s24
       end if                                                           3d9s22
       do it=1,nsopt                                                     10d8s21
        if(isopt(1,it).eq.isymc.and.isopt(3,it).eq.idels.and.           3d29s24
     $       irif*(isopt(2,it)+1-irori).eq.0)then                       3d29s24
         itu=it                                                         3d19s24
        end if                                                          3d19s24
       end do                                                            10d8s21
       itt=itu                                                          2d17s22
      end if                                                            2d17s22
      if(itu.lt.0)then
       do i=1,nrootk                                                    11d3s21
        do j=1,nrootb                                                   11d3s21
         xout(j,i)=0d0                                                  11d3s21
        end do                                                          11d3s21
       end do                                                           11d3s21
       return                                                           11d3s21
      end if
      ibcoffo=ibcoff                                                    8d3s21
      if(iprtr(26).ne.0.and.lpr)then                                            3d2s22
       ixout=ibcoff                                                      11d12s21
       ibcoff=ixout+nrootb*nrootk                                        11d12s21
       call enough('psisopsi.  1',bc,ibc)
       call testsome(iwbra,iwket,i2smb,i2smk,nsymb,irefo,nvirt,         2d8s23
     $     nbasdws,multh,mdon,ism,irel,norb,ih0n(1,itt),nh0,            2d17s22
     $     iall(1,1,1,itu),iifmx(1,itu),ntype(itu),isopt(1,itt),irori,  2d23s22
     $     irw0,irw1,irw2,isymc,bc(ixout),l2e,bc,ibc)                   11d9s22
      end if                                                            3d2s22
      call dws_synca                                                    8d11s22
      call ddi_iaccword0                                                8d11s22
c     1 is spin 2 is sym 6 is packed 2=nlzz,3=lamdaci
      ipack8=iwbra(6)                                                   3d14s24
      ipack8=iwket(6)                                                   3d14s24
      nwiacc=0                                                          8d11s22
      icsf=1                                                            7d27s21
      icsf2=1                                                           7d27s21
      i2sk=iwket(1)-1                                                   8d30s21
      i2sb=iwbra(1)-1                                                   8d30s21
      nff0k=iwket(4)+iwket(9)                                           7d27s21
      jff0k=iwket(4)+iwket(10)                                          7d27s21
      nec=iwket(7)                                                      7d27s21
      mdookp=iwket(21)                                                  8d31s21
      mdoobp=iwbra(21)                                                  8d31s21
      ilook=iwket(4)+iwket(13)                                          7d27s21
      ipack8=ibc(ilook)                                                 5d12s21
      ncsftk=ipack4(1)
      ilook=ilook+1
      imff1=ilook+nrootk*(ncsftk+1)
      if(lpr)write(6,*)('for ket, imff1 = '),imff1,ibc(imff1),loop
      mff1=ibc(imff1)                                                   7d9s21
      inext=imff1+1                                                     7d21s21
      if(mff1.gt.0)then                                                 7d9s21
       ihsdiag=imff1+1                                                  7d9s21
       nff1=ihsdiag+nsymb*mdookp*2                                       7d9s21
       iff1=nff1+nsymb*mdookp                                            7d9s21
       icsfkq=iff1+ibc(imff1)                                             7d9s21
       inext=icsfkq+mdookp-mdon                                            7d21s21
      end if
      mff2=ibc(inext)                                                   7d21s21
      if(lpr)write(6,*)('for ket, mff2 = '),mff2,loop
      ihddiag=1                                                         1d13s23
      nff2=1                                                            1d13s23
      iff2=1                                                            1d13s23
      icsf=1                                                            1d13s23
      icsf2=1                                                           1d13s23
      nfdat=1                                                           1d13s23
      ivdk=1                                                            1d13s23
      if(mff2.gt.0)then                                                 7d27s21
       ihddiag=inext+1                                                  7d21s21
       nff2=ihddiag+nsymb*mdookp*2                                       7d21s21
       iff2=nff2+nsymb*mdookp                                            7d21s21
       icsf=iff2+mff2                                                   7d21s21
       icsf2=icsf+mdookp-mdon                                            7d21s21
       nfdat=1                                                          8d3s21
       ivdk=1                                                           8d3s21
      else if(mff2.lt.0)then                                            8d3s21
       mdoubstore=inext+1
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoub=ipack4(1)                                                  8d12s21
       mdoub=ipack4(2)                                                  8d3s21
       mff2a=iabs(mff2)                                                 8d3s21
       iff2=mdoubstore+1                                                8d3s21
       nff2=iff2+mff2a                                                  8d3s21
       nfdatk=nff2+mdookp*nsymb                                         11d15s21
       ivdk=nfdatk+10*nsymb                                             11d15s21
       ivdknon=ivdk+nrootk*ndoub                                        8d12s21
       icsfk=ivdknon+nrootk*mdoub                                        8d12s21
       icsfk2=icsfk+mdookp-mdon                                         11d15s21
      else                                                              8d3s21
      end if                                                            7d27s21
      nff0b=iwbra(4)+iwbra(9)                                           7d27s21
      jff0b=iwbra(4)+iwbra(10)                                          7d27s21
      nec=iwbra(7)                                                      7d27s21
      iloob=iwbra(4)+iwbra(13)                                          7d27s21
      ipack8=ibc(iloob)                                                 5d12s21
      ncsftb=ipack4(1)
      iloob=iloob+1
      imff1b=iloob+nrootb*(ncsftb+1)
      if(lpr)write(6,*)('for bra, imff1b '),imff1b,ibc(imff1b),loop
      mff1b=ibc(imff1b)                                                   7d9s21
      inext=imff1b+1                                                     7d21s21
      icsfbq=1                                                          10d22s21
      igsdiag=ibcoff                                                    2d9s23
      ihsdiagb=ibcoff                                                   2d9s23
      nff1b=ibcoff                                                      2d9s23
      iff1b=ibcoff                                                      2d9s23
      if(mff1b.gt.0)then                                                 7d9s21
       ihsdiagb=imff1b+1                                                  7d9s21
       nff1b=ihsdiagb+nsymb*mdoobp*2                                       7d9s21
       iff1b=nff1b+nsymb*mdoobp                                            7d9s21
       icsfbq=iff1b+ibc(imff1b)                                            8d18s21
       inext=icsfbq+mdoobp-mdon                                            7d21s21
       igsdiag=ibcoff                                                   8d21s21
       ibcoff=igsdiag+nsymb*mdoobp                                       8d21s21
       call enough('psisopsi.  2',bc,ibc)
       call makegs(ibc(igsdiag),ibc(nff1b),mdon,mdoobp,nsymb,nvirt,       8d21s21
     $      iwbra(2),multh,ibc(icsfbq),nrootk,maxbxq,npadddi,bc,ibc)    11d10s22
      end if
      mff2b=ibc(inext)                                                   7d21s21
      if(lpr)write(6,*)('for bra, mff2b = '),mff2b,loop
      ilout=ibcoff                                                      7d27s21
      ibcoff=ilout+ncsftb*nrootk                                        7d27s21
      nwds=ncsftb*nrootk-1                                              8d18s21
      do i=ilout,ilout+nwds                                             8d18s21
       bc(i)=0d0                                                        7d27s21
      end do                                                            7d27s21
      call enough('psisopsi.  3',bc,ibc)
      ivdb=1                                                            8d3s21
      igdb=1                                                            8d3s21
      ihddiagb=1                                                        1d13s23
      igddiag=1                                                         1d13s23
      nff2b=1                                                           1d17s23
      nfdatb=1                                                          1d17s23
      iff2b=1                                                           2d9s23
      if(mff2b.gt.0)then                                                 7d27s21
       ihddiagb=inext+1                                                  7d21s21
       nff2b=ihddiagb+nsymb*mdoobp*2                                       7d21s21
       iff2b=nff2b+nsymb*mdoobp                                            7d21s21
       icsfbd=iff2b+mff2b                                                   7d21s21
       icsf2=icsfbd+mdoobp-mdon                                         10d22s21
       igddiag=ibcoff                                                   8d21s21
       ibcoff=igddiag+nsymb*mdoobp                                       8d21s21
       call enough('psisopsi.  4',bc,ibc)
       call makegd(ibc(igddiag),ibc(nff2b),mdon,mdoobp,nsymb,nvirt,      8d24s21
     $     iwbra(2),multh,ibc(icsfbd),ibc(icsf2),nrootk,maxbxdq,npadddi,11d15s22
     $      bc,ibc)                                                     11d15s22
       ivdb=1                                                           8d3s21
       igdb=1                                                           8d3s21
      else if(mff2b.lt.0)then
       mdoubstore=inext+1
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoubb=ipack4(1)                                                  8d3s21
       mdoubb=ipack4(2)                                                  8d3s21
       mff2a=iabs(mff2b)                                                 8d3s21
       iff2b=mdoubstore+1                                                8d3s21
       nff2b=iff2b+mff2a                                                  8d3s21
       nfdatb=nff2b+mdoobp*nsymb                                           8d3s21
       ivdb=nfdatb+10*nsymb                                               8d3s21
       ivdbnon=ivdb+ndoubb*nrootb                                        8d21s21
       icsfbq=ivdbnon+nrootb*mdoubb                                     11d18s21
       icsfb2=icsfbq+mdoobp-mdon                                        11d18s21
        igdb=ibcoff                                                      8d3s21
        ibcoff=igdb+mdoubb*nrootk                                         8d3s21
        nwds=mdoubb*nrootk-1                                             8d18s21
        call enough('psisopsi.  5',bc,ibc)
        do iz=igdb,igdb+nwds                                             8d18s21
         bc(iz)=0d0                                                      8d3s21
        end do                                                           8d3s21
      else
      end if                                                            7d27s21
      if(lpr)write(6,*)('calling hccsfbk4 '),ibcoff,ilout,loop
      idorbb=ibcoff                                                     8d31s21
      isorbb=idorbb+norb                                                8d31s21
      idorbk=isorbb+norb                                                8d31s21
      isorbk=idorbk+norb                                                8d31s21
      icode=isorbk+norb                                                 8d31s21
      imap=icode+norb                                                   8d31s21
      ibcoff=imap+norb                                                  8d31s21
      if(lpr)write(6,*)('ibcoff b4 hccsfbk4 '),ibcoff,loop
      call enough('psisopsi.  6',bc,ibc)
      call second(time1)
      call hccsfbk4(bc(ilout),ncsftb,ibc(nff0b),ibc(jff0b),mdoobp,      10d20s21
     $      i2sb,i2smb,bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),mdookp,   10d8s21
     $      i2sk,i2smk,nrootk,mdon,nec,ism,irel,irefo,norb,multh,nsymb, 10d8s21
     $      ih0,i4o,isymc,ibc(idorbb),ibc(isorbb),ibc(idorbk),          10d8s21
     $      ibc(isorbk),ibc(icode),ibc(imap),irori,irw0,irw1,irw2,      10d8s21
     $      ih0n(1,itt),nh0,i4or(1,1,1,itu),iifmx(1,itu),ntype(itu),    2d17s22
     $      isopt(1,itt),l2e,bc,ibc)                                    11d9s22
      call second(time2)                                                5d4s22
      telap=time2-time1-tovr                                            5d4s22
      tso(1)=tso(1)+telap                                               5d4s22
      time1=time2                                                       5d4s22
      if(mff1b.gt.0)then                                                 7d9s21
      if(lpr)write(6,*)('calling hcsibk4'),loop
       call hcsibk4(                                                    10d20s21
     $ ibc(igsdiag),ibc(nff1b),ibc(iff1b),i2sb,i2smb,mdoobp,ibc(icsfbq),10d20s21
     $     ibc(nff0k),ibc(jff0k),bc(ilook),ncsftk,i2sk,i2smk,mdookp,nec,10d20s21
     $     mdon,nsymb,multh,irw1,irw2,nvirt,isymc,irori,isopt(1,itt),   2d23s22
     $     maxbx,ism,irel,irefo,norb,nrootk,ih0n(1,itt),nh0,            2d17s22
     $     ionexr(1,1,1,itu),iifmx(1,itu),ntype(itu),iwbra(2),l2e,      1d26s23
     $     nwiacc,bc,ibc)                                               1d26s23
       call second(time2)                                               5d4s22
       telap=time2-time1-tovr                                           5d4s22
       tso(2)=tso(2)+telap                                              5d4s22
       time1=time2                                                      5d4s22
      if(lpr)write(6,*)('calling hcisbk4'),loop
       call hcisbk4(ibc(ihsdiag),ibc(nff1),ibc(iff1),i2sk,i2smk,mdookp, 10d26s21
     $    ibc(icsfkq),ibc(nff0b),ibc(jff0b),bc(ilout),ncsftb,i2sb,i2smb,10d25s21
     $     mdoobp,nec,mdon,nsymb,multh,irw1,irw2,nvirt,nrootk,          10d25s21
     $      .false.,isymc,irori,isopt(1,itt),ism,irel,irefo,norb,        2d23s22
     $      ih0n(1,itt),nh0,ionexr(1,1,1,itu),iifmx(1,itu),ntype(itu),  2d17s22
     $     iwket(2),l2e,bc,ibc)                                         11d9s22
       call second(time2)                                               5d4s22
       telap=time2-time1-tovr                                           5d4s22
       tso(3)=tso(3)+telap                                              5d4s22
       time1=time2                                                      5d4s22
      if(lpr)write(6,*)('calling hcssbk4'),loop
       call hcssbk4(ibc(igsdiag),ibc(nff1b),ibc(iff1b),i2sb,i2smb,      11d5s21
     $      mdoobp,ibc(icsfbq),iwbra(2),                                11d5s21
     $      ibc(ihsdiag),ibc(nff1),ibc(iff1),i2sk,i2smk,mdookp,         11d5s21
     $      ibc(icsfkq),iwket(2),                                       11d5s21
     $      nec,mdon,nsymb,multh,irw0,irw1,irw2,nvirt,nrootk,           11d5s21
     $      .false.,isymc,irori,isopt(1,itt),ism,irel,irefo,norb,        2d23s22
     $      ih0n(1,itt),nh0,i4or(1,1,1,itu),jmatsr(1,1,1,itu),          2d17s22
     $      kmatsr(1,1,1,itu),iifmx(1,itu),ntype(itu),maxbx,l2e,        1d26s23
     $      nwiacc,bc,ibc)                                              1d26s23
       call second(time2)                                               5d4s22
       telap=time2-time1-tovr                                           5d4s22
       tso(4)=tso(4)+telap                                              5d4s22
       time1=time2                                                      5d4s22
      end if                                                            7d9s21
      if(mff2b.lt.0)then                                                11d15s21
       if(l2e.eq.0)then                                                 2d17s22
        if(lpr)write(6,*)('calling hcdibk4'),loop
        call hcdibk4(ibc(nff2b),ibc(iff2b),ibc(iff2b),ibc(nfdatb),       11d15s21
     $     ibc(icsfbq),ibc(icsfb2),bc(igdb),i2sb,i2smb,mdoobp,iwbra(2), 11d18s21
     $     ibc(nff0k),ibc(jff0k),bc(ilook),ncsftk,i2sk,i2smk,mdookp,    11d16s21
     $     nec,mdon,nsymb,multh,irw1,irw2,nvirt,nrootk,.false.,isymc,    11d15s21
     $     irori,isopt(1,itt),ism,irel,irefo,norb,kmatsrb(1,1,1,itu),   2d23s22
     $     iifmx(1,itu),ntype(itu),srh,bc,ibc)                          11d9s22
        call second(time2)                                               5d4s22
        telap=time2-time1-tovr                                           5d4s22
        tso(5)=tso(5)+telap                                              5d4s22
        time1=time2                                                      5d4s22
        if(lpr)write(6,*)('calling hcidbk4'),loop
        call hcidbk4(ibc(nff0b),ibc(jff0b),bc(ilout),ncsftb,i2sb,i2smb,  11d17s21
     $      mdoobp,ibc(nff2),ibc(iff2),ibc(iff2),ibc(nfdatk),           11d17s21
     $     ibc(icsfk),ibc(icsfk2),bc(ivdknon),i2sk,i2smk,mdookp,        11d17s21
     $      iwket(2),nec,mdon,nsymb,multh,irw1,irw2,nvirt,nrootk,       2d3s22
     $      .false.,isymc,irori,isopt(1,itt),ism,irel,irefo,norb,       2d23s22
     $      kmatsrb(1,1,1,itu),iifmx(1,itu),ntype(itu),srh,bc,ibc)      11d9s22
        call second(time2)                                               5d4s22
        telap=time2-time1-tovr                                           5d4s22
        tso(6)=tso(6)+telap                                              5d4s22
        time1=time2                                                      5d4s22
       end if                                                           2d17s22
       ib4=ibcoff                                                       12d3s21
       iuall=0
      if(lpr)write(6,*)('calling hcddjkbk4'),loop
       call hcddjkbk4(ibc(nff2b),ibc(iff2b),ibc(iff2b),ibc(nfdatb),     11d19s21
     $     ibc(icsfbq),ibc(icsfb2),bc(igdb),i2sb,i2smb,mdoobp,iwbra(2), 11d19s21
     $      ibc(nff2),ibc(iff2),ibc(iff2),ibc(nfdatk),ibc(icsfk),       11d19s21
     $      ibc(icsfk2),bc(ivdknon),i2sk,i2smk,mdookp,iwket(2),         11d19s21
     $      nec,mdon,nsymb,multh,irw0,irw1,irw2,nvirt,nrootk,.false.,   2d3s22
     $      isymc,irori,isopt(1,itt),ism,irel,irefo,norb,               2d23s22
     $      ih0n(1,itt),nh0,i4or(1,1,1,itu),jmatsr(1,1,1,itu),          2d17s22
     $      kmatsr(1,1,1,itu),iifmx(1,itu),ntype(itu),sr2,srh,          11d19s21
     $      iall(1,1,1,itu),idv4,ndv4,iuall,l2e,bc,ibc,n4vso)           2d8s23
       call second(time2)                                               5d4s22
       telap=time2-time1-tovr                                           5d4s22
       tso(7)=tso(7)+telap                                              5d4s22
       time1=time2                                                      5d4s22
       if(iuall.eq.0.and.l2e.eq.0.and.n4vso.ne.0)then                   2d8s23
        if(lpr)write(6,*)('calling vdmo2sobk4'),loop
        call vdmo2sobk4(idv4,iorb,multh,ibc(nfdatb),iwket(2),ntype(itu), 12d27s21
     $      nbaspc,iaddr,nfcn,nrootk,nbasp,srh,idorel,bc,ibc)           11d14s22
        call paraerikbk4(natom,ngaus,ibdat,isym,iapair,ibstor,isstor,    12d20s20
     $      multh,ascale,nbasp,iaddr,nfcn,nbaspc,iorb,nrootk,           12d20s20
     $      ibc(nfdatb),iwbra(2),iwket(2),bc(igdb),idorel,itu,          12d30s21
     $      ntype(itu),srh,bc,ibc)                                      11d14s22
        call second(time2)                                               5d4s22
        telap=time2-time1-tovr                                           5d4s22
        tso(8)=tso(8)+telap                                              5d4s22
        time1=time2                                                      5d4s22
       end if
       ibcoff=ib4                                                       12d3s21
       if(mff1b.gt.0.and.mff2b.ne.0)then                                2d2s22
      if(lpr)write(6,*)('calling hcdsbk4'),loop
        call hcdsbk4(ibc(nff2b),ibc(iff2b),ibc(iff2b),ibc(nfdatb),      12d6s21
     $     ibc(icsfbq),ibc(icsfb2),bc(igdb),i2sb,i2smb,mdoobp,iwbra(2), 12d6s21
     $      ibc(ihsdiag),ibc(nff1),ibc(iff1),i2sk,i2smk,mdookp,         12d6s21
     $    ibc(icsfkq),iwket(2),nec,mdon,nsymb,multh,irw1,irw2,nvirt,    12d6s21
     $      nrootk,irori,isopt(1,itt),ism,irel,irefo,norb,ih0n(1,itt),  2d23s22
     $      nh0,ionexr(1,1,1,itu),i3xr(1,1,1,itu),iifmx(1,itu),         12d6s21
     $      ntype(itu),maxbx,sr2,l2e,bc,ibc)                            11d14s22
        call second(time2)                                               5d4s22
        telap=time2-time1-tovr                                           5d4s22
        tso(9)=tso(9)+telap                                              5d4s22
        time1=time2                                                      5d4s22
      if(lpr)write(6,*)('calling hcsdbk4'),loop
        call hcsdbk4(ibc(igsdiag),ibc(nff1b),ibc(iff1b),i2sb,i2smb,     12d9s21
     $       mdoobp,ibc(icsfbq),iwbra(2),                               12d9s21
     $       ibc(nff2),ibc(iff2),ibc(iff2),ibc(nfdatk),                 12d9s21
     $     ibc(icsfk),ibc(icsfk2),bc(ivdknon),i2sk,i2smk,mdookp,        12d9s21
     $       iwket(2),nec,mdon,nsymb,multh,irw1,irw2,nvirt,             12d9s21
     $      nrootk,irori,isopt(1,itt),ism,irel,irefo,norb,ih0n(1,itt),  2d23s22
     $      nh0,ionexr(1,1,1,itu),i3xr(1,1,1,itu),iifmx(1,itu),         12d6s21
     $      ntype(itu),maxbx,sr2,l2e,nwiacc,bc,ibc)                     1d26s23
        call second(time2)                                               5d4s22
        telap=time2-time1-tovr                                           5d4s22
        tso(10)=tso(10)+telap                                              5d4s22
        time1=time2                                                      5d4s22
       end if                                                           12d6s21
      end if                                                            11d15s21
      if(lpr)write(6,*)('gsum '),ncsftb,nrootk,loop
      call dws_gsumf(bc(ilout),ncsftb*nrootk)                           8d18s21
      call ddi_iaccword2(nwiacc)                                        8d11s22
      if(lpr)write(6,*)('dotvall'),loop
      call dotvall(bc(ilout),ibc(igsdiag),ibc(igddiag),nrootk,          8d24s21
     $     bc(iloob),                                                   8d19s21
     $    ibc(ihsdiagb),ibc(ihddiagb),nrootb,ncsftb,mff1,mff2,
     $     ibc(nff1b),                                                  8d24s21
     $     ibc(nff2b),mdon,mdoobp,iwbra(2),multh,nsymb,ibc(icsfbq),     10d22s21
     $      ibc(icsf2),nvirt,xout,ibc(iff1b),ibc(iff2b),norb,irefo,     8d24s21
     $      ibc(nff0b),ibc(jff0b),ibc(nfdatb),bc(ivdb),bc(igdb),mdoubb,   8d3s21
     $      ndoubb,bc(ivdk),0,dum,i2eop,ixmt,phase2,sr2,                12d3s21
     $     idumb,nbasdws,ioverwrite,1,idv4,ndv4,iwket(2),itu,           2d8s23
     $     ntype(itu),iifmx(1,itu),bc,ibc)                              11d10s22
      if(lpr)write(6,*)('back from dotvall'),ibcoff,loop
       call second(time2)                                               5d4s22
       telap=time2-time1-tovr                                           5d4s22
       tso(11)=tso(11)+telap                                            5d4s22
      if(iprtr(26).ne.0.and.lpr)then                                            3d2s22
       write(6,*)('xout from dotvall ')
       call prntm2(xout,nrootb,nrootk,nrootb)
       rms=0d0
       jxout=ixout-1
       do ik=1,nrootk                                                    11d12s21
        do ib=1,nrootb                                                   11d12s21
         rms=rms+(xout(ib,ik)-bc(jxout+ib))**2                           11d12s21
        end do                                                           11d12s21
        jxout=jxout+nrootb                                               11d12s21
       end do                                                            11d12s21
       rms=sqrt(rms/dfloat(nrootk*nrootb))                               11d12s21
       write(6,*)('rms difference = '),rms
       write(6,*)('icall = '),icall
       if(.not.(rms.lt.1d-10))then                                       11d12s21
        write(6,*)('rms difference is kind of large ...')                11d12s21
        write(6,*)('xout from testsome: ')
        call prntm2(bc(ixout),nrootb,nrootk,nrootb)
        call dws_synca
        call dws_finalize
        stop
       end if                                                            11d12s21
      end if                                                            3d2s22
      if(lpr)write(6,*)('cleanup'),loop
      if(mff2b.gt.0)then                                                8d21s21
       call unmakegd(ibc(igddiag),ibc(ihddiagb),ibc(nff2),mdon,mdoobp,  8d31s21
     $     nsymb,nvirt,iwbra(2),multh,ncsf,ibc(icsf2),nrootk,0,bc,ibc)  11d10s22
      end if                                                            8d18s21
      if(mff1b.gt.0)then                                                8d21s21
       call unmakegs(ibc(igsdiag),ibc(ihsdiagb),ibc(nff1),mdon,mdoobp,  8d31s21
     $      nsymb,nvirt,iwbra(2),multh,ncsf,nrootk,0,bc,ibc)            11d10s22
      end if                                                            8d21s21
      ibcoff=ilout                                                      8d30s21
      if(lpr)then
       write(6,*)('returning'),loop
      end if
      return                                                            7d27s21
      end                                                               7d27s21
      subroutine lookatgb(gb,m,n)
      implicit real*8 (a-h,o-z)
      dimension gb(*)
      do i=1,n*m
       if(abs(gb(i)).gt.1d3)then
        write(6,*)('big gb! '),i,gb(i)
       end if
      end do
      return
      end
