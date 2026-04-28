c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine psipsiden(iwbra,nrootbk,iden,iwket,mdon,mdoop,         3d10s25
     $     ixw1,ixw2,ism,irel,                                          3d10s25
     $     irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,multh,  7d27s21
     $     nsymb,ncsf,bc,ibc)                                           3d10s25
c
c     if ioverwrite is nonzero, overwrite ivb with igb.                 8d18s21
      implicit real*8 (a-h,o-z)                                         7d27s21
      external second                                                   5d4s22
      integer*4 ipack4(2)                                               7d27s21
      integer*8 ipack8                                                  7d27s21
      logical lpr,lpr2                                                  1d5s23
      equivalence (ipack8,ipack4)                                       7d27s21
      dimension iwbra(*),iden(8),iwket(*),ism(*),irel(*),irefo(*),      3d10s25
     $     nbasdws(*),idoubo(*),nvirt(*),multh(8,8),ncsf(*)             3d10s25
      include "common.store"                                            7d27s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data icall/0/                                                     8d23s21
      common/cpucom/tovr,top(10),tso(11)                                5d4s22
      save icall                                                        8d23s21
      integer*8 ibcoffo,igdb,igddiag,igsdiag,ihddiagb,ilout
      write(6,*)('Hi, my name is psipsiden '),loc(bc),loc(ibc)
      icall=icall+1                                                     8d23s21
      lpr=icall.eq.-1
      lpr2=icall.eq.-1
      if(lpr)then                                                       8d25s21
       pp=phase1*2d0
       write(6,*)('testme would have been next ...')
       if(bc(132).ne.-132d0)then
        call dws_synca
        call dws_finalize
        stop 'psipsiden'
       end if
      end if                                                            8d25s21
      ibcoffo=ibcoff                                                    8d3s21
      icsf=1                                                            7d27s21
      icsf2=1                                                           7d27s21
      nff0k=iwket(4)+iwket(9)                                           7d27s21
      jff0k=iwket(4)+iwket(10)                                          7d27s21
      nec=iwket(7)                                                      7d27s21
      ilook=iwket(4)+iwket(13)                                          7d27s21
      ipack8=ibc(ilook)                                                 5d12s21
      ncsftk=ipack4(1)
      ilook=ilook+1
      imff1=ilook+nrootbk*(ncsftk+1)                                    3d10s25
      if(lpr)write(6,*)('for ket, imff1 = '),imff1
      mff1=ibc(imff1)                                                   7d9s21
      inext=imff1+1                                                     7d21s21
      ihsdiag=1                                                         1d17s23
      nff1=1                                                            1d17s23
      iff1=1                                                            1d17s23
      if(mff1.gt.0)then                                                 7d9s21
       ihsdiag=imff1+1                                                  7d9s21
       nff1=ihsdiag+nsymb*mdoop*2                                       7d9s21
       iff1=nff1+nsymb*mdoop                                            7d9s21
       icsf=iff1+ibc(imff1)                                             7d9s21
       inext=icsf+mdoop-mdon                                            7d21s21
      end if
      mff2=ibc(inext)                                                   7d21s21
      ivdknon=1                                                         1d17s23
      if(lpr)write(6,*)('for ket, mff2 = '),mff2
       ihddiag=1                                                        8d3s21
       nff2=1                                                           8d3s21
       iff2=1                                                           8d3s21
       icsf=1                                                           8d3s21
       icsf2=1                                                          8d3s21
       nfdat=1                                                          8d3s21
       ivdk=1                                                           8d3s21
      if(mff2.gt.0)then                                                 7d27s21
       ihddiag=inext+1                                                  7d21s21
       nff2=ihddiag+nsymb*mdoop*2                                       7d21s21
       iff2=nff2+nsymb*mdoop                                            7d21s21
       icsf=iff2+mff2                                                   7d21s21
       icsf2=icsf+mdoop-mdon                                            7d21s21
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
       nfdat=nff2+mdoop*nsymb                                           8d3s21
       ivdk=nfdat+10*nsymb                                               8d3s21
       ivdknon=ivdk+nrootbk*ndoub                                       3d10s25
       icsf=ivdknon+nrootbk*mdoub                                       3d10s25
       icsf2=icsf+mdoop-mdon                                            8d3s21
      else                                                              8d3s21
      end if                                                            7d27s21
      nff0b=iwbra(4)+iwbra(9)                                           7d27s21
      jff0b=iwbra(4)+iwbra(10)                                          7d27s21
      nec=iwbra(7)                                                      7d27s21
      iloob=iwbra(4)+iwbra(13)                                          7d27s21
      ipack8=ibc(iloob)                                                 5d12s21
      ncsftb=ipack4(1)
      iloob=iloob+1
      imff1b=iloob+nrootbk*(ncsftb+1)                                   3d10s25
      if(lpr)write(6,*)('for bra, imff1b '),imff1b,nrootbk,ncsftb
      mff1b=ibc(imff1b)                                                   7d9s21
      inext=imff1b+1                                                     7d21s21
      ihsdiagb=1                                                        1d17s23
      nff1b=1                                                           1d17s23
      iff1b=1                                                           1d17s23
      igsdiag=1                                                         1d17s23
      if(mff1b.gt.0)then                                                 7d9s21
       ihsdiagb=imff1b+1                                                  7d9s21
       nff1b=ihsdiagb+nsymb*mdoop*2                                       7d9s21
       iff1b=nff1b+nsymb*mdoop                                            7d9s21
       icsfq=iff1b+ibc(imff1b)                                            8d18s21
       inext=icsfq+mdoop-mdon                                            7d21s21
       igsdiag=ibcoff                                                   8d21s21
       ibcoff=igsdiag+nsymb*mdoop                                       8d21s21
       call enough('psipsiden.  1',bc,ibc)
       call makegs(ibc(igsdiag),ibc(nff1b),mdon,mdoop,nsymb,nvirt,       8d21s21
     $      iwbra(2),multh,ncsf,nrootbk,maxbxq,npadddi,bc,ibc)           11d10s22
      end if
      mff2b=ibc(inext)                                                   7d21s21
      if(lpr)write(6,*)('for bra, mff2b = '),mff2b
       ilout=ibcoff                                                      7d27s21
       ibcoff=ilout+ncsftb*nrootbk                                      3d10s25
       nwds=ncsftb*nrootbk-1                                            3d10s25
       do i=ilout,ilout+nwds                                             8d18s21
        bc(i)=0d0                                                        7d27s21
       end do                                                            7d27s21
       call enough('psipsiden.  2',bc,ibc)
       ivdb=1                                                           8d3s21
       igdb=1                                                           8d3s21
       igddiag=ibcoff                                                    1d13s23
       ihddiagb=ibcoff                                                  1d13s23
       nff2b=1                                                          1d17s23
       iff2b=1                                                          1d17s23
      nfdatb=1                                                          1d17s23
      if(mff2b.gt.0)then                                                 7d27s21
       ihddiagb=inext+1                                                  7d21s21
       nff2b=ihddiagb+nsymb*mdoop*2                                       7d21s21
       iff2b=nff2b+nsymb*mdoop                                            7d21s21
       icsf=iff2b+mff2b                                                   7d21s21
       icsf2=icsf+mdoop-mdon                                            7d21s21
       igddiag=ibcoff                                                   8d21s21
       ibcoff=igddiag+nsymb*mdoop                                       8d21s21
       call enough('psipsiden.  3',bc,ibc)
       call makegd(ibc(igddiag),ibc(nff2b),mdon,mdoop,nsymb,nvirt,      8d24s21
     $      iwbra(2),multh,ncsf,ibc(icsf2),nrootbk,maxbxdq,npadddi,bc,   11d15s22
     $      ibc)                                                        11d15s22
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
       nfdatb=nff2b+mdoop*nsymb                                           8d3s21
      ivdb=nfdatb+10*nsymb                                               8d3s21
       ivdbnon=ivdb+ndoubb*nrootbk                                      3d10s25
        igdb=ibcoff                                                      8d3s21
        ibcoff=igdb+mdoubb*nrootbk                                      3d10s25
        nwds=mdoubb*nrootbk-1                                           3d10s25
        call enough('psipsiden.  4',bc,ibc)
        do iz=igdb,igdb+nwds                                             8d18s21
         bc(iz)=0d0                                                      8d3s21
        end do                                                           8d3s21
      else
      end if                                                            7d27s21
      write(6,*)('hccsfbkden is next '),iloob,ncsftb,ilook,ncsftk,nff0b,
     $     nff0k,jff0b,jff0k,nrootbk
      call hccsfbkden(bc(iloob),ncsftb,bc(ilook),ncsftk,ibc(nff0b),     3d10s25
     $     ibc(nff0k),ncsf,ibc(jff0b),ibc(jff0k),nrootbk,mdon,mdoop,    3d10s25
     $     ixw1,ixw2,nec,ism,irel,irefo,norb,iden,multh,nbasdws,idoubo, 3d10s25
     $     nsymb,bc,ibc)                                                3d10s25
      igoal=iden(1)+8-1+28*(3-1)
      write(6,*)('goal after hccsfbkden '),igoal,bc(igoal)
      write(6,*)('back from hccsfbkden')
      return                                                            7d27s21
      end                                                               7d27s21
