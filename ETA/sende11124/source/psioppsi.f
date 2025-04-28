c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine psioppsi(iwbra,nrootb,lxmt,phase1,ixmt,iosym,n2e,      7d27s21
     $     i2eop,phase2,iwket,nrootk,xout,mdon,mdoop,ixw1,ixw2,ism,irel,7d27s21
     $     irefo,norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,multh,  7d27s21
     $     nsymb,ilout,ncsf,ioverwrite,npadddi,bc,ibc)                  11d10s22
c
c     if ioverwrite is nonzero, overwrite ivb with igb.               8d18s21
      implicit real*8 (a-h,o-z)                                         7d27s21
      external second                                                   5d4s22
      integer*4 ipack4(2)                                               7d27s21
      integer*8 ipack8                                                  7d27s21
      logical lpr,lpr2                                                  1d5s23
      equivalence (ipack8,ipack4)                                       7d27s21
      dimension iwbra(*),ixmt(8,*),iosym(*),i2eop(2,*),iwket(*),        7d27s21
     $     xout(nrootb,*),ism(*),irel(*),irefo(*),nbasdws(*),idoubo(*), 7d27s21
     $     nvirt(*),multh(8,8),ixmtf(8),ncsf(*)                         8d18s21
      include "common.store"                                            7d27s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data icall/0/                                                     8d23s21
      common/cpucom/tovr,top(10),tso(11)                                5d4s22
      save icall                                                        8d23s21
      icall=icall+1                                                     8d23s21
      lpr=icall.eq.-1
      lpr2=icall.eq.-1
      if(lpr)then                                                       8d25s21
       pp=phase1*2d0
       call testme(iwbra,iwket,nsymb,idoubo,irefo,nvirt,nbasdws,        8d25s21
     $     multh,mdon,mdoop,ism,irel,norb,ncsf,lxmt,ixmt,iosym,i2eop,   8d25s21
     $      n2e,phase1,phase2,ioverwrite,ixw1,ixw2,maxbx,bc,ibc)        11d10s22
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
      imff1=ilook+nrootk*(ncsftk+1)
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
       ivdknon=ivdk+nrootk*ndoub                                        8d12s21
       icsf=ivdknon+nrootk*mdoub                                        8d12s21
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
      imff1b=iloob+nrootb*(ncsftb+1)
      if(lpr)write(6,*)('for bra, imff1b '),imff1b,nrootb,ncsftb
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
       call enough('psioppsi.  1',bc,ibc)
       call makegs(ibc(igsdiag),ibc(nff1b),mdon,mdoop,nsymb,nvirt,       8d21s21
     $      iwbra(2),multh,ncsf,nrootk,maxbxq,npadddi,bc,ibc)           11d10s22
      end if
      mff2b=ibc(inext)                                                   7d21s21
      if(lpr)write(6,*)('for bra, mff2b = '),mff2b
       ilout=ibcoff                                                      7d27s21
       ibcoff=ilout+ncsftb*nrootk                                        7d27s21
       nwds=ncsftb*nrootk-1                                             8d18s21
       do i=ilout,ilout+nwds                                             8d18s21
        bc(i)=0d0                                                        7d27s21
       end do                                                            7d27s21
       call enough('psioppsi.  2',bc,ibc)
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
       call enough('psioppsi.  3',bc,ibc)
       call makegd(ibc(igddiag),ibc(nff2b),mdon,mdoop,nsymb,nvirt,      8d24s21
     $      iwbra(2),multh,ncsf,ibc(icsf2),nrootk,maxbxdq,npadddi,bc,   11d15s22
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
       ivdbnon=ivdb+ndoubb*nrootb                                        8d21s21
        igdb=ibcoff                                                      8d3s21
        ibcoff=igdb+mdoubb*nrootk                                         8d3s21
        nwds=mdoubb*nrootk-1                                             8d18s21
        call enough('psioppsi.  4',bc,ibc)
        do iz=igdb,igdb+nwds                                             8d18s21
         bc(iz)=0d0                                                      8d3s21
        end do                                                           8d3s21
      else
      end if                                                            7d27s21
      if(lpr)write(6,*)('calling hcbk '),ibcoff,ilout,bc(ilout)
      call hcbk(bc(ilout),ncsftb,bc(ilook),ncsftk,ibc(nff0b),ibc(nff0k),7d27s21
     $     ncsf,ibc(icsf2),ibc(jff0b),ibc(jff0k),nrootk,mdon,           8d26s21
     $     mdoop,ixw1,ixw2,nec,ism,irel,irefo,norb,lxmt,ixmt,multh,     7d27s21
     $     phase1,nbasdws,idoubo,iosym,n2e,i2eop,phase2,nsymb,shift,    7d27s21
     $     mff1,mff2,ibc(igsdiag),ibc(nff1b),ibc(iff1b),ibc(igddiag),     8d21s21
     $     ibc(nff2b),ibc(iff2b),iwbra(2),ibc(ihsdiag),ibc(nff1),       8d24s21
     $     ibc(iff1),ibc(ihddiag),ibc(nff2),ibc(iff2),iwket(2),         7d27s21
     $     nvirt,maxbx,maxbxd,srh,sr2,ibc(nfdatb),ibc(nfdat),bc(igdb),
     $     bc(ivdknon),lpr2,bc,ibc)                                     1d5s23
      if(lpr)write(6,*)('gsum '),ncsftb,nrootk,ibcoff
      call second(time1)
      call dws_gsumf(bc(ilout),ncsftb*nrootk)                           8d18s21
      if(ioverwrite.ne.0)then                                           8d27s21
       do i=0,ncsftb*nrootk-1                                           8d27s21
        bc(iloob+i)=bc(iloob+i)+bc(ilout+i)                             8d27s21
       end do                                                           8d27s21
      end if                                                            8d27s21
      if(lpr)write(6,*)('dotvall')
      call dotvall(bc(ilout),ibc(igsdiag),ibc(igddiag),nrootk,          8d24s21
     $     bc(iloob),                                                   8d19s21
     $    ibc(ihsdiagb),ibc(ihddiagb),nrootb,ncsftb,mff1,mff2,
     $     ibc(nff1b),                                                  8d24s21
     $     ibc(nff2b),mdon,mdoop,iwbra(2),multh,nsymb,ncsf,             8d24s21
     $      ibc(icsf2),nvirt,xout,ibc(iff1b),ibc(iff2b),norb,irefo,     8d24s21
     $      ibc(nff0b),ibc(jff0b),ibc(nfdatb),bc(ivdb),bc(igdb),mdoubb,   8d3s21
     $      ndoubb,bc(ivdk),n2e,iosym,i2eop,ixmt,phase2,sr2,             8d12s21
     $     idoubo,nbasdws,ioverwrite,0,idum,idum,idum,idum,idum,idum,bc,11d10s22
     $     ibc)                                                         11d10s22
      call second(time2)                                                5d4s22
      telap=time2-time1-tovr                                            5d4s22
      top(10)=top(10)+telap                                             11d17s22
      time1=time2                                                       5d4s22
      if(lpr)write(6,*)('cleanup'),ibcoff
      if(lpr)call prntm2(xout,nrootb,nrootk,nrootb)
      if(ioverwrite.ne.0.and.mff2b.lt.0)then                            8d23s21
       nwds=mdoubb*nrootk-1                                             8d18s21
       do iz=0,nwds                                                     8d27s21
        bc(ivdb+iz)=bc(ivdb+iz)+bc(igdb+iz)                             8d27s21
       end do                                                           8d27s21
      if(lpr)write(6,*)('tofrob'),ndoubb,mdoubb
       call tofrob(bc(ivdb),bc(ivdbnon),nrootk,bc(nfdatb),nvirt,nsymb,  8d26s21
     $     multh,iwbra(2),1,ndoubb,mdoubb,bc(iff2b),bc,ibc)             11d10s22
      end if                                                            8d21s21
      if(mff2b.gt.0)then                                                8d21s21
       call unmakegd(ibc(igddiag),ibc(ihddiagb),ibc(nff2),mdon,mdoop,   8d21s21
     $     nsymb,nvirt,iwbra(2),multh,ncsf,ibc(icsf2),nrootk,ioverwrite,11d10s22
     $     bc,ibc)                                                      11d10s22
      end if                                                            8d18s21
      if(mff1b.gt.0)then                                                8d21s21
      if(lpr)write(6,*)('unmakegs'),igsdiag,ihsdiagb,nff1,
     $     ibc(igsdiag),ibcoff
       call unmakegs(ibc(igsdiag),ibc(ihsdiagb),ibc(nff1),mdon,mdoop,   8d21s21
     $      nsymb,nvirt,iwbra(2),multh,ncsf,nrootk,ioverwrite,bc,ibc)   11d10s22
      end if                                                            8d21s21
      if(lpr)write(6,*)('shift = '),shift
      if(abs(shift).gt.1d-10)then                                       7d27s21
       if(nrootb.ne.nrootk)then
        write(6,*)('we have shift, but roots are different!! '),nrootb,
     $       nrootk
        call dws_synca
        call dws_finalize
        stop
       end if
       do i=1,nrootb                                                    7d27s21
        xout(i,i)=xout(i,i)+shift                                       7d27s21
       end do                                                           7d27s21
      end if
      if(ioverwrite.eq.0)ibcoff=ilout                                   8d19s21
      return                                                            7d27s21
      end                                                               7d27s21
