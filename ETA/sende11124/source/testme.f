c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine testme(iwbra,iwket,nsymb,idoubo,irefo,nvirt,nbasdws,   8d25s21
     $     multh,mdon,mdoop,ism,irel,norb,ncsf,llzz,ixmt,iosym,i2eop,   8d25s21
     $     n2e,phase1,phase2,ioverwrite,ixw1,ixw2,maxbx,bc,ibc)         11d10s22
      implicit real*8 (a-h,o-z)
      integer*8 ipack8
      equivalence (ipack8,ipack4)
      dimension iwbra(*),iwket(*),ipack4(2),idoubo(*),irefo(*),nvirt(*),8d25s21
     $     nbasdws(*),multh(8,8),ism(*),irel(*),ncsf(*),ixmt(8,*),      8d25s21
     $     iosym(*),ixmtf(8),i2eop(2,3),iden(8),nh0av(8)                3d18s22
      include "common.store"
      write(6,*)('hi, my name is testme'),phase1,loc(phase1)
      ibcoffo=ibcoff
      nrootb=iwbra(3)                                                   8d25s21
      nec=iwbra(7)                                                      7d27s21
      iloob=iwbra(4)+iwbra(13)                                          8d25s21
      ipack8=ibc(iloob)                                                 5d12s21
      ncsftb=ipack4(1)                                                  5d12s21
      write(6,*)('ncsftb = '),ncsftb
      iloob=iloob+1                                                     5d12s21
      nff0b=iwbra(4)+iwbra(9)                                           8d25s21
      jff0b=iwbra(4)+iwbra(10)                                          8d25s21
      nff0k=iwket(4)+iwket(9)                                           8d25s21
      jff0k=iwket(4)+iwket(10)                                          8d25s21
      nrootk=iwket(3)                                                   8d25s21
      ilook=iwket(4)+iwket(13)                                          8d25s21
      ipack8=ibc(ilook)                                                 5d12s21
      ncsftk=ipack4(1)                                                  5d12s21
      write(6,*)('ncsftk = '),ncsftk
      ilook=ilook+1                                                     5d12s21
      imff1b=iloob+nrootb*(ncsftb+1)
      mff1b=ibc(imff1b)                                                 8d25s21
      write(6,*)('imff1b '),imff1b,mff1b
      write(6,*)('going for '),iwbra(2)
      imff1k=ilook+nrootk*(ncsftk+1)                                    8d19s21
      mff1k=ibc(imff1k)                                                  7d15s21
      write(6,*)('imff1k '),imff1k,mff1k
      inextk=imff1k+1
      inextb=imff1b+1
      ntotb=ncsftb
      ntotk=ncsftk
      write(6,*)('setting ntotk to ncsftk = '),ntotk
      nntot=0
      do isb=1,nsymb
       nntot=nntot+nbasdws(isb)-idoubo(isb)
      end do
      write(6,*)('nntot = '),nntot
      itest=ibcoff
      ibcoff=itest+max(nrootb,nrootk)*nrootk
      call enough('testme.  1',bc,ibc)
      call dgemm('t','n',nrootk,nrootk,ncsftk,1d0,bc(ilook),ncsftk,
     $      bc(ilook),ncsftk,0d0,bc(itest),nrootk,
     d' testme.  1')
      write(6,*)('ortho test for ket after internals: ')
      call prntm2(bc(itest),nrootk,nrootk,nrootk)
      if(mff1k.gt.0)then
       ihsdiagk=imff1k+1                                                  7d9s21
       nff1k=ihsdiagk+nsymb*mdoop*2                                       7d9s21
       iff1k=nff1k+nsymb*mdoop                                            7d9s21
       icsfk=iff1k+ibc(imff1k)                                          10d25s21
       inextk=icsfk+mdoop-mdon                                          10d25s21
       ivx=ibcoff                                                       7d16s21
       ibcoff=ivx+max(nrootb,nrootk)                                                 7d16s21
       call enough('testme.  2',bc,ibc)
       write(6,*)('sotest for ket')
       call sotest(ibc(ihsdiagk),ibc(nff1k),ibc(icsfk),iwket(2),        10d25s21
     $     nvirt,bc(itest),nrootk,mdon,mdoop,nsymb,multh,nsingk,bc(ivx),11d10s22
     $      bc,ibc)                                                     11d10s22
       write(6,*)('nsing from sotest: '),nsingk
       ihsdiagb=imff1b+1
       write(6,*)('ihsdiagb: '),ihsdiagb
       nff1b=ihsdiagb+nsymb*mdoop*2
       iff1b=nff1b+nsymb*mdoop
       icsfb=iff1b+ibc(imff1b)
       inextb=icsfb+mdoop-mdon
       ntotk=ntotk+nsingk
       write(6,*)('adding nsingk = '),nsingk,('to ntotk to get '),
     $      ntotk
       itestb=ibcoff
       ivxb=itestb+nrootb*nrootb
       ibcoff=ivxb+nrootb
       call enough('testme.  3',bc,ibc)
       do iz=itestb,ibcoff-1                                            8d25s21
        bc(iz)=0d0                                                      8d25s21
       end do                                                           8d25s21
       write(6,*)('sotest for bra')
       call sotest(ibc(ihsdiagb),ibc(nff1b),ibc(icsfb),iwbra(2),        8d25s21
     $     nvirt,bc(itestb),nrootb,mdon,mdoop,nsymb,multh,nsingb,       8d25s21
     $     bc(ivxb),bc,ibc)                                             11d10s22
       ntotb=ntotb+nsingb
       ibcoff=itestb
      else                                                              7d15s21
       ihsdiagk=1                                                        7d9s21
       nff1k=1                                                           7d9s21
       iff1k=1                                                           7d9s21
       ihsdiagb=1                                                       8d25s21
       nff1b=1
       iff1k=1
      end if                                                            7d9s21
      mff2b=ibc(inextb)
      mff2k=ibc(inextk)
      write(6,*)('what we have for mff2: '),mff2b,mff2k
      if(mff2k.gt.0)then
       ihddiagk=inextk+1                                                  7d21s21
       nff2k=ihddiagk+nsymb*mdoop*2                                       7d21s21
       iff2k=nff2k+nsymb*mdoop                                            7d21s21
       icsfk=iff2k+mff2k                                                   7d21s21
       icsf2=icsfk+mdoop-mdon                                            7d21s21
       ivx2=ibcoff
       ibcoff=ivx2+max(nrootb,nrootk)                                                8d3s21
       call enough('testme.  4',bc,ibc)
       write(6,*)('sotest2 for ket')
       call prntm2(bc(itest),nrootk,nrootk,nrootk)
       call sotest2(ibc(ihddiagk),ibc(nff2k),ibc(icsfk),ibc(icsf2),        7d21s21
     $      iwket(2),nvirt,bc(itest),nrootk,mdon,mdoop,nsymb,multh,     7d21s21
     $      ndoubk,bc(ivx2),bc,ibc)                                     11d10s22
       ihddiagb=inextb+1
       nff2b=ihddiagb+nsymb*mdoop*2                                       7d21s21
       iff2b=nff2b+nsymb*mdoop                                            7d21s21
       icsfb=iff2b+mff2b                                                8d27s21
       ntotk=ntotk+ndoubk
       write(6,*)('adding ndoubk = '),ndoubk,(' to ntotk to get '),
     $      ntotk
       itestb=ibcoff
       ivx2b=itestb+nrootb*nrootb                                       8d25s21
       ibcoff=ivx2b+nrootb
       call enough('testme.  5',bc,ibc)
       do iz=itestb,ibcoff-1                                            8d25s21
        bc(iz)=0d0                                                      8d25s21
       end do                                                           8d25s21
       write(6,*)('sotest2 for bra')
       call sotest2(ibc(ihddiagb),ibc(nff2b),ibc(icsfb),ibc(icsf2),        7d21s21
     $      iwbra(2),nvirt,bc(itestb),nrootb,mdon,mdoop,nsymb,multh,     7d21s21
     $      ndoubb,bc(ivx2b),bc,ibc)                                    11d10s22
       ibcoff=itestb
       ntotb=ntotb+ndoubb
      else if(mff2k.lt.0)then
       mff2a=-mff2k                                                      8d3s21
       mdoubstore=inextk+1                                               8d3s21
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoubk=ipack4(1)                                                  8d3s21
       mdoubk=ipack4(2)                                                  8d3s21
       write(6,*)('ndoubk,mdoubk: '),ndoubk,mdoubk
       iff22k=mdoubstore+1                                               8d3s21
       nff22k=iff22k+mff2a                                                8d3s21
       nfdatk=nff22k+mdoop*nsymb                                          8d3s21
       ivdk=nfdatk+10*nsymb                                               8d3s21
       icsfk=ivdk+nrootk*(ndoubk+mdoubk)                                     8d12s21
       icsf2=icsfk+mdoop-mdon                                            8d3s21
       ivx2=ibcoff
       ibcoff=ivx2+max(nrootk,nootb)                                                8d3s21
       call enough('testme.  6',bc,ibc)
       mff2a=-mff2b                                                     8d25s21
       ivdb=inextb+1+1+mff2a+mdoop*nsymb+10*nsymb
       ivdknon=ivdk+nrootk*ndoubk                                           8d12s21
       write(6,*)('sotest2c for ket')
       call sotest2c(bc(iff22k),ibc(nff22k),ibc(nfdatk),bc(ivdk),
     $       ibc(icsfk),ibc(icsf2),iwket(2),nvirt,bc(itest),nrootk,mdon,8d25s21
     $       mdoop,nsymb,multh,ndoubxk,bc(ivx2),bc(iff22k),bc,ibc)      11d10s22
       call tofrob(bc(ivdk),bc(ivdknon),nrootk,ibc(nfdatk),nvirt,nsymb, 8d26s21
     $      multh,iwket(2),1,ndoubk,mdoubk,bc(iff22k),bc,ibc)           11d10s22
       mdoubstore=inextb+1                                              8d25s21
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoubb=ipack4(1)                                                  8d3s21
       mdoubb=ipack4(2)                                                  8d3s21
       iff22b=mdoubstore+1                                               8d3s21
       nff22b=iff22b+mff2a                                                8d3s21
       nfdatb=nff22b+mdoop*nsymb                                          8d3s21
       ivdbnon=ivdb+nrootb*ndoubb
       icsfb=ivdb+nrootb*(ndoubb+mdoubb)                                8d25s21
       call tofrob(bc(ivdk),bc(ivdknon),nrootk,ibc(nfdatk),nvirt,
     $       nsymb,multh,iwket(2),1,ndoubk,mdoubk,bc(iff22k),bc,ibc)    11d10s22
       ntotk=ntotk+ndoubxk
       write(6,*)('adding ndoubxk = '),ndoubxk,('to ntotk to get '),
     $      ntotk
       itestb=ibcoff
       ivx2b=itestb+nrootb*nrootb
       ibcoff=ivx2b+nrootb
       call enough('testme.  7',bc,ibc)
       do iz=itestb,ibcoff-1                                            8d25s21
        bc(iz)=0d0                                                      8d25s21
       end do                                                           8d25s21
       call sotest2c(bc(iff22b),ibc(nff22b),ibc(nfdatb),bc(ivdb),
     $       ibc(icsfb),ibc(icsf2),iwbra(2),nvirt,bc(itestb),nrootb,
     $      mdon,mdoop,nsymb,multh,ndoubxb,bc(ivx2b),bc(iff22b),bc,ibc) 11d10s22
       ibcoff=itestb
       ntotb=ntotb+ndoubxb
      else
       ihddiagk=1                                                        7d22s21
       ihddiagb=1
       nff2k=1                                                           7d22s21
       iff2k=1                                                           7d22s21
       nff2b=1                                                           7d22s21
       iff2b=1                                                           7d22s21
       icsf2=1                                                          7d22s21
      end if
      irsum=ibcoff
      ibcoff=irsum+max(nrootb,nrootk)**2                                8d25s21
      do i=irsum,ibcoff-1
       bc(i)=0d0
      end do
      write(6,*)('going after vnew: '),ibcoff,ntotk
      ivnew=ibcoff
      nff10k=ivnew+ntotk*nrootk
      ismv=nff10k+mdoop*3
      irelv=ismv+nntot
      iff10k=irelv+nntot
      write(6,*)('iff10k = '),iff10k,loc(bc(iff10k))
      if(mff1k.ne.0.and.mff2k.eq.0)then                                   7d27s21
       if(bc(132).ne.-132d0)then
         write(6,*)('b4 convert1tor ')
        call dws_synca
        call dws_finalize
        stop
       end if
       call convert1tor(ibc(ihsdiagk),ibc(nff1k),ibc(iff1k),
     $      ibc(icsfk),iwket(2),nvirt,nrootk,mdon,mdoop,nsymb,multh,
     $      ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),
     $      iff10k,norb,irefo,ibc(iff10k),bc(ivx),bc(ilook),ncsftk,
     $      ibc(nff0k),ibc(jff0k),bc,ibc)                               11d14s22
       if(bc(132).ne.-132d0)then
         write(6,*)('back from convert1tor ')
        call dws_synca
        call dws_finalize
        stop
       end if
      else if(mff1k.eq.0.and.mff2k.ne.0)then                              7d27s21
       if(mff2k.gt.0)then                                                8d3s21
        call convert2tor(ibc(ihddiagk),ibc(nff2k),ibc(iff2k),           8d25s21
     $       ibc(icsfk),ibc(icsf2),iwket(2),nvirt,nrootk,mdon,mdoop,    8d25s21
     $       nsymb,multh,ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,8d25s21
     $       ibc(irelv),iff10k,norb,irefo,ibc(iff10k),bc(ivx2),         8d25s21
     $       bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),bc,ibc)             11d14s22
       else                                                             8d3s21
        call convert2torc(bc(iff22k),ibc(nff22k),ibc(nfdatk),
     $        bc(ivdknon),ibc(icsfk),ibc(icsf2),iwket(2),nvirt,nrootk,  8d25s21
     $        mdon,mdoop,nsymb,multh,bc(ivx2),ntotk,bc(ivnew),
     $        ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),iff10k,norb,
     $        irefo,ibc(iff10k),bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),
     $        bc(iff22k),bc(irsum),bc,ibc)                              11d14s22
       end if                                                           8d3s21
      else if(mff1k.ne.0.and.mff2k.ne.0)then                              7d27s21
       if(mff2k.gt.0)then                                                8d3s21
        write(6,*)('b4 convert12tor, ntotk = '),ntotk
        call convert12tor(ibc(ihsdiagk),ibc(ihddiagk),ibc(nff1k),       8d25s21
     $      ibc(iff1k),ibc(nff2k),ibc(iff2k),ibc(icsfk),                8d25s21
     $      ibc(icsf2),iwket(2),nvirt,nrootk,mdon,mdoop,nsymb,multh,    8d25s21
     $      ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),  8d25s21
     $      iff10k,norb,irefo,ibc(iff10k),bc(ivx2),bc(ilook),ncsftk,       7d22s21
     $      ibc(nff0k),ibc(jff0k),bc,ibc)                               11d14s22
        write(6,*)('after convert12tor, ntotk = '),ntotk
       else
        write(6,*)('calling convert12torc ')
        call convert12torc(bc(iff22k),ibc(nff22k),ibc(nfdatk),
     $        bc(ivdknon),ibc(icsfk),ibc(icsf2),iwket(2),nvirt,nrootk,  8d25s21
     $        mdon,mdoop,nsymb,multh,bc(ivx2),ntotk,bc(ivnew),
     $        ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),iff10k,norb,
     $        irefo,ibc(iff10k),bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),
     $        bc(iff22k),bc(irsum),ibc(ihsdiagk),ibc(nff1k),ibc(iff1k), 11d14s22
     $       bc,ibc)                                                    11d14s22
        write(6,*)('back')
       end if
      end if                                                            7d27s21
      write(6,*)('ibcoff after conversion: '),ibcoff
      call dgemm('t','n',nrootk,nrootk,ntotk,1d0,bc(ivnew),ntotk,
     $      bc(ivnew),ntotk,0d0,bc(ibcoff),nrootk,
     d' testme.  2')
      call prntm2(bc(ibcoff),nrootk,nrootk,nrootk)
      write(6,*)('going after gnew ')
      ignew=ibcoff                                                      8d19s21
      ivold=ignew+ntotb*nrootk                                          8d25s21
      nff10b=ivold+ntotb*nrootb                                         8d25s21
      ismv=nff10b+mdoop*3
      irelv=ismv+nntot                                                  8d19s21
      iff10b=irelv+nntot                                                8d19s21
      if(mff1b.ne.0.and.mff2b.eq.0)then                                   7d27s21
       write(6,*)('convert1tor ')
       call convert1tor(ibc(ihsdiagb),ibc(nff1b),ibc(iff1b),
     $      ibc(icsfb),iwbra(2),nvirt,nrootb,mdon,mdoop,nsymb,multh,    8d25s21
     $      ntotb,bc(ivold),ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),
     $      iff10b,norb,irefo,ibc(iff10b),bc(ivx),bc(iloob),ncsftb,
     $      ibc(nff0b),ibc(jff0b),bc,ibc)                               11d14s22
      else if(mff1b.eq.0.and.mff2b.ne.0)then                              7d27s21
       if(mff2b.gt.0)then                                                8d3s21
        write(6,*)('convert2tor')
        call convert2tor(ibc(ihddiagb),ibc(nff2b),ibc(iff2b),
     $        ibc(icsfb),ibc(icsf2),iwbra(2),nvirt,nrootb,mdon,mdoop,
     $        nsymb,multh,ntotb,bc(ivold),ibc(nff10b),ism,ibc(ismv),
     $        irel,ibc(irelv),iff10b,norb,irefo,ibc(iff10b),bc(ivx),
     $        bc(iloob),ncsftb,ibc(nff0b),ibc(jff0b),bc,ibc)            11d14s22
       else                                                             8d3s21
        write(6,*)('convert2torc')
        call convert2torc(bc(iff22b),ibc(nff22b),ibc(nfdatb),
     $        bc(ivdbnon),ibc(icsfb),ibc(icsf2),iwbra(2),nvirt,nrootb,
     $        mdon,mdoop,nsymb,multh,bc(ivx2),ntotb,bc(ivold),
     $        ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),iff10b,norb,
     $        irefo,ibc(iff10b),bc(iloob),ncsftb,ibc(nff0b),ibc(jff0b),
     $        bc(iff22b),bc(irsum),bc,ibc)                              11d14s22
       end if                                                           8d3s21
      else if(mff1b.ne.0.and.mff2b.ne.0)then                              7d27s21
       if(mff2b.gt.0)then                                                8d3s21
        write(6,*)('for convert12tor: '),ihsdiagb,ihddiagb,nff1b,iff1b,
     $      nff2b,iff2b,icsfb,icsf2,nrootb,ntotb,ivold,nff10b,ismv,
     $      irelv,iff10b,ivx,iloob,nff0b,jff0b
        call convert12tor(ibc(ihsdiagb),ibc(ihddiagb),ibc(nff1b),       8d25s21
     $      ibc(iff1b),ibc(nff2b),ibc(iff2b),ibc(icsfb),                8d25s21
     $      ibc(icsf2),iwbra(2),nvirt,nrootb,mdon,mdoop,nsymb,multh,     8d25s21
     $      ntotb,bc(ivold),ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),  8d25s21
     $      iff10b,norb,irefo,ibc(iff10b),bc(ivx),bc(iloob),ncsftb,     8d25s21
     $      ibc(nff0b),ibc(jff0b),bc,ibc)                               11d14s22
        write(6,*)('back ')
       else
        write(6,*)('convert12torc')
        call convert12torc(bc(iff22b),ibc(nff22b),ibc(nfdatb),
     $        bc(ivdbnon),ibc(icsfb),ibc(icsf2),iwbra(2),nvirt,nrootb,
     $        mdon,mdoop,nsymb,multh,bc(ivx2),ntotb,bc(ivold),
     $        ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),iff10b,norb,
     $        irefo,ibc(iff10b),bc(iloob),ncsftb,ibc(nff0b),ibc(jff0b),
     $        bc(iff22b),bc(irsum),ibc(ihsdiagb),ibc(nff1b),ibc(iff1b), 11d14s22
     $       bc,ibc)                                                    11d14s22
       end if
      end if                                                            7d27s21
      if(max(iabs(mff2b),mff1b).gt.0)then                                                7d19s21
       phasez=0d0
       write(6,*)('zeroing ignew '),ignew,ntotb,nrootk
       do iz=0,ntotb*nrootk-1                                            8d21s21
        bc(ignew+iz)=0d0                                                8d21s21
       end do                                                           8d21s21
       write(6,*)('calling hccsfbk ')
       call hccsfbk(bc(ignew),ntotb,bc(ivnew),ntotk,ibc(nff10b),           7d19s21
     $     ibc(nff10k),ncsf,ibc(iff10b),ibc(iff10k),nrootk,mdon,mdoop,
     $       ixw1,ixw2,nec,ibc(ismv),ibc(irelv),irefo,nntot,llzz,ixmt,
     $       multh,phase1,nbasdws,idoubo,iosym,n2e,i2eop,phase2,nsymb,
     $       shift,ixmtf,bc,ibc)                                        11d10s22
       write(6,*)('back from hccsfbk')
       call dws_gsumf(bc(ignew),ntotb*nrootk)                             7d19s21
       if(ioverwrite.ne.0)then
        write(6,*)('g dot g')
        call dgemm('t','n',nrootk,nrootk,ntotb,1d0,bc(ignew),ntotb,
     $      bc(ignew),ntotb,0d0,bc(ibcoff),nrootk,
     d' testme.  3')
       else
        write(6,*)('v dot g')
        call dgemm('t','n',nrootb,nrootk,ntotb,1d0,bc(ivold),ntotb,
     $      bc(ignew),ntotb,0d0,bc(ibcoff),nrootb,
     d' testme.  4')
        write(6,*)('ntotb '),ntotb
        do i=0,ntotb-1
         jgnew=ignew+i
         sz=0d0
         do j=0,nrootk-1
          sz=sz+bc(jgnew)**2
          jgnew=jgnew+ntotb
         end do
         sz=sqrt(sz/dfloat(nrootk))
         jvold=ivold+i
         sv=0d0
         do j=0,nrootb-1
          sv=sv+bc(jvold)**2
          jvold=jvold+ntotb
         end do
         sv=sqrt(sv/dfloat(nrootb))
         if(min(sv,sz).gt.1d-10)then
          write(6,*)i+1
          write(6,*)(bc(ignew+i+ntotb*j),j=0,nrootk-1),
     $         loc(bc(ignew+i+ntotb*2))
          write(6,*)(bc(ivold+i+ntotb*j),j=0,nrootb-1)
         end if
        end do
       end if
       call prntm2(bc(ibcoff),nrootb,nrootk,nrootb)
      end if
      ibcoff=ibcoffo
      return
      end
