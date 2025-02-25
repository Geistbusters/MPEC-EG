c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine testsome(iwbra,iwket,i2smb,i2smk,nsymb,irefo,          2d8s23
     $     nvirt,nbasdws,multh,mdon,ism,irel,norb,ih0n,nh0,iall,        10d20s21
     $     iifmx,ntype,isopt,irori,irw0,irw1,irw2,isymc,xout,l2e,bc,ibc)11d9s22
      implicit real*8 (a-h,o-z)
      integer*8 ipack8
      equivalence (ipack8,ipack4)
      dimension iwbra(*),iwket(*),ipack4(2),irefo(*),nvirt(*),          2d8s23
     $     nbasdws(*),multh(8,8),ism(*),irel(*),ih0n(*),nh0(*),
     $     iall(8,8,8),iifmx(*),isopt(*),xout(*)                        11d12s21
      include "common.store"
      common/n0virtcm/n0virt                                            10d25s21
      data icall/0/                                                     2d28s22
      save icall                                                        2d28s22
      icall=icall+1                                                     2d28s22
      write(6,*)('testsome for call '),icall
       ihddiagk=1                                                        7d22s21
       ihddiagb=1
       nff2k=1                                                           7d22s21
       iff2k=1                                                           7d22s21
       nff2b=1                                                           7d22s21
       iff2b=1                                                           7d22s21
       icsfb=1                                                          11d16s21
       icsfb2=1                                                         11d16s21
       icsfk=1                                                          11d16s21
       icsfbk=1                                                         11d16s21
      n0virt=norb+1                                                     10d25s21
      ibcoffo=ibcoff
      mdookp=iwket(21)                                                  8d31s21
      mdoobp=iwbra(21)                                                  8d31s21
      i2sk=iwket(1)-1                                                   8d30s21
      i2sb=iwbra(1)-1                                                   8d30s21
      nrootb=iwbra(3)                                                   8d25s21
      nec=iwbra(7)                                                      7d27s21
      iloob=iwbra(4)+iwbra(13)                                          8d25s21
      ipack8=ibc(iloob)                                                 5d12s21
      ncsftb=ipack4(1)                                                  5d12s21
      iloob=iloob+1                                                     5d12s21
      nff0b=iwbra(4)+iwbra(9)                                           8d25s21
      jff0b=iwbra(4)+iwbra(10)                                          8d25s21
      nff0k=iwket(4)+iwket(9)                                           8d25s21
      jff0k=iwket(4)+iwket(10)                                          8d25s21
      nrootk=iwket(3)                                                   8d25s21
      ilook=iwket(4)+iwket(13)                                          8d25s21
      ipack8=ibc(ilook)                                                 5d12s21
      ncsftk=ipack4(1)                                                  5d12s21
      ilook=ilook+1                                                     5d12s21
      imff1b=iloob+nrootb*(ncsftb+1)
      mff1b=ibc(imff1b)                                                 8d25s21
      imff1k=ilook+nrootk*(ncsftk+1)                                    8d19s21
      mff1k=ibc(imff1k)                                                  7d15s21
      inextk=imff1k+1
      inextb=imff1b+1
      ntotb=ncsftb
      ntotk=ncsftk
      nntot=0
      do isb=1,nsymb
       nntot=nntot+nh0(isb)                                             10d20s21
      end do
      itest=ibcoff
      ibcoff=itest+max(nrootb,nrootk)*nrootk
      call enough('testsome.  1',bc,ibc)
      call dgemm('t','n',nrootk,nrootk,ncsftk,1d0,bc(ilook),ncsftk,
     $      bc(ilook),ncsftk,0d0,bc(itest),nrootk,
     d' testsome.  1')
      if(mff1k.gt.0)then
       ihsdiagk=imff1k+1                                                  7d9s21
       nff1k=ihsdiagk+nsymb*mdookp*2                                    10d20s21
       iff1k=nff1k+nsymb*mdookp                                         10d20s21
       icsfk=iff1k+ibc(imff1k)                                             7d9s21
       inextk=icsfk+mdookp-mdon                                         10d20s21
       ivx=ibcoff                                                       7d16s21
       ibcoff=ivx+max(nrootb,nrootk)                                                 7d16s21
       call enough('testsome.  2',bc,ibc)
       call sotest(ibc(ihsdiagk),ibc(nff1k),ibc(icsfk),iwket(2),
     $    nvirt,bc(itest),nrootk,mdon,mdookp,nsymb,multh,nsingk,bc(ivx),11d10s22
     $      bc,ibc)                                                     11d10s22
       ihsdiagb=imff1b+1
       nff1b=ihsdiagb+nsymb*mdoobp*2                                    10d20s21
       iff1b=nff1b+nsymb*mdoobp                                         10d20s21
       icsfb=iff1b+ibc(imff1b)
       inextb=icsfb+mdoobp-mdon                                         10d20s21
       ntotk=ntotk+nsingk
       itestb=ibcoff
       ivxb=itestb+nrootb*nrootb
       ibcoff=ivxb+nrootb
       call enough('testsome.  3',bc,ibc)
       do iz=itestb,ibcoff-1                                            8d25s21
        bc(iz)=0d0                                                      8d25s21
       end do                                                           8d25s21
       call sotest(ibc(ihsdiagb),ibc(nff1b),ibc(icsfb),iwbra(2),        8d25s21
     $     nvirt,bc(itestb),nrootb,mdon,mdoobp,nsymb,multh,nsingb,      10d20s21
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
      if(mff2k.gt.0)then
       ihddiagk=inextk+1                                                  7d21s21
       nff2k=ihddiagk+nsymb*mdookp*2                                    10d20s21
       iff2k=nff2k+nsymb*mdookp                                         10d20s21
       icsfk=iff2k+mff2k                                                   7d21s21
       icsfk2=icsfk+mdookp-mdon                                          10d20s21
       ivx2=ibcoff
       ibcoff=ivx2+max(nrootb,nrootk)                                                8d3s21
       call enough('testsome.  4',bc,ibc)
       call sotest2(ibc(ihddiagk),ibc(nff2k),ibc(icsfk),ibc(icsfk2),        7d21s21
     $      iwket(2),nvirt,bc(itest),nrootk,mdon,mdookp,nsymb,multh,    10d20s21
     $      ndoubk,bc(ivx2),bc,ibc)                                     11d10s22
       ihddiagb=inextb+1
       nff2b=ihddiagb+nsymb*mdoobp*2                                    10d20s21
       iff2b=nff2b+nsymb*mdoobp                                         10d20s21
       icsfb=iff2b+mff2b                                                8d27s21
       icsfb2=icsfb+mdoobp-mdon
       ntotk=ntotk+ndoubk
       itestb=ibcoff
       ivx2b=itestb+nrootb*nrootb                                       8d25s21
       ibcoff=ivx2b+nrootb
       call enough('testsome.  5',bc,ibc)
       do iz=itestb,ibcoff-1                                            8d25s21
        bc(iz)=0d0                                                      8d25s21
       end do                                                           8d25s21
       call sotest2(ibc(ihddiagb),ibc(nff2b),ibc(icsfb),ibc(icsfb2),        7d21s21
     $      iwbra(2),nvirt,bc(itestb),nrootb,mdon,mdoobp,nsymb,multh,   10d20s21
     $      ndoubb,bc(ivx2b),bc,ibc)                                    11d10s22
       ibcoff=itestb
       ntotb=ntotb+ndoubb
      else if(mff2k.lt.0)then
       mff2a=-mff2k                                                      8d3s21
       mdoubstore=inextk+1                                               8d3s21
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoubk=ipack4(1)                                                  8d3s21
       mdoubk=ipack4(2)                                                  8d3s21
       iff22k=mdoubstore+1                                               8d3s21
       nff22k=iff22k+mff2a                                                8d3s21
       nfdatk=nff22k+mdookp*nsymb                                       10d20s21
       ivdk=nfdatk+10*nsymb                                               8d3s21
       icsfk=ivdk+nrootk*(ndoubk+mdoubk)                                     8d12s21
       icsfk2=icsfk+mdookp-mdon                                          10d20s21
       ivx2=ibcoff
       ibcoff=ivx2+max(nrootk,nootb)                                                8d3s21
       call enough('testsome.  6',bc,ibc)
       mff2a=-mff2b                                                     8d25s21
       ivdb=inextb+1+1+mff2a+mdoobp*nsymb+10*nsymb                      10d20s21
       ivdknon=ivdk+nrootk*ndoubk                                           8d12s21
       call sotest2c(ibc(iff22k),ibc(nff22k),ibc(nfdatk),bc(ivdk),
     $      ibc(icsfk),ibc(icsfk2),iwket(2),nvirt,bc(itest),nrootk,mdon,8d25s21
     $       mdookp,nsymb,multh,ndoubxk,bc(ivx2),ibc(iff22k),bc,ibc)    11d10s22
       call tofrob(bc(ivdk),bc(ivdknon),nrootk,ibc(nfdatk),nvirt,nsymb, 8d26s21
     $      multh,iwket(2),1,ndoubk,mdoubk,ibc(iff22k),bc,ibc)          11d10s22
       mdoubstore=inextb+1                                              8d25s21
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoubb=ipack4(1)                                                  8d3s21
       mdoubb=ipack4(2)                                                  8d3s21
       iff22b=mdoubstore+1                                               8d3s21
       nff22b=iff22b+mff2a                                                8d3s21
       nfdatb=nff22b+mdoobp*nsymb                                       10d20s21
       ivdbnon=ivdb+nrootb*ndoubb
       icsfb=ivdb+nrootb*(ndoubb+mdoubb)                                8d25s21
       icsfb2=icsfb+mdoobp-mdon                                          11d16s21
       call tofrob(bc(ivdk),bc(ivdknon),nrootk,ibc(nfdatk),nvirt,
     $       nsymb,multh,iwket(2),1,ndoubk,mdoubk,ibc(iff22k),bc,ibc)   11d10s22
       ntotk=ntotk+ndoubxk
       itestb=ibcoff
       ivx2b=itestb+nrootb*nrootb
       ibcoff=ivx2b+nrootb
       call enough('testsome.  7',bc,ibc)
       do iz=itestb,ibcoff-1                                            8d25s21
        bc(iz)=0d0                                                      8d25s21
       end do                                                           8d25s21
       call sotest2c(bc(iff22b),ibc(nff22b),ibc(nfdatb),bc(ivdb),
     $       ibc(icsfb),ibc(icsfb2),iwbra(2),nvirt,bc(itestb),nrootb,   11d16s21
     $      mdon,mdoobp,nsymb,multh,ndoubxb,bc(ivx2b),bc(iff22b),bc,ibc)11d10s22
       ibcoff=itestb
       ntotb=ntotb+ndoubxb
      else
      end if
      irsum=ibcoff
      ibcoff=irsum+max(nrootb,nrootk)**2                                8d25s21
      do i=irsum,ibcoff-1
       bc(i)=0d0
      end do
      ivnew=ibcoff
      nff10k=ivnew+ntotk*nrootk
      ismv=nff10k+mdookp*3                                              10d20s21
      irelv=ismv+nntot
      iff10k=irelv+nntot
      if(mff1k.ne.0.and.mff2k.eq.0)then                                   7d27s21
       call convert1tor(ibc(ihsdiagk),ibc(nff1k),ibc(iff1k),
     $      ibc(icsfk),iwket(2),nvirt,nrootk,mdon,mdookp,nsymb,multh,   10d20s21
     $      ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),
     $      iff10k,norb,irefo,ibc(iff10k),bc(ivx),bc(ilook),ncsftk,
     $      ibc(nff0k),ibc(jff0k),bc,ibc)                               11d14s22
      else if(mff1k.eq.0.and.mff2k.ne.0)then                              7d27s21
       if(mff2k.gt.0)then                                                8d3s21
        call convert2tor(ibc(ihddiagk),ibc(nff2k),ibc(iff2k),           8d25s21
     $       ibc(icsfk),ibc(icsfk2),iwket(2),nvirt,nrootk,mdon,mdookp,  11d16s21
     $       nsymb,multh,ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,8d25s21
     $       ibc(irelv),iff10k,norb,irefo,ibc(iff10k),bc(ivx2),         8d25s21
     $       bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),bc,ibc)             11d14s22
       else                                                             8d3s21
        call convert2torc(ibc(iff22k),ibc(nff22k),ibc(nfdatk),
     $        bc(ivdknon),ibc(icsfk),ibc(icsfk2),iwket(2),nvirt,nrootk, 11d16s21
     $        mdon,mdookp,nsymb,multh,bc(ivx2),ntotk,bc(ivnew),         10d20s21
     $        ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),iff10k,norb,
     $        irefo,ibc(iff10k),bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),
     $        ibc(iff22k),bc(irsum))
       end if                                                           8d3s21
      else if(mff1k.ne.0.and.mff2k.ne.0)then                              7d27s21
       if(mff2k.gt.0)then                                                8d3s21
        call convert12tor(ibc(ihsdiagk),ibc(ihddiagk),ibc(nff1k),       8d25s21
     $      ibc(iff1k),ibc(nff2k),ibc(iff2k),ibc(icsfk),                8d25s21
     $      ibc(icsfk2),iwket(2),nvirt,nrootk,mdon,mdookp,nsymb,multh,  11d16s21
     $      ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),  8d25s21
     $      iff10k,norb,irefo,ibc(iff10k),bc(ivx2),bc(ilook),ncsftk,       7d22s21
     $      ibc(nff0k),ibc(jff0k),bc,ibc)                               11d14s22
       else
        call convert12torc(ibc(iff22k),ibc(nff22k),ibc(nfdatk),
     $        bc(ivdknon),ibc(icsfk),ibc(icsfk2),iwket(2),nvirt,nrootk, 11d16s21
     $        mdon,mdookp,nsymb,multh,bc(ivx2),ntotk,bc(ivnew),         10d20s21
     $        ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),iff10k,norb,
     $        irefo,ibc(iff10k),bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),
     $        ibc(iff22k),bc(irsum),ibc(ihsdiagk),ibc(nff1k),ibc(iff1k),11d14s22
     $       bc,ibc)                                                    11d14s22
       end if
      end if                                                            7d27s21
      ignew=ibcoff                                                      8d19s21
      ivold=ignew+ntotb*nrootk                                          8d25s21
      nff10b=ivold+ntotb*nrootb                                         8d25s21
      ismv=nff10b+mdoobp*3                                              10d20s21
      irelv=ismv+nntot                                                  8d19s21
      iff10b=irelv+nntot                                                8d19s21
      if(mff1b.ne.0.and.mff2b.eq.0)then                                   7d27s21
       call convert1tor(ibc(ihsdiagb),ibc(nff1b),ibc(iff1b),
     $      ibc(icsfb),iwbra(2),nvirt,nrootb,mdon,mdoobp,nsymb,multh,   10d20s21
     $      ntotb,bc(ivold),ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),
     $      iff10b,norb,irefo,ibc(iff10b),bc(ivx),bc(iloob),ncsftb,
     $      ibc(nff0b),ibc(jff0b),bc,ibc)                               11d14s22
      else if(mff1b.eq.0.and.mff2b.ne.0)then                              7d27s21
       if(mff2b.gt.0)then                                                8d3s21
        call convert2tor(ibc(ihddiagb),ibc(nff2b),ibc(iff2b),
     $        ibc(icsfb),ibc(icsfb2),iwbra(2),nvirt,nrootb,mdon,mdoobp, 11d16s21
     $        nsymb,multh,ntotb,bc(ivold),ibc(nff10b),ism,ibc(ismv),
     $        irel,ibc(irelv),iff10b,norb,irefo,ibc(iff10b),bc(ivx),
     $        bc(iloob),ncsftb,ibc(nff0b),ibc(jff0b),bc,ibc)            11d14s22
       else                                                             8d3s21
        call convert2torc(bc(iff22b),ibc(nff22b),ibc(nfdatb),
     $        bc(ivdbnon),ibc(icsfb),ibc(icsfb2),iwbra(2),nvirt,nrootb, 11d16s21
     $        mdon,mdoobp,nsymb,multh,bc(ivx2),ntotb,bc(ivold),         10d20s21
     $        ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),iff10b,norb,
     $        irefo,ibc(iff10b),bc(iloob),ncsftb,ibc(nff0b),ibc(jff0b),
     $        bc(iff22b),bc(irsum))                                     8d25s21
       end if                                                           8d3s21
      else if(mff1b.ne.0.and.mff2b.ne.0)then                              7d27s21
       if(mff2b.gt.0)then                                                8d3s21
        call convert12tor(ibc(ihsdiagb),ibc(ihddiagb),ibc(nff1b),       8d25s21
     $      ibc(iff1b),ibc(nff2b),ibc(iff2b),ibc(icsfb),                8d25s21
     $      ibc(icsfb2),iwbra(2),nvirt,nrootb,mdon,mdoobp,nsymb,multh,  11d16s21
     $      ntotb,bc(ivold),ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),  8d25s21
     $      iff10b,norb,irefo,ibc(iff10b),bc(ivx),bc(iloob),ncsftb,     8d25s21
     $      ibc(nff0b),ibc(jff0b),bc,ibc)                               11d14s22
       else
        call convert12torc(bc(iff22b),ibc(nff22b),ibc(nfdatb),
     $        bc(ivdbnon),ibc(icsfb),ibc(icsfb2),iwbra(2),nvirt,nrootb, 11d16s21
     $        mdon,mdoobp,nsymb,multh,bc(ivx2),ntotb,bc(ivold),         10d20s21
     $        ibc(nff10b),ism,ibc(ismv),irel,ibc(irelv),iff10b,norb,
     $        irefo,ibc(iff10b),bc(iloob),ncsftb,ibc(nff0b),ibc(jff0b),
     $        bc(iff22b),bc(irsum),ibc(ihsdiagb),ibc(nff1b),ibc(iff1b), 11d14s22
     $       bc,ibc)                                                    11d14s22
       end if
      end if                                                            7d27s21
      if(max(iabs(mff2b),mff1b).gt.0)then                                                7d19s21
       phasez=0d0
       do iz=0,ntotb*nrootk-1                                            8d21s21
        bc(ignew+iz)=0d0                                                8d21s21
       end do                                                           8d21s21
       idorbb=ibcoff                                                     8d31s21
       isorbb=idorbb+nntot                                              10d20s21
       idorbk=isorbb+nntot                                              10d20s21
       isorbk=idorbk+nntot                                              10d20s21
       icode=isorbk+nntot                                               10d20s21
       imap=icode+nntot                                                 10d20s21
       ibcoff=imap+nntot                                                10d20s21
       call enough('testsome.  8',bc,ibc)
       call hccsfbk4(                                                   10d20s21
     $      bc(ignew),ntotb,ibc(nff10b),ibc(iff10b),mdoobp,i2sb,i2smb,  10d20s21
     $      bc(ivnew),ntotk,ibc(nff10k),ibc(iff10k),mdookp,i2sk,i2smk,  10d20s21
     $      nrootk,mdon,nec,ibc(ismv),ibc(irelv),nh0,nntot,multh,       10d20s21
     $      nsymb,idum,idum,isymc,ibc(idorbb),ibc(isorbb),ibc(idorbk),  10d20s21
     $      ibc(isorbk),ibc(icode),ibc(imap),irori,irw0,irw1,irw2,ih0n, 10d20s21
     $      nh0,iall,iifmx,ntype,isopt,l2e,bc,ibc)                      11d9s22
       call dws_gsumf(bc(ignew),ntotb*nrootk)                             7d19s21
        call dgemm('t','n',nrootb,nrootk,ntotb,1d0,bc(ivold),ntotb,
     $      bc(ignew),ntotb,0d0,xout,nrootb,                            11d12s21
     d' testsome.  2')
       if(mff2b.lt.0)then                                               8d9s21
        if(mff1b.eq.0)then                                              8d9s21
         igc=ibcoff                                                     8d9s21
         ibcoff=igc+mdoubb*nrootk                                        8d9s21
         call enough('testsome.  9',bc,ibc)
         call deconvert2torc(bc(iff22b),ibc(nff22b),ibc(nfdatb),
     $         ibc(icsfb),ibc(icsfb2),                                  11d16s21
     $        iwbra(2),nvirt,nrootk,mdon,mdoobp,nsymb,multh,ntotb,      10d20s21
     $        bc(ignew),bc(iff22b),mdoubb,ibc(nff0b),bc(igc),ndoubb,bc, 11d10s22
     $        ibc)                                                      11d10s22
         ibcoff=igc                                                     8d9s21
        else                                                            8d9s21
         igc=ibcoff                                                     8d9s21
         ibcoff=igc+mdoubb*nrootk                                        8d9s21
         call enough('testsome. 10',bc,ibc)
         call deconvert12torc(bc(iff22b),ibc(nff22b),ibc(nfdatb),
     $         ibc(icsfb),ibc(icsfb2),                                  11d16s21
     $        iwbra(2),nvirt,nrootk,mdon,mdoobp,nsymb,multh,ntotb,      10d20s21
     $        bc(ignew),bc(iff22b),mdoubb,ibc(nff0b),bc(igc),ndoubb,       8d13s21
     $         ibc(nff1b),bc,ibc)                                       11d10s22
         ibcoff=igc                                                     8d9s21
        end if                                                          8d9s21
       end if                                                           8d9s21
      end if
      ibcoff=ibcoffo
      return
      end
