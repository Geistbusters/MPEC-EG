c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cal1int(ih0,iorb,idorel,ascale,nbasdwsc,natom,ngaus,   5d25s18
     $     ibdat,isym,iapair,ibstor,isstor,idwsdeb,nowpro,              4d22s21
     $     iapairg,nbasisp,ipassr,nextrad,iextrad,nwavef,nwavrec,nder,  8d11s22
     $     nrestart,iftype,fstgth,ndfld,lambdaci,bc,ibc,idarot)         12d6s23
      implicit real*8 (a-h,o-z)                                         5d25s18
      logical lget                                                      11d6s19
      include "common.hf"                                               5d25s18
      include "common.store"                                            5d25s18
      integer*2 ipack2(2)                                               12d5s22
      equivalence (ipack4,ipack2)                                       12d5s22
      dimension nbasdwsc(8),iso(8),isym(3,8),iapair(3,*),ibstor(*),     5d25s18
     $     isstor(*),morb(8),isymg(3,8),iapairg(3,*),ibcode(8),iorb(8), 2d15s19
     $     morbc(8),morbp(8),nbasisp(*),idum4(64),idumsym(6),istinfo(11)1d2s20
     $     ,ipassr(*),iread(8),makegbas(8),icanog(8),iftype(*),fstgth(*)12d6s23
      common/lowersymcm/nsymbgx,iptno(8),ipts(8),nhsz(8),ipao(8)        4d25s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      myguessi=1                                                        4d22s21
      ierr=0                                                            1d11s23
      if(nsymb.eq.0)then                                                11d6s19
       lget=.true.                                                      11d6s19
      else                                                              11d6s19
       lget=.false.                                                     11d6s19
      end if                                                            11d6s19
      if(.not.lget)then                                                 11d6s19
       ncomp=1                                                           2d15s19
       if(idorel.ne.0)ncomp=2                                            2d15s19
       nsqbas=0                                                          8d8s14
       do isb=1,nsymb                                                    8d8s14
        iso(isb)=nsqbas                                                  5d25s18
        nsqbas=nsqbas+nbasdwsc(isb)**2                                   8d31s15
        iptno(isb)=ibcoff                                                5d9s19
        ipts(isb)=iptno(isb)+nbasisp(isb)                                5d9s19
        ipao(isb)=ipts(isb)+nbasisp(isb)                                 5d9s19
        ibcoff=ipao(isb)+3*nbasisp(isb)                                  5d9s19
       end do                                                            8d8s14
       ih0=ibcoff                                                        5d25s18
       idwss=ih0+nsqbas*2                                                  5d25s18
       ibcoff=idwss+nsqbas*2                                               5d25s18
       call enough('cal1int.  1',bc,ibc)
       call parah0(natom,ngaus,ibdat,ndum,bc(ih0),bc(idwss),isym,        5d25s18
     $            iapair,ibstor,isstor,iso,nsqbas,idwsdeb,idorel,ascale,8d29s22
     $      iftype,fstgth,ndfld,bc,ibc)                                 12d6s23
       ivecr=ibcoff                                                      5d25s18
       ibcoff=ivecr+nsqbas                                               5d25s18
       newbasdws=ibcoff                                                  2d15s19
       ibcoff=newbasdws+nsymb                                            2d15s19
       ibcoff=ibcoff+1                                                  8d22s22
       nsqbas=nsqbas+nsymb                                               2d15s19
       call enough('cal1int.  2',bc,ibc)
       do i=ivecr,ibcoff-1                                              10d25s20
        bc(i)=0d0                                                       10d25s20
       end do                                                           10d25s20
      end if                                                            11d6s19
      if(nowpro.eq.0)then                                               3d19s07
       open(unit=1,file='forbs',status='old',form='unformatted',        10d25s20
     $      iostat=ios)                                                 10d25s20
       if(ios.ne.0)then                                                 10d25s20
        write(6,*)('failed to open file forbs ... aborting')            10d25s20
        bc(ivecr)=132d0                                                 10d25s20
        go to 222                                                       10d25s20
       end if                                                           10d25s20
       rewind(unit=1)                                                   11d6s19
       if(iabs(myguessi).eq.1)then                                      3d14s16
        write(6,*)('orbitals from ao basis are used')                   2d18s16
       else                                                             2d18s16
        write(6,*)('orbitals from orthogonal basis are used ')          2d18s16
       end if                                                           2d18s16
       read(1)nsymbg,idorelg,ngausg,natomg,nwcontg,idum1,idum2,idum3,   11d21s19
     $      idum4,dum1,dum2,idumsym,istinfo,nextrad                     5d25s21
       if(nwavef.ne.0)then                                              3d23s21
        nwavrec=nwavrec+1
        if(nrestart.lt.0)then                                           8d12s22
         read(2)nsymbgr,idorelgr,ngausgr,natomgr,nwcontgr               8d11s22
         ier=0                                                          8d11s22
         if(nsymbg.ne.nsymbgr)ier=ier+1                                 8d11s22
         if(idorelg.ne.idorelgr)ier=ier+1                               8d11s22
         if(ngausg.ne.ngausgr)ier=ier+1                                 8d11s22
         if(natomg.ne.natomgr)ier=ier+1                                 8d11s22
         if(nwcontg.ne.nwcontgr)ier=ier+1                               8d11s22
         if(ier.ne.0)then                                                8d11s22
          write(6,*)('restart data differs for record '),nwavrec         8d11s22
          write(6,*)('gotA  '),nsymbgr,idorelgr,ngausgr,natomgr,nwcontgr 8d11s22
          write(6,*)('want '),nsymbg,idorelg,ngausg,natomg,nwcontg      8d11s22
          stop 'restart'                                                 8d11s22
         end if                                                          8d11s22
        else                                                            8d11s22
         write(2)nsymbg,idorelg,ngausg,natomg,nwcontg,idum1,idum2,idum3,  11d21s19
     $      idum4,dum1,dum2,idumsym,istinfo,nextrad                     5d25s21
        end if                                                          8d11s22
       end if                                                           3d23s21
       idorelg=iabs(idorelg)                                            6d24s16
       if(lget)then                                                     11d6s19
        nsymb=nsymbg                                                    11d6s19
        idorel=idorelg                                                  11d6s19
        ngaus=ngausg                                                    11d6s19
        natom=natomg                                                    11d6s19
        nwcont=nwcontg                                                  11d6s19
       end if                                                           11d6s19
       if(nsymbg.ne.nsymb)then                                          2d18s16
        write(6,*)('point group doesn''t match current one: '),nsymbg   2d18s16
       end if                                                           2d18s16
       if(idorelg.ne.idorel)then                                        2d18s16
        if(idorelg.eq.0)then                                            2d18s16
         write(6,*)('using non-relativistic orbitals as guess ')        2d18s16
        else if(idorel.eq.0)then                                        2d18s16
         write(6,*)('using relativistic orbitals as guess ')            2d18s16
        end if                                                          2d18s16
        if(abs(myguessi).ne.1)then                                      4d27s16
         write(6,*)('need to use ao rather than ob guess in this case') 4d27s16
         call dws_sync                                                  4d27s16
         call dws_finalize                                              4d27s16
         stop                                                           4d27s16
        end if                                                          4d27s16
       end if                                                           2d18s16
       if(ngausg.ne.ngaus)then                                          2d18s16
        write(6,*)('number of guassians does not match ')               2d18s16
       end if                                                           2d18s16
       nextrad=nextrad+1                                                1d10s19
       iextrad=ibcoff                                                   1d10s19
       ibcoff=iextrad+nextrad                                           1d10s19
       call enough('cal1int.  3',bc,ibc)
       inbasg=ibcoff                                                     2d18s16
       inbasgp=inbasg+nsymbg                                            1d28s19
       inbasgc=inbasgp+nsymbg                                           1d28s19
       ibdatg=inbasgc+nsymbg                                            2d15s19
       nbdat=9*ngausg+nwcontg                                           2d15s19
       ibcoff=ibdatg+nbdat                                              2d15s19
       icang=ibcoff                                                     5d3s23
       ibcoff=icang+nsymbg                                              5d3s23
       idelt1=loc(nbasdwsc(2))-loc(nbasdwsc(1))
       idelt2=loc(ibc(inbasg+1))-loc(ibc(inbasg))
       if(idelt1.ne.idelt2)then
        read(1)(morb(isb),isb=1,nsymbg)
        if(nwavef.ne.0)then                                             5d4s21
         nwavrec=nwavrec+1                                              5d4s21
         if(nrestart.lt.0)then                                          8d12s22
          read(2)(iread(isb),isb=1,nsymbg)                              8d11s22
          do isb=1,nsymbg                                               8d11s22
           if(morb(isb).ne.iread(isb))ier=ier+1                         8d11s22
          end do                                                        8d11s22
          if(ier.ne.0)then                                                8d11s22
           write(6,*)('restart data differs for record '),nwavrec         8d11s22
           write(6,*)('gotB  '),(iread(isb),isb=1,nsymbg)                8d11s22
           write(6,*)('want '),(morb(isb),isb=1,nsymbg)                 8d11s22
           stop 'restart'                                                 8d11s22
          end if                                                          8d11s22
         else                                                           8d11s22
          write(2)(morb(isb),isb=1,nsymbg)                               5d4s21
         end if                                                         8d11s22
        end if                                                          5d4s21
        do isb=1,nsymbg
         ibc(inbasg+isb-1)=morb(isb)
        end do
        read(1)(morbc(isb),isb=1,nsymbg)                                2d15s19
        if(nwavef.ne.0)then                                             5d4s21
         nwavrec=nwavrec+1                                              5d4s21
         if(nrestart.lt.0)then                                          8d12s22
          read(2)(iread(isb),isb=1,nsymbg)                              8d11s22
          do isb=1,nsymbg                                               8d11s22
           if(morbc(isb).ne.iread(isb))ier=ier+1                        8d11s22
          end do                                                        8d11s22
          if(ier.ne.0)then                                                8d11s22
           write(6,*)('restart data differs for record '),nwavrec         8d11s22
           write(6,*)('gotC  '),(iread(isb),isb=1,nsymbg)                8d11s22
           write(6,*)('want '),(morbc(isb),isb=1,nsymbg)                8d11s22
           stop 'restart'                                                 8d11s22
          end if                                                          8d11s22
         else                                                           8d11s22
          write(2)(morbc(isb),isb=1,nsymbg)                              5d4s21
         end if                                                         8d11s22
        end if                                                          5d4s21
        do isb=1,nsymbg                                                 2d15s19
         ibc(inbasgc+isb-1)=morbc(isb)                                  2d15s19
        end do                                                          2d15s19
        read(1)(morbp(isb),isb=1,nsymbg)
        if(nwavef.ne.0)then                                             5d4s21
         nwavrec=nwavrec+1                                              5d4s21
         if(nrestart.lt.0)then                                          8d12s22
          read(2)(iread(isb),isb=1,nsymbg)                              8d11s22
          do isb=1,nsymbg                                               8d11s22
           if(morbp(isb).ne.iread(isb))ier=ier+1                        8d11s22
          end do                                                        8d11s22
          if(ier.ne.0)then                                                8d11s22
           write(6,*)('restart data differs for record '),nwavrec         8d11s22
           write(6,*)('gotD  '),(iread(isb),isb=1,nsymbg)                8d11s22
           write(6,*)('want '),(morbp(isb),isb=1,nsymbg)                8d11s22
           stop 'restart'                                                 8d11s22
          end if                                                          8d11s22
         else                                                           8d11s22
          write(2)(morbp(isb),isb=1,nsymbg)                              5d4s21
         end if                                                         8d11s22
        end if                                                          5d4s21
        do isb=1,nsymbg
         ibc(inbasgp+isb-1)=morbp(isb)
        end do
       else
        read(1)(ibc(inbasg+isb),isb=0,nsymbg-1)                          2d18s16
        if(nwavef.ne.0)then
         nwavrec=nwavrec+1                                              5d4s21
         if(nrestart.lt.0)then                                          8d12s22
          ireadr=ibcoff                                                 8d11s22
          ibcoff=ireadr+nsymbg                                          8d11s22
          read(2)(ibc(ireadr+isb),isb=0,nsymbg-1)                       8d11s22
          do isb=0,nsymbg-1                                             8d11s22
           if(ibc(ireadr+isb).ne.ibc(inbasg+isb))ier=ier+1              8d11s22
          end do                                                        8d11s22
          ibcoff=ireadr                                                 8d11s22
          if(ier.ne.0)then                                                8d11s22
           write(6,*)('restart data differs for record '),nwavrec         8d11s22
           write(6,*)('gotE  '),(ibc(ireadr+isb),isb=0,nsymbg-1)         8d11s22
           write(6,*)('want '),(ibc(inbasg+isb),isb=0,nsymbg-1)         8d11s22
           stop 'restart'                                                 8d11s22
          end if                                                          8d11s22
         else                                                           8d11s22
          write(2)(ibc(inbasg+isb),isb=0,nsymbg-1)                       5d4s21
         end if                                                         8d11s22
        end if                                                          5d4s21
        read(1)(ibc(inbasgc+isb),isb=0,nsymbg-1)                        1d28s19
        if(nwavef.ne.0)then
         nwavrec=nwavrec+1                                              5d4s21
         if(nrestart.lt.0)then                                          8d12s22
          ireadr=ibcoff                                                 8d11s22
          ibcoff=ireadr+nsymbg                                          8d11s22
          call enough('cal1int.  4',bc,ibc)
          read(2)(ibc(ireadr+isb),isb=0,nsymbg-1)                       8d11s22
          do isb=0,nsymbg-1                                             8d11s22
           if(ibc(ireadr+isb).ne.ibc(inbasgc+isb))ier=ier+1             8d11s22
          end do                                                        8d11s22
          ibcoff=ireadr                                                 8d11s22
          if(ier.ne.0)then                                                8d11s22
           write(6,*)('restart data differs for record '),nwavrec         8d11s22
           write(6,*)('gotF  '),(ibc(ireadr+isb),isb=0,nsymbg-1)         8d11s22
           write(6,*)('want '),(ibc(inbasgc+isb),isb=0,nsymbg-1)        8d11s22
           stop 'restart'                                                 8d11s22
          end if                                                          8d11s22
         else                                                           8d11s22
          write(2)(ibc(inbasgc+isb),isb=0,nsymbg-1)                      5d4s21
         end if                                                         8d11s22
        end if
        read(1)(ibc(inbasgp+isb),isb=0,nsymbg-1)                        1d28s19
        if(nwavef.ne.0)then
         nwavrec=nwavrec+1                                              5d4s21
         if(nrestart.lt.0)then                                          8d12s22
          ireadr=ibcoff                                                 8d11s22
          ibcoff=ireadr+nsymbg                                          8d11s22
          call enough('cal1int.  5',bc,ibc)
          read(2)(ibc(ireadr+isb),isb=0,nsymbg-1)                       8d11s22
          do isb=0,nsymbg-1                                             8d11s22
           if(ibc(inbasgp+isb).ne.ibc(ireadr+isb))ier=ier+1             8d11s22
          end do                                                        8d11s22
          ibcoff=ireadr                                                 8d11s22
          if(ier.ne.0)then                                                8d11s22
           write(6,*)('restart data differs for record '),nwavrec         8d11s22
           write(6,*)('gotG  '),(ibc(ireadr+isb),isb=0,nsymbg-1)         8d11s22
           write(6,*)('want '),(ibc(inbasgp+isb),isb=0,nsymbg-1)        8d11s22
           stop 'restart'                                                 8d11s22
          end if                                                          8d11s22
         else                                                           8d11s22
          write(2)(ibc(inbasgp+isb),isb=0,nsymbg-1)                      5d4s21
         end if                                                         8d11s22
        end if                                                          5d4s21
        do isb=1,nsymbg
         morb(isb)=ibc(inbasg+isb-1)
         morbc(isb)=ibc(inbasgc+isb-1)                                  1d28s19
         morbp(isb)=ibc(inbasgp+isb-1)                                  1d28s19
        end do
       end if
       nallg=0                                                          10d23s17
       do isb=1,nsymbg                                                  10d23s17
        nallg=nallg+ibc(inbasgp+isb-1)                                  5d9s19
       end do                                                           10d23s17
       idum1=ibcoff                                                     11d21s19
       ibcoff=idum1+natomg*6                                            11d21s19
       write(6,*)('reading isymg '),12+natomg*6,nextrad,nsymb
       read(1)isymg,(bc(idum1+i),i=0,natomg*6-1),
     $      (bc(iextrad+i),i=0,nextrad-1),                              1d10s19
     $      (ipassr(i),i=1,nsymb*2)                                     1d10s19
     $      ,nlzzq                                                      1d2s20
       ii=1
       do i=1,nsymb*2,2                                                 5d4s23
        ipack4=ipassr(i)                                                5d4s23
        write(6,*)ii,ipack2(2),ipack2(1),ipassr(i)                      5d4s23
        icanog(ii)=ipack2(2)                                            5d4s23
        ii=ii+1                                                         5d4s23
       end do                                                           5d4s23
       ipack2(1)=nlzzq                                                  12d5s22
       ipack2(2)=lambdaci                                               12d5s22
       nlzzq=ipack4                                                     12d5s22
       if(nwavef.ne.0)then                                              5d4s21
        nwavrec=nwavrec+1                                               5d4s21
        if(nrestart.lt.0)then                                           8d12s22
         read(2)idum                                                    8d11s22
        else                                                            8d11s22
         write(6,*)('writing '),12+natomg*6,nextrad,nsymb
         write(2)isymg,(bc(idum1+i),i=0,natomg*6-1),                     5d4s21
     $      (bc(iextrad+i),i=0,nextrad-1),                              1d10s19
     $      (ipassr(i),i=1,nsymb*2),nlzzq                               1d2s20
        end if                                                          8d11s22
       end if                                                           5d4s21
       nvec=natom-1                                                     5d27s21
       nexpdata=0                                                       5d27s21
       if(nvec.eq.1)then                                                 5d26s21
        nexpdata=1                                                       5d26s21
       else if(nvec.gt.1)then                                            5d26s21
        nexpdata=3*nvec-3                                                5d26s21
       end if                                                            5d26s21
       nexpdata=nexpdata+1                                              5d27s21
       ntag=6+natom                                                     12d22s22
       nexpdata2=nexpdata+ntag                                          12d22s22
       if(nexpdata+2.eq.nextrad.or.nexpdata2+2.eq.nextrad)then          12d22s22
        nadd=2                                                          12d22s22
       else                                                             12d22s22
        nadd=0                                                          12d22s22
       end if                                                           12d22s22
       nexpdata=nexpdata+nadd                                           12d22s22
       if(nexpdata.le.nextrad)then                                      5d27s21
        bc(iextrad+nexpdata-1)=bc(iextrad+nextrad-1)                    5d27s21
        nextrad=nexpdata                                                5d27s21
       end if                                                           5d27s21
       read(1)((iapairg(j,i),j=1,3),i=1,natomg)                         10d27s17
       if(nwavef.ne.0)then                                              5d4s21
        nwavrec=nwavrec+1                                               5d4s21
        if(nrestart.lt.0)then                                           8d12s22
         read(2)idum                                                    8d11s22
        else                                                            8d11s22
         write(2)((iapairg(j,i),j=1,3),i=1,natomg)                       5d4s21
        end if                                                          8d11s22
       end if                                                           5d4s21
       isstorg=ibcoff                                                    10d23s17
       ibstorg=isstorg+nallg*3                                          11d22s19
       ibcoff=ibstorg+nallg*3                                           11d22s19
       call enough('cal1int.  6',bc,ibc)
       read(1)(bc(ibdatg+i),i=0,nbdat-1)                                2d15s19
       if(nwavef.ne.0)then                                              5d4s21
        nwavrec=nwavrec+1
        if(nrestart.lt.0)then                                           8d12s22
         ireadr=ibcoff                                                  8d11s22
         ibcoff=ireadr+nbdat                                            8d11s22
         call enough('cal1int.  7',bc,ibc)
         read(2)(bc(ireadr+i),i=0,nbdat-1)                              8d11s22
         do i=0,nbdat-1                                                 8d11s22
          if(ibc(ireadr+i).ne.ibc(ibdatg+i))then                          8d16s22
           ier=ier+1                                                    8d16s22
          end if                                                        8d16s22
         end do                                                         8d11s22
         ibcoff=ireadr                                                  8d11s22
         if(ier.ne.0)then                                                8d11s22
          write(6,*)('restart data differs for record '),nwavrec         8d11s22
          write(6,*)('gotH  '),(ibc(ireadr+i),i=0,nbdat-1)                8d11s22
          write(6,*)('want '),(ibc(ibdatg+i),i=0,nbdat-1)                8d11s22
          stop 'restart'                                                 8d11s22
         end if                                                          8d11s22
        else                                                            8d11s22
         write(2)(bc(ibdatg+i),i=0,nbdat-1)                              5d4s21
        end if                                                          8d11s22
       end if                                                           5d4s21
       nallg3=nallg*3                                                   1d2s20
       read(1)(ibc(ibstorg+i),ibc(isstorg+i),i=0,nallg3-1)              1d2s20
       if(nwavef.ne.0)then
        nwavrec=nwavrec+1
        if(nrestart.lt.0)then                                           8d11s22
         iread1=ibcoff                                                  8d11s22
         iread2=iread1+nallg3                                           8d16s22
         ibcoff=iread2+nallg3                                           8d16s22
         call enough('cal1int.  8',bc,ibc)
         read(2)(ibc(iread1+i),ibc(iread2+i),i=0,nallg3-1)              8d11s22
         do i=0,nallg3-1                                                8d11s22
          if(ibc(iread1+i).ne.ibc(ibstorg+i))ier=ier+1                  8d11s22
          if(ibc(iread2+i).ne.ibc(isstorg+i))ier=ier+1                  8d11s22
         end do                                                         8d11s22
         ibcoff=iread1                                                  8d11s22
         if(ier.ne.0)then                                                8d11s22
          write(6,*)('restart data differs for record '),nwavrec         8d11s22
          write(6,*)('gotI  '),
     $         (ibc(iread1+i),ibc(iread2+i),i=0,nallg3-1)               8d16s22
          write(6,*)('want '),(ibc(ibstorg+i),ibc(isstorg+i),           8d11s22
     $      i=0,nallg3-1)                                               8d11s22
          stop 'restart'                                                 8d11s22
         end if                                                          8d11s22
        else                                                            8d11s22
         write(2)(ibc(ibstorg+i),ibc(isstorg+i),                         5d4s21
     $      i=0,nallg3-1)                                               3d23s21
        end if                                                          8d11s22
       end if                                                           5d4s21
       do isb=1,nsymbg                                                  10d23s17
        ibcode(isb)=ibcoff                                              10d23s17
        ibcoff=ibcoff+ibc(inbasg+isb-1)                                 10d23s17
        call enough('cal1int.  9',bc,ibc)
        read(1)(ibc(ibcode(isb)+i),i=0,ibc(inbasg+isb-1)-1)             10d23s17
       end do                                                           10d23s17
       nn=0                                                             2d18s16
       do isb=0,nsymbg-1                                                2d18s16
        nn=nn+ibc(inbasg+isb)*ibc(inbasgp+isb)                          1d28s19
       end do                                                           2d18s16
       nn=nn*ncomp                                                      2d15s19
       ivguess=ibcoff                                                   2d18s16
       ibcoff=ivguess+nn                                                2d18s16
       ivdum=ibcoff                                                     2d18s16
       ibcoff=ivdum+nn                                                  2d18s16
       call enough('cal1int. 10',bc,ibc)
       jvguess=ivguess                                                  2d18s16
       do isb=0,nsymbg-1                                                2d18s16
        nh=ncomp*ibc(inbasgp+isb)*ibc(inbasg+isb)                       2d15s19
        if(iabs(myguessi).eq.1)then                                     3d14s16
         read(1)(bc(jvguess+i),i=0,nh-1)                                2d18s16
         if(nwavef.ne.0.and.nh.gt.0)then
          nwavrec=nwavrec+1
          if(nrestart.lt.0)then                                         8d12s22
           ireadr=ibcoff                                                8d11s22
           ibcoff=ireadr+nh                                             8d11s22
           call enough('cal1int. 11',bc,ibc)
           read(2)(bc(ireadr+i),i=0,nh-1)                               8d11s22
           do i=0,nh-1                                                  8d11s22
            if(bc(ireadr+i).ne.bc(jvguess+i))ier=ier+1                  8d11s22
           end do                                                       8d11s22
           ibcoff=ireadr                                                8d11s22
          else                                                          8d11s22
           write(2)(bc(jvguess+i),i=0,nh-1)                              5d4s21
          end if                                                        8d11s22
          nrow=ncomp*ibc(inbasgp+isb)
          ncol=ibc(inbasg+isb)
         end if                                                         5d4s21
         jvguess=jvguess+nh                                             2d15s19
        else                                                            2d18s16
         read(1)(bc(ivdum+i),i=0,nh-1)                                  2d18s16
         if(nwavef.ne.0.and.nh.gt.0)then
          nwavrec=nwavrec+1
          if(nrestart.lt.0)then                                         8d12s22
           read(2)dum                                                   8d11s22
          else                                                          8d11s22
           write(2)(bc(ivdum+i),i=0,nh-1)                                5d4s21
          end if                                                        8d11s22
         end if                                                         5d4s21
        end if                                                          2d18s16
        nh=ibc(inbasg+isb)*ibc(inbasgc+isb)                             5d4s23
        if(myguessi.eq.2)then                                           2d18s16
         read(1)(bc(jvguess+i),i=0,nh-1)                                2d18s16
         if(nwavef.ne.0.and.nh.gt.0)then
          nwavrec=nwavrec+1
          if(nrestart.lt.0)then                                         8d12s22
           ireadr=ibcoff                                                8d11s22
           ibcoff=ireadr+nh                                             8d11s22
           call enough('cal1int. 12',bc,ibc)
           read(2)(bc(ireadr+i),i=0,nh-1)                               8d11s22
           do i=0,nh-1                                                  8d11s22
            if(bc(ireadr+i).ne.bc(jvguess+i))ier=ier+1                  8d11s22
           end do                                                       8d11s22
           ibcoff=ireadr                                                8d11s22
           if(ier.ne.0)then                                                8d11s22
            write(6,*)('restart data differs for record '),nwavrec         8d11s22
            write(6,*)('gotK  '),(bc(ireadr+i),i=0,nh-1)                 8d11s22
            write(6,*)('want '),(bc(jvguess+i),i=0,nh-1)                8d11s22
            stop 'restart'                                                 8d11s22
           end if                                                          8d11s22
          else                                                          8d11s22
           write(2)(bc(jvguess+i),i=0,nh-1)                              5d4s21
          end if                                                        8d11s22
          nrow=ibc(inbasg+isb)
         end if                                                         5d4s21
         jvguess=jvguess+nh                                             2d15s19
        else                                                            2d18s16
         read(1)(bc(ivdum+i),i=0,nh-1)                                  2d18s16
         if(nwavef.ne.0.and.nh.gt.0)then                                5d4s21
          nwavrec=nwavrec+1
          if(nrestart.lt.0)then                                         8d12s22
           read(2)dum                                                   8d11s22
          else                                                          8d11s22
           write(2)(bc(ivdum+i),i=0,nh-1)                                5d4s21
          end if                                                        8d11s22
         end if                                                         5d4s21
        end if                                                          2d18s16
        nb=ibc(inbasg+isb)                                              10d23s17
       end do                                                           2d18s16
       ibcoff=ivdum                                                     2d18s16
       nsymbgx=nsymbg                                                   10d24s17
       idwsdeb=000
       call makeguess(idorel,idorelg,nsymb,nsymbg,ngaus,ngausg,         2d18s16
     $      bc(ibdat),bc(ibdatg),nbasdws,morb,bc(ivecr),                2d15s19
     $      bc(ivguess),bc(idwss),myguessi,ierr,idwsdeb,ibstor,isstor,  10d23s17
     $      ibc(ibstorg),ibc(isstorg),ibcode,iptno,ipts,iapair,isym,    10d27s17
     $      iapairg,isymg,nhsz,ipao,morbp,nbasisp,1,0,dum,idum,idum,    3d27s20
     $      idum,idum,idum,idum,idum,idum,idum,idum,idum,idum,          5d7s21
     $      bc(idwss),0,nvguess,ircode,bc,ibc,makegbas,icanog)          5d3s23
c
c     check for derivative information ...
c
       write(6,*)('check for derivative information '),ibcoff,iextrad
       nder=0                                                           9d12s22
       read(1,end=28)nstate                                             9d12s22
       if(nstate.eq.0)go to 28                                          12d2s22
       do i=1,nstate                                                    7d28s22
        read(1,end=28)idum                                              10d12s22
        read(1,end=28)dum                                               10d12s22
       end do                                                           7d28s22
       idarot=ibcoff                                                    7d11s23
   27  continue                                                         7d28s22
        read(1,end=28)ixyz,ia,ipass,ia2                                 7d28s22
        write(6,*)('for derivative '),ixyz,ia,ipass,ia2,ibcoff
        if(ipass.lt.0)then                                              7d13s23
         if(ixyz.le.3)then                                              8d3s23
          write(6,*)('for '),char(ichar('w')+ixyz),                      7d13s23
     $       (' component of dipole ')                                  7d13s23
         else if(ixyz.eq.4)then                                         8d3s23
          write(6,*)('for Q0')                                          8d3s23
         else if(ixyz.eq.5)then                                         8d3s23
          write(6,*)('for ReQ2')                                          8d3s23
         else if(ixyz.eq.6)then                                         8d3s23
          write(6,*)('for ImQ2')                                          8d3s23
         else if(ixyz.eq.7)then                                         8d3s23
          write(6,*)('for ReQ1')                                          8d3s23
         else                                                           8d3s23
          write(6,*)('for ImQ1')                                        8d3s23
         end if                                                         8d3s23
        else                                                            7d13s23
         if(ia2.eq.0)then                                               7d13s23
          write(6,*)('for d/d'),char(ichar('w')+ixyz),ia                7d13s23
         else                                                           7d13s23
          if(ipass.eq.1)then                                            7d13s23
           write(6,*)('for d/d'),char(ichar('w')+ixyz),ia,('+'),ia2     7d13s23
          else                                                          7d13s23
           write(6,*)('for d/d'),char(ichar('w')+ixyz),ia,('-'),ia2     7d13s23
          end if                                                        7d13s23
         end if                                                         7d13s23
        end if                                                          7d13s23
        ibc(ibcoff)=ixyz                                                7d13s23
        ibc(ibcoff+1)=ia                                                7d13s23
        ibc(ibcoff+2)=ipass                                             7d13s23
        ibc(ibcoff+3)=ia2                                               7d13s23
        ibcoff=ibcoff+4                                                 7d13s23
        do isb=1,nsymb                                                  7d28s22
         if(nbasdws(isb).gt.0)then                                      7d28s22
          ireadr=ibcoff                                                  8d15s23
          nn=nbasdws(isb)**2                                            7d28s22
          ibcoff=ireadr+nn                                               8d15s23
          call enough('calint.idarot',bc,ibc)                           7d11s23
          read(1)(bc(ireadr+i),i=0,nn-1)                                 8d15s23
          if(ipass.ge.0)then                                            8d15s23
           ireadr=ibcoff                                                  8d15s23
           ibcoff=ireadr+nn                                               8d15s23
           call enough('calint.idarot',bc,ibc)                           7d11s23
           read(1)(bc(ireadr+i),i=0,nn-1)                                 8d15s23
          end if                                                        8d15s23
         end if                                                         7d28s22
        end do
        nder=nder+1                                                     7d28s22
        go to 27                                                        7d28s22
   28  continue                                                         7d28s22
       close(unit=1)                                                     3d19s07
       do isb=1,nsymb                                                   2d15s19
        im=isb-1                                                        2d15s19
        ibc(newbasdws+im)=nbasdws(isb)                                  2d15s19
       end do                                                           2d15s19
      end if                                                             3d19s07
  222 continue                                                          10d25s20
      iarg1=nsqbas+1                                                    7d28s22
      ibc(ivecr+nsqbas)=nder                                            7d28s22
      if(ierr.ne.0)then                                                 2d18s16
       bc(ivecr)=132d0                                                  2d18s16
      end if                                                            2d18s16
      call dws_bcast(bc(ivecr),iarg1)                                   2d24s10
      if(abs(bc(ivecr)-132d0).lt.1d-8)then                              2d18s16
       write(6,*)('error from orbital guess routine ')                  2d18s16
       call dws_sync                                                    2d18s16
       call dws_finalize                                                2d18s16
       stop                                                             2d18s16
      end if                                                            2d18s16
      nder=ibc(ivecr+nsqbas)                                            7d28s22
      mder=0                                                            7d11s23
      do isb=1,nsymb                                                    2d15s19
       im=isb-1                                                         2d15s19
       nbasdws(isb)=ibc(newbasdws+im)                                   2d15s19
       if(nder.ne.0)then                                                7d11s23
        mder=mder+nbasdws(isb)**2                                       7d11s23
       end if                                                           7d11s23
      end do                                                            2d15s19
      if(nder.ne.0)then                                                 7d11s23
       mder=mder+4                                                      7d13s23
       mder=mder*nder                                                   7d11s23
       mder=mder*2                                                      8d15s23
       if(mynowprog.ne.0)then                                           7d11s23
        idarot=ibcoff                                                   7d11s23
        ibcoff=idarot+mder                                              7d11s23
        call enough('cal1int.idarot',bc,ibc)                             7d11s23
       end if                                                           7d11s23
       call dws_bcast(bc(idarot),mder)                                  7d11s23
      end if                                                            7d11s23
      jvecr=ivecr                                                       5d29s18
      jh0=ih0                                                           5d29s18
      do isb=1,nsymb
       iorb(isb)=jvecr                                                  5d29s18
       nbas=nbasisp(isb)*ncomp                                          2d15s19
       call square(bc(jh0),nbas)
       itmp=ibcoff
       ibcoff=itmp+nbas**2                                              2d15s19
       call enough('cal1int. 13',bc,ibc)
       nrow=nbas                                                        2d15s19
       do ipass=1,2
        if(nrow.gt.0)then                                               2d22s19
        call dgemm('n','n',nrow,nbasdws(isb),nbas,                      2d15s19
     $       1d0,bc(jh0),nrow,bc(jvecr),nbas,0d0,                       2d15s19
     $       bc(itmp),nrow,                                             6d6s18
     d' cal1int.  1')
        end if                                                          2d22s19
        do i=0,nbasdws(isb)-1                                           6d6s18
         do j=0,nrow-1                                                  6d6s18
          ji=itmp+j+nrow*i                                              6d6s18
          ij=jh0+i+nbasdws(isb)*j                                       6d6s18
          bc(ij)=bc(ji)                                                 5d29s18
         end do                                                         5d29s18
        end do                                                          5d29s18
        nrow=nbasdws(isb)                                               6d6s18
       end do                                                           5d29s18
       jh0=jh0+nbas*nbas                                                2d15s19
       jvecr=jvecr+nbas*nbasdws(isb)                                    2d15s19
      end do
      return
      end
