c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine runit(bc,ibc)                                          11d9s22
c mpec2.1 version zeta copyright u.s. government
      implicit real*8 (a-h,o-z)
      common/ddidwscm/dws_me,dws_np                                     10d23s14
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/drsigncm/drsign                                            8d20s24
      include "common.store"
      include "common.hf"
      include "common.cas"                                              8d7s14
      include "common.basis"
      include "common.spher"
      include "common.input"                                            3d1s10
      include "common.rys"                                              6d22s12
      include "common.mrci"                                             5d24s18
      include "common.print"                                            1d3s20
      include "common.mympi"                                            1d29s21
      integer dws_me,dws_np                                             5d24s18
      common/xcom/iadd1,iadd2,iadd3                                     7d14s15
      dimension inv(2,8,8,8)                                            4d9s18
      common/fnd2cm/inv                                                 4d9s18
      dimension pass(6),noc(8),iapairg(3,ida)                           1d19s23
      dimension jmats(idbk),kmats(idbk),ioooo(idbk),ionex(idbk)         3d23s12
     $     ,i3x(idbk)                                                   11d9s22
      character*3 symlab(8,8)                                           4d28s10
      character*50 pwd,cfile                                            12d13s22
      integer*1 idogrado1(4)                                            6d18s22
      equivalence (idogrado,idogrado1)                                  6d18s22
      data symlab(1,1)/'A  '/
      data (symlab(i,2),i=1,8)/'Ag ','B3u','B2u','B1g','B1u','B2g',     4d28s10
     $     'B3g','Au '/                                                 4d28s10
      data (symlab(i,3),i=1,2)/'A'' ','A" '/                            4d28s10
      data (symlab(i,4),i=1,2)/'A  ','B  '/                             4d28s10
      data (symlab(i,5),i=1,2)/'Ag ','Au '/                             4d28s10
      data (symlab(i,6),i=1,4)/'A1 ','B1 ','B2 ','A2 '/                 4d28s10
      data (symlab(i,7),i=1,4)/'Ag ','Au ','Bu ','Bg '/                 4d28s10
      data (symlab(i,8),i=1,4)/'A  ','B3 ','B2 ','B1 '/                 4d28s10
      common/unitcm/iunit                                               11d9s17
      common/ilimcm/ilimcode                                            5d7s19
      real*16 fl                                                        5d27s19
      COMMON/FACT16/FL(922),NCALL                                         9d23s99
      save
      drsign=1d0                                                        8d20s24
      iagrp(1)=0                                                        5d26s21
      ioffsx=1                                                          11d9s22
      ihighwtr=0                                                        3d2s09
      ienough=0                                                         3d2s09
      lmax=2
      ibcoff=1
      intmul=1
      iadd1=0                                                           7d14s15
      iadd2=0                                                           7d14s15
      iadd3=0                                                           7d14s15
      if(mynowprog.eq.0)then                                            1d31s21
       write(6,57)
   57  format(//,'Welcome to MPEC',
     $       //,'Massively',/'Parallel',/,'Electron',/,'Correlation',
     $       //'version 2.1 eta')                                       6d29s23
       write(6,59)
       include 'mpec.date'                                              1d31s21
       is=1
       call delim(pwd,is,ie)
       inewe=ie                                                         2d27s19
       do ix=1,ie                                                       2d27s19
        ixp=ix+5                                                        2d27s19
        ixm=ix-1                                                        2d27s19
        if(pwd(ix:ixp).eq.'source'.or.pwd(ix:ix+2).eq.'obj')then        10d31s22
         inewe=ix+4                                                     2d27s19
         pwd(ix:inewe)='data/'                                          2d27s19
         go to 1066                                                     2d27s19
        end if                                                          2d27s19
       end do
 1066  continue                                                         2d27s19
       ie=inewe                                                         2d27s19
       call basisz(ibdat,ibstor,isstor,pwd(is:ie),ifname,ncore,nmn,bc,  11d9s22
     $      ibc,cfile,ncfile)                                           12d13s22
       pass(1)=dfloat(ncore)                                             1d31s21
       pass(2)=dfloat(maxbct)                                            1d31s21
       pass(3)=dfloat(nmn)                                              8d12s22
      else                                                              5d25s21
       is=1                                                             5d25s21
       ie=1                                                             5d25s21
      end if                                                            1d31s21
      call dws_bcasta(pass,3)                                           8d12s22
      ncore=nint(pass(1))                                                1d31s21
      maxbct=nint(pass(2))                                               1d31s21
      nmn=nint(pass(3))                                                 8d12s22
      call dws_init(ncore,nmn,bc,ibc)                                   11d15s22
      nbytes=(loc(first(2))-loc(first(1)))
      nmove=(loc(last)-loc(first(1)))/nbytes
      nmove=nmove+1                                                     3d1s10
      call dws_bcast(first,nmove)                                       3d1s10
      if(irtrn.eq.2)then                                                2d6s23
       if(mynowprog.eq.0)then                                           2d6s23
        call addcomma4(ihighwtr,pwd,0)                                  2d6s23
        write(6,*)('high water mark '),pwd                              2d6s23
        write(6,*)('double precision words on this core')               6d7s23
       end if                                                           2d6s23
       return                                                           2d6s23
      end if                                                            2d6s23
      if(irtrn.ne.0)then                                                11d4s19
       if(mynowprog.eq.0)then                                           8d12s22
        write(6,*)('stopping because of input problem!!')               8d12s22
       end if                                                           8d12s22
       call dws_synca                                                   11d18s20
       call dws_finalize                                                11d4s19
       stop                                                             11d4s19
      end if                                                            11d4s19
      dws_me=mynowprog                                                  10d23s14
      dws_np=mynprocg                                                   10d23s14
      ncall=0                                                           5d27s19
      call f3j(0,0,0,0,0,0,2)                                           5d27s19
      ilimcode=1                                                        5d7s19
      inv(1,1,1,1)=1                                                    4d12s18
      inv(2,1,1,1)=1                                                    4d12s18
      iunit=6                                                           11d9s17
c
c     mem_dws is millions of bytes summed over all procs
c     columbia has 121.23 Gb/64 proc = 1940 Mb/proc=242 MW/proc
c     maxbc is words in common.store
c     mcode is words of code
c     ivy bridge nodes: 20 core, 64 GB per node = 3.2 GB per core = 429 MW per core
c     skylake nodes: 40 core, 192 GB per node = 4.8 GB per core = 644 MW per core
c
      if(mynowprog.eq.0)then
       write(6,*)('going after "'),pwd(is:ie)//'/rystabc.ufort',('"')
       open(unit=1,file=pwd(is:ie)//'/rystabc.ufort',form='unformatted')              6d22s12
       read(1)nqx,nwrys,mxrysi                                          7d5s12
       write(6,*)('maximum order of rys quadrature '),nqx               6d22s12
       ibcrys=ibcoff                                                    6d22s12
       ibcoff=ibcoff+nwrys+2*nqx                                        6d22s12
       call enough('start.  1',bc,ibc)
       read(1)(ibc(ibcrys+j),j=0,2*nqx-1)                               6d22s12
       read(1)(bc(ibcrys+2*nqx+j),j=0,nwrys-1)                          6d22s12
       close(unit=1)                                                    6d22s12
       pass(1)=dfloat(nqx)                                              6d22s12
       pass(2)=dfloat(nwrys)                                            7d5s12
       pass(3)=dfloat(ibcrys)                                           7d5s12
       pass(4)=dfloat(mxrysi)                                           7d5s12
       write(6,*)('going after "'),pwd(is:ie)//'/ghtab.ufort',('"')     1d19s23
       open(unit=1,file=pwd(is:ie)//'/ghtab.ufort',form='unformatted')  1d19s23
       read(1)nghx,nmqh                                                 1d19s23
       write(6,*)('maximum order of gauss hermite quadrature: '),nghx   1d19s23
       write(6,*)('number of words of storage: '),nmqh                  1d19s23
       ibcgh=ibcoff                                                     1d19s23
       ibcoff=ibcgh+nmqh                                                1d19s23
       call enough('start. 1gh',bc,ibc)                                 1d19s23
       read(1)(bc(ibcgh+j),j=0,nmqh-1)                                  1d19s23
       pass(5)=dfloat(nghx)                                             1d19s23
       pass(6)=dfloat(nmqh)                                             1d19s23
      end if
      call dws_bcast(pass,6)                                            6d22s12
      nqx=nint(pass(1))                                                 6d22s12
      nwrys=nint(pass(2))                                               6d22s12
      mxrysi=nint(pass(4))                                              7d5s12
      ntot=nwrys+2*nqx                                                  6d22s12
      nghx=nint(pass(5))                                                1d19s23
      nmqh=nint(pass(6))                                                1d19s23
      if(mynowprog.ne.0)then                                            11d1s22
       ibcrys=ibcoff                                                    11d1s22
       ibcgh=ibcrys+ntot                                                1d19s23
       ibcoff=ibcgh+nmqh                                                1d19s23
       call enough('runit.1',bc,ibc)
      end if                                                            11d1s22
      ntot=ntot+nmqh                                                    1d19s23
      call dws_bcast(bc(ibcrys),ntot)                                   6d22s12
      do i=1,nqx                                                        11d1s22
       iorg=ibc(ibcrys+(i-1)*2)                                         11d1s22
       ito=iorg+ibcrys+2*nqx-1                                          11d1s22
       ibc(ibcrys+(i-1)*2)=ito                                          11d9s22
      end do                                                            11d1s22
      nproc=mynprocg
      if(mynowprog.eq.0)then
       write(6,446)mynprocg                                             12d2s22
  446  format(' number of compute processors = ',i8)                    12d2s22
       gbyte=dfloat(maxbct)/1024d0                                      1d31s21
       gbyte=gbyte/1024d0
       gbyte=gbyte/1024d0
       gbyte=gbyte*8d0
       write(6,447)maxbc,gbyte                                          12d2s22
  447  format(' scratch memory space: ',i10,(' words or '),f5.1,x,       12d2s22
     $      'Gbytes.')                                                  12d2s22
      end if                                                            2d16s10
      if(nextradatx.gt.0)then                                           1d10s19
       if(mynowprog.ne.0)then                                           1d10s19
        iextradatad=ibcoff                                                1d10s19
        ibcoff=iextradatad+nextradatx                                   1d10s19
        call enough('start.  2',bc,ibc)
       else                                                             1d10s19
       end if                                                           1d10s19
       call dws_bcast(bc(iextradatad),nextradatx)                       1d10s19
      else                                                              1d10s19
       iextradatad=ibcoff                                               1d10s19
      end if                                                            1d10s19
      nmove=(loc(casspadd(2))-loc(iacto(1)))/nbytes                     8d7s14
      call dws_bcast(iacto,nmove)                                       8d11s14
      nmove=idprt/2                                                     1d3s20
      call dws_bcast(iprtr,nmove)                                       1d3s20
      nmove=(loc(xmrci(1))-loc(nhole(1)))/nbytes                        5d24s18
      nmove=nmove+1                                                     5d24s18
      call dws_bcast(nhole(1),nmove)
      potdws=first(1)                                                   2d9s12
      nbasall=nbasis                                                    3d24s10
      nsymb=nsymbb                                                      3d3s10
      call dws_bcast(nbasdws,8)                                         3d29s21
      nbasdwst=nbasdws(1)                                               1d17s20
      do i=2,nsymb                                                      1d17s20
       nbasdwst=nbasdwst+nbasdws(i)                                     1d17s20
      end do                                                            1d17s20
      if(nbasdwst.eq.0)then                                             1d17s20
       do i=1,nsymb                                                      3d3s10
        nbasdws(i)=nbasb(i)                                              3d3s10
       end do                                                            3d3s10
      end if                                                            1d17s20
      call spher(lmax,bc,ibc)                                           11d9s22
      iusecart=iusecartr                                                2d21s20
      call dws_bcast(xcart,3*natom)                                     2d19s10
      call dws_bcast(atnum,3*natom)                                     8d17s15
      nbdat=ngaus*9+nwcont                                              9d12s17
      if(mynowprog.ne.0)then                                            2d19s10
       ibdat=ibcoff                                                     2d19s10
       ibcoff=ibdat+nbdat                                               1d28s19
       ibstor=ibcoff                                                    5d3s10
c
c     see note in basisz                                                2d14s19
c
       isstor=ibstor+nbasis*3                                           2d14s19
       ibcoff=isstor+nbasis*3                                           2d14s19
       call enough('start.  3',bc,ibc)
      end if                                                            2d19s10
c
c     see note in basisz                                                2d14s19
c
      nbdat=nbdat+6*nbasis                                              2d14s19
      call dws_bcast(bc(ibdat),nbdat)                                   1d28s19
      if(idogrado1(4).ne.0)then                                         7d28s22
       igopt=idogrado1(4)                                               7d28s22
       idata=0                                                          7d28s22
       call gopt(igopt,idata,nunique,nder,bc,ibc)                       3d13s23
      end if                                                            7d28s22
      idenergy=ibcoff                                                   7d28s22
      ibcoff=idenergy+3*natom                                           7d28s22
      medws=mynowprog
      itoomo=0                                                          4d16s08
      isend1=ibcoff                                                     6d1s18
      isend2=isend1+mynprocg                                            6d1s18
      isend3=isend2+mynprocg                                            6d1s18
      isend4=isend3+mynprocg                                            6d1s18
      ibcoff=isend4+mynprocg                                            6d1s18
      if(idohf.ne.0)then                                                8d8s14
       call parahf(medws,nproc,jmats,kmats,mgmat,                        3d19s12
     $     nocc,i3x,itmax,ieigr,ivecr,ehf,ih0ao,maxdws,ispt1,icore,     7d26s16
     $     nbasis,isymdws,myguess,natom,ngaus,ibdat,isym,iapair,        5d3s10
     $     multh,ibc(ibstor),ibc(isstor),potdws,ioooo,ionex,iwerner,    7d9s13
     $     idavopt,0,nlzz,idorel,ascale,idogrado,ipropsym,inocan,       10d27s17
     $     iapairg,dynw,nbasb,nbasc,nwcont,lmax,xcart,numeminus,atnum,  1d10s19
     $     bc(iextradatad),nextradatx,idorel4c,smallest,bc(idenergy),   8d29s22
     $     iftype,fstgth,ndfld,bc,ibc,0)                                12d6s23
       if(idomp2.ne.0)then                                               8d8s14
        call paramp2(kmats,nocc,ieigr,mynowprog,ehf,nproc,icore,nbasb,  8d9s24
     $      bc,ibc)                                                     8d9s24
       end if                                                            8d8s14
      end if                                                            8d8s14
      if(idocas.ne.0)then                                               8d8s14
       ibcprior=ibcoff                                                  11d6s19
       call parahf(medws,nproc,jmats,kmats,mgmat,                        3d19s12
     $     nocc,i3x,itmax,ieigr,ivecr,ehf,ih0ao,maxdws,ispt1,icore,     7d26s16
     $     nbasis,isymdws,myguess,natom,ngaus,ibdat,isym,iapair,        5d3s10
     $     multh,ibc(ibstor),ibc(isstor),potdws,ioooo,ionex,iwerner,    7d9s13
     $     idavopt,1,nlzz,idorel,ascale,idogrado,ipropsym,inocan,       10d27s17
     $     iapairg,dynw,nbasb,nbasc,nwcont,lmax,xcart,numeminus,atnum,  1d10s19
     $     bc(iextradatad),nextradatx,idorel4c,smallest,bc(idenergy),   8d29s22
     $      iftype,fstgth,ndfld,bc,ibc,0)                               12d6s23
       ibcoff=ibcprior                                                  11d6s19
      end if                                                            8d8s14
      if(idomrci.ne.0)then                                              5d24s18
       call mrci(ioooo,ionex,jmats,kmats,i3x,ascale,idorel,natom,ngaus, 5d25s18
     $     ibdat,isym,iapair,multh,ibc(ibstor),ibc(isstor),             4d22s21
     $     iapairg,potdws,ibc(isend1),ibc(isend2),ibc(isend3),          6d1s18
     $     ibc(isend4),nbasb,nbasc,iftype,fstgth,ndfld,bc,ibc,cfile,    12d6s23
     $     ncfile,ipropsym)                                             10d30s24
      end if                                                            5d24s18
      if(idogrado1(4).ne.0)then                                         3d13s23
       icasguess=3                                                      3d13s23
       call move2min(bc,ibc,idata,nunique,cfile,ehf,i3x,iapairg,ibdat,  3d13s23
     $     ibstor,idenergy,ieigr,ih0ao,ionex,ioooo,ispt1,isstor,isymdws,3d13s23
     $     ivecr,jmats,kmats,maxdws,mgmat,ncfile,potdws,nder)           3d13s23
      end if                                                            3d13s23
      if(nfname.gt.0)then                                               11d22s19
       nspc=2*11                                                        8d31s21
       nwavedat=nfneed*nspc                                             12d5s22
       iwavedat=ibcoff                                                  11d24s19
       ibcoff=iwavedat+nwavedat/2                                       5d12s21
       call enough('start.  4',bc,ibc)
       call prop(nfname,ifname,ibc(iwavedat),nspc,nsymb,nbasc,natom,    12d16s19
     $      ngaus,ibdat,isym,iapair,ibc(ibstor),ibc(isstor),idorel,     12d16s19
     $      ascale,multh,nbasb,nwavedat,pwd(is:ie),iprop,bc,ibc)        11d9s22
      end if                                                            11d22s19
      if(modecc.ne.0)then                                               3d22s12
      end if                                                            3d22s12
      if(mynowprog.eq.0)then
       call addcomma4(ihighwtr,pwd,0)                                   2d11s21
       write(6,*)('high water mark '),pwd                               6d7s23
       write(6,*)('double precision words on this core')                6d7s23
      end if
      return                                                            11d15s22
      end
