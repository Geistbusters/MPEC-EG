c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mrci(ioooo,ionex,jmats,kmats,i3x,ascale,idorel,natom,  5d25s18
     $     ngaus,ibdat,isym,iapair,multh,ibstor,isstor,iapairg,         4d22s21
     $     potdws,isend1,isend2,isend3,isend4,nbasisp,nbasisc,iftype,   8d29s22
     $     fstgth,ndfld,bc,ibc,cfile,ncfile,ipropsym)                   10d30s24
      implicit real*8 (a-h,o-z)                                         5d24s18
      external second                                                   12d14s18
      logical lprint                                                    5d24s18
c
c
c     mrci calculations                                            5d24s18
c
      parameter (idroot=100)                                            2d26s19
      parameter (idsln=1000)                                            11d21s22
      character*19 ostng(idroot)                                        2d26s19
      character*6 termsym                                               1d5s20
      character*2 spinm                                                 1d5s20
      character*1 ccode(3)                                              12d7s21
      character*50 numbers,numbers2                                     2d3s21
      character*(idsln) sline                                           11d21s22
      character*(*) cfile                                               12d13s22
      include "common.store"                                            6d1s18
      include "common.mrci"                                             5d24s18
      include "common.hf"                                               5d24s18
      include "common.print"                                            1d3s20
      integer*8 i0a,i0b,ndoubx,i18,i28,in8,ncupc8,ipack8,ihi,igxx,      4d14s21
     $     icsfxv,itop8,ibot8                                           4d14s21
      integer*4 ipack4(2)                                               10d7s24
      integer*2 ipack2(2)                                               10d7s24
      integer isend1(*),isend2(*),isend3(*),isend4(*)                   10d7s24
      logical l4,l3,l2,l1,ldebug,bintvec,ldebugf,ldavid,lconv           6d20s24
      equivalence (ipack8,ipack4)                                       12d8s20
      equivalence (ipack4,ipack2)                                       10d7s24
      dimension ioooo(*),ionex(*),jmats(*),kmats(*),i3x(*),isym(3,8),   5d25s18
     $     iapair(3,*),ibstor(*),isstor(*),iapairg(*),nbasdwsc(8),      5d29s18
     $     iorb(8),noc(8),iptoh(8,8,8),ivirtc(8),ih0a(8),nrules(6,8),   3d4s19
     $     multh(8,8),ncs(2,8),csfs(2,2),ieigm(2,8),ivm(2,2,8),         6d28s18
     $     ndetm(2,2,8),nosing(2,8),ih0av(8),nh0av(8),nra(6),nrb(6),    3d4s19
     $     ngxs(2),ncd(3,8,2),nvxv(2,8),csfd(3,3),ieigmm(3,8,2),        10d1s18
     $     ivmm(3,8,2),ndetmm(3,3,8,2),idoit(6),nvcul(8),timex(9),      3d9s21
     $     ncdx(3,8,2),ncsx(2,8),nbasisp(*),nbasisc(*),avgsptot(4),     5d7s19
     $     wgtsptot(4),ivmt(2,2,8),nrule(8),ivintref(8),irwt(8),        5d28s19
     $     isorts(2,8),ivmcsf(2,8),iphs(3,8,2),nfcn(8),ibasis(8),       7d15s19
     $     ibasis1(8),nfcn1(8),ibasis2(8),nfcn2(8),nct(8),nct1(8),      7d15s19
     $     nct2(8),nct2c(4,8),i3xb(512),ixlzz(8,6),islz(3),ionexb(512), 4d8s20
     $     ih0ae(8),ibasisf(8),nctf(8),nfcnf(8),ih0ave(8),ih0vve(8),    4d16s20
     $     ih0aek(8),nbasis(8),nbasis1(8),nbasisf(8),nbasis2(8),        4d24s20
     $     jmatt(512),kmatt(512),ionexc(512),ibasis3(8),nbasis3(8),     7d20s20
     $    nfcn3(8),nct3(8),kmatd(8,8),ibasisfi(8),nfcnfi(8),nbasisfi(8),10d8s20
     $    ibasis1i(8),nfcn1i(8),nctfi(8),nct1i(8),nbasis1i(8),nfill0(2),10d8s20
     $     nhole0(3),ibasis2i(8),nfcn2i(8),nct2i(8),nbasis2i(8),        11d18s20
     $     nfdat(5,4,8),ieigvspace(8),i3x3(8),tdendd(7),tdends(8),      3d23s21
     $     ionexbt(512),jvintref(8),jeigvspace(8),nameciu(6),           6d11s21
     $     i4x(512),i4xb(512),iden1e(4,8),idenhvv(4,8,8),idenj(4,8,8,8),11d7s22
     $     idenk(4,4,8,8,8),iff(8),idvi(8),idviu(8),fstgth(*),iftype(*),8d7s24
     $     idervcv(2,8)                                                 8d7s24
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/kmfind/invk1(2,8,8,8,2)
      common/fnd2cm/inv(2,8,8,8)                                        9d2s20
      common/fnd4xcm/inv4x(2,8,8,8)                                     7d4s21
      common/tpartcm/tpart(10)
      data ccode/' ','c','s'/                                           12d7s21
      data idoit/6*1/                                                   10d18s18
      data nfill0/2*0/                                                  10d8s20
      data nhole0/3*0/                                                  10d8s20
      common/ilimcm/ilimcode                                            5d7s19
      common/hcdscm/timds(12)
      common/hcsicm/timsi(10)
      common/wxcm/ixw1,ixw2                                             5d28s20
      igoal=1840251
      sr2=sqrt(2d0)                                                     1d12s21
      srh=sqrt(0.5d0)                                                   1d14s21
      navg=10                                                           3d1s21
      do iavg=0,navg                                                    3d1s21
       call second(bc(ibcoff+iavg))                                     3d1s21
      end do                                                            3d1s21
      call avgr(bc(ibcoff),navg,tovr)                                   3d1s21
      timds=0d0
      timds(1)=tovr                                                     3d18s21
      timsi=0d0                                                         4d16s20
      ibcoffo=ibcoff                                                    11d23s19
      idvmt=ibcoff                                                      9d1s23
      ibcoff=idvmt+nsymb*2                                              9d1s23
      if(iprtr(4).eq.0)then                                             1d3s20
       ldebug=.false.                                                   1d3s20
      else                                                              1d3s20
       ldebug=.true.                                                    1d3s20
      end if                                                            1d3s20
      ldebugf=.false.                                                   5d5s21
      bintvec=.true.                                                    1d17s20
      intvec=0                                                          1d17s20
      lprint=mynowprog.eq.0                                             5d24s18
      ilimcode=2                                                        5d7s19
      do i=1,4                                                          5d7s19
       avgsptot(i)=0d0                                                  5d7s19
       wgtsptot(i)=0d0                                                  5d7s19
      end do                                                            5d7s19
      if(idorel.eq.0)then                                               8d20s15
       ncomp=1                                                          8d20s15
      else                                                              8d20s15
       ncomp=2                                                          8d20s15
      end if                                                            8d20s15
      do isb=1,nsymb                                                    5d25s18
       nbasdws(isb)=nbasisp(isb)                                        11d6s19
       nbasdwsc(isb)=nbasdws(isb)*ncomp                                 5d25s18
      end do                                                            5d25s18
      idwsdeb=0                                                         2d19s19
      iwopen=0                                                          8d11s22
      if(lprint.and.nwavef.ne.0)then                                    3d23s21
       iwopen=1                                                         8d11s22
       open(unit=2,file='wavef',form='unformatted')
       nwavrec=1                                                        5d3s21
       if(nrestart.lt.0)then                                            8d12s22
        read(2)ismultr,isymmrcir,nrootr,norbr,nsymbr,nameciu            8d11s22
        ier=0                                                           8d11s22
        if(ismultr.ne.ismult)ier=ier+1                                  8d11s22
        if(isymmrcir.ne.isymmrci)ier=ier+1                              8d11s22
        if(nrootr.ne.nroot)ier=ier+1                                    8d11s22
        if(norbr.ne.norb)ier=ier+1                                      8d11s22
        if(nsymbr.ne.nsymb)ier=ier+1                                    8d11s22
        do i=1,6                                                        8d11s22
         if(nameciu(i).ne.nameci(i))ier=ier+1                           8d11s22
        end do                                                          8d11s22
        if(ier.gt.0)then                                                8d11s22
         write(6,*)('restart data differs for record '),nwavrec         8d16s22
         write(6,*)('gota  '),ismultr,isymmrcir,nrootr,norbr,nsymbr,     8d11s22
     $        nameciu                                                   8d11s22
         write(6,*)('want '),ismult,isymmrci,nroot,norb,nsymb,nameci    8d11s22
         stop 'restart'                                                 8d11s22
        end if                                                          8d11s22
       else                                                             8d11s22
        write(2)ismult,isymmrci,nroot,norb,nsymb,nameci                  12d27s19
       end if                                                           8d11s22
       nwavrec=nwavrec+1
       if(nrestart.lt.0)then                                            8d12s22
        read(2)(ibasisfi(i),nfcnfi(i),nbasisfi(i),i=1,nsymb)            8d11s22
        do i=1,nsymb                                                    8d11s22
         if(idoubo(i).ne.ibasisfi(i))ier=ier+1                          8d11s22
         if(irefo(i).ne.nfcnfi(i))ier=ier+1                             8d11s22
         if(nbasdws(i).ne.nbasisfi(i))ier=ier+1                         8d11s22
        end do                                                          8d11s22
        if(ier.ne.0)then                                                8d11s22
         write(6,*)('restart data differs for record '),nwavrec         8d16s22
         write(6,*)('gotb  '),(ibasisfi(i),nfcnfi(i),nbasisfi(i),        8d11s22
     $        i=1,nsymb)                                                8d11s22
         write(6,*)('want '),(idoubo(i),irefo(i),nbasdws(i),i=1,nsymb)  8d11s22
         stop 'restart'                                                 8d11s22
        end if                                                          8d11s22
       else                                                             8d11s22
        write(2)(idoubo(i),irefo(i),nbasdws(i),i=1,nsymb)                11d25s19
       end if                                                           8d11s22
       nwavrec=nwavrec+1
       if(nrestart.lt.0)then                                            8d12s22
        read(2)(i4x(i),i4xb(i),i=1,norb)                                8d11s22
        do i=1,norb                                                     8d11s22
         if(irel(i).ne.i4x(i))ier=ier+1                                 8d11s22
         if(ism(i).ne.i4xb(i))ier=ier+1                                 8d11s22
        end do                                                          8d11s22
        if(ier.ne.0)then                                                8d11s22
         write(6,*)('restart data differs for record '),nwavrec         8d16s22
         write(6,*)('gotc  '),(i4x(i),i4xb(i),i=1,norb)                  8d11s22
         write(6,*)('want '),(irel(i),ism(i),i=1,norb)                  8d11s22
         stop 'restart'                                                 8d11s22
        end if                                                          8d11s22
       else                                                             8d11s22
        write(2)(irel(i),ism(i),i=1,norb)                                11d25s19
       end if                                                           8d11s22
      end if                                                            2d22s19
      ipassr=ibcoff                                                     1d2s20
      ibcoff=ipassr+nsymb                                               1d2s20
      call cal1int(ih0,iorb,idorel,ascale,nbasdwsc,natom,ngaus,ibdat,    2d15s19
     $     isym,iapair,ibstor,isstor,idwsdeb,mynowprog,iapairg,         4d22s21
     $     nbasisp,ibc(ipassr),nextrad,iextrad,nwavef,nwavrec,nder,     8d11s22
     $     nrestart,iftype,fstgth,ndfld,lambdaci,bc,ibc,idarot)         12d6s23
      if(nder.ne.0)then                                                 7d11s23
       do isb=1,nsymb                                                    7d11s23
        iff(isb)=ibcoff                                                  7d11s23
        ibcoff=iff(isb)+nbasdws(isb)*nbasdws(isb)                        7d11s23
       end do                                                            7d11s23
       call enough('mrci.iff',bc,ibc)                                   7d11s23
       do iz=iff(1),ibcoff-1                                             7d11s23
        bc(iz)=0d0                                                       7d11s23
       end do                                                            7d11s23
      end if                                                            7d11s23
      ih0copy=ibcoff                                                    1d31s21
      jh0copy=ih0copy                                                   1d31s21
      jh0=ih0                                                           1d31s21
      do isb=1,nsymb                                                    1d31s21
       nbasdwsf=nbasisp(isb)*ncomp                                      1d31s21
       do i=0,nbasdwsf*nbasdwsf-1                                       1d31s21
        bc(jh0copy+i)=bc(jh0+i)                                         1d31s21
       end do                                                           1d31s21
       jh0copy=jh0copy+nbasdwsf*nbasdwsf                                1d31s21
       jh0=jh0+nbasdwsf*nbasdwsf                                        1d31s21
      end do                                                            1d31s21
      ibcoff=jh0copy                                                    1d31s21
      if(itestmrci.eq.1.or.itestmrci.eq.2)then                          5d15s19
       ii=0                                                             5d15s19
       do isb=1,nsymb                                                   5d15s19
        nn=nbasdws(isb)-idoubo(isb)                                     5d15s19
        do i=1,nn                                                       5d15s19
         ii=ii+1                                                        5d15s19
         if(ii.gt.ido)then                                              5d15s19
          write(6,*)('no. orbs in reference space too large ')          5d15s19
          write(6,*)('limit is ido = '),ido                             5d15s19
          stop                                                          5d15s19
         end if                                                         5d15s19
         ism(ii)=isb                                                    5d15s19
         irel(ii)=i                                                     5d15s19
        end do                                                          5d15s19
       end do                                                           5d15s19
       norb=ii                                                          5d15s19
       if(lprint.and.itestmrci.eq.2)then                                5d15s19
        open(unit=7,file='hdata')                                       5d15s19
        write(7,*)ismult
        write(7,*)'csfhami'                                             5d15s19
        write(7,*)norb                                                  5d15s19
        do i=1,norb                                                     5d15s19
         write(7,*)i,irel(i),ism(i)                                     5d15s19
        end do                                                          5d15s19
       end if                                                           5d15s19
       ioffo=1                                                          5d15s19
       do isb=1,nsymb                                                   5d15s19
        no=idoubo(isb)+irefo(isb)                                       5d15s19
        nvcul(isb)=irefo(isb)                                           5d15s19
        joffo=ioffo+irefo(isb)                                          5d15s19
        do io=no,nbasdws(isb)-1                                         5d15s19
         irefo(isb)=irefo(isb)+1                                        5d15s19
c
c singles and doubles
c
         nfill(2)=nfill(2)+1                                            5d15s19
         ifill(nfill(2),2)=joffo                                        5d15s19
c
c singles only
c
         joffo=joffo+1                                                  5d15s19
        end do                                                          5d15s19
        ioffo=ioffo+nbasdws(isb)-idoubo(isb)                            5d15s19
       end do                                                           5d15s19
      end if                                                            5d15s19
      ndoub=0                                                           5d24s18
      do isb=1,nsymb                                                    5d24s18
       noc(isb)=idoubo(isb)+irefo(isb)                                  5d29s18
       ndoub=ndoub+idoubo(isb)                                          5d24s18
       nvirt(isb)=nbasdws(isb)-idoubo(isb)-irefo(isb)                   5d24s18
      end do                                                            5d24s18
      nec=ne-ndoub*2                                                    6d1s18
      if(lprint)then                                                    5d24s18
       write(6,*)(' ')
       write(6,*)('MRCI')                                               5d24s18
       if(ndoub.gt.0)write(6,860)(idoubo(isb),isb=1,nsymb)              10d16s20
  860  format('   doubly occupied orbitals: ',8i5)                      10d16s20
       write(6,861)(irefo(isb),isb=1,nsymb)                             10d16s20
  861  format('orbitals in reference space: ',8i5)                      10d16s20
       write(6,862)(nvirt(isb),isb=1,nsymb)                             10d16s20
  862  format('           virtual orbitals: ',8i5)                      10d16s20
       do if=1,3                                                        5d24s18
        if(nhole(if).gt.0)then                                          5d24s18
         ifm=if-1                                                       5d24s18
         write(6,1)ifm,(irel(ihole(i,if)),ism(ihole(i,if)),             5d25s18
     $        i=1,nhole(if))                                            5d25s18
   1     format(' orbs with ',i1,' holes in ref. space: ',12(i4,'s',i1)) 5d25s18
        end if                                                          5d24s18
       end do                                                           5d24s18
       do if=1,4                                                        4d19s23
        if(nfill(if).gt.0)then                                          5d24s18
         write(6,*)(ifill(i,if),i=1,nfill(if))
         write(6,2)if,(irel(ifill(i,if)),ism(ifill(i,if)),              7d2s18
     $       i=1,nfill(if))                                             5d25s18
    2    format(' orbs with at most ',i1,' e- in ref. space: ',          5d25s18
     $       12(i4,'s',i1))                                             5d25s18
        end if                                                          5d24s18
       end do                                                           5d24s18
       write(6,*)('total number of electrons: '),ne                     5d24s18
       write(6,*)('number of electrons to correlate '),nec              5d24s18
       write(6,*)('spin multiplicity: '),ismult                         5d24s18
       ntop=0                                                           12d24s19
       do i=1,6                                                         12d24s19
        if(nameci(i).ne.0)then                                          1d5s20
         ntop=i                                                         1d5s20
         termsym(i:i)=char(nameci(i))                                   1d5s20
        else                                                            1d5s20
         termsym(i:i)=' '                                               1d5s20
        end if                                                          1d5s20
       end do                                                           12d24s19
       write(spinm,1066)ismult                                          1d5s20
 1066  format(i2)                                                       1d5s20
       write(6,*)('symmetry of state: '),isymmrci,                      12d24s19
     $      ('('),spinm,(char(nameci(i)),i=1,ntop),(')')                1d5s20
       write(6,*)('number of roots to extract: '),nroot                 5d24s18
       write(6,*)(' ')                                                  3d9s21
       write(6,305)econvci                                              3d9s21
  305  format('energy convergence criterion for external iterations: ', 3d9s21
     $      es8.1)                                                      3d9s21
       if(nroot.gt.1)write(6,306)econvcii                               5d23s23
  306  format('energy convergence criterion for internal',              5d23s23
     $      ' iterations: ',es8.1)                                      3d11s21
      end if                                                            5d24s18
      idwsdeb=0                                                         8d24s18
      nsingx=0
      nrxx=0                                                            5d7s21
      do isb=1,nsymb                                                    5d7s21
       nrxx=max(nrxx,irxinfo(isb))                                      5d7s21
      end do                                                            5d7s21
      if(nrxx.gt.0.and.nrxinfos.eq.0)then                               5d7s21
       write(6,*)('adding to rxinfos ')
       do isb=1,nsymb                                                   5d7s21
        if(irxinfo(isb).gt.0)then                                       5d7s21
         nrxinfos=nrxinfos+1                                            5d7s21
         irxinfos(1,nrxinfos)=isb                                       5d7s21
         irxinfos(2,nrxinfos)=irxinfo(isb)*nroot                        5d7s21
         irxinfos(3,nrxinfos)=nlzzci                                    5d7s21
         irxinfos(4,nrxinfos)=lambdaci                                  5d7s21
         write(6,*)isb,(irxinfos(j,nrxinfos),j=1,4)
        end if                                                          5d7s21
       end do                                                           5d7s21
      end if                                                            5d7s21
      jdenpt=ibcoff                                                     6d12s19
      ibcoff=jdenpt+idbk                                                6d12s19
      call make1dm(nsymb,noc,nbasdws,lprint,maxddi,jmats,kmats,jmats,   2d15s19
     $   ibc(jdenpt),isblk,isblkk,isblkh,nsdlk,nsdlkk,nsdlkh,multh,idbk,6d12s19
     $     isblk1,nsdlk1,isblkd,nsdlkd,1)                               7d5s21
      if(iunc(2).ne.0.or.nder.ne.0)then                                 7d18s23
       call make4xdm(nsymb,nvirt,isblk4x,n4x,multh)                     7d4s21
      end if                                                            7d4s21
      ibcoff=jdenpt                                                     6d12s19
      nnewj=0                                                           9d2s20
       nnewj=nsdlk                                                       9d2s20
       do i=1,nsdlk                                                      9d2s20
        if(isblk(3,i).ne.isblk(4,i))then                                 9d2s20
         nnewj=nnewj+1                                                   9d2s20
         if(nnewj.gt.idbk)then                                           9d2s20
          write(6,*)('trying to add '),nnewj,                           10d12s20
     $        ('Jmat but dimensions are '),idbk,(' (in common.hf)')     10d12s20
          write(6,*)('we are at i = '),i,(' out of '),nsdlk              9d2s20
          call dws_sync                                                  9d1d20
          call dws_finalize                                              9d1d20
          stop                                                           9d1d20
         end if                                                          9d1d20
         isblk(1,nnewj)=isblk(1,i)                                        9d2s20
         isblk(2,nnewj)=isblk(2,i)                                        9d2s20
         isblk(3,nnewj)=isblk(4,i)                                        9d2s20
         isblk(4,nnewj)=isblk(3,i)                                        9d2s20
        end if                                                           9d2s20
       end do                                                            9d2s20
       nsdlk=nnewj                                                       9d2s20
c
c     inv is set in first call to ifind2
c     call it now so inv set below does not get overwritten.
c
      idum=ifind2(1,1,1,1,icase)                                        9d2s20
      nnew=0
      nnew=nsdlkk                                                       9d1d20
      do i=1,nsdlkk                                                     9d1d20
       if(isblkk(3,i).ne.isblkk(4,i))then                               7d5s18
        nnew=nnew+1                                                     9d1d20
        if(nnew.gt.idbk)then                                            9d1d20
         write(6,*)('trying to add '),nnew,                             10d12s20
     $       (' Kmat but dimensions are '),idbk,(' (in common.hf)')     10d12s20
         write(6,*)('we are at i = '),i,(' out of '),nsdlkk             9d1d20
         call dws_sync                                                  9d1d20
         call dws_finalize                                              9d1d20
         stop                                                           9d1d20
        end if                                                          9d1d20
        isblkk(1,nnew)=isblkk(2,i)                                      9d1d20
        isblkk(2,nnew)=isblkk(1,i)                                      9d1d20
        isblkk(3,nnew)=isblkk(4,i)                                      9d1d20
        isblkk(4,nnew)=isblkk(3,i)                                      9d1d20
       end if                                                           6d30s18
      end do                                                            6d30s18
      nsdlkk=nnew                                                       9d1d20
      do i3=1,nsymb
       do i2=1,nsymb
        do i1=1,nsymb
         invk1(1,i1,i2,i3,1)=1                                          1d12s23
         invk1(1,i1,i2,i3,2)=1                                          1d12s23
         do i=1,nsdlkk
          if(isblkk(1,i).eq.i1.and.isblkk(2,i).eq.i2.and.
     $         isblkk(3,i).eq.i3)then                                   6d28s18
           invk1(1,i1,i2,i3,1)=i                                        6d30s18
           invk1(2,i1,i2,i3,1)=1                                        6d30s18
          else if(isblkk(2,i).eq.i1.and.isblkk(1,i).eq.i2.and.          6d28s18
     $          isblkk(4,i).eq.i3.and.nnew.eq.0)then                    9d1d20
           invk1(1,i1,i2,i3,1)=i                                        6d30s18
           invk1(2,i1,i2,i3,1)=2                                        6d30s18
          end if                                                        6d28s18
         end do                                                         6d28s18
         do i=1,nsdlk1                                                  6d30s18
          if(isblk1(3,i).eq.i3)then                                     6d30s18
           if(isblk1(1,i).eq.i1.and.isblk1(2,i).eq.i2)then              6d30s18
            invk1(1,i1,i2,i3,2)=i                                       6d30s18
            invk1(2,i1,i2,i3,2)=1                                       6d30s18
           else if(isblk1(1,i).eq.i2.and.isblk1(2,i).eq.i1)then         6d30s18
            invk1(1,i1,i2,i3,2)=i                                       6d30s18
            invk1(2,i1,i2,i3,2)=2                                       6d30s18
           end if                                                       6d30s18
          end if                                                        6d30s18
         end do                                                         6d30s18
         inv(1,i1,i2,i3)=1                                              1d19s23
         inv(2,i1,i2,i3)=0                                              9d2s20
         do i=1,nsdlk                                                   9d2s20
          if(isblk(1,i).eq.i1.and.isblk(2,i).eq.i2.and.                 9d2s20
     $        isblk(3,i).eq.i3)then                                     9d2s20
           inv(1,i1,i2,i3)=i                                            9d2s20
           inv(2,i1,i2,i3)=1                                            9d2s20
          else if(isblk(1,i).eq.i2.and.isblk(2,i).eq.i1.and.            9d2s20
     $         isblk(3,i).eq.i3)then                                    9d2s20
           inv(1,i1,i2,i3)=i                                            9d2s20
           inv(2,i1,i2,i3)=4                                            9d2s20
          end if                                                        9d2s20
         end do                                                         9d2s20
        end do
       end do                                                           6d28s18
      end do                                                            6d28s18
      do i1=1,nsymb                                                     7d22s14
       do i2=1,nsymb                                                    7d22s14
        do i3=1,nsymb                                                   7d22s14
         iptoh(i3,i2,i1)=0                                              7d22s14
        end do                                                          7d22s14
       end do                                                           7d22s14
      end do                                                            7d22s14
      do i=1,nsdlkh                                                     5d12s10
       iptoh(isblkh(1,i),isblkh(2,i),isblkh(3,i))=i                     5d12s10
      end do                                                            5d12s10
      if((itestmrci.eq.1.or.itestmrci.eq.2).and.bc(132).ne.-132d0)then  11d4s19
      write(6,*)('zeroing 2e ints with 4 virts: '),nvcul,nsymb
      do isb1=1,nsymb
       do isb2=1,nsymb
        isb12=multh(isb1,isb2)
        do isb3=1,nsymb
         isb4=multh(isb3,isb12)
          i2eu=ifind2(isb4,isb3,isb2,isb1,icase)
          write(6,*)isb4,isb3,isb2,isb1,i2eu,icase
          if(i2eu.gt.0.and.icase.eq.1)then
           i2eu=ioooo(i2eu)
           if(isb4.eq.isb3)then
            ncol=(irefo(isb1)*(irefo(isb2)+1))/2                        10d23s18
            nrow=(irefo(isb3)*(irefo(isb4)+1))/2
            write(6,*)('starting matrix: ')
            call prntm2(bc(i2eu),nrow,ncol,nrow)
            do i1=nvcul(isb1),irefo(isb1)-1
             iad1=i2eu+((i1*(i1+1))/2)*nrow                             10d23s18
             do i2=nvcul(isb2),i1                                       10d23s18
              iad2=iad1+i2*nrow
              do i3=nvcul(isb3),irefo(isb3)-1
               iad3=((i3*(i3+1))/2)+iad2
               do i4=nvcul(isb4),i3
                bc(iad3+i4)=0d0
               end do
              end do
             end do
            end do
            write(6,*)('ending matrix: ')
            call prntm2(bc(i2eu),nrow,ncol,nrow)
           else
            ncol=irefo(isb1)*irefo(isb2)                                 10d23s18
            write(6,*)('starting matrix: ')
            call prntm2(bc(i2eu),nrow,ncol,nrow)
            nrow=irefo(isb3)*irefo(isb4)
            do i1=nvcul(isb1),irefo(isb1)-1
             iad1=i2eu+i1*irefo(isb2)*irefo(isb3)*irefo(isb4)
             do i2=nvcul(isb2),irefo(isb2)-1
              iad2=iad1+i2*irefo(isb4)*irefo(isb3)
              do i3=nvcul(isb3),irefo(isb3)-1
               iad3=iad2+i3*irefo(isb4)
               do i4=nvcul(isb4),irefo(isb4)-1
                bc(iad3+i4)=0d0
               end do
              end do
             end do
            end do
            write(6,*)('ending matrix: ')
            call prntm2(bc(i2eu),nrow,ncol,nrow)
           end if
          end if
        end do
       end do
      end do
      end if                                                            8d5s19
c
c     count space for pointers for det maps
c
      norb=0                                                            6d1s18
      do isb=1,nsymb                                                    6d1s18
       norb=norb+irefo(isb)                                             6d1s18
      end do                                                            6d1s18
      n2s=ismult-1                                                      6d1s18
      nbeta=(nec-n2s)/2                                                 6d1s18
      nalpha=nbeta+n2s                                                  6d1s18
      ntest=nalpha+nbeta                                                6d1s18
      if(ntest.ne.nec)then                                              6d1s18
       write(6,*)('n2s inconsistent with nec '),n2s,nalpha,nbeta,ntest  11d11s18
       call dws_sync                                                    6d1s18
       call dws_finalize                                                6d1s18
       stop                                                             6d1s18
      end if                                                            6d1s18
      nbetap=nbeta+1                                                    6d1s18
      ndetub=0                                                          6d1s18
      natot=0                                                           6d1s18
      idoe=max(nalpha,nbetap)                                           6d4s18
      mxh=0                                                             6d4s18
      do ie=1,idoe                                                      6d4s18
       if1=ie
       if2=norb-ie
       ifx=max(if1,if2)
       ifn=min(if1,if2)
       if(ifn.eq.0)then
        nhere=1
       else
        itop8=ifx+1                                                     2d26s21
        in=itop8                                                        2d26s21
        ibot8=1                                                         2d26s21
        do j=ifx+2,norb
         itop8=itop8*j                                                  2d26s21
        end do
        do j=1,ifn
         ibot8=ibot8*j                                                  2d26s21
        end do
        nhere=itop8/ibot8                                               2d26s21
        ipack8=nhere*ibot8                                                 2d26s21
        irem=itop8-ipack8                                               2d26s21
        if(irem.ne.0)then
         write(6,*)('itop,ibot '),itop8,ibot8,nhere,ipack8              2d26s21
         write(6,*)('remainder is not zero! '),irem
         write(6,*)norb,ie
         call dws_sync
         call dws_finalize
         stop
        end if
       end if                                                           5d17s18
       natot=natot+nhere*ie                                             6d1s18
       ndetub=ndetub+nhere                                              6d1s18
       mxh=max(mxh,nhere)                                               6d4s18
      end do                                                            6d1s18
      spin=0.5d0*(dfloat(ismult)-1d0)                                   2d19s19
      isref=nint(2d0*spin)                                              2d19s19
      icsf=ibcoff                                                       6d5s19
      idet=icsf+norb+1                                                  11d19s20
      ibcoff=idet+norb+1                                                11d19s20
      call second(time1)
      call gencsf1(nec,spin,norb,nos,nod,mdon,mdoo,ibc(icsf),maxopens,  6d12s19
     $     ixsoo,3,ibc(idet),nwavef,nfill,ifill,ido)                    8d28s24
      nod=max(nod,1)                                                    1d3s20
      mdoop=mdoo+1                                                      5d7s20
      call puttoiprop(ismult,nec,mdon,mdoop)                            8d16s22
      xms=0.5d0*(nalpha-nbeta)                                          6d12s19
      call second(time1)
      nfcnpx=ibcoff                                                     9d20s20
      ibcoff=nfcnpx+4*nsymb*mdoop                                       10d19s20
      idorb=ibcoff
      nod8=2*nsymb*nod/8                                                1d2s20
      if(nod8*8.lt.nod*nsymb)nod8=nod8+1                                7d16s19
      isorb=idorb+nod8                                                  6d12s19
      nos=nsymb*2*nos*8                                                 7d15s19
      if(nos.gt.10 000 000)then                                            8d29s24
       if(lprint)write(6,*)('nos = '),nos,('is too large ')             8d28s24
       nos=10 000 000                                                      8d29s24
       if(lprint)write(6,*)('reset it to '),nos                         8d28s24
      end if                                                            8d28s24
      nos8=nos/8                                                        6d12s19
      if(nos8*8.lt.nos)nos8=nos8+1                                      6d12s19
      iptr=isorb+nos8                                                   6d12s19
      ibasis(1)=iptr+mdoop*4*nsymb                                      5d7s20
      call enough('mrci.  1',bc,ibc)
      call second(time1)
      maxopens4=maxopens+4                                              10d26s20
      maxopens2=maxopens+2                                              10d26s20
      icsf2=ibcoff                                                      1d22s21
      call gencsf2(nec,spin,norb,ibc(idorb),ibc(isorb),ibc(iptr),       6d4s19
     $     nsymb,mdon,mdoo,ibc(ibasis(1)),multh,ism,ibc(icsf),nfcn,nct, 7d15s19
     $     nod,nos,maxopens,nhole,ihole,nfill,ifill,ido,0,itestmrci,    9d5s19
     $     ismult,idum,idum,idum,idum,idum,idum,idum,idum,interacts,    11d24s19
     $     nodu,nosu,nsymb,idum,idum,idum,idum,idum,idum,ibc(nfcnpx),   10d18s20
     $     idum,idum,ibc(icsf2),1,bc,ibc)                               11d14s22
      need=0                                                            7d15s19
      jbasis=ibasis(1)                                                  7d15s19
      nbss=0                                                            4d17s20
      do isb=1,nsymb                                                    7d15s19
       nhere=(nfcn(isb)*3)+mod(nfcn(isb),2)                             7d15s19
       nbasis(isb)=nbss                                                 4d17s20
       nbss=nbss+nhere                                                  4d17s20
       need=need+nhere                                                  7d15s19
       ibasis(isb)=jbasis                                               7d15s19
       jbasis=jbasis+(nhere/2)                                          7d15s19
      end do                                                            7d15s19
      need8=need/2                                                      7d15s19
      ibcoff=ibasis(1)+need8                                            7d15s19
      iptrbit=ibcoff                                                    5d7s20
      ibcoff=iptrbit+nsymb*mdoop                                        5d7s20
      call int1tobit(ibc(iptr),ibc(iptrbit),ibc(idorb),ibc(isorb),      5d7s20
     $     idorbbit,isorbbit,mdoop,nsymb,nec,bc,ibc)                    11d10s22
      if(mynowprog.eq.0.and.iunc(1).eq.0.and.iunc(2).eq.0.and.          11d24s19
     $     ethred.eq.0d0.and.nwavef.ne.0)then                           3d23s21
       jptr=iptr+2*mdoop*(isymmrci-1)                                   11d24s19
       write(6,*)('dumpbas1 '),nwavrec
       call dumpbas(nfcn(isymmrci),ibc(jptr),ibc(ibasis(isymmrci)),     11d24s19
     $      ibc(idorb),ibc(isorb),nodu,nosu,mdoop,nec,1,nct(isymmrci),  5d12s21
     $      ibc(icsf),mdon,nwavrec,nrestart,bc,ibc)                     11d14s22
       nwavrec=nwavrec+1
       if(nrestart.lt.0)then                                            8d12s22
        read(2)(ibasisfi(i),i=1,2),ethredr                              8d12s22
        do i=1,2                                                        8d11s22
         if(ibasisfi(i).ne.iunc(i))ier=ier+1                            8d12s22
        end do                                                          8d11s22
        if(ethredr.ne.ethred)ier=ier+1                                  8d11s22
        if(ier.ne.0)then                                                8d11s22
         write(6,*)('restart data differs for record '),nwavrec         8d16s22
         write(6,*)('gotd  '),(ibasisfi(i),i=1,2),ethredr                8d12s22
         write(6,*)('want '),iunc,ethred                                8d11s22
         stop 'restart'                                                 8d11s22
        end if                                                          8d11s22
       else                                                             8d11s22
        write(2)iunc,ethred                                              11d24s19
       end if                                                           8d11s22
      end if                                                            11d24s19
      nffr=ibcoff                                                       11d25s20
      iffr=nffr+mdoop*2                                                 11d25s20
      ibcoff=iffr+nfcn(isymmrci)*2                                      11d25s20
      call enough('mrci.  2',bc,ibc)
      jptrbit=iptrbit+mdoop*(isymmrci-1)                                11d26s20
      call genf0(ibc(jptrbit),ibc(ibasis(isymmrci)),nfcn(isymmrci),     11d25s20
     $     mdon,mdoo,ibc(nffr),ibc(iffr),ibc(icsf),bc,ibc)              11d11s22
      nff2=ibcoff                                                       11d12s20
      nff22=nff2+nsymb*mdoop*4                                          12d18s20
      ibcoff=nff22+nsymb*mdoop                                          12d18s20
      call enough('mrci.  3',bc,ibc)
      nneed=0                                                           11d17s20
c
c     irxinfo(3,i) is isb, nlzz, nroots, i=1,nrxinfos                   4d22s21
      do isb=1,nsymb                                                    11d17s20
       do i=1,nrxinfos                                                  8d24s21
        if(irxinfos(1,i).eq.isb)then                                    8d24s21
         nneed=nneed+1                                                   11d17s20
         nbasis2(nneed)=isb                                              11d17s20
         go to 241                                                      8d24s21
        end if                                                           11d17s20
       end do                                                           8d24s21
  241  continue                                                         8d24s21
      end do                                                            11d17s20
      call genfcn2(ibc(ibasis(1)),nbasis,ibc(iptr),ibc(idorb),          11d12s20
     $     ibc(isorb),nfcn,nec,mdon,mdoo,multh,nbasis2,nneed,ibc(nff2), 11d17s20
     $     nsymb,iff2,iptrto,nct2,ibc(icsf),bc,ibc)                     11d11s22
      ibcoff=iptrto                                                     12d7s20
      nff1=ibcoff                                                       11d12s20
      ibcoff=nff1+nsymb*mdoop*2                                         12d11s20
      call enough('mrci.  4',bc,ibc)
      call genfcnp(nec,mdon,mdoo,multh,ibc(nff2),ibc(iff2),ibc(nff1),   11d12s20
     $     iff1,nsymb,0,nct1,ibc(icsf),ntotf,bc,ibc)                    11d11s22
      nff0=ibcoff                                                       11d12s20
      ibcoff=nff0+nsymb*mdoop*2                                         12d11s20
      call enough('mrci.  5',bc,ibc)
      call genfcnp(nec,mdon,mdoo,multh,ibc(nff1),ibc(iff1),ibc(nff0),   11d12s20
     $     iff0,nsymb,isymmrci,nctf(isymmrci),ibc(icsf),ntotf,bc,ibc)   11d11s22
      nfcnf(isymmrci)=ntotf                                             1d21s21
      jptrf=ibcoff                                                      1d21s21
      ibcoff=jptrf+2*(mdoo+1)                                           1d21s21
      ibasisf(isymmrci)=ibcoff                                          1d21s21
      nww=ntotf*3                                                       1d21s21
      nww2=nww/2                                                        1d21s21
      if(nww2*2.ne.nww)nww2=nww2+1                                      1d21s21
      idorbf=ibasisf(isymmrci)+nww2                                     1d21s21
      nww=ntotf*norb                                                    1d21s21
      nww8=nww/8                                                        1d21s21
      if(nww8*8.ne.nww)nww8=nww8+1                                      1d21s21
      isorbf=idorbf+nww8                                                1d21s21
      ibcoff=isorbf+nww8                                                1d21s21
      call maptoold(mdon,mdoo,ibc(nff0),ibc(iff0),nfcnf(isymmrci),      1d21s21
     $    ibc(ibasisf(isymmrci)),ibc(idorbf),ibc(isorbf),ibc(jptrf),    1d21s21
     $     norb,nec,bc,ibc)                                             11d10s22
      iptrfbit=ibcoff                                                   1d21s21
      ibcoff=iptrfbit+nsymb*mdoop                                       1d21s21
      call int1tobit(ibc(jptrf),ibc(iptrfbit),ibc(idorbf),ibc(isorbf),  1d21s21
     $     idorbfbit,isorbfbit,mdoop,1,nec,bc,ibc)                      11d10s22
      ivint=ibcoff                                                      7d24s19
      ibcoff=ivint+nctf(isymmrci)*nroot                                 4d8s20
      ivintlast=ibcoff                                                  6d20s24
      ibcoff=ivintlast+nctf(isymmrci)*nroot                             6d20s24
      call enough('mrci.  6',bc,ibc)
      icsfpd=ibcoff                                                     6d11s19
      ibcoff=icsfpd+mdoop-mdon                                          5d7s20
      icsf2=ibcoff                                                      8d2s19
      ibcoff=icsf2+4*(mdoop-mdon)                                       5d7s20
      do isb=1,nsymb                                                    12d6s20
       ieigvspace(isb)=ibcoff                                           12d6s20
       if(isb.eq.isymmrci)then                                          12d6s20
        nkeep=max(nroot,irxinfo(isb))                                   12d7s20
       else                                                             12d6s20
        nkeep=irxinfo(isb)                                              12d6s20
       end if                                                           12d6s20
       irwt(isb)=ieigvspace(isb)+nkeep*(nct(isb)+1)                     12d6s20
       ibcoff=irwt(isb)+nkeep                                           12d6s20
      end do                                                            12d6s20
      call enough('mrci.  7',bc,ibc)
      if(nlzzci.ne.0)then                                               1d5s20
       call propcas(natom,ngaus,ibdat,isym,iapair,ibstor,isstor,idorel, 12d20s19
     $       ascale,multh,nbasisp,iorb,ixlzz,noc,islz,nlzzci,0,idum,0,  4d19s21
     $     idum,bc,ibc,epssymci)                                        9d19s23
       if(nlzzci.gt.2)then                                              1d5s20
        do isb=1,nsymb                                                  1d2s20
         nn=idoubo(isb)+irefo(isb)                                      1d2s20
c     4 is lz^2 5 is Lx^2 and 6 is ly^2.                                5d14s21
c     replace 4 with l^2 and 5 with lz^2.                               5d14s21
         do i=0,nn*nn-1                                                 1d2s20
          xlz=bc(ixlzz(isb,4)+i)                                        5d14s21
          xll=bc(ixlzz(isb,4)+i)+bc(ixlzz(isb,5)+i)+bc(ixlzz(isb,6)+i)  5d14s21
          bc(ixlzz(isb,4)+i)=xll                                        5d14s21
          bc(ixlzz(isb,5)+i)=xlz                                        5d14s21
         end do                                                         1d2s20
        end do                                                          1d2s20
       end if                                                           1d2s20
      end if                                                            12d24s19
      call second(time1)
      ibasisc=ibcoff                                                    3d24s21
      iptrcb=ibcoff                                                     3d24s21
      nfcnpxfi=ibcoff                                                   1d12s23
       if(iuse2den.ne.1)then                                             10d8s20
c
c     if using resolution of idenity, we mustn't skip fcns
c
        nfcnpxfi=ibcoff                                                 10d8s20
        ibcoff=nfcnpxfi+4*nsymb*mdoop                                   10d19s20
        idorbfi=ibcoff                                                  10d8s20
        isorbfi=idorbfi+nod8                                            10d8s20
        iptrfi=isorbfi+nos8                                             10d8s20
        ibasisfi(1)=iptrfi+(mdoo+1)*4*nsymb                             10d8s20
        nfcnx=1 000 000                                                   10d2s19
        ibcoff=ibasisfi(1)+3*nfcnx                                      10d8s20
        call enough('mrci.  8',bc,ibc)
        nfcnxx=nfcnx*6                                                  8d29s24
        call gencsf2b(nec,spin,norb,ibc(idorbfi),ibc(isorbfi),          11d19s20
     $       ibc(iptrfi),nsymb,mdon,mdoo,ibc(ibasisfi(1)),multh,ism,    10d8s20
     $     ibc(icsf),nfcnfi,nctfi,nod,nos,maxopens4,nhole0,ihole,nfill0, 9d5s24
     $       ifill,ido,ismult,0,nodu,nosu,nsymb,ibc(nfcnpxfi),          11d19s20
     $       ibc(icsf2),0,bc,ibc,nfcnxx)                                8d29s24
        need=0                                                           4d8s20
        jbasisfi=ibasisfi(1)                                            10d8s20
        nbss=0                                                           4d17s20
        do isb=1,nsymb                                                   4d8s20
         nhere=(nfcnfi(isb)*3)+mod(nfcnfi(isb),2)                       10d8s20
         nbasisfi(isb)=nbss                                             10d8s20
         nbss=nbss+nhere                                                 4d17s20
         need=need+nhere                                                 4d8s20
         ibasisfi(isb)=jbasisfi                                         10d8s20
         jbasisfi=jbasisfi+(nhere/2)                                    10d8s20
        end do                                                           4d8s20
        need8=need/2                                                     4d8s20
        ibcoff=ibasisfi(1)+need8                                        10d8s20
        ibasisc=ibcoff                                                  3d24s21
        iptrcb=ibasisc+need8                                            3d24s21
        ibcoff=iptrcb+mdoop                                             3d24s21
        call enough('mrci.  9',bc,ibc)
        call consolidate1(ibasisfi,ibc(iptrfi),mdoop,nfcnfi,            3d24s21
     $       ibc(idorbfi),ibc(isorbfi),nsymb,ibc(ibasisc),ibc(iptrcb),  3d24s21
     $       nfcnc,nec,bc,ibc)                                          11d14s22
        iptrfibit=ibcoff                                                10d8s20
        ibcoff=iptrfibit+nsymb*mdoop                                    10d8s20
        call int1tobit(ibc(iptrfi),ibc(iptrfibit),ibc(idorbfi),         10d8s20
     $       ibc(isorbfi),idorbfibit,isorbfibit,mdoop,nsymb,nec,bc,ibc) 11d10s22
       end if                                                           10d8s20
      if(lprint)then                                                    2d2s21
       call addcomma8(ibcoff+1-ioffsx,numbers,0)                        11d1s22
       write(6,*)('ibcoff b4 stepwisecsf '),numbers                     2d2s21
       call addcomma8(ihighwtr,numbers,0)                               2d2s21
       write(6,*)('high water mark: '),numbers                          2d2s21
      end if                                                            2d2s21
      if(nder.gt.0)then                                                 8d30s23
       do isb=1,nsymb                                                    8d30s23
        idvi(isb)=0                                                     8d30s23
       end do                                                            8d30s23
       do i=1,nrxinfos                                                  8d30s23
        isb=irxinfos(1,i)                                               8d30s23
        nrth=irxinfos(2,i)                                              8d30s23
        idvi(isb)=idvi(isb)+nrth*(1+nct(isb))                           8d30s23
       end do                                                           8d30s23
       do isb=1,nsymb                                                   8d30s23
        ineed=idvi(isb)*nder                                            8d30s23
        idvi(isb)=ibcoff                                                8d30s23
        ibcoff=idvi(isb)+ineed                                          8d30s23
       end do                                                           8d30s23
       call enough('mrci.dvi',bc,ibc)                                   8d30s23
      end if                                                            8d30s23
      ib4=ibcoff                                                        7d10s20
      call stepwisecsf(ixsoo,spin,xms,dum,idum,idum,1,idata,bc,ibc)     11d9s22
      call addcomma8(ihighwtr,numbers,0)                                2d2s21
      if(lprint)write(6,*)('ihighwtr after stepwisecsf: '),numbers      2d2s21
      call second(time1bx)
      telap=time1bx-time1-tovr
      if(lprint)write(6,*)('time for stepwisecsf = '),telap             1d25s21
      call gencsf3(nec,spin,mdon,mdoo,norb,ibc(icsfpd),ismult,          8d2s19
     $     ibc(icsf2),.false.,ixw1,ixw2,ibc(icsf),1,ib4,icsfxv,bc,ibc)  11d11s22
      ib4gencsf3=ib4                                                    12d6s20
      iafgencsf3=ibcoff                                                 12d6s20
      call addcomma8(ihighwtr,numbers,0)                                2d2s21
      if(lprint)write(6,*)('ihighwtr after gencsf3 '),numbers           2d2s21
      call second(timaa)
      ibcb4=ibcoff                                                      7d7s21
      ifull4x=1
      ipair=0                                                           2d26s24
      call paraeri(natom,ngaus,ibdat,idum,ihmat,iorb,noc,               4d24s20
     $             ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,  11d29s12
     $     iter1,idwsdeb,idorel,ascale,nbasisp,ifull4x,0,idum,0,idum,   5d20s24
     $     idum,bc,ibc)                                                 5d20s24
      icol=ibcoff
      ibcoff=icol+mynprocg
      imsg=ibcoff                                                       6d3s10
      ibcoff=imsg+mynnode                                               6d3s10
      iter1=0
      idwsdeb=0
      call parajkfromh(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,         11d29s12
     $     iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      5d30s18
     $     ioooo,ionex,irefo,iter1,idwsdeb,ncomp,nvirt,i3x,1,ih0,       5d31s18
     $     idoubo,shift,nbasisp,idorel,ifull4x,i4x,i4xb,ibcb4,bc,ibc,   8d21s23
     $     isend1,isend2,isend3,isend4,nder,iff)                        7d11s23
      call second(timaaa)
      telap=timaaa-timaa-tovr
      shift=shift+potdws                                                8d21s18
      shift0=0d0                                                        6d1s21
      jh0=ih0                                                           6d5s18
      do isb=1,nsymb                                                    6d5s18
       ih0a(isb)=ibcoff                                                 6d5s18
       ibcoff=ih0a(isb)+irefo(isb)*irefo(isb)                           6d5s18
       call enough('mrci. 10',bc,ibc)
       ndim=nbasdws(isb)-idoubo(isb)                                    6d5s18
       sym=0d0
       do i=0,irefo(isb)-1                                              6d5s18
        do j=0,irefo(isb)-1                                             6d5s18
         iad1=jh0+j+ndim*i                                              6d5s18
         iad2=ih0a(isb)+j+irefo(isb)*i                                  6d5s18
         bc(iad2)=bc(iad1)                                              6d5s18
         if(iad2.eq.igoal)write(6,*)('getting from ih0 '),bc(iad2),
     $        iad1
         if(j.ne.i)then
          iad1t=jh0+i+ndim*j                                              6d5s18
          diff=bc(iad1t)-bc(iad1)                                       10d28s20
          sym=sym+diff**2                                               10d28s20
         end if
        end do                                                          6d5s18
       end do                                                           6d5s18
       ih0av(isb)=jh0
       nh0av(isb)=ndim                                                  6d28s18
       jh0=jh0+ndim*ndim                                                6d5s18
      end do                                                            6d5s18
c
c     replicate oooo on all procs
c
      do i=1,mynprocg                                                   6d1s18
       isend1(i)=0                                                      6d1s18
      end do                                                            6d1s18
      ntot=0                                                            6d1s18
      do is=1,nsdlk                                                     6d1s18
       if(isblk(1,is).eq.isblk(2,is))then                               6d1s18
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              6d1s18
       else                                                             6d1s18
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      6d1s18
       end if                                                           6d1s18
       if(isblk(3,is).eq.isblk(4,is))then                               6d1s18
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              6d1s18
       else                                                             6d1s18
        ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       6d1s18
       end if                                                           6d1s18
       if(nrow.gt.0)then                                                6d1s18
        ntot=ntot+ncol*nrow                                              6d1s18
        do ip=0,mynprocg-1                                               6d1s18
         ipp=ip+1                                                        6d1s18
         call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,ip,  6d1s18
     $       il,ih,i1s,i1e,i2s,i2e)                                     6d1s18
         i10=i1s                                                        6d1s18
         i1n=irefo(isblk(3,is))                                         6d1s18
         nhere=0                                                        6d1s18
         do i2=i2s,i2e                                                  6d1s18
          if(i2.eq.i2e)i1n=i1e                                          6d1s18
          do i1=i10,i1n                                                 6d1s18
           if(i1.ge.i2.or.isblk(3,is).ne.isblk(4,is))nhere=nhere+1      6d1s18
          end do                                                        6d1s18
          i10=1                                                         6d1s18
         end do                                                         6d1s18
         isend1(ipp)=isend1(ipp)+nhere*nrow                              6d1s18
        end do                                                           6d1s18
       end if                                                           6d1s18
      end do                                                            6d1s18
      isum=0
      do i=1,mynprocg
       isum=isum+isend1(i)
      end do
      idest=ibcoff                                                      6d1s18
      ibcoff=idest+ntot                                                 6d1s18
      ibufr=ibcoff                                                      6d1s18
      ibcoff=ibufr+ntot                                                 6d1s18
      call enough('mrci. 11',bc,ibc)
      jdest=idest                                                       6d1s18
      do is=1,nsdlk                                                     6d1s18
       if(isblk(1,is).eq.isblk(2,is))then                               6d1s18
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              6d1s18
       else                                                             6d1s18
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      6d1s18
       end if                                                           6d1s18
       if(nrow.gt.0)then                                                6d1s18
        call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,      6d1s18
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            6d1s18
        i10=i1s                                                         6d1s18
        i1n=irefo(isblk(3,is))                                          6d1s18
        i4o=ioooo(is)                                                   6d1s18
        do i2=i2s,i2e                                                   6d1s18
         if(i2.eq.i2e)i1n=i1e                                           6d1s18
         do i1=i10,i1n                                                  6d1s18
          if(i1.ge.i2.or.isblk(3,is).ne.isblk(4,is))then                6d1s18
           do irow=0,nrow-1                                             6d1s18
            bc(jdest+irow)=bc(i4o+irow)                                 6d1s18
           end do                                                       6d1s18
           jdest=jdest+nrow                                             6d1s18
          end if                                                        6d1s18
          i4o=i4o+nrow                                                  6d1s18
         end do                                                         6d1s18
         i10=1                                                          6d1s18
        end do                                                          6d1s18
       end if                                                           6d1s18
      end do                                                            6d1s18
      do i=1,mynprocg                                                   6d1s18
       isend2(i)=0                                                      6d1s18
       isend3(i)=isend1(mynowprog+1)                                    6d1s18
      end do                                                            6d1s18
      isend4(1)=0                                                       6d1s18
      do i=1,mynprocg-1                                                 6d1s18
       im=i-1                                                           6d1s18
       isend4(i+1)=isend4(i)+isend1(i)                                  6d1s18
      end do                                                            6d1s18
      call dws_all2allvb(bc(idest),isend3,isend2,bc(ibufr),             6d1s18
     $     isend1,isend4)                                               6d1s18
      do i=0,mynprocg-1                                                 6d1s18
       iarg=isend1(i+1)                                                  6d1s18
       jbufr=ibufr+isend4(i+1)
       isend4(i+1)=jbufr                                                6d1s18
      end do                                                            6d1s18
      do is=1,nsdlk                                                     6d1s18
       ioooo(is)=idest                                                  6d1s18
       if(isblk(1,is).eq.isblk(2,is))then                               6d1s18
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              6d1s18
       else                                                             6d1s18
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      6d1s18
       end if                                                           6d1s18
       if(isblk(3,is).eq.isblk(4,is))then                               6d1s18
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              6d1s18
       else                                                             6d1s18
        ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       6d1s18
       end if                                                           6d1s18
       if(nrow*ncol.gt.0)then                                           6d1s18
        idest=idest+nrow*ncol                                            6d1s18
        do i=0,mynprocg-1                                                6d1s18
         ip=i+1                                                          6d1s18
         call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,i,   6d1s18
     $       il,ih,i1s,i1e,i2s,i2e)                                     6d1s18
         i10=i1s                                                        6d1s18
         i1n=irefo(isblk(3,is))                                         6d1s18
         do i2=i2s,i2e                                                  6d1s18
          if(i2.eq.i2e)i1n=i1e                                          6d1s18
          do i1=i10,i1n                                                 6d1s18
           if(i1.ge.i2.or.isblk(3,is).ne.isblk(4,is))then               6d1s18
            if(isblk(3,is).eq.isblk(4,is))then                          6d1s18
             ito=ioooo(is)+nrow*(((i1*(i1-1))/2)+i2-1)                   6d1s18
            else                                                        6d1s18
             ito=ioooo(is)+nrow*(i1-1+irefo(isblk(3,is))*(i2-1))        6d1s18
            end if                                                      6d1s18
            do irow=1,nrow                                              6d1s18
             bc(ito)=bc(isend4(ip))                                     6d1s18
             ito=ito+1                                                     6d1s18
             isend4(ip)=isend4(ip)+1                                       6d1s18
            end do                                                      6d1s18
           end if                                                       6d1s18
          end do                                                        6d1s18
          i10=1                                                         6d1s18
         end do                                                         6d1s18
        end do                                                           6d1s18
       end if                                                           6d1s18
      end do                                                            6d1s18
      ibcoff=ibufr                                                      6d1s18
c
c     fold -0.5*(ij|jl) into h0a(il) for resolution of identity         4d8s20
c
      do isb=1,nsymb                                                    4d8s20
       ih0ae(isb)=ibcoff                                                4d8s20
       ibcoff=ih0ae(isb)+irefo(isb)*irefo(isb)                          4d24s20
       call enough('mrci. 12',bc,ibc)
       do i=0,irefo(isb)*irefo(isb)-1                                   4d8s20
        bc(ih0ae(isb)+i)=bc(ih0a(isb)+i)                                4d8s20
       end do                                                           4d8s20
       do i=0,irefo(isb)-1                                              4d8s20
        do l=0,irefo(isb)-1                                             4d8s20
         iad=ih0ae(isb)+l+irefo(isb)*i                                  4d8s20
         do jsb=1,nsymb                                                   4d8s20
          do j=0,irefo(jsb)-1                                             4d8s20
           term=getint(ioooo,isb,jsb,jsb,isb,i,j,j,l,bc,ibc)            11d15s22
           orig=bc(iad)
           bc(iad)=bc(iad)-0.5d0*term                                   4d29s20
          end do                                                        4d8s20
         end do                                                         4d8s20
        end do                                                          4d8s20
       end do                                                           4d8s20
      end do
      if(iunc(1).eq.0.and.iunc(2).eq.0.and.ethred.eq.0d0)then           6d20s23
       igoon=0                                                           1d10s19
      else                                                              6d20s23
       igoon=1                                                          6d20s23
      end if                                                            6d20s23
      call second(time1)
      pthrs0=pthrsci                                                    6d6s18
      if(lprint)write(6,*)('spin = '),spin
      s2=spin*(spin+1d0)                                                6d7s18
      if(lprint)write(6,*)('desired S2 eigenvalue '),s2                           6d7s18
c
c     for n electron determinants.
c
      s20=s2                                                            8d24s18
c
c     weights for internal contraction
c
      do isb=1,nsymb                                                    7d26s19
       ivintref(isb)=0                                                  7d26s19
      end do                                                            7d26s19
       nrootu=nroot                                                     10d16s20
      call enough('mrci. 13',bc,ibc)
      do i=0,nrootu-1                                                   5d27s19
       bc(irwt(isymmrci)+i)=1d0                                         5d23s19
      end do                                                            3d1s19
      call second(time1)
      jptr=iptr+(mdoo+1)*2*(isymmrci-1)                                 7d15s19
      nlzzx=0                                                           5d1s20
      nlzzx=nlzzci                                                      2d25s21
      if(iuse2den.eq.1)then                                             10d18s20
      call intcsf(ibc(idorb),ibc(isorb),mdon,mdoo,ibc(ibasis(isymmrci)),7d15s19
     $     ibc(jptr),ibc(icsf),nfcn(isymmrci),ih0a,ioooo,shift,pthrs0,  7d15s19
     $     nec,ibc(icsfpd),multh,isymmrci,0,dum,dum,nrootu,0d0,0.1d0,
     $     lprint,ieigint,ivintref(isymmrci),maxpsci,ncsfint,nvcul,     12d24s19
     $     nlzzx,islz,ixlzz,igoon,bintvec,intvec,nfcn,ibc(iptr),ibasis, 5d1s20
     $     ih0ae,ibc(iptrbit),ixw1,ixw2,nbasis,ibc(nfcnpx),nsymb,       10d15s20
     $     ibc(iptrbit),ibc(idorb),ibc(isorb),idum,idum,idum,vidotvi,   3d31s21
     $     vihvi,nameci,0,bc,ibc,isaved,ostng,idroot)                   7d13s23
      else                                                              10d18s20
      call intcsf(ibc(idorb),ibc(isorb),mdon,mdoo,ibc(ibasis(isymmrci)),7d15s19
     $     ibc(jptr),ibc(icsf),nfcn(isymmrci),ih0a,ioooo,shift,pthrs0,  7d15s19
     $     nec,ibc(icsfpd),multh,isymmrci,0,dum,dum,nrootu,0d0,0.1d0,
     $     lprint,ieigint,ivintref(isymmrci),maxpsci,ncsfint,nvcul,     12d24s19
     $     nlzzx,islz,ixlzz,igoon,bintvec,intvec,nfcnfi,ibc(iptrfi),    10d15s20
     $     ibasisfi,ih0ae,ibc(iptrbit),ixw1,ixw2,nbasisfi,ibc(nfcnpx),  10d15s20
     $     nsymb,ibc(iptrfibit),ibc(idorbfi),ibc(isorbfi),nfcnc,        3d24s21
     $      ibc(ibasisc),ibc(iptrcb),vidotvi,vihvi,nameci,0,bc,ibc,     2d3s23
     $      isaved,ostng,idroot)                                        7d13s23
      end if                                                            10d18s20
      do i=0,nroot-1                                                    12d7s20
       bc(ieigvspace(isymmrci)+i)=bc(ieigint+i)                         12d7s20
      end do                                                            12d7s20
      ivnew=ieigvspace(isymmrci)+max(nroot,irxinfo(isymmrci))           4d28s21
      do i=0,nroot*nct(isymmrci)-1                                      12d7s20
       bc(ivnew+i)=bc(ivintref(isymmrci)+i)                             12d7s20
      end do                                                            12d7s20
      ibcreset=ieigint                                                  3d12s21
      ieigint=ieigvspace(isymmrci)                                      12d7s20
      ivintref(isymmrci)=ivnew                                          12d7s20
      if(iunc(1).eq.0.and.iunc(2).eq.0.and.ethred.eq.0d0)then           11d24s19
       if(mynowprog.eq.0)then                                           11d24s19
        open(unit=42,file='summary')                                    1d10s19
        write(42,*)(bc(iextrad+i),i=0,nextrad-1),ncsfint,spinm,termsym, 1d10s19
     $       (bc(ieigint+i)+shift,i=0,nroot-1)                          5d27s21
        close(unit=42)                                                  1d10s19
        if(nwavef.ne.0)then                                             3d23s21
         ntot=0                                                          11d24s19
         nwavrec=nwavrec+1
         write(2)ncsfint,ntot                                            11d24s19
         nwavrec=nwavrec+1
         write(2)(bc(ivintref(isymmrci)+i),i=0,ncsfint*nroot-1)          11d24s19
         nwavrec=nwavrec+1
         write(2)(bc(ieigint+i)+shift,i=0,nroot-1)                      11d17s22
         ipack2(1)=nlzzci                                               10d7s24
         ipack2(2)=nder                                                 10d7s24
         write(2)ipack4(1),lambdaci                                     10d7s24
         izero=0                                                        7d19s21
         write(2)izero                                                  7d19s21
         write(2)izero                                                  7d19s21
         if(nlzzci.ne.0)then                                            5d12s21
          write(6,*)('now go on and dumpbas for equivalent symmetries')
          if(lambdaci.ne.0)then                                         5d12s21
           if(nlzzci.eq.2)then                                          5d12s21
            if(isymmrci.eq.1)then                                        5d12s21
             isymo=4                                                     5d12s21
            else if(isymmrci.eq.2)then                                   5d12s21
             isymo=3                                                     5d12s21
            else if(isymmrci.eq.3)then                                   5d12s21
             isymo=2                                                     5d12s21
            else if(isymmrci.eq.4)then                                   5d12s21
             isymo=1                                                     5d12s21
            else if(isymmrci.eq.5)then                                   5d12s21
             isymo=8                                                     5d12s21
            else if(isymmrci.eq.6)then                                   5d12s21
             isymo=7                                                     5d12s21
            else if(isymmrci.eq.7)then                                   5d12s21
             isymo=6                                                     5d12s21
            else                                                         5d12s21
             isymo=5                                                     5d12s21
            end if                                                       5d12s21
            write(2)isymo                                                5d12s21
            jptr=iptr+2*mdoop*(isymo-1)                                  5d12s21
            call dumpbas(nfcn(isymo),ibc(jptr),ibc(ibasis(isymo)),       5d12s21
     $      ibc(idorb),ibc(isorb),nodu,nosu,mdoop,nec,1,nct(isymo),     5d12s21
     $      ibc(icsf),mdon,nwavrec,0,bc,ibc)                            11d14s22
            nwavrec=nwavrec+1
           else                                                         5d12s21
            nl=2*lambdaci                                               5d12s21
            isymo=ibcoff                                                5d12s21
            ibcoff=isymo+nl                                             5d12s21
            call enough('mrci. 14',bc,ibc)
            write(6,*)('isymmrci = '),isymmrci
            if(isymmrci.eq.1.or.isymmrci.eq.4)then                      5d14s21
             iodd1=6                                                    5d14s21
             iodd2=7                                                    5d14s21
             ieven1=1                                                   5d14s21
             ieven2=4                                                   5d14s21
            else                                                        5d14s21
             iodd1=2                                                    5d14s21
             iodd2=3                                                    5d14s21
             ieven1=5                                                   5d14s21
             ieven2=8                                                   5d14s21
            end if                                                      5d14s21
            jsymo=isymo                                                 5d14s21
            do ll=1,lambdaci                                            5d14s21
             if(mod(ll,2).eq.0)then                                     5d14s21
              ibc(jsymo)=ieven1                                         5d14s21
              jsymo=jsymo+1                                             5d14s21
              ibc(jsymo)=ieven2                                         5d14s21
              jsymo=jsymo+1                                             5d14s21
             else                                                       5d14s21
              ibc(jsymo)=iodd1                                          5d14s21
              jsymo=jsymo+1                                             5d14s21
              ibc(jsymo)=iodd2                                          5d14s21
              jsymo=jsymo+1                                             5d14s21
             end if                                                     5d14s21
            end do                                                      5d14s21
            do is=0,nl-1                                                5d12s21
             mysym=ibc(isymo+is)                                        5d12s21
             write(2)mysym                                              5d12s21
             jptr=iptr+2*mdoop*(mysym-1)                                  5d12s21
             write(6,*)('dumpbas3 '),nwavrec
             call dumpbas(nfcn(mysym),ibc(jptr),ibc(ibasis(mysym)),       5d12s21
     $      ibc(idorb),ibc(isorb),nodu,nosu,mdoop,nec,1,nct(mysym),     5d12s21
     $      ibc(icsf),mdon,nwavrec,0,bc,ibc)                            11d14s22
             nwavrec=nwavrec+1
            end do                                                      5d12s21
           end if                                                       5d12s21
          end if                                                        5d12s21
         end if                                                         5d12s21
         close(unit=2)                                                   11d24s19
        end if                                                          3d23s21
       end if                                                           11d24s19
       jptrbit=iptrbit+mdoop*(isymmrci-1)                               3d12s21
       if(nder.eq.0)then                                                7d29s22
        call denmx(mdon,mdoo,ibc(ibasis(isymmrci)),ibc(jptrbit),         3d12s21
     $      ibc(icsf),nfcn(isymmrci),bc(ivintref(isymmrci)),ncsfint,    3d12s21
     $      isymmrci,irel,ism,irefo,nvirt,nsymb,multh,ixw1,ixw2,nh0av,  3d12s21
     $      nroot,idum,idum,idum,idum,idum,0,0,1,norb,idum,idum,        3d12s21
     $      dum,nec,ibc(icsf2),ih0copy,nbasdws,nbasisp,ncomp,lprint,sr2,6d22s21
     $      srh,norbci,iorb,idoubo,bc,ibc)                              11d10s22
       else                                                             7d29s22
        call denmx12(mdon,mdoo,ibc(ibasis(isymmrci)),ibc(jptrbit),         3d12s21
     $      ibc(icsf),nfcn(isymmrci),bc(ivintref(isymmrci)),ncsfint,    3d12s21
     $      isymmrci,irel,ism,irefo,multh,ixw1,ixw2,nh0av,              7d29s22
     $      nroot,idum,idum,idum,idum,idum,0,0,1,norb,idum,idum,        3d12s21
     $      dum,nec,ibc(icsf2),ih0av,nbasisp,ncomp,lprint,sr2,srh,      7d29s22
     $      ioooo,shift,bc,ibc,ionexbt,jmats,kmats,idoubo,iff,idarot,   7d11s23
     $       nder,natom,ngaus,ibdat,ihmat,iorb,noc,ipair,isym,iapair,   7d17s23
     $       ibstor,isstor,iptoh,idorel,ascale,i3x,i4x,i4xb,            7d17s23
     $       isend1,isend2,isend3,isend4,ih0copy,idum,idum,dum,idum,dum,9d28s23
     $       idum,idum,idum,idum,idum,idum,idum,idum,tovr,idum,idum,    5d2s24
     $       potdws,idervcv,ipropsym)                                   10d30s24
       end if                                                           7d29s22
       ibcoff=ibcoffo                                                   11d24s19
       return                                                           11d24s19
      end if                                                            11d24s19
      ibcoff=ibcreset                                                   3d12s21
      elow=bc(ieigint)                                                  5d23s19
      wsum=0d0                                                          5d27s19
      iwtmp=0
      dynwcii=1d0/dynwci(1)                                             9d20s21
      if(lprint)write(6,*)('dynwci = '),dynwci                          5d25s21
      dynwci(1)=dynwcii                                                 8d29s23
      do isb=1,nsymb                                                    4d25s21
       jeigvspace(isb)=ieigvspace(isb)                                  4d25s21
       nrootu=irxinfo(isb)                                              4d30s21
       if(isb.eq.isymmrci)nrootu=max(nrootu,nroot)                      4d30s21
       jvintref(isb)=jeigvspace(isb)+nrootu                             4d30s21
      end do                                                            4d25s21
      if(nreadc.ne.0)then                                               8d10s22
       if(mynowprog.eq.0)then                                           8d10s22
        write(6,3338)                                                   8d10s22
 3338   format('rather then using internal ci vectors to generate',
     $        ' contraction coefficients, read them in from ')          12d13s22
        write(6,*)cfile(1:ncfile)                                       12d13s22
        open(unit=1,file=cfile(1:ncfile),form='unformatted')            12d13s22
       end if                                                           8d10s22
      end if                                                            8d10s22
      do i=1,nrxinfos                                                   4d25s21
       if(nlzzci.eq.0)then                                              8d24s21
        if(lprint)write(6,237)                                          8d24s21
  237   format(4x,'#',9x,'sym',x,'roots')                               8d24s21
        if(lprint)write(6,240)i,(irxinfos(j,i),j=1,2)                   8d24s21
  240   format(i5,5x,5i6)                                               8d24s21
       else if(nlzzci.eq.2)then                                         8d24s21
        if(lprint)write(6,238)                                          8d24s21
  238   format(4x,'#',9x,'sym',x,'roots',2x,'nlzz',x,'Lambda')          8d24s21
        if(lprint)write(6,240)i,(irxinfos(j,i),j=1,4)                   8d24s21
       else if(nlzzci.eq.6)then                                         8d24s21
        if(lprint)write(6,239)                                          8d24s21
  239   format(4x,'#',9x,'sym',x,'roots',2x,'nlzz',4x,'L',4x,'Lz')      8d24s21
        if(lprint)write(6,240)i,(irxinfos(j,i),j=1,5)                   8d24s21
       end if                                                           8d23s21
       isb=irxinfos(1,i)                                                4d25s21
       nrth=irxinfos(2,i)                                               4d25s21
       nlzzh=irxinfos(3,i)                                              4d25s21
       lamh=irxinfos(4,i)                                               4d25s21
       lz2=irxinfos(5,i)                                                8d23s21
       jptr=iptr+(mdoo+1)*2*(isb-1)                                     4d25s21
       lamsav=lambdaci                                                  4d25s21
       lambdaci=lamh                                                    4d25s21
       do j=1,6                                                         4d25s21
        nameciu(j)=0                                                    4d25s21
       end do                                                           4d25s21
       if(nlzzh.eq.0)then                                               4d25s21
        nameciu(2)=ichar('S')                                           4d25s21
        write(numbers(1:1),65)isb                                       4d25s21
   65   format(i1)                                                      4d25s21
        nameciu(3)=ichar(numbers(1:1))                                  4d25s21
       else if(nlzzh.eq.2)then                                          4d25s21
        if(lamh.eq.0)then                                               4d25s21
         numbers(1:5)='Sig  '                                             4d25s21
         if(isb.eq.1.or.isb.eq.5)then                                   4d25s21
          numbers(4:4)='+'                                              4d25s21
         else                                                           4d25s21
          numbers(4:4)='-'                                              4d25s21
         end if                                                         4d25s21
         if(nsymb.eq.8)then                                             4d25s21
          if(isb.le.4)then                                               4d25s21
           numbers(5:5)='g'                                              4d25s21
          else                                                           4d25s21
           numbers(5:5)='u'                                              4d25s21
          end if                                                         4d25s21
         end if                                                         4d27s21
        else if(lamh.eq.1)then                                          4d25s21
         numbers(1:5)='Pi   '                                           4d25s21
         if(nsymb.eq.8)then                                             4d25s21
          if(isb.le.4)then                                              4d25s21
           numbers(3:3)='u'                                             4d25s21
          else                                                          4d25s21
           numbers(3:3)='g'                                             4d25s21
          end if                                                        4d25s21
         end if                                                         4d25s21
        else if(lamh.eq.2)then                                          4d25s21
         numbers(1:5)='Del  '                                           4d25s21
         if(nsymb.eq.8)then                                             4d25s21
          if(isb.le.4)then                                              4d25s21
           numbers(4:4)='g'                                             4d25s21
          else                                                          4d25s21
           numbers(4:4)='u'                                             4d25s21
          end if                                                        4d25s21
         end if                                                         4d25s21
        else if(lamh.eq.3)then                                          4d25s21
         numbers(1:5)='Phi  '                                           4d25s21
         if(nsymb.eq.8)then                                             4d25s21
          if(isb.le.4)then                                              4d25s21
           numbers(4:4)='u'                                             4d25s21
          else                                                          4d25s21
           numbers(4:4)='g'                                             4d25s21
          end if                                                        4d25s21
         end if                                                         4d25s21
        else if(lamh.eq.4)then                                          4d25s21
         numbers(1:5)='Gam  '                                           4d25s21
         if(nsymb.eq.8)then                                             4d25s21
          if(isb.le.4)then                                              4d25s21
           numbers(4:4)='g'                                             4d25s21
          else                                                          4d25s21
           numbers(4:4)='u'                                             4d25s21
          end if                                                        4d25s21
         end if                                                         4d25s21
        end if                                                          4d25s21
       else                                                             4d25s21
        if(lamh.eq.0)then                                               4d25s21
         numbers(1:5)='S    '                                           4d25s21
        else if(lamh.eq.1)then                                          4d25s21
         numbers(1:5)='P    '                                           4d25s21
        else if(lamh.eq.2)then                                          4d25s21
         numbers(1:5)='D    '                                           4d25s21
        else if(lamh.eq.3)then                                          4d25s21
         numbers(1:5)='F    '                                           4d25s21
        else if(lamh.eq.4)then                                          4d25s21
         numbers(1:5)='G    '                                           4d25s21
        end if                                                          4d25s21
        if(isb.eq.8.or.isb.eq.2.or.isb.eq.3.or.isb.eq.5)numbers(2:2)='o'4d25s21
       end if
       do j=1,5                                                         4d25s21
        nameciu(j)=ichar(numbers(j:j))                                  4d27s21
       end do                                                           4d25s21
       if(nreadc.ne.0)then                                              8d10s22
        call cfileget(ibc(idorb),ibc(isorb),mdon,mdoo,ibc(ibasis(isb)), 8d10s22
     $      ibc(iptrbit),ibc(icsf),nfcn(isb),multh,ieiginth,ivintinth,  8d10s22
     $      isb,nrth,nlzzh,lamh,lz2,nameciu,nbasisp,nbasdws,nvirt,nsymb,8d10s22
     $      nct(isb),norb,shift,lprint,idorel,bc,ibc,cfile(1:ncfile))   12d19s22
       else
        call intcsf(ibc(idorb),ibc(isorb),mdon,mdoo,ibc(ibasis(isb)),    4d25s21
     $     ibc(jptr),ibc(icsf),nfcn(isb),ih0a,ioooo,shift,pthrs0,       7d15s19
     $     nec,ibc(icsfpd),multh,isb,0,dum,dum,nrth,0d0,0.1d0,          4d25s21
     $     lprint,ieiginth,ivintinth,maxpsci,ncsfinth,nvcul,            4d25s21
     $      nlzzh,islz,ixlzz,1,bintvec,intvec,nfcnfi,ibc(iptrfi),
     $      ibasisfi,ih0ae,                                             4d25s21
     $     ibc(iptrbit),ixw1,ixw2,nbasisfi,ibc(nfcnpx),nsymb,           10d15s20
     $     ibc(iptrfibit),ibc(idorbfi),ibc(isorbfi),nfcnc,              3d24s21
     $      ibc(ibasisc),ibc(iptrcb),vidotvi,vihvi,nameciu,lz2,bc,ibc,  2d3s23
     $       isaved,ostng,idroot)                                       7d13s23
       end if                                                           8d10s22
       lambdaci=lamsav                                                  4d25s21
       do j=0,nrth-1                                                    4d25s21
        bc(jeigvspace(isb)+j)=bc(ieiginth+j)                            4d25s21
       end do                                                           4d25s21
       jeigvspace(isb)=jeigvspace(isb)+nrth                             4d25s21
       do j=0,nrth*nct(isb)-1                                           4d25s21
        bc(jvintref(isb)+j)=bc(ivintinth+j)                             4d25s21
       end do                                                           4d25s21
       l4=irxinfos(1,i).eq.isymmrci                                     5d24s23
       if(nlzzci.eq.2)then                                              5d24s23
        l4=l4.and.irxinfos(4,i).eq.lambdaci                             5d24s23
       else if(nlzzci.eq.6)then                                         5d24s23
        l4=l4.and.irxinfos(4,i).eq.lambdaci                             5d24s23
        l4=l4.and.irxinfos(5,i).eq.0                                    5d24s23
       end if                                                           5d24s23
       if(l4)lvintref=jvintref(isb)                                     5d24s23
       jvintref(isb)=jvintref(isb)+nrth*nct(isb)                        4d25s21
       ibcoff=ieiginth                                                  4d25s21
      end do                                                            4d25s21
      if(nder.eq.0)then                                                 8d29s23
       do isb=1,nsymb                                                    5d23s19
        ieiginth=ieigvspace(isb)                                        12d7s20
        iewgt=ibcoff                                                    8d29s23
        ibcoff=iewgt+irxinfo(isb)                                       8d29s23
        call enough('mrci.ewgt',bc,ibc)                                 8d29s23
        do i1=iewgt,ibcoff-1                                            8d29s23
         bc(i1)=1d0                                                     8d29s23
        end do                                                          8d29s23
        call dynwtr(irxinfo(isb),bc(ieiginth),bc(irwt(isb)),elow,dum,   8d29s23
     $       dynwci,bc(iewgt),wsum,dum,shift,dum,0)                     8d29s23
        ibcoff=iewgt                                                    8d29s23
       end do                                                           5d23s19
       wsumi=1d0/wsum                                                    5d27s19
      else                                                              8d29s23
       do isb=1,nsymb                                                    4d25s21
        jeigvspace(isb)=ieigvspace(isb)                                  4d25s21
        nrootu=irxinfo(isb)                                              4d30s21
        if(isb.eq.isymmrci)nrootu=max(nrootu,nroot)                      4d30s21
        jvintref(isb)=jeigvspace(isb)+nrootu                             4d30s21
       end do                                                            4d25s21
       call genivder(nder,idarot,nrxinfos,irxinfos,nlzzci,mdoo,         8d29s23
     $       iptr,ioooo,nsdlk,isblk,ionex,nsdlk1,isblk1,idoubo,irefo,   8d17s23
     $       noc,nvirt,nbasdws,nsymb,natom,ngaus,ibdat,isym,iapair,     8d17s23
     $       ibstor,isstor,idorel,ascale,multh,nbasisp,iorb,ih0copy,    8d17s23
     $       bc,ibc,isend1,isend2,isend3,isend4,iptoh,                  8d23s23
     $       mdon,ibc(idorb),ibc(isorb),ibasis,ibc(icsf),nfcn,pthrs0,   8d23s23
     $       nec,jvintref,nroot,isymmrci,nct,iptrbit,ixw1,ixw2,ih0ae,   8d25s23
     $       nfcnc,ibc(ibasisc),ibc(iptrcb),norb,irel,ism,ih0a,idvi,    7d22s24
     $       toldv,maxdiisci,maxrestci)                                 10d16s24
       do isb=1,nsymb                                                   8d29s23
        idviu(isb)=idvi(isb)                                            8d29s23
        do i=0,irxinfo(isb)-1                                           8d29s23
         if(abs(bc(ieigvspace(isb)+i)-elow).lt.1d-10)then               8d29s23
          ielows=isb                                                    8d29s23
          ielowr=i                                                      8d29s23
         end if                                                         8d29s23
        end do                                                          8d29s23
       end do                                                           8d29s23
       do ider=1,nder                                                   8d29s23
        wsum=0d0                                                        8d29s23
        dwsum=0d0                                                       8d29s23
        delow=bc(idviu(ielows)+ielowr)                                  8d29s23
        do isb=1,nsymb                                                  8d29s23
         if(irxinfo(isb).gt.0)then                                      8d29s23
          ieigd=ibcoff                                                  8d29s23
          iewgt=ieigd+2*irxinfo(isb)                                    8d29s23
          iwgtd=iewgt+irxinfo(isb)                                      8d29s23
          ibcoff=iwgtd+irxinfo(isb)*2                                   8d29s23
          call enough('mrci.eigd',bc,ibc)                               8d29s23
          jeigd=ieigd+irxinfo(isb)                                      8d29s23
          do i=0,irxinfo(isb)-1                                         8d29s23
           bc(iewgt+i)=1d0                                              8d29s23
           bc(ieigd+i)=bc(ieigvspace(isb)+i)                              8d29s23
           bc(jeigd+i)=bc(idviu(isb)+i)                                   8d29s23
          end do                                                        8d29s23
          call dynwtr(irxinfo(isb),bc(ieigd),bc(iwgtd),elow,delow,      8d29s23
     $           dynwci,bc(iewgt),wsum,dwsum,shift,0d0,1)               8d29s23
          jwgtd=iwgtd+irxinfo(isb)                                      8d29s23
          do i=0,irxinfo(isb)-1                                         8d29s23
           bc(irwt(isb)+i)=bc(iwgtd+i)                                  8d29s23
           bc(idviu(isb)+i)=bc(jwgtd+i)                                 8d29s23
          end do                                                        8d29s23
          ibcoff=ieigd                                                  8d29s23
         end if                                                         8d29s23
        end do                                                          8d29s23
        wsumi=1d0/wsum                                                  8d29s23
        dwsum=dwsum*wsumi                                               8d29s23
c
c     wi=wi/wsum
c     dwi=dwi/wsum-wi*dwsum/wsum**2
c
        do isb=1,nsymb                                                  8d29s23
         do i=0,irxinfo(isb)-1                                          8d29s23
          orig=bc(idviu(isb)+i)
          bc(idviu(isb)+i)=(bc(idviu(isb)+i)-bc(irwt(isb)+i)*dwsum)     8d29s23
     $         *wsumi                                                   8d29s23
         end do                                                         8d29s23
         idviu(isb)=idviu(isb)+irxinfo(isb)*(1+nct(isb))                8d29s23
        end do                                                          8d29s23
       end do                                                            5d23s19
      end if                                                            5d23s19
      if(lprint)write(6,*)('renormalized weights ')                     6d4s19
      do isb=1,nsymb                                                    5d27s19
       if(isb.eq.-isymmrci)then                                         12d29s21
        nrootu=max(nroot,irxinfo(isb))                                  5d27s19
       else                                                             5d27s19
        nrootu=irxinfo(isb)
       end if                                                           5d27s19
       do j=0,nrootu-1                                                  5d27s19
        bc(irwt(isb)+j)=bc(irwt(isb)+j)*wsumi                           5d27s19
       end do                                                           5d27s19
       if(lprint)write(6,*)isb,(bc(irwt(isb)+j),j=0,nrootu-1)
      end do                                                            5d27s19
      ieigintref=ieigint                                                2d26s19
      do i=0,ncsfint-1                                                  4d2s19
       do j=0,nroot-1                                                   4d2s19
        ij=ivintref(isymmrci)+i+ncsfint*j                               7d15s19
        ji2=ivint+j+nroot*i                                             5d27s19
        bc(ji2)=bc(ij)                                                  8d2s19
       end do                                                           4d2s19
      end do                                                            4d2s19
      ivintt=ivintref(isymmrci)                                         7d15s19
      nsing=0                                                           6d26s18
      ndoub=0
      if(ethres.ne.0d0.or.iabs(iunc(1)).ne.0)then                       7d15s19
      else                                                              11d4s19
       do isb=1,nsymb                                                   11d4s19
        nct1(isb)=0                                                     11d4s19
       end do                                                           11d4s19
      end if                                                            7d15s19
c
c     just in case variable doesn't get set later
c
      thrs=1d-6                                                         5d15s19
c
c     for n-1 electron fcns
c
      if(iunc(1).ne.0)then                                              11d24s19
       if(lprint)then                                                    2d19s19
        if(itestmrci.eq.3)write(7,*)('end of n-1')                        5d15s19
        write(6,*)(' ')
        write(6,*)('for single excitations ')
        write(6,3)                                                        6d26s18
       end if                                                            2d19s19
    3  format('sym       total           *nvirt')                       10d17s20
       do isb=1,nsymb                                                    6d26s18
        nlh=nct1(isb)                                                    7d16s19
        nsv=multh(isb,isymmrci)                                          6d26s18
        nhere=nlh*nvirt(nsv)                                             7d15s19
        nsing=nsing+nhere                                                7d15s19
        nsingx=nsingx+nhere                                              7d15s19
        call addcomma4(nhere,numbers,12)                                2d2s21
        call addcomma4(nct1(isb),numbers2,9)                            2d3s21
        if(lprint)write(6,4)isb,numbers2(1:9),numbers                   2d3s21
    4   format(i3,3x,a9,5x,a50)                                         2d3s21
       end do                                                            6d26s18
       if(lprint)then                                                    2d19s19
        write(6,*)(' ')                                                  3d18s19
        call addcomma4(nsing,numbers,0)                                 2d2s21
        write(6,2292)numbers                                            3d12s21
 2292   format(' total no. of single excitations: ',a50)                3d12s21
        write(6,*)                                                      3d12s21
       end if                                                            2d19s19
      end if                                                            11d24s19
      ndoub=0                                                           8d24s18
      mdoub=0                                                           12d21s20
      if(ethred.ne.0d0.or.iabs(iunc(2)).ne.0)then                       5d12s20
       if(ethred.ne.0d0)then                                            8d1s19
c
c     free up memory by putting coupling coefficients into distributed  12d6s20
c     memory.                                                           12d6s20
c
        ncupc=iafgencsf3-ib4gencsf3                                     12d6s20
        i18=1                                                           12d6s20
        in8=ncupc                                                       12d6s20
        call ddi_create(bc,ibc,i18,in8,ncupc8)                          11d15s22
        call ilimts(1,ncupc,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   12d6s20
        istrt=ib4gencsf3+il-1                                           12d6s20
        i28=il                                                          12d6s20
        in8=ih                                                          12d6s20
        call ddi_put(bc,ibc,ncupc8,i18,i18,i28,in8,bc(istrt))           11d15s22
        ibcoff=ib4gencsf3                                               12d7s20
        call second(time1)
        ib4=ibcoff                                                        7d10s20
        call stepwisecsf(ixsoo,spin,xms,dum,idum,idum,1,idata,bc,ibc)   11d9s22
        call addcomma8(ihighwtr,numbers,0)                              2d2s21
        if(lprint)write(6,*)('ihighwtr after stepwisecsf: '),numbers    2d2s21
        call second(time1bx)
        telap=time1bx-time1-tovr
        nff2x=ibcoff                                                    1d20s21
        ibcoff=nff2x+nsymb*mdoop*4                                      1d20s21
        call enough('mrci. 15',bc,ibc)
        if(lprint)write(6,*)('time for stepwisecsf = '),telap
        call genfcn2(ibc(ibasis(1)),nbasis,ibc(iptr),ibc(idorb),          11d12s20
     $     ibc(isorb),nfcn,nec,mdon,mdoo,multh,nbasis2,nneed,ibc(nff2x),12d7s20
     $     nsymb,iff2x,iptrto,nct2,ibc(icsf),bc,ibc)                    11d11s22
        thrs=1d1**ethred                                                8d1s19
        if(lprint)write(6,*)('thrs,ethred: '),thrs,ethred               8d26s22
        do isb=1,nsymb                                                    4d25s21
         jeigvspace(isb)=ieigvspace(isb)                                  4d25s21
         nrootu=irxinfo(isb)                                              4d30s21
         if(isb.eq.isymmrci)nrootu=max(nrootu,nroot)                      4d30s21
         jvintref(isb)=jeigvspace(isb)+nrootu                             4d30s21
        end do                                                            4d25s21
        call cont2ecsf(ibc(idorb),ibc(isorb),mdon,mdoo,ibasis,iptr,      7d26s19
     $      ibc(icsf),ibc(icsf2),nfcn,nct,nec,ibc(icsfpd),multh,        8d14s19
     $       jvintref,irwt,nsymb,thrs,                                  5d10s21
     $       nct2,nct2c,lprint,iff2,iptrto,ibc(nff2),nfdat,ibc(idet),   12d18s20
     $       ibc(nff22),icsfxv,ipaddout,bc,ibc,nder,idvi,ibc(idvmt),    8d7s24
     $       idervcv)                                                   8d7s24
        do isb=1,nsymb                                                  12d8s20
         if(nfdat(5,1,isb).ne.0)then
          ipack8=ibc(nfdat(5,1,isb))                                    12d8s20
         end if
        end do                                                          12d8s20
        nvectr=0                                                        11d18s20
        iorg=ibcoff                                                     11d19s20
        do isb=1,nsymb                                                  11d18s20
         ncd(1,isb,1)=nfdat(3,1,isb)                                    11d18s20
         nvectr=max(nvectr,nfdat(4,1,isb)+nfdat(3,1,isb)*nfdat(2,1,isb))11d18s20
         if(nfdat(5,1,isb).ne.0)then                                    11d19s20
          iorg=min(iorg,nfdat(5,1,isb))                                 11d19s20
         end if                                                         11d19s20
         do l=1,3                                                       11d18s20
          lp=l+1                                                        11d18s20
          nvectr=                                                        11d18s20
     $       max(nvectr,nfdat(4,lp,isb)+nfdat(3,lp,isb)*nfdat(2,lp,isb))11d18s20
          ncd(l,isb,2)=nfdat(3,lp,isb)                                  11d18s20
         end do                                                         11d18s20
         nvectr=nvectr+idervcv(2,isb)                                   8d7s24
        end do                                                          11d18s20
        nvectr=nvectr-iorg                                              11d19s20
        if(lprint)then                                                  1d25s21
         write(6,*)('memory mapping so far:')
         write(6,*)('space we need for coupling coefficients: '),
     $       ib4gencsf3,('to '),iafgencsf3
         write(6,*)('space we have for contraction coefficients: '),
     $       iorg,('to '),iorg+nvectr-1
        end if                                                          1d25s21
        ioff=iafgencsf3-iorg                                            12d7s20
        if(iorg.gt.iafgencsf3)then                                      12d7s20
         do i=0,nvectr-1                                                12d7s20
          bc(iafgencsf3+i)=bc(iorg+i)                                   12d7s20
         end do                                                         12d7s20
        else                                                            12d7s20
         do i=nvectr-1,0,-1                                             12d7s20
          bc(iafgencsf3+i)=bc(iorg+i)                                   12d7s20
         end do                                                         12d7s20
        end if                                                          12d7s20
        ibcoff=iafgencsf3+nvectr                                        12d7s20
        do isb=1,nsymb                                                  12d7s20
         if(nfdat(5,1,isb).ne.0)then                                    12d7s20
          nfdat(5,1,isb)=nfdat(5,1,isb)+ioff                            12d7s20
          idervcv(1,isb)=idervcv(1,isb)+ioff                            8d7s24
          ipack8=ibc(nfdat(5,1,isb))                                    12d8s20
          do l=1,4                                                      12d7s20
           if(nfdat(4,l,isb).ne.0)then                                  12d7s20
            nfdat(4,l,isb)=nfdat(4,l,isb)+ioff                           12d7s20
           end if                                                       12d7s20
          end do                                                        12d7s20
         end if                                                         12d7s20
        end do                                                          12d7s20
        if(lprint)write(6,*)('now reload coupling coefficients')        1d25s21
        i18=1                                                           12d6s20
        if(nder.gt.0)then                                               9d5s23
         ibufvcvds=ibcoff                                               9d5s23
        end if                                                          9d5s23
        in8=ncupc                                                       12d6s20
        call ddi_get(bc,ibc,ncupc8,i18,i18,i18,in8,bc(ib4gencsf3))      11d15s22
        call dws_synca                                                  12d7s20
        call ddi_destroy(ncupc8)                                        12d7s20
        call dws_synca                                                  12d7s20
       else                                                             9d18s19
        if(itestmrci.eq.3.and.lprint)then                               9d18s19
         do isb=1,nsymb                                                 9d18s19
          write(7,*)isb,0,0,0,0                                         9d18s19
         end do                                                         9d18s19
        end if                                                          9d18s19
       end if                                                           8d1s19
        call countem2(ibc(nff2),nsymb,ncdx,mdon,mdoo,ibc(icsf2),multh,  11d19s20
     $       nvirt,nvxv,bc,ibc)                                         11d14s22
       if(ethred.eq.0d0)then                                            8d2s19
        do isb=1,nsymb                                                  8d2s19
         ncd(1,isb,1)=ncdx(1,isb,1)                                     8d2s19
         do i=1,3                                                       8d2s19
          ncd(i,isb,2)=ncdx(i,isb,2)                                    8d2s19
         end do                                                         8d2s19
        end do                                                          8d2s19
       end if                                                           8d2s19
       if(lprint)then                                                   8d2s19
        write(6,*)(' ')
        write(6,*)('for double excitations ')                             6d26s18
        write(6,*)('     singlet     triplet coupled')                  1d19s23
        write(6,*)                                                      10d17s20
     $       ('sym  coupled lows mids highs   total    *virt pairs')    10d17s20
       end if                                                           8d2s19
       ndoubx=0
       do isb=1,nsymb                                                    6d26s18
        nlh=ncd(1,isb,2)+ncd(2,isb,2)+ncd(3,isb,2)                      8d2s19
        mlh=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)                12d21s20
        nsv=isb                                                         8d27s19
        ndoub=ndoub+ncd(1,isb,1)*nvxv(1,nsv)                            9d25s18
        mdoub=mdoub+nfdat(2,1,isb)*nvxv(1,nsv)+mlh*nvxv(2,nsv)          12d21s20
        i0a=ncdx(1,isb,1)                                               10d18s20
        i0b=nvxv(1,nsv)                                                 10d18s20
        ndoubx=ndoubx+i0a*i0b                                           10d18s20
        nhere=ncd(1,isb,1)*nvxv(1,nsv)                                  9d25s18
        i0a=ncdx(1,isb,2)+ncdx(2,isb,2)+ncdx(3,isb,2)                   10d18s20
        i0b=nvxv(2,nsv)                                                 10d18s20
        ndoub=ndoub+(ncd(1,isb,2)+ncd(2,isb,2)+ncd(3,isb,2))*nvxv(2,nsv)             9d25s18
        ndoubx=ndoubx+i0a*i0b                                           10d18s20
        nhere=nhere+(ncd(1,isb,2)+ncd(2,isb,2)+ncd(3,isb,2))*nvxv(2,nsv)             9d25s18
        call addcomma4(nhere,numbers,8)                                 2d2s21
        if(lprint)write(6,7)isb,ncd(1,isb,1),(ncd(j,isb,2),j=1,3),nlh,  2d19s19
     $       numbers                                                    2d2s21
    7   format(i3,3x,i5,2x,3i5,i8,5x,a50)                               2d2s21
       end do                                                            6d26s18
       if(lprint)then
        write(6,*)(' ')
        call addcomma4(ndoub,numbers,0)                                 2d2s21
        write(6,707)numbers                                             2d26s21
  707   format('total number of double excitations: ',a50)              6d7s23
        if(iunc(2).eq.0)then                                            10d22s21
         call addcomma4(mdoub,numbers,0)                                 2d2s21
         write(6,708)numbers                                            10d22s21
  708    format('space required for non-orthogonal doubles: ',a50)      12d2s22
         call addcomma8(ndoubx,numbers,0)                                2d2s21
         write(6,*)(' ')                                                 5d7s21
         write(6,709)numbers                                             2d26s21
  709    format('total number of uncontracted double excitations: ',      2d26s21
     $       a50)                                                       2d26s21
         ndoubx=ndoubx+nsingx+nctf(isymmrci)                             2d2s21
         call addcomma8(ndoubx,numbers,0)                                2d2s21
         write(6,710)numbers                                             2d26s21
  710    format('total number of uncontracted csfs: ',a50)                2d26s21
        end if                                                          10d22s21
       end if                                                           2d19s19
      end if                                                            7d15s19
      ntot=nsing+ndoub
      ntotr=ntot+nctf(isymmrci)                                         4d8s20
      if(lprint)then                                                    1d3s20
       write(6,*)                                                       2d2s21
       call addcomma4(ntot,numbers,0)                                   2d2s21
       write(6,711)numbers                                              2d26s21
  711  format('total number of single and double excitations: ',a50)     2d26s21
       call addcomma4(ntotr,numbers,0)                                  2d2s21
       write(6,712)numbers                                              2d26s21
  712  format('with reference part: ',a50)                              2d26s21
      end if                                                            1d3s20
c
c     diagonals
c
      if(ntot.ne.0)then                                                 9d18s18
      call second(timaa)
      jh0copy=ih0copy                                                   1d31s21
      jh0=ih0                                                           1d31s21
      do isb=1,nsymb                                                    1d31s21
       nbasdwsf=nbasisp(isb)*ncomp                                      1d31s21
       do i=0,nbasdwsf*nbasdwsf-1                                       1d31s21
        bc(jh0+i)=bc(jh0copy+i)                                         1d31s21
       end do                                                           1d31s21
       jh0copy=jh0copy+nbasdwsf*nbasdwsf                                1d31s21
       jh0=jh0+nbasdwsf*nbasdwsf                                        1d31s21
      end do                                                            1d31s21
      ibcb4=ibcoff                                                      7d7s21
      if(iunc(2).eq.0)then                                              8d18s23
       ifull4x=2                                                        8d18s23
      else                                                              8d18s23
       ifull4x=3                                                        8d18s23
      end if                                                            8d18s23
      ipair=0                                                           2d26s24
      call paraeri(natom,ngaus,ibdat,idum,ihmat,iorb,noc,               4d24s20
     $             ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,  11d29s12
     $     iter1,idwsdeb,idorel,ascale,nbasisp,ifull4x,0,idum,0,idum,   5d20s24
     $     idum,bc,ibc)                                                 5d20s24
      icol=ibcoff
      ibcoff=icol+mynprocg
      imsg=ibcoff                                                       6d3s10
      ibcoff=imsg+mynnode                                               6d3s10
      iter1=0
      idwsdeb=0
      call parajkfromh(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,         11d29s12
     $     iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      5d30s18
     $     ioooo,ionex,irefo,iter1,idwsdeb,ncomp,nvirt,i3x,1,ih0,       5d31s18
     $     idoubo,shiftdup,nbasisp,idorel,ifull4x,i4x,i4xb,ibcb4,bc,ibc,8d18s23
     $     isend1,isend2,isend3,isend4,0,iff)                           7d11s23
      call second(timaaa)
      telap=timaaa-timaa-tovr
      call trans2e(jmats,kmats,ionex,i3x,i3xb,irefo,iunc(2),ionexb,     4d24s20
     $     jmatt,kmatt,ionexc,kmatd,i3x3,ionexbt,bc,ibc)                11d10s22
      jh0=ih0                                                           6d5s18
      do isb=1,nsymb                                                    6d5s18
        ih0a(isb)=ibcoff                                                6d5s18
        ibcoff=ih0a(isb)+irefo(isb)*irefo(isb)                          6d5s18
        call enough('mrci. 16',bc,ibc)
        ndim=nbasdws(isb)-idoubo(isb)                                   6d5s18
        sym=0d0
        do i=0,irefo(isb)-1                                             6d5s18
         do j=0,irefo(isb)-1                                            6d5s18
          iad1=jh0+j+ndim*i                                             6d5s18
          iad2=ih0a(isb)+j+irefo(isb)*i                                 6d5s18
          bc(iad2)=bc(iad1)                                             6d5s18
          if(j.ne.i)then
           iad1t=jh0+i+ndim*j                                             6d5s18
           diff=bc(iad1t)-bc(iad1)                                      10d28s20
           sym=sym+diff**2                                              10d28s20
          end if
         end do                                                         6d5s18
        end do                                                          6d5s18
        ih0av(isb)=jh0
        nh0av(isb)=ndim                                                 6d28s18
        jh0=jh0+ndim*ndim                                               6d5s18
      end do                                                            6d5s18
c
c     replicate oooo on all procs
c
      do i=1,mynprocg                                                   6d1s18
       isend1(i)=0                                                      6d1s18
      end do                                                            6d1s18
      ntot=0                                                            6d1s18
      do is=1,nsdlk                                                     6d1s18
       if(isblk(1,is).eq.isblk(2,is))then                               6d1s18
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              6d1s18
       else                                                             6d1s18
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      6d1s18
       end if                                                           6d1s18
       if(isblk(3,is).eq.isblk(4,is))then                               6d1s18
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              6d1s18
       else                                                             6d1s18
        ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       6d1s18
       end if                                                           6d1s18
       if(nrow.gt.0)then                                                6d1s18
        ntot=ntot+ncol*nrow                                              6d1s18
        do ip=0,mynprocg-1                                               6d1s18
         ipp=ip+1                                                        6d1s18
         call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,ip,  6d1s18
     $       il,ih,i1s,i1e,i2s,i2e)                                     6d1s18
         i10=i1s                                                        6d1s18
         i1n=irefo(isblk(3,is))                                         6d1s18
         nhere=0                                                        6d1s18
         do i2=i2s,i2e                                                  6d1s18
          if(i2.eq.i2e)i1n=i1e                                          6d1s18
          do i1=i10,i1n                                                 6d1s18
           if(i1.ge.i2.or.isblk(3,is).ne.isblk(4,is))nhere=nhere+1      6d1s18
          end do                                                        6d1s18
          i10=1                                                         6d1s18
         end do                                                         6d1s18
         isend1(ipp)=isend1(ipp)+nhere*nrow                              6d1s18
        end do                                                           6d1s18
       end if                                                           6d1s18
      end do                                                            6d1s18
      isum=0
      do i=1,mynprocg
       isum=isum+isend1(i)
      end do
      idest=ibcoff                                                      6d1s18
      ibcoff=idest+ntot                                                 6d1s18
      ibufr=ibcoff                                                      6d1s18
      ibcoff=ibufr+ntot                                                 6d1s18
      call enough('mrci. 17',bc,ibc)
      jdest=idest                                                       6d1s18
      do is=1,nsdlk                                                     6d1s18
       if(isblk(1,is).eq.isblk(2,is))then                               6d1s18
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              6d1s18
       else                                                             6d1s18
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      6d1s18
       end if                                                           6d1s18
       if(nrow.gt.0)then                                                6d1s18
        call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,      6d1s18
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            6d1s18
        i10=i1s                                                         6d1s18
        i1n=irefo(isblk(3,is))                                          6d1s18
        i4o=ioooo(is)                                                   6d1s18
        do i2=i2s,i2e                                                   6d1s18
         if(i2.eq.i2e)i1n=i1e                                           6d1s18
         do i1=i10,i1n                                                  6d1s18
          if(i1.ge.i2.or.isblk(3,is).ne.isblk(4,is))then                6d1s18
           do irow=0,nrow-1                                             6d1s18
            bc(jdest+irow)=bc(i4o+irow)                                 6d1s18
           end do                                                       6d1s18
           jdest=jdest+nrow                                             6d1s18
          end if                                                        6d1s18
          i4o=i4o+nrow                                                  6d1s18
         end do                                                         6d1s18
         i10=1                                                          6d1s18
        end do                                                          6d1s18
       end if                                                           6d1s18
      end do                                                            6d1s18
      do i=1,mynprocg                                                   6d1s18
       isend2(i)=0                                                      6d1s18
       isend3(i)=isend1(mynowprog+1)                                    6d1s18
      end do                                                            6d1s18
      isend4(1)=0                                                       6d1s18
      do i=1,mynprocg-1                                                 6d1s18
       im=i-1                                                           6d1s18
       isend4(i+1)=isend4(i)+isend1(i)                                  6d1s18
      end do                                                            6d1s18
      call dws_all2allvb(bc(idest),isend3,isend2,bc(ibufr),             6d1s18
     $     isend1,isend4)                                               6d1s18
      do i=0,mynprocg-1                                                 6d1s18
       iarg=isend1(i+1)                                                  6d1s18
       jbufr=ibufr+isend4(i+1)
       isend4(i+1)=jbufr                                                6d1s18
      end do                                                            6d1s18
      do is=1,nsdlk                                                     6d1s18
       ioooo(is)=idest                                                  6d1s18
       if(isblk(1,is).eq.isblk(2,is))then                               6d1s18
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              6d1s18
       else                                                             6d1s18
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      6d1s18
       end if                                                           6d1s18
       if(isblk(3,is).eq.isblk(4,is))then                               6d1s18
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              6d1s18
       else                                                             6d1s18
        ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       6d1s18
       end if                                                           6d1s18
       if(nrow*ncol.gt.0)then                                           6d1s18
        idest=idest+nrow*ncol                                            6d1s18
        do i=0,mynprocg-1                                                6d1s18
         ip=i+1                                                          6d1s18
         call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,i,   6d1s18
     $       il,ih,i1s,i1e,i2s,i2e)                                     6d1s18
         i10=i1s                                                        6d1s18
         i1n=irefo(isblk(3,is))                                         6d1s18
         do i2=i2s,i2e                                                  6d1s18
          if(i2.eq.i2e)i1n=i1e                                          6d1s18
          do i1=i10,i1n                                                 6d1s18
           if(i1.ge.i2.or.isblk(3,is).ne.isblk(4,is))then               6d1s18
            if(isblk(3,is).eq.isblk(4,is))then                          6d1s18
             ito=ioooo(is)+nrow*(((i1*(i1-1))/2)+i2-1)                   6d1s18
            else                                                        6d1s18
             ito=ioooo(is)+nrow*(i1-1+irefo(isblk(3,is))*(i2-1))        6d1s18
            end if                                                      6d1s18
            do irow=1,nrow                                              6d1s18
             bc(ito)=bc(isend4(ip))                                     6d1s18
             ito=ito+1                                                     6d1s18
             isend4(ip)=isend4(ip)+1                                       6d1s18
            end do                                                      6d1s18
           end if                                                       6d1s18
          end do                                                        6d1s18
          i10=1                                                         6d1s18
         end do                                                         6d1s18
        end do                                                           6d1s18
       end if                                                           6d1s18
      end do                                                            6d1s18
      ibcoff=ibufr                                                      6d1s18
c
c     fold -0.5*(ij|jl) into h0a(il) for resolution of identity         4d8s20
c
      do isb=1,nsymb                                                    4d8s20
       ih0ae(isb)=ibcoff                                                4d8s20
       ibcoff=ih0ae(isb)+irefo(isb)*irefo(isb)                          4d24s20
       call enough('mrci. 18',bc,ibc)
       do i=0,irefo(isb)*irefo(isb)-1                                   4d8s20
        bc(ih0ae(isb)+i)=bc(ih0a(isb)+i)                                4d8s20
       end do                                                           4d8s20
       do i=0,irefo(isb)-1                                              4d8s20
        do l=0,irefo(isb)-1                                             4d8s20
         iad=ih0ae(isb)+l+irefo(isb)*i                                  4d8s20
         do jsb=1,nsymb                                                   4d8s20
          do j=0,irefo(jsb)-1                                             4d8s20
           term=getint(ioooo,isb,jsb,jsb,isb,i,j,j,l,bc,ibc)            11d15s22
           orig=bc(iad)
           bc(iad)=bc(iad)-0.5d0*term                                   4d29s20
          end do                                                        4d8s20
         end do                                                         4d8s20
        end do                                                          4d8s20
       end do                                                           4d8s20
      end do
       nlocals=0                                                        3d19s21
c                                                                       6d17s21
c     npadddi is number of extra words require for ddi_ixxx buffer      6d17s21
c                                                                       6d17s21
       npadddi=4*mynprocg+mynprocg                                               6d17s21
       ihsdiag=ibcoff                                                   2d10s23
       if(nsing.gt.0)then                                               9d18s18
        ibcoff=ihsdiag+nsymb*mdoop*2                                    11d19s20
        call enough('mrci. 19',bc,ibc)
        call prephss(ibc(ihsdiag),nsymb,mdon,mdoo,ibc(nff1),ibc(icsf),  11d19s20
     $       nvirt,isymmrci,multh,nroot,nlocals,maxbx,bc,ibc)           11d10s22
        maxbx=maxbx+npadddi                                             6d17s21
       end if                                                           11d26s20
       nlocald=0                                                        7d8s21
       iddoub=ibcoff                                                    2d1s23
       ihddiag=ibcoff                                                   2d10s23
       if(ndoub.gt.0)then                                               11d26s20
        if(iunc(2).eq.0)then                                            11d26s20
         iddoub=ibcoff                                                     11d19s20
         ibcoff=iddoub+ndoub                                            1d22s21
         call enough('mrci. 20',bc,ibc)
         do i=iddoub,ibcoff-1                                              11d19s20
          bc(i)=0d0                                                        11d19s20
         end do                                                            6d28s18
         nlocald=ndoub*nroot                                            6d9s21
        else                                                            11d26s20
         ihddiag=ibcoff                                                  11d19s20
         ibcoff=ihddiag+nsymb*nsymb*mdoop*2                             11d26s20
         call enough('mrci. 21',bc,ibc)
         call prephdd(ibc(ihddiag),nsymb,mdon,mdoo,ibc(nff2),ibc(icsf), 5d28s21
     $        ibc(icsf2),nvirt,isymmrci,multh,nroot,maxbxd,nlocald,bc,  11d10s22
     $        ibc)                                                      11d10s22
         maxbxd=maxbxd+npadddi                                          6d17s21
         mdoub=0                                                        6d9s21
        end if                                                          11d26s20
       end if                                                           11d26s20
       call dws_synca                                                   11d26s20
       call ddi_iaccword0                                               8d11s22
       nwiacc=0                                                         8d11s22
       issdig=ibcoff                                                    4d19s23
       if(nsing.gt.0)then                                               11d26s20
        if(ldebug)write(6,*)('calling hsscsf'),ixw2,ibc(32800)                          1d3s20
        call hsscsf(1,ibc(ihsdiag),nsing,ih0av,nh0av,ibc(nff1),         11d19s20
     $       ibc(iff1),ibc(icsf),nct1,                                  11d19s20
     $      mdon,mdoo,jmats,kmats,ntot,ioooo,multh,nec,shift0,          2d9s21
     $       nsymb,idum,tdenss,tovr,ixw1,ixw2,issdig,ldebugf,bc,ibc)    11d10s22
        if(ldebug)then                                                  4d8s20
         write(6,*)('back from hss')                                    4d8s20
        end if
       end if                                                           9d18s18
       idddig=ibcoff                                                    2d10s23
       if(ndoub.ne.0)then                                                6d28s18
        if(ldebug)write(6,*)('calling hcddcsfjk for diagonals')         1d3s20
        call second(timedd)                                             9d26s20
        timeddd=0d0                                                     1d19s23
        timeddde=0d0                                                    1d19s23
        if(iunc(2).eq.0)then                                            10d19s20
         call enough('mrci. 22',bc,ibc)
         do i=0,ndoub-1                                                 10d19s20
          bc(iddoub+i)=0d0                                              10d19s20
         end do                                                         10d19s20
         if(ldebug)write(6,*)('calling hcdddiag '),ibcoff
         call second(timeddd)                                           10d28s20
         call hcdddiag(ibc(nff2),                                       11d25s20
     $       nsymb,mdon,mdoo,multh,isymmrci,nvirt,ibc(icsf),nec,        11d19s20
     $       ibc(icsf2),nfdat,irel,ism,ih0av,nh0av,ioooo,shift0,kmats,  2d9s21
     $       jmats,irefo,ismult,bc(iddoub),ixw1,ixw2,norb,              1d25s21
     $        ldebugf,bc,ibc)                                           11d10s22
         call second(timeddde)                                           10d28s20
         if(ldebug)write(6,*)('back from hcddding'),ibcoff,
     $        timeddde-timeddd-tovr
        else                                                            11d26s20
         call hddcsf(ibc(ihddiag),ndoub,ih0av,nh0av,ibc(nff2),          4d9s21
     $       ibc(iff2),mdon,mdoo,jmats,kmats,ioooo,multh,nec,shift0,    5d28s21
     $       nsymb,idddig,ldebugf,bc,ibc)                                11d10s22
        end if                                                          10d19s20
        call second(timedde)                                            9d26s20
        telap=timedde-timedd-tovr                                       9d26s20
        telapd=timeddde-timeddd-tovr                                    10d28s20
        if(ldebug)write(6,*)('back')                                    1d3s20
       end if                                                            6d28s18
      end if                                                            6d28s18
       if(maxds.eq.1)then                                               12d29s21
        ldavid=.false.                                                  12d29s21
       else                                                             12d29s21
        ldavid=.true.                                                    12d27s21
       end if                                                           12d29s21
       maxdsr=maxds*nroot                                               12d27s21
      ihdi=ibcoff
      ibcoff=ihdi+nroot*ndoub                                           11d25s20
      call enough('mrci. 23',bc,ibc)
      do i=ihdi,ibcoff-1                                                4d17s20
       bc(i)=0d0                                                        4d17s20
      end do                                                            4d17s20
      if(nsing.ne.0)then
       if(itestmrci.eq.3.and.lprint)then                                10d1s19
        write(6,*)('ss diagonal after gs: ')
        jdsing=idsing                                                   7d15s19
        do isb=1,nsymb
         nsv=multh(isb,isymmrci)
         if(nct1(isb).gt.0)then
          write(6,*)('for symmetry block: '),isb,nsv
          call prntm2(bc(jdsing),nvirt(nsv),nct1(isb),nvirt(nsv))       7d15s19
         end if
         jdsing=jdsing+nvirt(nsv)*nct1(isb)                             7d15s19
        end do
       end if                                                           5d17s19
       if(ldebug)write(6,*)('calling hcsicsf')                          1d3s20
       call hcsi(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(nffr),            11d25s20
     $      ibc(iffr),bc(ivintt),ncsfint,ibc(icsf),nec,mdon,mdoo,       11d25s20
     $      nsymb,multh,ixw1,ixw2,ih0av,nh0av,ionexb,nvirt,1,           1d25s21
     $      ldebugf,maxbx,ionexbt,nwiacc,bc,ibc)                        11d10s22
       if(ldebug)write(6,*)('back ')                                    1d3s20
      end if
      if(ndoub.ne.0)then                                                9d18s18
       if(ldebug)write(6,*)('calling hcdi'),ihdi,ibcoff                 1d3s20
       if(iunc(2).eq.0)then                                             6d7s21
        call hcdi(ibc(nff2),nsymb,mdon,mdoo,multh,isymmrci,nvirt,
     $      ibc(icsf),nec,ibc(icsf2),nfdat,irel,ism,kmats,irefo,ismult, 12d8s20
     $      bc(ihdi),ixw1,ixw2,norb,ibc(nffr),ibc(iffr),bc(ivintt),     12d8s20
     $      ncsfint,nroot,ldebugf,srh,bc,ibc)                           11d10s22
        if(ldebug)write(6,*)('calling prehcddjk'),ibcoff                11d7s22
        idenput=ibcoff                                                  5d17s23
        ibcoff=idenput+nsymb*nsymb*2                                    5d17s23
        call enough('mrci.denput',bc,ibc)                               5d17s23
        call prehcddjk(ibc(nff22),nfdat,nsymb,mdon,mdoo,nec,multh,      11d7s22
     $       ibc(icsf),ibc(icsf2),irel,ism,irefo,ismult,ixw1,ixw2,norb, 11d7s22
     $       ih0av,nh0av,ioooo,shift0,iden1e,idenhvv,idenj,idenk,bc,ibc,5d17s23
     $       ibc(idenput),idenorg)                                      5d17s23
       else                                                             6d7s21
        call hcdiuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ibc(nffr),            11d25s20
     $      ibc(iffr),bc(ivintt),ncsfint,ibc(icsf),ibc(icsf2),nec,mdon, 6d7s21
     $       mdoo,nsymb,multh,ixw1,ixw2,kmats,nvirt,1,                  6d7s21
     $      ldebugf,maxbxd,srh,nwiacc,bc,ibc)                           11d10s22
       end if                                                           6d7s21
       if(ldebug)write(6,*)('back '),ibcoff                             12d13s21
      end if                                                            9d18s18
      if(ntot.gt.0)then                                                 9d13s18
       ieigintxyz=ibcoff                                                   10d7s21
       ibcoff=ieigintxyz+nroot                                             10d7s21
       ivintt=ibcoff                                                    10d7s21
       ibcoff=ivintt+nroot*nctf(isymmrci)                                10d7s21
       ivints=ibcoff                                                    10d7s21
       ibcoff=ivints+nroot*(nctf(isymmrci)+nroot+maxdsr)                12d29s21
       call enough('mrci. 24',bc,ibc)
       ieold=ibcoff                                                     11d13s18
       ibcoff=ieold+nroot                                               11d13s18
       bintvec=intvec.eq.0
       if(bintvec)then                                                   1d17s20
        intvec=ibcoff                                                   1d17s20
        nneed=nctf(isymmrci)*nroot                                      4d8s20
        call ilimts(nneed,1,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s,m2e)   1d17s20
        nhneed=mh+1-ml                                                  1d17s20
        call ilimts(nneed,1,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s,m2e)   1d17s20
        nhneed=nhneed+mh+1-ml                                           1d21s20
        nhneed=nhneed+nroot*nroot*2                                     1d17s20
        ibcoff=intvec+nhneed                                            1d17s20
       end if                                                            1d17s20
       ivdinout=ibcoff                                                      9d13s18
       ivdold=ivdinout+nlocald                                          12d8s21
       igdold=ivdold+nlocald*maxds                                      12d27s21
       ivsold=igdold+nlocald*maxds                                      12d27s21
       igsold=ivsold+nlocals*maxds                                      12d27s21
       igd=igsold+nlocals*maxds                                         12d27s21
       ivdnono=igd+mdoub*nroot                                          12d21s20
       ihixold=ivdnono+mdoub*nroot                                      1d23s21
       ibcoff=ihixold+nctf(isymmrci)*nroot*maxds                        12d27s21
       igss=ibcoff                                                      2d17s21
       igsso=igss+nlocals                                               2d17s21
       ihis=igsso+nlocals*maxds                                         1d24s23
       ihiso=ihis+nctf(isymmrci)*nroot                                  2d19s21
       idotsing=ihiso+nctf(isymmrci)*nroot*(maxds+1)                    1d24s23
       ibcoff=idotsing+nroot*2                                          2d18s21
       call enough('mrci. 25',bc,ibc)
       do i=ivdinout,ibcoff-1                                           1d24s21
        bc(i)=0d0                                                       1d24s21
       end do                                                           1d24s21
       ilastde=ibcoff                                                   10d22s21
       ibcoff=ilastde+nroot                                             10d22s21
       iamconverged=ibcoff                                              4d21s21
       ibcoff=iamconverged+nroot                                        4d21s21
       nlastde=ibcoff                                                   10d22s21
       icombo=nlastde+nroot                                             12d7s21
       ibcoff=icombo+nroot                                              12d7s21
       call enough('mrci. 26',bc,ibc)
       do i=iamconverged,ibcoff-1                                       4d21s21
        ibc(i)=0                                                        4d21s21
       end do                                                           4d21s21
       icit=0                                                           11d13s18
       if(nroot.eq.1.and.nrestart.ge.0)econvcii=1d-3                    7d17s24
       ibc1000=ibcoff                                                   11d13s18
       do i=1,9                                                         3d9s21
        timex(i)=0d0                                                    12d14s18
       end do                                                           12d14s18
       do i=1,10                                                        5d22s19
        tpart(i)=0d0                                                    5d22s19
       end do                                                           5d22s19
       i18=1                                                            12d12s20
       in8=nctf(isymmrci)*nroot                                         12d12s20
       call ddi_create(bc,ibc,i18,in8,ihi)                              11d15s22
       if(iunc(2).eq.0)then                                             7d7s21
        in8=mdoub*nroot                                                  1d18s21
        call ddi_create(bc,ibc,i18,in8,igxx)                            11d15s22
       end if                                                           7d7s21
       iepart=ibcoff                                                    7d12s21
       nroote=maxdsr                                                    1d12s23
       ibcoff=iepart+nroote*4*nsymb                                     1d12s23
       ividotvi=ibcoff                                                  3d31s21
       ivihvi=ividotvi+nroot                                            3d31s21
       ibcoff=ivihvi+nroot                                              3d31s21
       ipairs=ibcoff                                                    8d19s21
       ibcoff=ipairs+4*nroot                                            8d19s21
       ivstry=ibcoff                                                    11d10s20
       ibcoff=ivstry+nlocals                                            11d10s20
       call enough('mrci. 27',bc,ibc)
       ngots=0                                                          4d14s21
       ngotd=0                                                          4d14s21
       tdenss=0d0                                                        5d17s19
       tdends=0d0                                                       5d17s19
       tdendd=0d0                                                       5d17s19
       echange=0d0                                                      1d24s21
       nok2conv=0                                                       4d15s21
       ndoubs=ndoub                                                     4d14s21
       if(lprint)then                                                   9d28s20
        write(6,*)(' ')                                                 12d13s19
        write(6,717)                                                    12d13s19
  717   format('iteration  best energies(rounded!)')                    10d28s21
       end if                                                           9d28s20
       if(mynowprog.eq.0.and.nwavef.ne.0)then                           8d19s22
        call dumpbas(nfcnf(isymmrci),ibc(jptrf),ibc(ibasisf(isymmrci)), 8d16s22
     $      ibc(idorbf),ibc(isorbf),nodu,nosu,mdoop,nec,1,              7d8s21
     $       nctf(isymmrci),ibc(icsf),mdon,nwavrec,nrestart,bc,ibc)     11d14s22
       end if                                                           8d16s22
 1000  continue                                                         11d13s18
       call dws_synca                                                   1d25s21
       icit=icit+1                                                      11d13s18
       if(icit.eq.itabort)then                                          12d9s22
        if(lprint)write(6,*)('aborting calculation ')                   12d9s22
        call dws_synca                                                  12d9s22
        call dws_finalize                                               12d9s22
        stop                                                            12d9s22
       end if                                                           12d9s22
       ncors=0                                                          12d8s21
       if(icit.gt.2.and.nok2conv.eq.0.and.ndoub.eq.ndoubs)then          12d13s21
        nstag=0
       end if                                                           12d7s21
       if(icit.gt.2)then                                                6d20s24
        lconv=echange.lt.econvci                                         6d20s24
        if(wconvci.ne.0d0)then                                           6d20s24
         if(diffr.gt.wconvci)lconv=.false.                              6d20s24
        end if                                                          6d20s24
       else                                                             6d20s24
        lconv=.false.                                                   6d20s24
       end if                                                           6d20s24
       if((icit.gt.2.and.nok2conv.eq.0.and.(lconv.or.                   7d18s24
     $      ncors.eq.nroot).and.ndoub.eq.ndoubs).or.icit.gt.maxmit)then 7d18s24
        if(lprint)then
         write(6,*)(' ')                                                1d3s20
         if(nstag.eq.0.and.icit.lt.maxmit)then                          3d21s22
          write(6,*)('calculations converged ')                           12d6s21
         else if(nstag.lt.nroot.and.icit.lt.maxmit)then                 3d21s22
          write(6,*)('calculations converged for '),nroot-nstag,        12d7s21
     $         ('roots')                                                12d7s21
          write(6,*)('a total of '),nstag,('roots stagnated')           12d7s21
         else if(icit.gt.maxmit)then                                    3d21s22
          write(6,*)('maximum number of iterations exceeded!!')         3d21s22
         else                                                           12d7s21
          write(6,*)('all roots stagnated!!')                           12d7s21
         end if                                                         12d7s21
         write(6,*)('internal part ')                                   8d10s22
         call intsum(bc(ivint),nctf(isymmrci),nroot,ibc(idorbf),         1d21s22
     $       ibc(isorbf),ibc(ibasisf(isymmrci)),ibc(jptrf),              1d15s22
     $       nfcnf(isymmrci),mdon,nec,ibc(icsf),norb,irefo,nsymb)       2d10s22
         nrootu=max(nroot,irxinfo(isymmrci))                             5d27s19
         if(ndoub.eq.0)then                                              9d5s23
          write(6,221)                                                   9d5s23
  221     format(/18x,'mrci-s')                                          9d5s23
         else                                                            9d5s23
          if(iunc(2).eq.0)then                                           9d5s23
           write(6,222)                                                  9d5s23
  222      format(/15x,'i2cmrci',9x,'i2cmrci+q(f)',4x,'i2cmrci+q(r)',    9d5s23
     $       4x,'i2cmrci+q(sx)'                                         1d5s24
     $        7x,'ref',13x,'ecorr',8x,'Vref*Vi2cmrci',3x,               9d5s23
     $        'Frac Vi2cmrci in ref')                                   9d5s23
          else                                                           9d5s23
           write(6,22)                                                     10d1s19
   22      format(/18x,'mrci',12x,'mrci+q(f)',7x,'mrci+q(r)',              1d5s20
     $       6x,'mrci+q(sx)'                                            5d25s21
     $        7x,'ref',13x,'ecorr',11x,'Vref*Vmrci',6x,                 10d1s19
     $        'Frac Vmrci in ref')                                      10d1s19
          end if                                                          9d5s23
         end if                                                         9d6s23
        end if                                                          1d18s20
        if(mynowprog.eq.0)then                                          7d8s21
         if(iwopen.eq.0.and.nwavef.ne.0)then                            8d24s22
          open(unit=2,file='wavef',form='unformatted',status='old')     8d16s22
          rewind(2)                                                     8d16s22
          do i=1,nwavrec                                                8d16s22
           read(2)                                                      8d16s22
          end do                                                        8d16s22
         end if                                                         8d16s22
        end if                                                          7d8s21
        itimer=ibcoff                                                   3d23s21
        ibcoff=itimer+18                                                3d23s21
        call enough('mrci. 28',bc,ibc)
        jtimer=itimer-1                                                 3d23s21
        ktimer=jtimer+9                                                 3d23s21
        do i=1,9                                                        3d23s21
         bc(jtimer+i)=timex(i)                                          3d23s21
         bc(ktimer+i)=timex(i)**2                                       3d23s21
        end do                                                          3d23s21
        call dws_gsumf(bc(itimer),18)                                   3d23s21
        call dws_gsumf(bc(ivihvi),nroot)                                3d31s21
        if(iunc(2).ne.0)then                                            7d12s21
         nwds=nroote*4*nsymb                                            1d12s23
         call dws_gsumf(bc(iepart),nwds)                                7d12s21
        end if                                                          7d12s21
        if(lprint)then                                                  7d21s21
        do j=1,nroot                                                    3d18s19
         jm=j-1                                                         3d18s19
         llvintref=lvintref+ncsfint*jm                                  5d24s23
         jvint=ivint+nctf(isymmrci)*jm                                  4d8s20
         dot1=0d0
         dot2=0d0
         dotq=0d0                                                       1d18s20
         if(ndoub.gt.0)then                                             10d25s21
          dot2=bc(ividotvi+jm)                                           3d31s21
          dotq=bc(ivihvi+jm)                                             3d31s21
          if(ndecomp.ne.0)write(6,*)('dotq '),dotq,dotq/dot2
          jdotsing=idotsing+2*jm                                         3d31s21
          dot2s=dot2+bc(jdotsing)                                        2d18s21
          doubfrac=1d0-dot2s                                             7d27s21
          if(ndecomp.ne.0)                                              8d10s22
     $         write(6,*)('doubfrac '),dot2,bc(jdotsing),bc(jdotsing+1),1d23s23
     $         doubfrac                                                 1d23s23
          dotqs=dotq+bc(jdotsing+1)                                      2d18s21
          dototq=dotqs                                                   7d12s21
          partsum=0d0                                                    7d27s21
          do l=1,4                                                      3d24s23
           ipp=ipairs+l-1+nroot*jm                                      3d24s23
           if(ndecomp.ne.0)write(6,*)l,bc(ipp)                          3d24s23
           partsum=partsum+bc(ipp)                                      3d24s23
          end do                                                        3d24s23
          do isb=1,nsymb                                                1d15s22
           do l=1,4                                                       7d12s21
            iad=iepart+jm+nroote*(l-1+4*(isb-1))                        1d12s23
            dototq=dototq+bc(iad)
            if(ndecomp.ne.0)                                            8d10s22
     $           write(6,*)l,isb,bc(iad),bc(iad)*doubfrac               3d24s23
           end do                                                       1d15s22
          end do                                                        1d15s22
          if(ndecomp.ne.0)write(6,*)('dototq,dotqs '),dototq,dotqs                  8d10s22
          dotqs=dotqs/dot2s                                              2d18s21
          if(ndecomp.ne.0)                                              8d10s22
     $         write(6,*)('eis: '),dotqs+shift,('fracis '),dot2s
          davidsx=bc(ieigint+jm)-dotqs                                   2d18s21
          if(ndecomp.ne.0)                                              8d10s22
     $         write(6,*)('doub e? '),partsum,davidsx,partsum/davidsx   8d10s22
          val=1d0
          davidsx=dotqs+(davidsx/dot2s)+shift
          dotq=dotq/dot2                                                 1d21s20
          jptr=iptr+(mdoo+1)*2*(isymmrci-1)                              10d2s19
          call dotdvecs(bc(jvint),ibc(idorbf),ibc(isorbf),               4d8s20
     $        ibc(ibasisf(isymmrci)),ibc(jptrf),nfcnf(isymmrci),        4d8s20
     $        bc(llvintref),ibc(idorb),ibc(isorb),                      5d24s23
     $        ibc(ibasis(isymmrci)),ibc(jptr),nfcn(isymmrci),mdon,nec,  10d2s19
     $        ibc(icsf),dot1)                                           11d15s21
          ecorr=bc(ieigint+jm)-bc(ieigintref+jm)
          davidq=bc(ieigintref+jm)+ecorr/(dot1**2)                       3d18s19
          davidr=bc(ieigintref+jm)+ecorr/(dot2)                          3d18s19
          write(6,1676)spinm,termsym,                                    1d5s20
     $        j,bc(ieigint+jm)+shift,davidq+shift,davidr+shift,davidsx, 5d25s21
     $        bc(ieigintref+jm)+shift,                                  2d9s21
     $        ecorr,dot1,dot2,dotq+shift                                10d7s21
 1676    format(a2,a6,'>',i3,10f16.10)                                  5d25s21
        else
          write(6,1676)spinm,termsym,                                    1d5s20
     $        j,bc(ieigint+jm)+shift                                    10d25s21
        end if                                                          10d25s21
        end do                                                          3d18s19
        open(unit=42,file='summary')                                    1d10s19
        notok=0                                                         11d21s22
        is=1                                                            11d21s22
        do i=0,nextrad-1                                                11d21s22
         ie=is+14                                                       11d21s22
         if(ie.gt.idsln)then                                            11d21s22
          write(6,*)('too many characters for sline!')                  11d21s22
          notok=1                                                       11d21s22
          go to 1524                                                    11d21s22
         else                                                           11d21s22
          write(sline(is:ie),1522)bc(iextrad+i)                          11d21s22
 1522     format(es15.7)                                                 11d21s22
          is=ie+1                                                        11d21s22
         end if                                                         11d21s22
        end do                                                          11d21s22
        ie=is+11                                                        11d21s22
        if(ie.gt.idsln)then                                             11d21s22
         notok=1                                                        11d21s22
         go to 1524                                                     11d21s22
        end if                                                          11d21s22
        write(sline(is:ie),1523)ndoub                                   11d21s22
 1523   format(i12)                                                     11d21s22
        is=ie+1                                                         11d21s22
        ie=is+1                                                         11d21s22
        if(ie.gt.idsln)then                                             11d21s22
         notok=1                                                        11d21s22
         go to 1524                                                     11d21s22
        end if                                                          11d21s22
        sline(is:ie)=spinm                                              11d21s22
        is=ie+1                                                         11d21s22
        ie=is+6                                                         11d21s22
        if(ie.gt.idsln)then                                             11d21s22
         notok=1                                                        11d21s22
         go to 1524                                                     11d21s22
        end if                                                          11d21s22
        sline(is:ie)=termsym                                            11d21s22
        is=ie+1                                                         11d21s22
        do i=1,nroot                                                    11d21s22
         ie=is+18                                                       11d21s22
         if(ie.gt.idsln)then                                             11d21s22
          notok=1                                                        11d21s22
          go to 1524                                                     11d21s22
         end if                                                          11d21s22
         sline(is:ie)=ostng(i)                                          11d21s22
         is=ie+1                                                        11d21s22
        end do                                                          11d21s22
 1524   continue                                                        11d21s22
        if(notok.ne.0)then                                              11d21s22
         write(42,*)(bc(iextrad+i),i=0,nextrad-1),ndoub,spinm,termsym,   1d10s19
     $       (bc(ieigint+i)+shift,i=0,nroot-1)
        else                                                            11d21s22
         write(42,1525)(sline(i:i),i=1,ie)                              11d21s22
 1525    format(100000a1)                                               11d21s22
        end if                                                          11d21s22
        close(unit=42)                                                  1d10s19
        write(6,*)('cpu times: ')
        xproc=1d0/dfloat(mynprocg)                                      3d23s21
        do i=1,9                                                        3d23s21
         bc(jtimer+i)=bc(jtimer+i)*xproc                                3d23s21
         bc(ktimer+i)=bc(ktimer+i)*xproc                                3d23s21
         bc(ktimer+i)=sqrt(abs(bc(ktimer+i)-bc(jtimer+i)**2))           3d23s21
        end do                                                          3d23s21
        write(6,*)('intci: '),bc(jtimer+7),bc(ktimer+7)                 3d23s21
        write(6,*)('hcsi: '),bc(jtimer+1),bc(ktimer+1)                  3d23s21
        write(6,*)('hcis: '),bc(jtimer+8),bc(ktimer+8)                  3d23s21
        write(6,*)('hss: '),bc(jtimer+2),bc(ktimer+2)                   3d23s21
        write(6,*)('hcdi: '),bc(jtimer+4),bc(ktimer+4)                  3d23s21
        write(6,*)('hcid: '),bc(jtimer+9),bc(ktimer+9)                  3d23s21
        write(6,*)('hcds: '),bc(jtimer+3),bc(ktimer+3)                  3d23s21
        write(6,*)('hcddjk: '),bc(jtimer+5),bc(ktimer+5)                3d23s21
        write(6,*)('hcdd4v: '),bc(jtimer+6),bc(ktimer+6)                3d23s21
        if(mynowprog.eq.0.and.nwavef.ne.0)then                          3d23s21
         write(2)iunc,ethred                                            7d8s21
         write(2)nctf(isymmrci),ntot                                    4d8s20
         write(2)(bc(ivint+i),i=0,nctf(isymmrci)*nroot-1)               4d8s20
         write(2)(bc(ieigint+i)+shift,i=0,nroot-1)                      2d9s21
         ipack2(1)=nlzzci                                               10d7s24
         ipack2(2)=nder                                                 10d7s24
         write(2)ipack4(1),lambdaci                                     10d7s24
         if(nsing.ne.0)then
          call dumpsing(ibc(ihsdiag),ibc(nff1),ibc(iff1),               7d8s21
     $        ibc(icsf),mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,     7d9s21
     $         bc(ivsold),0,bc,ibc)                                     11d14s22
         else                                                           7d8s21
          izero=0                                                       7d8s21
          write(2)izero                                                 7d8s21
         end if                                                         7d8s21
         if(ndoub.gt.0)then                                             7d21s21
          if(iunc(2).eq.0)then                                          7d21s21
           call dumpdoubc(ibc(nff22),nfdat,bc(ivdinout),nsymb,ndoub,    8d3s21
     $         mdoub,nroot,mdon,mdoo,bc(icsf),bc(icsf2),0,bc,ibc,nder)  9d18s24
          else                                                          7d21s21
           call dumpdoubuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ibc(icsf),  7d21s21
     $        ibc(icsf2),mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,0,4,11d14s22
     $          bc,ibc)                                                 11d14s22
          end if                                                        7d21s21
         else                                                           7d21s21
          izero=0                                                       7d8s21
          write(2)izero                                                 7d8s21
         end if                                                         7d21s21
         if(nsing.ne.0)then                                             7d27s21
          call dumpsing(ibc(ihsdiag),ibc(nff1),ibc(iff1),               7d8s21
     $        ibc(icsf),mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,     7d9s21
     $         bc(ivsold),1,bc,ibc)                                     11d14s22
         end if                                                         7d27s21
         if(ndoub.gt.0)then                                             7d21s21
          if(iunc(2).eq.0)then                                          7d21s21
          else                                                          7d21s21
           call dumpdoubuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ibc(icsf),  7d21s21
     $        ibc(icsf2),mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,1,4,11d14s22
     $          bc,ibc)                                                 11d14s22
          end if                                                        7d21s21
         end if                                                         7d21s21
        end if                                                          11d24s19
        end if
        if(ndoub.ne.0)then                                              3d12s21
         if(iunc(2).eq.0)then                                           6d9s21
          call tofro(bc(ivdinout),bc(ivdnono),nroot,nfdat,nvirt,nsymb,     12d21s20
     $       multh,isymmrci,1,ndoub,mdoub,bc,ibc)                       11d10s22
         end if                                                         6d9s21
        end if
        if(iunc(2).eq.0)then                                            7d19s24
         if(nder.eq.0)then                                              7d19s24
          call denmx(mdon,mdoo,ibc(ibasisf(isymmrci)),ibc(iptrfbit),      3d3s21
     $       ibc(icsf),nfcnf(isymmrci),bc(ivint),nctf(isymmrci),        3d3s21
     $       isymmrci,irel,ism,irefo,nvirt,nsymb,multh,ixw1,ixw2,nh0av, 3d3s21
     $       nroot,ibc(ihsdiag),ibc(nff0),ibc(iff0),ibc(nff1),ibc(iff1),3d4s21
     $       nsing,ndoub,maxbx,norb,ibc(nff22),nfdat,bc(ivdnono),nec,   3d9s21
     $       ibc(icsf2),ih0copy,nbasdws,nbasisp,ncomp,lprint,sr2,srh,   1d24s22
     $       norbci,iorb,idoubo,bc,ibc)                                 11d10s22
         else                                                            5d5s23
          if(ndoub.eq.0)then                                             5d8s23
           nff22=1                                                       5d8s23
           ivdnono=1                                                     5d8s23
          end if                                                         5d8s23
          if(nsing.eq.0)then                                             5d8s23
           nff1=1                                                        5d8s23
           iff1=1                                                        5d8s23
          end if                                                         5d8s23
          call denmx12(mdon,mdoo,ibc(ibasisf(isymmrci)),ibc(iptrfbit),   5d8s23
     $      ibc(icsf),nfcnf(isymmrci),bc(ivint),nctf(isymmrci),         5d8s23
     $      isymmrci,irel,ism,irefo,multh,ixw1,ixw2,nh0av,              7d29s22
     $      nroot,ibc(ihsdiag),ibc(nff0),ibc(iff0),ibc(nff1),ibc(iff1),
     $      nsing,ndoub,maxbx,norb,ibc(nff22),nfdat,bc(ivdnono),        5d8s23
     $      nec,ibc(icsf2),ih0av,nbasisp,ncomp,lprint,sr2,srh,          5d8s23
     $      ioooo,shift,bc,ibc,ionexbt,jmats,kmats,idoubo,iff,idarot,   7d11s23
     $        nder,natom,ngaus,ibdat,ihmat,iorb,noc,ipair,isym,iapair,   7d17s23
     $       ibstor,isstor,iptoh,idorel,ascale,i3x,i4x,i4xb,            7d17s23
     $        isend1,isend2,isend3,isend4,ih0copy,iunc,nff2,            9d5s23
     $       bc(ivdinout),mdoub,bc(issdig),iden1e,idenhvv,idenj,idenk,  9d28s23
     $       ibc(idenput),idenorg,ncd,nbasdwsc,tovr,ibufvcvds,idvmt,    5d2s24
     $       potdws,idervcv,ipropsym)                                   10d30s24
         end if                                                          6d9s21
        end if                                                          7d19s24
        if(mynowprog.eq.0.and.nwavef.ne.0)then                          2d13s23
         if(nlzzci.ne.0)then                                            5d12s21
          write(6,*)('now go on and dumpbas for equivalent symmetries')
          if(lambdaci.ne.0)then                                         5d12s21
           if(nlzzci.eq.2)then                                          5d12s21
            if(isymmrci.eq.1)then                                        5d12s21
             isymo=4                                                     5d12s21
            else if(isymmrci.eq.2)then                                   5d12s21
             isymo=3                                                     5d12s21
            else if(isymmrci.eq.3)then                                   5d12s21
             isymo=2                                                     5d12s21
            else if(isymmrci.eq.4)then                                   5d12s21
             isymo=1                                                     5d12s21
            else if(isymmrci.eq.5)then                                   5d12s21
             isymo=8                                                     5d12s21
            else if(isymmrci.eq.6)then                                   5d12s21
             isymo=7                                                     5d12s21
            else if(isymmrci.eq.7)then                                   5d12s21
             isymo=6                                                     5d12s21
            else                                                         5d12s21
             isymo=5                                                     5d12s21
            end if                                                       5d12s21
            write(2)isymo                                                5d12s21
            call genfcnp(nec,mdon,mdoo,multh,ibc(nff1),ibc(iff1),       7d8s21
     $           ibc(nff0),iff0,nsymb,isymo,nctff,ibc(icsf),ntotf,bc,   11d11s22
     $           ibc)                                                   11d11s22
            call dumpbas2(nec,mdon,mdoo,ibc(nff0),ibc(iff0),            7d19s21
     $           ibc(icsf))                                             7d19s21
            nwavrec=nwavrec+1
           else                                                         5d12s21
            nl=2*lambdaci                                               5d12s21
            isymo=ibcoff                                                5d12s21
            ibcoff=isymo+nl                                             5d12s21
            call enough('mrci. 29',bc,ibc)
            if(isymmrci.eq.1.or.isymmrci.eq.4)then                      5d14s21
             iodd1=6                                                    5d14s21
             iodd2=7                                                    5d14s21
             ieven1=1                                                   5d14s21
             ieven2=4                                                   5d14s21
            else                                                        5d14s21
             iodd1=2                                                    5d14s21
             iodd2=3                                                    5d14s21
             ieven1=5                                                   5d14s21
             ieven2=8                                                   5d14s21
            end if                                                      5d14s21
            jsymo=isymo                                                 5d14s21
            do ll=1,lambdaci                                            5d14s21
             if(mod(ll,2).eq.0)then                                     5d14s21
              ibc(jsymo)=ieven1                                         5d14s21
              jsymo=jsymo+1                                             5d14s21
              ibc(jsymo)=ieven2                                         5d14s21
              jsymo=jsymo+1                                             5d14s21
             else                                                       5d14s21
              ibc(jsymo)=iodd1                                          5d14s21
              jsymo=jsymo+1                                             5d14s21
              ibc(jsymo)=iodd2                                          5d14s21
              jsymo=jsymo+1                                             5d14s21
             end if                                                     5d14s21
            end do                                                      5d14s21
            do is=0,nl-1                                                5d12s21
             mysym=ibc(isymo+is)                                        5d12s21
             write(2)mysym                                              5d12s21
            call genfcnp(nec,mdon,mdoo,multh,ibc(nff1),ibc(iff1),       8d24s21
     $           ibc(nff0),iff0,nsymb,mysym,nctff,ibc(icsf),ntotf,bc,   11d11s22
     $            ibc)                                                  11d11s22
            call dumpbas2(nec,mdon,mdoo,ibc(nff0),ibc(iff0),            8d24s21
     $           ibc(icsf))                                             7d19s21
             nwavrec=nwavrec+1
            end do                                                      5d12s21
           end if                                                       5d12s21
          end if                                                        5d12s21
         end if                                                         5d12s21
         close(unit=2)                                                   11d22s19
        end if                                                          2d13s23
        call dws_synca                                                  8d19s21
        ibcoff=ibcoffo                                                  11d23s19
        return                                                          11d22s19
       end if                                                           11d13s18
       if(nrestart.gt.0.and.mod(icit,max(1,nrestart)).eq.0.and.         8d19s22
     $      mynowprog.eq.0.and.nwavef.ne.0)then                         8d19s22
        write(6,*)('saving restart data ')
        call savewf(iwopen,nwavrec,iunc,ethred,nctf(isymmrci),nroot,    8d16s22
     $      ntot,bc(ivint),bc(ieigint),shift,nlzzci,lambdaci,nsing,     8d19s22
     $      ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(icsf),mdon,mdoo,nsymb, 8d16s22
     $      multh,nvirt,isymmrci,bc(ivsold),ndoub,nff22,nfdat,ivdinout, 8d16s22
     $      mdoub,bc(icsf2),ihddiag,nff2,iff2,bc,ibc,nder)              10d7s24
       end if                                                           8d16s22
       if(ldebug)write(6,*)('calling updatex '),ibcoff                  7d15s20
       if(iunc(2).eq.0)then                                             6d9s21
        call updatexc(bc(ivdinout),bc(ihdi),bc(iddoub),                  1d22s21
     $       bc(ivdold),ndoub,ibc(ihsdiag),bc(ivsold),                  1d24s21
     $      bc(issdig),ibc(nff1),ibc(icsf),nsymb,mdon,mdoo,nvirt,       11d10s20
     $      isymmrci,bc(ieigint),nroot,isto,nfdat,bc(ivstry),nsing,     1d25s21
     $      ldebugf,ngots,ngotd,ibc(iamconverged),bc(ipairs),ldavid,bc, 11d10s22
     $      ibc)                                                        11d10s22
       else                                                             6d9s21
        call updatexuc(bc(idddig),bc(ivdold),ibc(ihddiag),              6d9s21
     $       ibc(ihsdiag),bc(ivsold),                                   6d9s21
     $      bc(issdig),ibc(nff1),ibc(nff2),ibc(icsf),ibc(icsf2),nsymb,  6d9s21
     $       mdon,mdoo,nvirt,                                           6d9s21
     $      isymmrci,bc(ieigint),nroot,isto,bc(ivstry),nsing,           6d9s21
     $      ldebugf,ngots,ngotd,ibc(iamconverged),nlocald,shift0,       8d19s21
     $       bc(ipairs),bc,ibc)                                         11d10s22
       end if                                                           6d9s21
       if(nrestart.lt.0)then                                            8d16s22
        call recoverwf(iunc,nroot,ntot,nlzzci,lambdaci,nsing,           8d16s22
     $      ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(icsf),mdon,mdoo,nsymb, 8d16s22
     $      multh,nvirt,isymmrci,bc(ivsold),bc(ivstry),ndoub,nff22,     8d19s22
     $      nfdat,ivdinout,mdoub,bc(icsf2),ibc(ihddiag),ibc(nff2),iff2, 8d19s22
     $      norb,bc(issdig),bc(idddig),bc,ibc)                          11d10s22
        iwopen=0                                                        8d19s22
        nrestart=-nrestart                                              8d16s22
       end if                                                           8d16s22
       nok2conv=0                                                       4d15s21
       nrootu=isto                                                      9d14s18
       if(icit.gt.2.and.echange.lt.1d2*econvci.and.ndoub.ne.ndoubs)then 4d14s21
        if(lprint)write(6,*)                                            4d14s21
     $      ('singles have converged pretty well: now turn on doubles') 4d14s21
        nok2conv=1                                                      4d15s21
        ngots=0                                                         4d14s21
        ngotd=0                                                         4d14s21
        ndoub=ndoubs                                                     4d14s21
        if(iunc(2).eq.0)then                                            6d9s21
         do i=0,ndoub-1                                                  4d14s21
          bc(ivdinout+i)=0d0                                             4d14s21
         end do                                                          4d14s21
        else                                                            6d9s21
         write(6,*)('we probably need a reloadvd here ')
         call dws_synca
         call dws_finalize
         stop
        end if
        call reloadvs(bc(ivsold),bc(ivstry),ibc(ihsdiag),nroot,         4d14s21
     $       ibc(nff1),ibc(icsf),nvirt,mdon,mdoo,nsymb,multh,isymmrci,  4d14s21
     $       bc(issdig),bc,ibc)                                         11d9s22
       end if                                                            4d14s21
       if(ndoub.gt.0.and.iunc(2).eq.0)then                              6d9s21
        if(ldebug)write(6,*)('calling tofro '),ibcoff                   11d17s21
        call tofro(bc(ivdinout),bc(ivdnono),nrootu,nfdat,nvirt,nsymb,   3d30s21
     $       multh,isymmrci,1,ndoub,mdoub,bc,ibc)                       11d10s22
         do i=0,nrootu*mdoub-1                                          3d30s21
          bc(igd+i)=0d0                                                  12d23s20
         end do                                                          12d23s20
        end if                                                          2d17s21
       if(ldebug)write(6,*)('back'),ibcoff                              7d15s20
       if(icit.le.2.and.itestmrci.eq.3.and.lprint)then                  5d15s19
        write(7,*)ntot
        do i=0,ntot-1                                                   5d15s19
         write(7,*)bc(ivnew+i)                                          5d15s19
        end do                                                          5d15s19
       end if                                                           5d15s19
      end if                                                            9d13s18
      call dws_synca                                                    8d11s22
      call ddi_iaccword0                                                8d11s22
      nwiacc=0                                                          8d11s22
      call ddi_zero(bc,ibc,ihi)                                         11d15s22
      if(ndoub.gt.0)then
       if(iunc(2).eq.0)then                                             7d7s21
        call ddi_zero(bc,ibc,igxx)                                      11d15s22
       end if                                                           7d7s21
      end if                                                            7d7s21
      if(nsing.gt.0)then                                                9d14s18
       call dws_synca                                                   1d25s21
       call second(time1)                                               12d14s18
       if(ldebug)write(6,*)('calling hcss'),ibcoff,icit                    7d15s20
       call hcss(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(icsf),nec,mdon,   12d12s20
     $      mdoo,nsymb,multh,ixw1,ixw2,ih0av,nh0av,ioooo,jmats,kmats,   12d14s20
     $     nvirt,nrootu,ism,irel,irefo,isymmrci,norb,ldebugf,maxbx,nct1,5d5s21
     $     tdenss,tovr,nwiacc,bc,ibc)                                   1d27s23
       if(ldebug)write(6,*)('back'),ibcoff                              7d15s20
       call second(time2)                                               12d14s18
       telap=time2-time1-tovr                                           12d14s18
       timex(2)=timex(2)+telap                                          12d14s18
       call second(time1)
       if(ldebug)write(6,*)('calling hcsicsf'),ibcoff,
     $      icit
       call hcis(bc(ivstry),ibc(nff1),ibc(iff1),ibc(nff0),              12d11s20
     $      ibc(iff0),ihi,nctf(isymmrci),ibc(icsf),nec,mdon,mdoo,       12d12s20
     $      nsymb,multh,ixw1,ixw2,ih0av,nh0av,ionexb,nvirt,nrootu,      1d25s21
     $      ldebugf,ionexbt,bc,ibc)                                     11d10s22
       call dws_synca
       i18=1                                                            1d8s21
       in8=nrootu*nctf(isymmrci)                                        1d8s21
       call ddi_get(bc,ibc,ihi,i18,i18,i18,in8,bc(ihis))                11d15s22
       if(ldebug)write(6,*)('back'),ibcoff                              7d15s20
       call second(time2)
       telap=time2-time1-tovr                                           12d14s18
       timex(8)=timex(8)+telap                                          3d9s21
       call getgss(igss,nsymb,isymmrci,mdon,mdoo,ibc(nff1),nvirt,multh, 1d27s23
     $      ibc(icsf),nrootu,ibc(ihsdiag),bc,ibc)                       1d27s23
       if(ndoub.gt.0)then                                               10d2s18
        call second(time1)                                              12d14s18
        if(ldebug)write(6,*)('calling hcdscsf'),ibcoff,
     $       icit
        if(iunc(2).eq.0)then                                            6d9s21
         call hcds(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(nff22),nfdat,    12d18s20
     $       bc(igd),bc(ivdnono),nsymb,mdon,mdoo,nec,multh,isymmrci,    12d21s20
     $       nvirt,ibc(icsf),ibc(icsf2),irel,ism,irefo,ismult,ixw1,ixw2,12d18s20
     $       norb,nrootu,ih0av,nh0av,ionexb,i3x,i3x3,bc(ivdinout),mdoub,1d6s21
     $       ndoub,ldebugf,maxbx,tdends,tovr,ionexbt,sr2,nwiacc,bc,ibc) 11d10s22
        else                                                            6d9s21
         call hcdsuc(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(ihddiag),     6d9s21
     $        ibc(nff2),ibc(iff2),nsymb,mdon,mdoo,nec,multh,isymmrci,   6d9s21
     $       nvirt,ibc(icsf),ibc(icsf2),irel,ism,irefo,ismult,ixw1,ixw2,12d18s20
     $       norb,nrootu,ih0av,nh0av,ionexb,i3x,i3x3,                   6d9s21
     $       nlocald,ldebugf,maxbx,maxbxd,tdends,tovr,ionexbt,sr2,      8d11s22
     $        nwiacc,bc,ibc)                                            11d10s22
        end if                                                          6d9s21
        if(ldebug)write(6,*)('back'),ibcoff
        call second(time2)                                              12d14s18
        telap=time2-time1-tovr                                          12d14s18
        timex(3)=timex(3)+telap                                         12d14s18
       end if                                                           10d2s18
      end if                                                            9d14s18
      if(ndoub.gt.0)then                                                10d31s18
       call second(time1)
       if(ldebug)write(6,*)('calling hcid'),ibcoff,
     $      icit
       if(iunc(2).eq.0)then                                             6d17s21
        call hcid(ibc(nff22),nsymb,mdon,mdoo,multh,isymmrci,nvirt,       1d7s21
     $      ibc(icsf),nec,ibc(icsf2),nfdat,irel,ism,kmats,irefo,ismult, 12d8s20
     $      bc(ivdnono),ixw1,ixw2,norb,ibc(nff0),ibc(iff0),ihi,         1d7s20
     $      nctf(isymmrci),nrootu,mdoub,ldebug,tovr,tpart,sr2,srh,bc,   11d10s22
     $      ibc)                                                        11d10s22
       else                                                             6d17s21
        if(ldebug)write(6,*)('calling hciduc')
        call hciduc(ibc(ihddiag),ibc(nff2),ibc(iff2),ibc(nff0),         6d19s21
     $      ibc(iff0),ihi,nctf(isymmrci),ibc(icsf),ibc(icsf2),nec,mdon, 6d19s21
     $       mdoo,nsymb,multh,ixw1,ixw2,kmats,nvirt,ldebugf,maxbxd,     6d19s21
     $      nrootu,srh,bc,ibc)                                          11d10s22
       end if                                                           6d17s21
       if(ldebug)write(6,*)('back'),ibcoff,ixw2,ibc(32800)                              7d15s20
       call second(time2)                                               12d14s18
       telap=time2-time1-tovr                                           12d14s18
       timex(9)=timex(9)+telap                                          3d9s21
       call second(time1)
       if(iunc(2).eq.0)then                                             6d21s21
        if(ldebug)write(6,*)('calling hcddcsfjk'),ibcoff,
     $      icit
        call hcddjk(ibc(nff22),nfdat,bc(igd),bc(ivdnono),nsymb,mdon,    6d21s21
     $       mdoo,                                                      6d21s21
     $      nec,multh,isymmrci,nvirt,ibc(icsf),ibc(icsf2),irel,ism,     1d8s21
     $      irefo,ismult,ixw1,ixw2,norb,nrootu,ih0av,nh0av,ioooo,       1d8s21
     $      jmats,kmats,ndoub,mdoub,shift0,tdendd,tovr,sr2,srh,iden1e,  11d7s22
     $       idenhvv,idenj,idenk,bc,ibc,ibc(idenput),idenorg)           5d17s23
        i18=1                                                            1d8s21
        in8=nrootu*mdoub                                                 1d18s21
        call ddi_acc(bc,ibc,igxx,i18,i18,i18,in8,bc(igd))               11d15s22
        call dws_synca                                                   1d18s21
        call ddi_get(bc,ibc,igxx,i18,i18,i18,in8,bc(igd))               11d15s22
        call tofro(bc(ihdi),bc(igd),nrootu,nfdat,nvirt,nsymb,multh,      1d18s21
     $      isymmrci,2,ndoub,mdoub,bc,ibc)                              11d10s22
        if(ldebugf)then                                                  5d5s21
         write(6,*)('ihdi from tofro ')
         call prtit(bc(ihdi),nrootu,nfdat,nvirt,nsymb,multh,isymmrci,     1d18s21
     $      ndoub,bc,ibc)                                               11d10s22
        end if                                                           1d25s21
        call reordergv(bc(ihdi),nrootu,nfdat,nvirt,nsymb,multh,isymmrci, 1d20s21
     $      ndoub,1,bc,ibc)                                             11d9s22
        ivtmp=ibcoff                                                     1d20s21
        ibcoff=ivtmp+ndoub*nrootu                                        1d20s21
        call enough('mrci. 31',bc,ibc)
        do i=0,ndoub*nrootu-1                                            1d20s21
         bc(ivtmp+i)=bc(ivdinout+i)                                      1d20s21
        end do                                                           1d20s21
        call reordergv(bc(ivtmp),nrootu,nfdat,nvirt,nsymb,multh,        6d21s21
     $       isymmrci,ndoub,1,bc,ibc)                                   11d9s22
        if(ldebug)write(6,*)('back'),ibcoff                              7d15s20
       else                                                             6d21s21
        if(ldebug)write(6,*)('calling hcdduc')
        call hcdduc(ibc(ihddiag),ibc(nff2),ibc(iff2),bc(icsf),bc(icsf2),6d22s21
     $      nec,mdon,mdoo,nsymb,multh,ixw1,ixw2,ih0av,nh0av,ioooo,jmats,6d22s21
     $       kmats,i4x,i4xb,nvirt,nrootu,ism,irel,irefo,isymmrci,norb,  7d6s21
     $       ldebugf,maxbxd,sr2,srh,timex(6),tovr,nwiacc,bc,ibc)        11d10s22
       end if                                                           6d21s21
       call second(time2)                                               12d14s18
       telap=time2-time1-tovr                                           12d14s18
       timex(5)=timex(5)+telap                                          12d14s18
       if(iunc(2).eq.0)then                                             3d24s23
        call second(time1)                                               8d19s22
        if(ldebug)write(6,*)('calling hcdd4v'),ibcoff,icit                   7d15s20
        if(iskip4v.eq.0)then                                            3d24s23
         call hcdd4v(bc(ivtmp),bc(ihdi),ncd,multh,iorb,nbasdwsc,natom,    1d20s21
     $      ngaus,ibdat,ipair,isym,iapair,ibstor,isstor,idorel,         11d2s20
     $      ascale,ndoub,nrootu,nbasisp,ndoub,sr2,srh,bc,ibc)           11d9s22
        end if                                                          3d24s23
        ibcoff=ivtmp                                                     1d20s21
        if(ldebug)write(6,*)('back'),ibcoff                              1d18s21
        call reordergv(bc(ihdi),nrootu,nfdat,nvirt,nsymb,multh,isymmrci, 1d20s21
     $      ndoub,-1,bc,ibc)                                            11d9s22
        if(ldebugf)then                                                  5d5s21
         write(6,*)('after hcdd4v')
         call prtit(bc(ihdi),nrootu,nfdat,nvirt,nsymb,multh,isymmrci,     1d18s21
     $      ndoub,bc,ibc)                                               11d10s22
        end if                                                           1d25s21
        call second(time2)                                               12d14s18
        telap=time2-time1-tovr                                           12d14s18
        timex(6)=timex(6)+telap                                          12d14s18
       end if                                                           6d21s21
      end if                                                            12d13s18
      ngot=max(ngots,ngotd)                                             4d14s21
      nvdim=ngot+nrootu                                                 1d21s21
      ihxy=ibcoff                                                       11d9s18
      ioxy=ihxy+nvdim*nvdim                                             1d21s21
      ihix=ioxy+nvdim*nvdim                                             1d21s21
      ibcoff=ihix+nctf(isymmrci)*nvdim                                  1d21s21
      call enough('mrci. 32',bc,ibc)
      do i=ihxy,ibcoff-1                                                1d21s21
       bc(i)=0d0                                                        1d21s21
      end do                                                            1d21s21
      call ddi_iaccword2(nwiacc)                                        8d11s22
      if(nsing.gt.0)then                                                1d21s21
       call hsingpart(bc(ihxy),bc(ioxy),nvdim,nrootu,ibc(ihsdiag),      1d21s21
     $     ibc(nff1),ibc(icsf),nvirt,mdon,mdoo,nsymb,multh,isymmrci,    1d23s21
     $     bc(ivsold),bc(igsold),ngots,bc,ibc)                          11d10s22
       if(ndoub.eq.0.or.iunc(2).eq.0)then                               7d7s21
        nwds=nvdim*nvdim*2                                               1d24s21
        call dws_gsumf(bc(ihxy),nwds)                                    1d24s21
       end if                                                           7d7s21
      end if                                                            1d21s21
      if(ndoub.gt.0)then                                                1d21s21
       if(iunc(2).eq.0)then                                             7d7s21
        call hdoubpart(bc(ihxy),bc(ioxy),nvdim,nrootu,bc(ivdinout),      1d21s21
     $     bc(ihdi),nfdat,nvirt,nsymb,multh,isymmrci,bc(ivdold),        1d23s21
     $     bc(igdold),ngotd,bc(iepart),nroote,bc,ibc)                   1d12s23
       else                                                             7d7s21
        call hdoubpartuc(bc(ihxy),bc(ioxy),nvdim,nrootu,ibc(ihddiag),      1d21s21
     $     ibc(nff2),ibc(icsf),ibc(icsf2),nvirt,mdon,mdoo,nsymb,multh,
     $       isymmrci,bc(ivdold),bc(igdold),ngotd,bc(iepart),nroote,bc, 1d12s23
     $       ibc)                                                       11d10s22
        nwds=nvdim*nvdim*2                                               1d24s21
        call dws_gsumf(bc(ihxy),nwds)                                    1d24s21
       end if                                                           7d7s21
      end if                                                            1d21s21
      if(ldebug)then                                                    5d5s21
       write(6,*)('what we have for hxy: ')
       call prntm2(bc(ihxy),nvdim,nvdim,nvdim)
       write(6,*)('what we have for oxy: ')
       call prntm2(bc(ioxy),nvdim,nvdim,nvdim)
      end if                                                            1d25s21
      rms=0d0                                                           1d25s21
      rmsh=0d0                                                          2d11s21
      do i=0,nvdim-2                                                    1d21s21
       do j=i+1,nvdim-1                                                 1d21s21
        ji=ihxy+j+nvdim*i                                               1d21s21
        ij=ihxy+i+nvdim*j                                               1d21s21
        rmsh=rmsh+(bc(ji)-bc(ij))**2
        bc(ji)=bc(ij)                                                   1d24s21
        ji=ioxy+j+nvdim*i                                               1d21s21
        ij=ioxy+i+nvdim*j                                               1d21s21
        rms=rms+bc(ji)**2+bc(ij)**2                                     1d25s21
       end do                                                           1d13s20
       ii=ioxy+i*(nvdim+1)                                              1d25s21
       rms=rms+(bc(ii)-1d0)**2                                          1d25s21
      end do                                                            1d13s20
      ii=ioxy+(nvdim-1)*(nvdim+1)                                       1d25s21
      rms=rms+(bc(ii)-1d0)**2                                           1d25s21
      rms=sqrt(rms/dfloat(nvdim*nvdim))                                 1d25s21
      rmsh=sqrt(rmsh/dfloat(nvdim*nvdim))                                 1d25s21
      oxy=bc(ioxy)                                                      1d25s21
      bc(ioxy)=rms                                                      1d25s21
      if(ldebugf)then                                                   5d5s21
       write(6,*)('symmetrized hxy: '),rmsh
       call prntm2(bc(ihxy),nvdim,nvdim,nvdim)
      end if                                                            1d25s21
      nwds=nvdim*nvdim+1                                                1d25s21
      call dws_bcast(bc(ihxy),nwds)                                     1d25s21
      if(bc(ioxy).gt.1d-10)then                                         1d25s21
       write(6,*)('external vectors have lost orthogonality!!!')
       write(6,*)('rms overlap test: '),bc(ioxy)
       bc(ioxy)=oxy                                                     1d25s21
       write(6,*)('overlap matrix: ')
       call prntm2(bc(ioxy),nvdim,nvdim,nvdim)                          1d25s21
       stop                                                             1d25s21
      end if                                                            1d25s21
      i18=1                                                             1d21s21
      in8=nctf(isymmrci)*nrootu                                         1d21s21
      call ddi_get(bc,ibc,ihi,i18,i18,i18,in8,bc(ihix))                 11d15s22
      if(ldebug )then                                                   5d5s21
       write(6,*)('what we have for hix')
       call prntm2(bc(ihix),nctf(isymmrci),nrootu,nctf(isymmrci))
      end if                                                            1d25s21
      do i=0,ngot-1                                                     1d23s21
       ip=i+nrootu                                                      1d23s21
       iad=ihix+nctf(isymmrci)*ip                                       1d23s21
       jad=ihixold+nctf(isymmrci)*i                                     1d23s21
       do j=0,nctf(isymmrci)-1                                          1d23s21
        bc(iad+j)=bc(jad+j)                                             1d23s21
       end do                                                           1d23s21
      end do                                                            1d23s21
      if(ngot.gt.0.and.ldebug )then                                     5d5s21
       write(6,*)('with saved part ')
       call prntm2(bc(ihix),nctf(isymmrci),nvdim,nctf(isymmrci))
      end if                                                            1d23s21
      if(ldavid)then                                                    12d27s21
       if(nvdim.le.maxdsr)then                                          12d27s21
        do i=0,nvdim*nctf(isymmrci)-1                                   12d27s21
         bc(ihixold+i)=bc(ihix+i)                                       12d27s21
        end do                                                          12d27s21
       end if                                                           12d27s21
      end if                                                            12d27s21
      if(itestmrci.ne.0.and.bc(132).ne.-132d0)then
       call dws_sync
       call dws_finalize
       stop
      end if
      call second(time1)                                                5d22s19
      if(ldebug)write(6,*)('calling intcsf'),ibcoff,icit                     7d15s20
      ibcb4i=ibcoff                                                     10d7s21
      call intcsf(ibc(idorbf),ibc(isorbf),mdon,mdoo,
     $     ibc(ibasisf(isymmrci)),                                      4d8s20
     $     ibc(jptrf),ibc(icsf),nfcnf(isymmrci),ih0a,ioooo,shift,pthrs0,4d8s20
     $     nec,ibc(icsfpd),multh,1,nvdim,bc(ihix),bc(ihxy),             1d25s21
     $    nroot,0d0,2.0d0,lprint,ieigint,ivint,maxpsci,nctf(isymmrci),  10d7s21
     $     nvcul,                                                       4d8s20
     $   nlzzci,islz,ixlzz,1,bintvec,intvec,nfcnfi,ibc(iptrfi),ibasisfi,10d14s20
     $    ih0ae,ibc(iptrfbit),ixw1,ixw2,nbasisfi,ibc(nfcnpxfi),nsymb,   1d12s23
     $    ibc(iptrfibit),ibc(idorbfi),ibc(isorbfi),nfcnc,               3d24s21
     $      ibc(ibasisc),ibc(iptrcb),bc(ividotvi),bc(ivihvi),nameci,0,  11d9s22
     $     bc,ibc,isaved,ostng,idroot)                                  7d13s23
      do i=0,nroot-1                                                    10d7s21
       bc(ieigintxyz+i)=bc(ieigint+i)                                   10d7s21
      end do                                                            10d7s21
      ieigint=ieigintxyz                                                10d7s21
      do i=0,nroot*(nctf(isymmrci)+nvdim)-1                             12d29s21
       bc(ivints+i)=bc(ivint+i)                                         10d7s21
      end do                                                            10d7s21
      ivint=ivints                                                      10d7s21
      if(wconvci.ne.0d0)then                                            6d20s24
       if(icit.ne.1.and.ngot.gt.0)then                                  7d12s24
        diffr=0d0                                                        6d20s24
        do i=0,nroot-1                                                   6d20s24
         jvintlast=ivintlast+nctf(isymmrci)*i                            6d20s24
         jvints=ivints+(nctf(isymmrci)+nvdim)*i                          6d20s24
         do j=0,nctf(isymmrci)-1                                         6d20s24
          diffr=diffr+(bc(jvintlast+j)-bc(jvints+j))**2                  6d20s24
         end do                                                          6d20s24
        end do                                                           6d20s24
        diffr=sqrt(diffr/dfloat(nroot*nctf(isymmrci)))                   6d20s24
       end if                                                            6d20s24
       do i=0,nroot-1                                                    6d20s24
        jvintlast=ivintlast+nctf(isymmrci)*i                             6d20s24
        jvints=ivints+(nctf(isymmrci)+nvdim)*i                           6d20s24
        do j=0,nctf(isymmrci)-1                                          6d20s24
         bc(jvintlast+j)=bc(jvints+j)                                    6d20s24
        end do                                                           6d20s24
       end do                                                            6d20s24
      end if                                                            6d20s24
      ibcoff=ibcb4i                                                     10d7s21
      if(isaved.ne.0)then                                               2d3s23
       bintvec=.false.                                                   1d17s20
      else                                                              2d3s23
       bintvec=.true.                                                   2d3s23
      end if                                                            2d3s23
      if(ldebug)write(6,*)('back'),ibcoff                               7d15s20
      call second(time2)                                                5d22s19
      telap=time2-time1-tovr                                            5d22s19
      timex(7)=timex(7)+telap                                           5d22s19
      do i=0,nctf(isymmrci)-1                                           4d8s20
       do j=0,nroot-1                                                   4d3s19
        ji=ivint+j+nroot*i                                              4d3s19
        ij=ivintt+i+nctf(isymmrci)*j                                    4d8s20
        bc(ij)=bc(ji)                                                   4d3s19
       end do                                                           4d3s19
      end do                                                            4d3s19
      if(icit.ne.1.and.lprint)then                                      2d22s19
       changem=1d10                                                     3d12s21
       if(nroot.le.idroot)then                                          2d26s19
        if(ngot.ne.0)then                                               7d17s24
         do i=1,nroot
          im=i-1
          write(ostng(i),1551)bc(ieigint+im)+shift                       2d9s21
 1551     format(f19.12)
          change=bc(ieigint+im)-bc(ieold+im)                             2d26s19
          if(change.ne.0d0)then                                            2d26s19
           changem=min(changem,abs(change))                              3d12s21
           ndig=nint(-log10(abs(change)))+1                              2d26s19
           ist=min(19,7+ndig)                                            1d12s23
           if(ostng(i)(ist:ist).ne.' ')then                              1d6s20
            if(ichar(ostng(i)(ist:ist)).lt.ichar('0').or.                1d21s20
     $        ichar(ostng(i)(ist:ist)).gt.ichar('9'))then               1d21s20
             if(ist.lt.20)then                                           10d20s20
              write(6,*)('bad last digit: '),ostng(i)(ist:ist)
              write(6,*)('ndig,ist '),ndig,ist
             end if                                                      10d20s20
             left=0                                                      1d21s20
            else                                                         1d21s20
             read(ostng(i)(ist:ist),*)left                                 11d17s19
            end if                                                       1d21s20
            if(left.gt.5)then                                             11d17s19
             iadd=1                                                       11d17s19
            else if(left.lt.5)then                                        11d17s19
             iadd=0                                                       11d17s19
            else                                                          11d17s19
             if(ist.lt.19)then                                            11d17s19
              istp=ist+1                                                  11d17s19
              read(ostng(i)(istp:istp),*)left2                            11d17s19
              if(left2.ge.1)then                                          11d17s19
               iadd=1                                                     11d17s19
              else                                                        11d17s19
               iadd=0                                                     11d17s19
              end if                                                      11d17s19
             else                                                         11d17s19
              iadd=1                                                      11d17s19
             end if                                                       11d17s19
            end if                                                        11d17s19
            if(iadd.eq.1)then                                             11d17s19
             fact=1d1/(1d1**ndig)                                         11d17s19
             write(ostng(i),1551)bc(ieigint+im)+shift-fact               2d9s21
            end if                                                        11d17s19
           end if                                                        1d6s20
           do k=ist,19                                                   2d26s19
            ostng(i)(k:k)=' '                                            2d26s19
           end do                                                        2d26s19
          end if                                                         2d26s19
          iccode=1                                                       12d7s21
          if(ibc(iamconverged+im).ne.0)iccode=2                          12d7s21
          if(ibc(nlastde+im).ne.0)iccode=3                               12d7s21
          ibc(icombo+im)=iccode                                          12d7s21
         end do                                                          2d26s19
        end if                                                          7d17s24
        if(ngot.eq.0)then                                               12d29s21
         if(wconvci.ne.0d0)then                                         6d20s24
          write(6,5251)diffr,wconvci                                    9d17s24
 5251     format(x,'contract Davidson vectors down to best set. ',        9d17s24
     $         'diffr vs. wconvci: ',2es9.1)                            9d17s24
         else                                                           6d20s24
          write(6,*)('contract Davidson vectors down to best set')       12d29s21
         end if                                                         6d20s24
        else                                                            12d29s21
        write(6,2551)icit,(ostng(i),                                    4d21s21
     $       ccode(ibc(icombo+i-1)),                                    12d8s21
     $       i=1,nroot)                                                 4d21s21
        end if                                                          12d29s21
 2551   format(i5,50(a19,a1))                                           4d21s21
 3551   format(5x,(es10.2,9x))                                          2d26s19
       else                                                             2d26s19
        write(6,2552)icit,(bc(ieigint+i)+shift,i=0,nroot-1)             2d9s21
 2552   format(i5,(f19.12))                                             2d26s19
        write(6,3551)(bc(ieigint+i)-bc(ieold+i),i=0,nroot-1)            2d26s19
        do i=0,nroot-1                                                  3d12s21
         delta=abs(bc(ieigint+i)-bc(ieold+i))                           3d12s21
         changem=min(delta,changem)                                     3d12s21
        end do                                                          3d12s21
       end if
       if(nroot.eq.1)econvcii=min(econvcii,changem*0.1d0)               3d31s21
      end if
      call dws_bcast(econvcii,1)                                        3d12s21
      if(ngot.gt.0.or.icit.eq.1)then                                    12d29s21
       echange=0d0
       nb4=0                                                             4d23s21
       do i=0,nroot-1                                                    11d13s18
        if(dynwconv(2).eq.0d0)then                                       9d20s21
         delta=bc(ieigint+i)-bc(ieigint)                                  10d28s20
         deltaw=delta/dynwconv(1)                                         9d20s21
         ex=max(1d-10,exp(-deltaw))                                     1d19s23
         exi=1d0/ex                                                       10d28s20
         fff=(2d0/(ex+exi))**2                                            10d28s20
        else                                                             9d20s21
         wde=(bc(ieigint+i)+shift-dynwconv(2))/dynwconv(1)               9d20s21
         fff=4d0/(dynwconv(3)+exp(2d0*wde))                              9d20s21
        end if                                                           9d20s21
        dele=fff*abs(bc(ieigint+i)-bc(ieold+i))                          12d7s21
        echange=max(echange,fff*abs(bc(ieigint+i)-bc(ieold+i)))          10d28s20
        if(icit.gt.1)then                                                12d7s21
         if(dele.gt.econvci)then                                         12d7s21
          ratio=dele/bc(ilastde+i)                                       12d7s21
          if(ratio.gt.0.9d0.and.ratio.lt.2d0)then                        12d13s21
           if(lprint)write(6,*)('energy lowering ratio = '),ratio       3d21s22
           ibc(nlastde+i)=1                                              12d7s21
          else                                                           12d7s21
           ibc(nlastde+i)=0                                              12d7s21
          end if                                                         12d7s21
         end if                                                          12d7s21
        else                                                             12d7s21
         ibc(nlastde+i)=0                                                12d7s21
        end if                                                           12d7s21
        bc(ilastde+i)=dele                                               12d7s21
        ibc(iamconverged+i)=0                                            4d23s21
        if(abs(bc(ieigint+i)-bc(ieold+i)).lt.econvci.and.ndoub.eq.ndoubs 4d21s21
     $      .and.max(ngotd,ngots).gt.0)then                             4d23s21
         if(nb4.eq.i)then                                                4d23s21
          nb4=nb4+1                                                       4d23s21
         end if                                                          4d23s21
        end if                                                           4d23s21
        bc(ieold+i)=bc(ieigint+i)                                        11d13s18
       end do                                                            11d13s18
      end if                                                            12d29s21
      if(wconvci.eq.0d0)then                                            7d15s24
       call dws_bcast(echange,1)                                         4d14s21
      else                                                              7d15s24
       bc(ibcoff)=echange                                                7d12s24
       bc(ibcoff+1)=diffr                                                7d12s24
       call dws_bcast(bc(ibcoff),2)                                     7d15s24
       echange=bc(ibcoff)                                               7d15s24
       diffr=bc(ibcoff+1)                                               7d15s24
      end if                                                            7d15s24
      call dws_bcast(ibc(iamconverged),nroot)                           4d23s21
      call dws_bcast(ibc(ilastde),nroot)                                12d7s21
      call dws_bcast(ibc(nlastde),nroot)                                12d7s21
      if(ndoub.gt.0)then                                                2d18s21
       if(iunc(2).eq.0)then                                             7d7s21
        call ddi_zero(bc,ibc,igxx)                                      11d15s22
       end if                                                           7d7s21
       call dws_synca                                                    2d8s21
      end if                                                            2d18s21
      ivv=ivint+nroot*nctf(isymmrci)                                    4d8s20
      ivvn=ibcoff                                                       1d22s21
      ibcoff=ivvn+nvdim*nroot                                           1d22s21
      call enough('mrci. 33',bc,ibc)
      call vnormx(bc(ivv),bc(ivvn),nvdim,nroot,ldebugf)                 5d5s21
      if(.not.ldavid)then                                               12d27s21
       call dgemm('n','n',nctf(isymmrci),nroot,nvdim,1d0,bc(ihix),       1d23s21
     $     nctf(isymmrci),bc(ivvn),nvdim,0d0,bc(ihixold),nctf(isymmrci),1d23s21
     d' mrci.  1')
      end if                                                            12d27s21
      if(ndoub.ne.0)then                                                1d22s21
       if(ldebug)write(6,*)('mkdoub')
       if(iunc(2).eq.0)then                                             7d7s21
        call mkdoub(bc(ivvn),bc(ivv),nvdim,nrootu,nroot,bc(ivdinout),    1d23s21
     $     bc(ihdi),bc(ivdold),bc(igdold),ngotd,nfdat,nvirt,nsymb,multh,4d14s21
     $     isymmrci,ndoub,ldavid,maxdsr,bc,ibc)                         11d10s22
       else                                                             7d7s21
        call mkdoubuc(bc(ivvn),bc(ivv),nvdim,nrootu,nroot,ibc(ihddiag), 7d7s21
     $     bc(ivdold),bc(igdold),ngotd,ibc(nff2),ibc(icsf),ibc(icsf2),  7d7s21
     $       nvirt,mdon,mdoo,nsymb,multh,isymmrci,ldavid,maxdsr,bc,ibc) 11d10s22
       end if                                                           7d7s21
      end if                                                            1d22s21
      if(nsing.ne.0)then                                                1d22s21
       if(ldebug)write(6,*)('mksing')
       do j=0,nroot*2-1                                                 2d18s21
        bc(idotsing+j)=0d0                                              2d18s21
       end do                                                           2d18s21
       ifhis=ibcoff                                                      2d19s21
       itmp1=ifhis+nctf(isymmrci)*nvdim                                  2d19s21
       ibcoff=itmp1+nctf(isymmrci)*nroot                                 2d19s21
       call enough('mrci. 34',bc,ibc)
       do i=0,nctf(isymmrci)*nrootu-1                                    2d19s21
        bc(ifhis+i)=bc(ihis+i)                                           2d19s21
       end do                                                            2d19s21
       jfhis=ifhis+nctf(isymmrci)*nrootu                                 2d19s21
       do i=0,nctf(isymmrci)*ngots-1                                    4d14s21
        bc(jfhis+i)=bc(ihiso+i)                                          2d19s21
       end do                                                            2d19s21
       if(ldavid)then                                                   1d24s23
        do i=0,nctf(isymmrci)*nvdim-1                                    1d24s23
         bc(ihiso+i)=bc(ifhis+i)                                         1d24s23
        end do                                                           1d24s23
       end if                                                           1d24s23
       call dgemm('n','n',nctf(isymmrci),nroot,nvdim,1d0,bc(ifhis),      2d19s21
     $     nctf(isymmrci),bc(ivv),nvdim,0d0,bc(itmp1),nctf(isymmrci),   2d19s21
     d' mrci.  2')
       do ir=0,nroot-1                                                   2d19s21
        jtmp1=itmp1+nctf(isymmrci)*ir                                    2d19s21
        jvint=ivint+nctf(isymmrci)*ir                                   7d18s24
        dtt=0d0                                                          2d19s21
        do i=mynowprog,nctf(isymmrci)-1,mynprocg                         2d19s21
         dtt=dtt+bc(jtmp1+i)*bc(jvint+i)                                 2d19s21
        end do                                                           2d19s21
        jdotsing=idotsing+1+2*ir                                         2d19s21
        bc(jdotsing)=bc(jdotsing)+2d0*dtt                                2d19s21
       end do                                                            2d19s21
       if(.not.ldavid)then                                              1d24s23
        call dgemm('n','n',nctf(isymmrci),nroot,nvdim,1d0,bc(ifhis),      2d19s21
     $     nctf(isymmrci),bc(ivvn),nvdim,0d0,bc(ihiso),nctf(isymmrci),  2d19s21
     d' mrci.  3')
       end if                                                           1d24s23
       ibcoff=ifhis                                                      2d19s21
       call mksing(bc(ivvn),bc(ivv),nvdim,nrootu,nroot,ibc(ihsdiag),    1d23s21
     $     bc(ivsold),bc(igsold),ngots,ibc(nff1),ibc(icsf),nvirt,       4d14s21
     $     mdon,mdoo,nsymb,multh,isymmrci,bc(igss),bc(igsso),           2d18s21
     $     bc(idotsing),ldavid,maxdsr,nlocals,bc,ibc)                   11d10s22
       call dws_synca                                                   1d25s21
      end if                                                            1d22s21
      if(ldebugf)then                                                   3d24s23
       do i=0,2*nvdim*nvdim-1                                           1d27s21
        bc(ihxy+i)=0d0                                                  1d27s21
       end do                                                           1d27s21
       if(nsing.gt.0)then                                                1d21s21
        call hsingpart(bc(ihxy),bc(ioxy),nvdim,0,ibc(ihsdiag),          1d27s21
     $     ibc(nff1),ibc(icsf),nvirt,mdon,mdoo,nsymb,multh,isymmrci,    1d23s21
     $     bc(ivsold),bc(igsold),nroot,bc,ibc)                          11d10s22
        if(ndoub.eq.0.or.iunc(2).eq.0)then                              7d7s21
         nwds=nvdim*nvdim*2                                               1d24s21
         call dws_gsumf(bc(ihxy),nwds)                                    1d24s21
        end if                                                          7d7s21
       end if                                                            1d21s21
       if(ndoub.gt.0)then                                                1d21s21
        if(iunc(2).eq.0)then                                            7d7s21
         call hdoubpart(bc(ihxy),bc(ioxy),nvdim,0,bc(ivdinout),          1d27s21
     $     bc(ihdi),nfdat,nvirt,nsymb,multh,isymmrci,bc(ivdold),        1d23s21
     $     bc(igdold),nroot,bc(iepart),nroote,bc,ibc)                   1d12s23
        else                                                            7d7s21
         call hdoubpartuc(bc(ihxy),bc(ioxy),nvdim,0,ibc(ihddiag),       7d7s21
     $     ibc(nff2),ibc(icsf),ibc(icsf2),nvirt,mdon,mdoo,nsymb,multh,  7d7s21
     $       isymmrci,bc(ivdold),bc(igdold),nroot,bc(iepart),nroote,bc, 1d12s23
     $        ibc)                                                      11d10s22
         nwds=nvdim*nvdim*2                                               1d24s21
         call dws_gsumf(bc(ihxy),nwds)                                    1d24s21
        end if                                                          7d7s21
       end if                                                            1d21s21
      else
       if(ndoub.gt.0)then                                                1d21s21
        if(iunc(2).eq.0)then                                            7d7s21
         call hdoubpart(bc(ihxy),bc(ioxy),nvdim,0,bc(ivdinout),          1d27s21
     $     bc(ihdi),nfdat,nvirt,nsymb,multh,isymmrci,bc(ivdold),        1d23s21
     $     bc(igdold),nroot,bc(iepart),nroote,bc,ibc)                   1d12s23
        else                                                            7d7s21
         call hdoubpartuc(bc(ihxy),bc(ioxy),nvdim,0,ibc(ihddiag),       7d7s21
     $     ibc(nff2),ibc(icsf),ibc(icsf2),nvirt,mdon,mdoo,nsymb,multh,  7d7s21
     $       isymmrci,bc(ivdold),bc(igdold),nroot,bc(iepart),nroote,bc, 1d12s23
     $        ibc)                                                      11d10s22
        end if                                                          7d21s21
       end if                                                           7d21s21
      end if                                                            1d27s21
      if(nsing.gt.0)then                                                11d13s18
       call second(time1)                                               12d14s18
       if(ldebug)write(6,*)('calling hcsicsf'),ibcoff,
     $      icit
       if(ldavid)then                                                   12d28s21
        if(ngots+nroot.le.maxdsr)then                                   7d18s24
         ngots=ngots+nroot                                              12d28s21
        else                                                            12d28s21
         ngots=0                                                        12d28s21
        end if                                                          12d28s21
       else                                                             12d28s21
        ngots=nroot                                                      4d14s21
       end if                                                           12d28s21
       call hcsi(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(nff0),            11d25s20
     $      ibc(iff0),bc(ivint),nctf(isymmrci),ibc(icsf),nec,mdon,mdoo, 11d26s20
     $      nsymb,multh,ixw1,ixw2,ih0av,nh0av,ionexb,nvirt,0,ldebugf,   5d5s21
     $      maxbx,ionexbt,nwiacc,bc,ibc)                                11d10s22
       if(ldebug)write(6,*)('back'),ibcoff                              7d15s20
       call second(time2)                                               12d14s18
       telap=time2-time1-tovr                                           12d14s18
       timex(1)=timex(1)+telap                                          12d14s18
      end if                                                            11d13s18
      if(ndoub.ne.0)then                                                9d18s18
       call second(time1)                                               12d14s18
       if(ldebug)
     $      write(6,*)('calling hcdicsf after intcsf'),ibcoff,icit
       if(ldavid)then                                                   12d28s21
        if(ngotd+nroot.le.maxdsr)then                                   7d18s24
         ngotd=ngotd+nroot                                              12d28s21
        else                                                            12d28s21
         ngotd=0                                                        12d28s21
        end if                                                          12d28s21
       else                                                             12d28s21
        ngotd=nroot                                                      4d14s21
       end if                                                           12d28s21
       if(iunc(2).eq.0)then                                             7d7s21
        call hcdi(ibc(nff2),nsymb,mdon,mdoo,multh,isymmrci,nvirt,
     $      ibc(icsf),nec,ibc(icsf2),nfdat,irel,ism,kmats,irefo,ismult, 12d8s20
     $      bc(ihdi),ixw1,ixw2,norb,ibc(nff0),ibc(iff0),bc(ivint),      1d23s21
     $      nctf(isymmrci),nroot,ldebugf,srh,bc,ibc)                    11d10s22
       else                                                             7d7s21
        call hcdiuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ibc(nff0),         7d7s21
     $      ibc(iff0),bc(ivint),nctf(isymmrci),ibc(icsf),ibc(icsf2),nec,7d7s21
     $       mdon,mdoo,nsymb,multh,ixw1,ixw2,kmats,nvirt,0,             7d7s21
     $      ldebugf,maxbxd,srh,nwiacc,bc,ibc)                           11d10s22
       end if                                                           7d7s21
       if(ldebug)write(6,*)('back'),ibcoff                              7d15s20
       call second(time2)                                               12d14s18
       telap=time2-time1-tovr                                           12d14s18
       timex(4)=timex(4)+telap                                          12d14s18
      end if                                                            9d18s18
c     ?
      ibcoff=ihxy                                                       3d22s23
      go to 1000                                                        11d13s18
      return                                                            5d24s18
      end                                                               5d24s18
