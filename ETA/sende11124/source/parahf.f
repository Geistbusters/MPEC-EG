c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parahf(nowpro,numpro,jmats,kmats,mgmat,noc,i3x,itmax,  7d26s16
     $                  ieigr,ivecr,ehf,ih0,maxddi,                     2d24s10
     $                  ispt1,icore,nbasis,isymdws,myguessi,natom,ngaus,2d24s10
     $                  ibdat,isym,iapair,multh,ibstor,isstor,potdws,   3d23s12
     $                  ioooo,ionex,iwerner,idavopt,icas,nlzz,idorel,   8d20s15
     $                  ascale,idogrado,ipropsym,inocan,iapairg,dynw,   1d28s19
     $                  nbasisp,nbasisc,nwcont,lmax,xcartx,numeminus,   5d26s21
     $                  atnumx,extradata,nextradata,idorel4c,smallest,  7d28s22
     $                  denergy,iftype,fstgth,ndfld,bc,ibc,iforbs)      12d6s23
      implicit real*8 (a-h,o-z)
      real*16 fl                                                        2d26s20
      external second
c
c     parahf returns
c     1: fock matrix eigenvalues
c     2: mo coefficients
c     3: hf energy
c     4: h0 matrix in ao basis
c     5: ispt1: pointers
c
c     nbasisp(isb)*ncomp: number of primitive basis functions
c     nbasisc(isb): number of contracted basis functions
c     nbasdws(isb): number of optimized basis functions.
c     if non-relativistic caln, nbasdws(isb)=nbasisc(isb), otherwise
c     nbasdws(isb) lt nbasisc(isb), for we only optimize in the space
c     of h0 electronic functions (as opposed to positronic functions)
c
      logical lprint,myguess,lprt,lbug                                  5d13s22
      dimension times(26,3),timb(8,3),xiicnt(2),ieigs(8),ivecs(8),      6d28s04
     $          ivect(8),noc(8),iorbx(8),iorbn(8),iden1(8),iden2a(8),   6d28s04
     $          iden2b(8,8),iden2c(8,8),ipkla(8),ipklb(8,8),ipklc(8,8), 6d28s04
     $          iqkla(8),iqklb(8,8),iqklc(8,8),igmat(8,8),mgmat(8,8),   9d15s04
     $          igdig(8),nv4(8),iyc(8),ixyc(8),ibig2(8),icore(8),       9d10s07
     $          isymdws(8),llprod(3),multh(8,8),isym(3,8),iapair(3,*),  10d27s17
     $          ibstor(1),isstor(1),iso(8),iptoh(8,8,8),igmato(8,8),    3d19s12
     $          ieigv(8),idbly(8),iact(8),nbasdwsc(8),nvirtc(8),        1d13s16
     $          ieighs(8),icanon(8),morb(8),ixinv(8),ihessa(8,8),       4d25s17
     $          ihessb(8,8),ihessc(8,8),ibcode(8),isymg(3,8),morbp(8),  1d28s19
     $          iapairg(3,*),ihessd(8),iamatx(8),ignew(8),iumatb(8),    2d16s18
     $     itmatb(8),iamatb(8),rmssz0(8),rmssz(8),nbasisp(8),nbasisc(8),2d15s19
     $     morbc(8),xcartx(*),atnumx(*),ixlzz(8,6),ibasis(8),nfcn(8),   5d26s21
     $     nct(8),islz(3),extradata(*),nryd(8),ismile(8),ifockeig(8),   7d21s21
     $     iorbsymc(8),ivecso(8),isavu(8),dptype(4,2),denergy(*),       2d2s23
     $     makegbas(8),nposs(8),icanog(8),idoubg(8),iactog(8),          12d6s23
     $     iftype(*),fstgth(*)                                          12d6s23
      data dptype/1d0,-1d0,1d0,-1d0,1d0,1d0,-1d0,-1d0/                  5d9s22
      integer*8 ibstor,isstor                                           5d6s10
      integer*2 ipack2(2)                                               5d3s23
      equivalence (ipack2,ipack4)                                       5d3s23
      character*8 timlab(26)
      character*19 ostng                                                4d20s18
      logical log1
      character*10 output
      character*1 cdum                                                  1d24s06
      dimension ind(2),mjmat(1),mkmat(1),i3x(1),ipropsym(6),iovr(8),    5d7s21
     $     ipropsymg(6),isinfog(11)                                     5d7s21
      common/lowersymcm/nsymbgx,iptno(8),ipts(8),nhsz(8),ipao(8)        4d25s18
      COMMON/FACT16/FL(922),NCALL                                       2d26s20
      include "common.hf"
      include "common.cas"                                              8d8s14
      include "common.rys"
      include "common.spher"                                            2d20s20
      include "common.store"
      include "common.print"                                            1d3s20
      include "common.basis"                                            5d26s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/timerocm/tovr,telapo(15)                                   4d26s18
      common/singcm/iuse,nff
      dimension jmats(idbk),jdenpt(idbk),ipk(idbk),iumat(8),itmat(8),   9d13s04
     $          itpmat(8),itqmat(8),j2den(idbk),itotm(8),itotmb(8)      4d5s18
      dimension kmats(idbk),iqk(idbk),ioooo(idbk),ionex(idbk),          4d19s13
     $     ioooo2(idbk),jmatd(idbk),kmatd(idbk),iooood(idbk),           2d23s16
     $     ionexd(idbk),noocc(8),irtyp(5),iorbsym(8),ibasym(8),         4d19s21
     $     iorbsymz(8),npc(8),iorbsymao(8),multhg(8,8),dynw(3),iobsym(8)5d4s23
      common/getlcom/nlindws,nblindws                                   2d2s04
      common/gencm/return(8)                                            1d5s18
      times=0d0                                                         4d19s23
      timb=0d0                                                          4d19s23
      npost=0                                                           2d22s23
      if(iprtr(7).eq.0)then                                             1d3s20
       lprt=.false.                                                     1d3s20
       idwsdeb=0                                                        1d22s20
      else                                                              1d3s20
       lprt=.true.                                                      1d3s20
       idwsdeb=99999
      end if                                                            1d3s20
      lprint=nowpro.eq.0.or.idwsdeb.gt.0                                10d13s04
      nbasallp=0                                                        5d9s19
      do isb=1,nsymb                                                    10d24s17
       nryd(isb)=0                                                      1d18s20
       nbasallp=nbasallp+nbasisp(isb)                                   1d2s20
       iptno(isb)=ibcoff                                                10d24s17
       ipts(isb)=iptno(isb)+nbasisp(isb)                                2d15s19
       ibcoff=ipts(isb)+nbasisp(isb)                                    2d15s19
       ipao(isb)=ibcoff                                                 4d25s18
       ibcoff=ipao(isb)+3*nbasisp(isb)                                  2d15s19
      end do                                                            10d24s17
      do i=1,15                                                         4d26s18
       telapo(i)=0d0                                                        5d4s12
      end do
      myguess=myguessi.ne.0                                             2d18s16
      enobreit=0d0                                                      10d19s15
      eiglastgood=-1d10                                                 10d19s15
      xlam=1d0                                                          4d2s18
      thrs=1d-5                                                         1d10s07
      thrsa=1d-14                                                       3d8s07
      ienough=0
      iconv=0                                                           2d19s14
      icasvec=0                                                         4d18s18
      timlab(1)='1e ints '
      timlab(2)='diag s,.'
      timlab(3)='2e ints '
      timlab(4)='ao to mo'
      timlab(5)='energy  '
      timlab(6)='amat    '
      timlab(7)='gmat    '
      timlab(8)='aughess '
      timlab(9)='hessyc  '
      timlab(10)='schmidt '
      timlab(11)='build4  '
      timlab(12)='total   '
      timlab(13)='cannon  '                                             2d20s14
      timlab(14)='preden  '                                             2d20s14
      timlab(15)='preeri  '                                             2d20s14
      timlab(16)='bmat    '                                             2d20s14
      timlab(17)='diag g  '                                             2d20s14
      if(idwsdeb.ne.0)then
       write(6,*)('idwsdeb= '),idwsdeb
      end if
      call second(time1)
      call second(time2)
      tovr=time2-time1
c
c     for returned quantities
c
      ibsstor=ibcoff                                                    1d11s23
      iblstor=ibcoff                                                    1d11s23
      ncomp=1                                                           1d17s23
      if(idorel.eq.0)then                                               8d20s15
      else                                                              8d20s15
       if(idorel4c.eq.0)then                                            2d20s20
        ncomp=2                                                          8d20s15
       else                                                             2d20s20
        pi=acos(-1d0)                                                     2d20s20
        srpi=sqrt(pi)                                                     2d20s20
        nbaslarge=0                                                     2d20s20
        nbassmall=0                                                     2d20s20
        ilsz=ibcoff                                                     2d20s20
        ibcoff=ilsz+lmax+1                                              2d20s20
        call enough('parahf.  1',bc,ibc)                                3d13s23
        do l=1,lmax+1                                                   2d20s20
         lm=l-1
         nl=0                                                           2d20s20
         do isb=1,8                                                     2d20s20
          if(ipt(isb,l).gt.0)then                                       2d20s20
           nl=nl+ibc(ipt(isb,l))                                        2d20s20
          end if                                                        2d20s20
         end do                                                         2d20s20
         ibc(ilsz+lm)=nl                                                2d20s20
        end do                                                          2d20s20
        isnorm=ibcoff                                                   2d20s20
        ibcoff=isnorm+ngaus                                             2d20s20
        ibsstor=ibcoff                                                  2d21s20
        iblstor=ibsstor+ngaus                                           2d21s20
        ibcoff=iblstor+ngaus                                            2d21s20
        call enough('parahf.  2',bc,ibc)                                3d13s23
c
c     for small component, use large component gausians, but            2d20s20
c     increase l by 1, but use cartesian rather than spherical fcns.    2d20s20
c                                                                       2d20s20
        do i=0,ngaus-1
         ip=i+1
         ll=ibc(ibdat+i)                                                2d20s20
         ls=ll+1                                                        2d20s20
         nl=ibc(ilsz+ll)                                                2d20s20
         ns=ibc(ilsz+ls)                                                2d20s20
         ibc(iblstor+i)=nbaslarge                                       2d21s20
         ibc(ibsstor+i)=nbassmall                                       2d21s20
         nbaslarge=nbaslarge+nl                                         2d20s20
         nbassmall=nbassmall+ns                                         2d20s20
       zeta=bc(ibdat+i+ngaus)                                           2d20s20
       xnorm=0.25d0*srpi/sqrt((2d0*zeta)**3)                            2d20s20
       do il=0,ibc(ibdat+i)                                             2d20s20
        ff=0.5d0*(dfloat(il)+1.5d0)                                     2d20s20
        xnorm=xnorm*ff/zeta                                             2d20s20
       end do                                                           2d20s20
       xnorm=1d0/sqrt(xnorm)                                            2d20s20
       bc(isnorm+i)=xnorm                                               2d20s20
        end do
        if(lprint)then                                                  2d20s20
         write(6,*)('total number of large component basis fcns is '),   2d20s20
     $       nbaslarge                                                  2d20s20
         write(6,*)('total number of small component basis fcns is '),   2d20s20
     $       nbassmall                                                  2d20s20
        end if                                                          2d20s20
       end if                                                           2d20s20
      end if                                                            8d20s15
      ispt1=ibcoff
      ibcoff=ispt1+nbasis*2
      ispt1=ispt1                                                       8d1s12
      ieigr=ibcoff
      if(ivecr.ne.0)then                                                12d13s04xX
       ivecold=ivecr                                                    12d13s04xX
      end if                                                            12d13s04xX
      if(idorel4c.ne.0)then                                             2d20s20
       nbtot=2*(nbaslarge+nbassmall)                                    2d20s20
       nsqbas=nbtot*nbtot*2                                             2d20s20
      else                                                              2d20s20
       nsqbas=0                                                          8d8s14
       do isb=1,nsymb                                                    8d8s14
        nsqbas=nsqbas+(ncomp*nbasisp(isb))**2                            2d15s19
       end do                                                            8d8s14
      end if                                                            2d20s20
      ivecr=ieigr+nbasis                                                2d24s10
      ibcoff=ivecr+nsqbas                                               8d8s14
      ih0=ibcoff
      ibcoff=ih0+nsqbas                                                 2d15s19
      if(ibcoff.lt.0)then
       write(6,*)('ibcoff !!! '),ieigr,ivecr,ih0,nsqbas,nbasis
       call dws_sync
       call dws_finalize
       stop
      end if
      itopbc=ibcoff
      tstart=time2
      if(lprint)then                                                    10d13s04
       write(6,*)('entering parahf ')                                    2d24s10
       write(6,*)('output of processor no. '),nowpro
       write(6,*)('idorel = '),idorel
       write(6,*)('idogrado = '),idogrado
         write(6,*)('icasguess, myguess '),icasguess,myguess
       if(nsumcas.ne.0)then                                             5d3s21
        open(unit=42,file='summarycas')                                 5d3s21
       end if                                                           5d3s21
      end if                                                            2d2s16
      if(idorel.ne.0)then                                               2d14s19
       eposcut=-37557.724836488218d0                                    5d4s20
c
c     because of round of error ...
c
       if(lprint)write(6,*)('cut off energy for positron states: '),    2d14s19
     $      eposcut                                                     2d14s19
      end if                                                            2d14s19
      if(icas.eq.0)then                                                 8d8s14
       if(lprint)
     $     write(6,*)('performing single configuration HF calculation')    8d8s14
       do isb=1,nsymb                                                   8d20s14
        idbly(isb)=noc(isb)                                             8d20s14
        iact(isb)=0                                                     8d20s14
       end do                                                           8d20s14
       norba=0                                                          6d5s20
      else                                                              8d8s14
       myguess=icasguess.eq.3                                           5d3s16
       ncasr=0                                                          4d20s18
       nveccnt=ibcoff                                                   4d17s23
       ibcoff=nveccnt+nstate                                            4d17s23
       call enough('parahf.veccnt',bc,ibc)                              4d18s23
       do i=1,nstate                                                    4d20s18
        ncasr=ncasr+isinfo(4,i)                                         4d20s18
       end do                                                           4d20s18
       if(lprint)then
        write(6,*)('total no. of roots '),ncasr                         4d20s18
        write(6,*)('performing CAS calculation')                        8d8s14
        write(6,*)('icasguess = '),icasguess
        write(6,*)('nsymb: '),nsymb
        write(6,1327)(idoub(isb),isb=1,nsymb)                           6d12s23
 1327   format('idoub:',8i5)                                            6d12s23
        write(6,1328)(iacto(isb),isb=1,nsymb)                           6d12s23
 1328   format('iacto:',8i5)                                            6d12s23
        write(6,1329)(nocc2(isb),isb=1,nsymb)                           6d12s23
 1329   format('nocc :',8i5)                                            6d12s23
        write(6,*)('dynw = '),dynw
        write(6,*)('nusecsf = '),nusecsf                                12d24s19
        write(6,*)('nstate: '),nstate
       end if
       if(ncasr.eq.0)then                                               4d20s18
        if(lprint)write(6,*)('there are no roots to be calculated ...')
        call dws_sync                                                   4d20s18
        call dws_finalize                                               4d20s18
        stop                                                            4d20s18
       end if                                                           4d20s18
       icallcas=1                                                       8d8s14
       norba=0                                                          12d28s19
       do isb=1,nsymb                                                   8d20s14
        idbly(isb)=idoub(isb)                                            8d20s14
        iact(isb)=iacto(isb)                                            8d20s14
        norba=norba+iacto(isb)                                          12d28s19
        noc(isb)=nocc2(isb)                                             8d20s14
        noocc(isb)=ibcoff                                               4d5s18
        ibcoff=ibcoff+iacto(isb)                                        4d5s18
       end do                                                           8d20s14
      end if                                                            8d11s14
      call enough('parahf.  3',bc,ibc)                                  3d13s23
      if(icasguess.eq.2.and.icas.ne.0)then                              5d3s16
        jvecr=ivecold
        do i=1,nsymb                                                    8d8s14
         write(6,*)('initial orbs for symmetry '),i                     8d21s15
         call prntm2(bc(jvecr),nbasisp(i)*ncomp,nbasisc(i),             2d15s19
     $        nbasisp(i)*ncomp)                                         2d15s19
         jvecr=jvecr+(nbasisp(i)*ncomp)**2                              2d15s19
        end do                                                          8d8s14
      end if                                                            10d13s04
      if(lprint)write(6,*)('number of basis fcns = '),nbasis            2d24s10
      nbb=0                                                             5d3s10
      ieigv(1)=ieigr                                                    3d19s12
      do isb=1,nsymb                                                    5d3s10
       if(isb.gt.1)then                                                 3d19s12
        ieigv(isb)=ieigv(isb-1)+nbasisp(isb-1)                          2d15s19
       end if                                                           3d19s12
       iso(isb)=nbb                                                     5d3s10
       nbb=nbb+(ncomp*nbasisp(isb))**2                                  2d8d19
      end do                                                            5d3s10
      if(idorel4c.ne.0)nbb=nsqbas                                       2d21s20
      irelx=ibcoff                                                      12d28s19
      ismx=irelx+norba                                                  12d28s19
      ibcoff=ismx+norba                                                 12d28s19
      call enough('parahf.  4',bc,ibc)                                  3d13s23
      if(nusecsf.ge.0)then                                              12d28s19
       call precsf(nusecsf,nsymb,ibc(irelx),                            12d28s19
     $      ibc(ismx),norba,multh,ibasisp,icsfp,nfcnp,nctp,ixw1p,       8d2s22
     $     ixw2p,iptrbitp,ipoint2p,bc,ibc)                              11d9s22
      else                                                              12d28s19
       idorb=ibcoff                                                     12d28s19
       isorb=ibcoff                                                     12d28s19
       icsf=ibcoff                                                      12d28s19
       icsfpd=ibcoff                                                    12d28s19
       iptr=ibcoff                                                      12d28s19
      end if                                                            12d28s19
      ibstorc=ibcoff                                                    9d26s17
      ibcoff=ibstorc+nbasis                                             9d26s17
      idwsa=ibcoff                                                      9d2s06
      ibcoff=idwsa+nbb*20                                               5d3s10
      idwss=ibcoff                                                      9d2s06
      if(idorel4c.ne.0)then                                             2d20s20
       ibcoff=idwss+nsqbas                                              2d20s20
      else                                                              2d20s20
       ibcoff=idwss+nbb*20                                               5d3s10
      end if                                                            2d20s20
      if(iwerner.eq.0)then                                              11d29s17
       c0=1d0                                                           11d29s17
       c1=0d0                                                           11d29s17
       cjk=0d0                                                          11d29s17
       c3=0d0                                                           11d29s17
      else if(iwerner.eq.1)then                                         11d29s17
       c0=-3d0                                                          11d29s17
       c1=1d0                                                           11d29s17
       cjk=0d0                                                          11d29s17
       c3=0d0                                                           11d29s17
      else if(iwerner.eq.2)then                                         11d29s17
       c0=3d0                                                           11d29s17
       c1=-2d0                                                          11d29s17
       cjk=1d0                                                          11d29s17
       c3=0d0                                                           11d29s17
      else                                                              11d29s17
       c0=-1d0                                                          11d29s17
       c1=1d0                                                           11d29s17
       cjk=-1d0                                                         11d29s17
       c3=1d0                                                           11d29s17
      end if                                                            11d29s17
      if(lprint)then                                                    10d13s04
       if(idorel4c.ne.0)then                                            10d5s22
        write(6,*)('DHF scf method ')                                   10d5s22
       else                                                             10d5s22
        if(iwerner.ge.1)then                                              9d13s04
         write(6,*)('using werner method of order '),iwerner             11d29s17
         write(6,1330)c0,c1,cjk,c3                                      6d12s23
 1330    format('c0 =',f5.1,' c1 =',f5.1,' cjk = ',f5.1,' c3 = ',f5.1)  6d12s23
        else                                                              9d13s04
         write(6,*)('using 2nd order augmented hessian method ')          9d13s04
        end if                                                            9d13s04
       end if                                                           10d5s22
      end if                                                            10d13s04
c
c     get one electron integrals
c
      call second(time1)
      if(idorel4c.ne.0)then                                             2d20s20
       call parah04c(natom,ngaus,ibdat,nbasis,bc(ih0),bc(idwss),        2d24s20
     $       ibc(iblstor),ibc(ibsstor),nsqbas,idwsdeb,ascale,           2d24s20
     $            nbaslarge,nbassmall,ilsz,isnorm,nbtot,bc,ibc)         11d9s22
       nnbb=nbtot*nbtot                                                 2d24s20
       iscopy=ibcoff                                                    2d24s20
       ibcoff=iscopy+nnbb                                               2d24s20
       call enough('parahf.  5',bc,ibc)                                 3d13s23
       do i=0,nnbb-1                                                    2d24s20
        bc(iscopy+i)=bc(idwss+i)                                        2d24s20
       end do                                                           2d24s20
      else                                                              2d20s20
       do ifi=1,ndfld                                                   12d6s23
        if(iftype(ifi).gt.0.and.iftype(ifi).le.6)then                   12d6s23
         if(ipropsym(iftype(ifi)).ne.1)then                             12d6s23
          if(lprint)write(6,*)('finite field would break symmetry! ')    8d29s22
          call dws_synca                                                 8d29s22
          call dws_finalize                                              8d29s22
          stop 'parahf'                                                 12d6s23
         end if                                                          8d29s22
        end if                                                           8d29s22
       end do                                                           12d6s23
       call parah0(natom,ngaus,ibdat,nbasis,bc(ih0),bc(idwss),isym,      5d3s10
     $            iapair,ibstor,isstor,iso,nbb,idwsdeb,idorel,ascale,   8d29s22
     $      iftype,fstgth,ndfld,bc,ibc)                                 12d6s23
       if(nlzz.ne.0)then                                                 4d3s21
        do isb=1,nsymb                                                   4d3s21
         iorbsym(isb)=ibcoff                                             4d3s21
         iorbsymao(isb)=iorbsym(isb)+nbasisp(isb)*ncomp                 4d23s21
         iorbsymc(isb)=iorbsymao(isb)+nbasisp(isb)*ncomp                       4d23s21
         iobsym(isb)=iorbsymc(isb)+nbasisc(isb)                         5d4s23
         ibcoff=iobsym(isb)+nbasisc(isb)                                5d4s23
        end do                                                           4d3s21
        if(nlzz.gt.2)then                                                4d19s21
         do isb=1,nsymb                                                  4d19s21
          iorbsymz(isb)=ibcoff                                           4d19s21
          ibcoff=iorbsymz(isb)+nbasisp(isb)*ncomp                        4d19s21
         end do                                                          4d19s21
        end if                                                           4d19s21
        call setbqn(ngaus,ibdat,ibstor,isstor,iorbsym,iorbsymz,         4d19s21
     $       idorel,nlzz,nsymb,nbasisp,nbasisc,iapair,bc,ibc)           11d9s22
        if(nlzz.eq.2)then                                               4d23s21
         do isb=1,nsymb                                                 4d23s21
 7171     format(500i1)
          do i=0,nbasisp(isb)*ncomp-1                                   4d23s21
           ibc(iorbsymao(isb)+i)=ibc(iorbsym(isb)+i)                    4d23s21
          end do                                                        4d23s21
         end do                                                         4d23s21
        else                                                            4d23s21
         do isb=1,nsymb                                                 4d23s21
          do i=0,nbasisp(isb)*ncomp-1                                   4d23s21
           ibc(iorbsymao(isb)+i)=ibc(iorbsym(isb)+i)                    4d23s21
     $          +100*ibc(iorbsymz(isb)+i)                               4d23s21
          end do                                                        4d23s21
         end do                                                         4d23s21
        end if                                                          4d23s21
       end if                                                            4d3s21
      end if                                                            2d20s20
       isend=ibcoff
       irecv=isend+numpro                                               2d22s19
       iptx=irecv+numpro                                                2d20s20
       ipf=iptx+numpro                                                  2d20s20
       ibcoff=ipf+numpro                                                12d24s19
       call enough('parahf.  6',bc,ibc)                                 3d13s23
      if(idorel4c.eq.0)then                                             2d21s20
       nwst=0                                                            5d2s18
       do isb=1,nsymb                                                    5d2s18
        nsze=nbasisp(isb)                                                1d30s18
        if(idorel.ne.0)nsze=nsze*2                                       1d30s18
        nwst=nwst+nsze**2                                                1d30s18
       end do                                                            5d2s18
       iscopy=ibcoff
       ibcoff=iscopy+nwst*2
       call enough('parahf.  7',bc,ibc)                                 3d13s23
       iscopyp=iscopy+nwst
       jscopy=iscopy
       jdwss=idwss
       do isb=1,nsymb
        nsze=nbasisp(isb)                                                1d30s18
        if(idorel.ne.0)nsze=nsze*2                                       1d30s18
        nn=nsze*nsze                                                     1d30s18
        do i=0,nn-1
         bc(jscopy+i)=bc(jdwss+i)
        end do
        call square(bc(jscopy),nsze)                                     1d30s18
        if(idwsdeb.gt.5)then
         write(6,*)('copy of overlap matrix for symmetry block '),isb,
     $      jscopy,jscopy-iscopy,loc(bc(jscopy))
         call prntm2(bc(jscopy),nsze,nsze,nsze)                           1d30s18
        end if                                                           5d31s22
        iovr(isb)=jscopy                                                 2d8d19
        jscopy=jscopy+nn
        jdwss=jdwss+nn
       end do
       call contractg(bc(idwss),ngaus,ibdat,ibstor,isstor,idorel,0,      9d26s17
     $     nbasisp,nbasisc,iapair,nlzz,iorbsym,iorbsymz,idum,0,bc,ibc)  11d9s22
       jscopyp=iscopyp
       jdwss=idwss
       do isb=1,nsymb
        nsze=nbasisc(isb)                                                1d30s18
        nn=nsze**2                                                       1d30s18
        do i=0,nn-1
         bc(jscopyp+i)=bc(jdwss+i)
        end do
        call square(bc(jscopyp),nsze)
        nsze=nbasisp(isb)                                                1d30s18
        if(idorel.ne.0)nsze=nsze*2                                       1d30s18
        nn=nsze**2                                                       1d30s18
        jscopyp=jscopyp+nn
        jdwss=jdwss+nn
       end do
c     not used?
       jspt1=ispt1
       do isb=1,nsymb
        do i=1,nbasisp(isb)                                              2d15s19
         ibc(jspt1)=i
         jspt1=jspt1+1
         ibc(jspt1)=isb
         jspt1=jspt1+1
        end do
       end do
       call second(time2)
       telap=time2-time1-tovr
       times(1,1)=times(1,1)+telap
       times(1,2)=times(1,2)+telap**2
       times(1,3)=times(1,3)+1d0
      else                                                              2d21s20
       ncomp=1                                                          2d21s20
       nbasisc(1)=nbtot                                                 2d21s20
       nbasisp(1)=nbtot                                                 2d21s20
      end if                                                            2d21s20
c
      ioff=idwss
      nsymt=0                                                           6d28s04
      call second(time1)                                                10d19s04
      if(noc(1).lt.0)then                                               6d28s04
       not=iabs(noc(1))                                                 6d28s04
       if(lprint)then                                                   10d13s04
       write(6,*)('program to figure out sym. orb. occupancy based')
       write(6,*)('on eigenvalues of h0. total number occ orbs = '),not 6d28s04
       end if                                                           10d13s04
       ieigh=ibcoff                                                     6d28s04
       isymh=ieigh+min(nsymb*not,nbasis*ncomp)                          8d20s15
       ibcoff=isymh+min(nsymb*not,nbasis*ncomp)                         8d20s15
       jsymh=isymh                                                      6d28s04
       jeigh=ieigh                                                      6d28s04
      end if                                                            6d28s04
      iorb=ibcoff                                                       6d28s04
      jorb=iorb                                                         6d28s04
      nworb=0                                                           3d19s07
      do isb=1,nsymb                                                    6d28s04
       nc=nbasisc(isb)                                                  2d15s19
       np=nbasisp(isb)*ncomp                                            5d8s18
       jorb=jorb+nc*np                                                  1d28s19
       nworb=nworb+nc*np                                                1d28s19
      end do                                                            6d28s04
      ibcoff=jorb                                                       6d28s04
      do isb=1,nsymb                                                    1d13s16
       ieighs(isb)=ibcoff                                               1d13s16
       ibcoff=ibcoff+nbasisc(isb)                                       2d15s19
      end do                                                            1d13s16
      nvguess=0                                                         1d11s23
      ircode=0                                                          1d11s23
      ierr=0                                                            1d11s23
      if(myguess)then
       call enough('parahf.  8',bc,ibc)                                 3d13s23
       if(nowpro.eq.0)then                                              3d19s07
        if(iforbs.eq.0)then                                             3d13s23
         write(6,*)('try and grab orbs from file orbs')                       3d19s07
         open(unit=1,file='orbs',status='old',form='unformatted')        5d27s06
        else                                                            3d13s23
         write(6,*)('try and grab orbs from file forbs')                3d13s23
         open(unit=1,file='forbs',status='old',form='unformatted')      3d13s23
        end if                                                          3d13s23
        if(iabs(myguessi).eq.1)then                                     3d14s16
         write(6,*)('orbitals from ao basis are used')                  2d18s16
        else                                                            2d18s16
         write(6,*)('orbitals from orthogonal basis are used ')         2d18s16
        end if                                                          2d18s16
        read(1)nsymbg,idorelg,ngausg,natomg,nwcontg,idum1,idum2,idum3,  5d7s21
     $       multhg,dum,dum,ipropsymg,isinfog,nextradatag               5d25s21
        write(6,*)('nsymbg etc '),nsymbg,idorelg,ngausg,natomg          10d27s17
        idorelg=iabs(idorelg)                                           6d24s16
        if(idorelg.ne.0.and.idorel.ne.0)then                            2d13s20
         idorelg=idorel                                                 2d13s20
        end if                                                          2d13s20
        if(nsymbg.ne.nsymb)then                                         2d18s16
         write(6,*)('point group doesn''t match current one: '),nsymbg  2d18s16
        end if                                                          2d18s16
        if(idorelg.ne.idorel)then                                       2d18s16
         write(6,*)('myguessi = '),myguessi
         if(idorelg.eq.0)then                                           2d18s16
          write(6,*)('using non-relativistic orbitals as guess ')       2d18s16
          write(6,*)('force myguessi to be unity ')                     6d1s21
          myguessi=1                                                    6d1s21
         else if(idorel.eq.0)then                                       2d18s16
          write(6,*)('using relativistic orbitals as guess ')           2d18s16
         end if                                                         2d18s16
         if(abs(myguessi).ne.1)then                                     4d27s16
          write(6,*)('need to use ao rather than ob guess in this case')4d27s16
          stop 'parahf'                                                 10d7s22
         end if                                                         4d27s16
        end if                                                          2d18s16
        if(ngausg.ne.ngaus)then                                         2d18s16
         write(6,*)('number of gaussians does not match '),ngausg,ngaus              2d18s16
        end if                                                          2d18s16
        if(idorelg.eq.0)then                                            3d14s20
         ncompg=1                                                       3d14s20
        else                                                            3d14s20
         ncompg=2                                                       3d14s20
        end if                                                          3d14s20
        inbasg=ibcoff                                                    2d18s16
        inbasgp=inbasg+nsymbg                                           1d28s19
        inbasgc=inbasgp+nsymbg                                           1d28s19
        ibdatg=inbasgc+nsymbg                                           2d15s19
        nbdatg=ngausg*9+nwcontg
        ibcoff=ibdatg+nbdatg                                            1d28s19
        icang=ibcoff                                                    5d3s23
        ibcoff=icang+nsymbg                                             5d3s23
        idelt1=loc(nbasdwsc(2))-loc(nbasdwsc(1))
        idelt2=loc(ibc(inbasg+1))-loc(ibc(inbasg))
        if(idelt1.ne.idelt2)then
         read(1)(morb(isb),isb=1,nsymbg)
         do isb=1,nsymbg
          ibc(inbasg+isb-1)=morb(isb)
         end do
         read(1)(morbc(isb),isb=1,nsymbg)                               2d15s19
         do isb=1,nsymbg                                                2d15s19
          ibc(inbasgc+isb-1)=morbc(isb)                                 2d15s19
         end do                                                         2d15s19
         read(1)(morbp(isb),isb=1,nsymbg)
         do isb=1,nsymbg
          ibc(inbasgp+isb-1)=morbp(isb)
         end do
        else
         read(1)(ibc(inbasg+isb),isb=0,nsymbg-1)                         2d18s16
         read(1)(ibc(inbasgc+isb),isb=0,nsymbg-1)                       1d28s19
         read(1)(ibc(inbasgp+isb),isb=0,nsymbg-1)                       1d28s19
         do isb=1,nsymbg
          morb(isb)=ibc(inbasg+isb-1)
          morbc(isb)=ibc(inbasgc+isb-1)                                 1d28s19
          morbp(isb)=ibc(inbasgp+isb-1)                                 1d28s19
         end do
        end if
        nallg=0                                                         10d23s17
        do isb=1,nsymbg                                                 10d23s17
         nallg=nallg+ibc(inbasgp+isb-1)                                 5d9s19
        end do                                                          10d23s17
        nddy=natomg*6+nextradatag+1                                     5d7s21
        read(1)isymg,(bc(ibcoff+j),j=0,nddy-1),(idoubg(isb),iactog(isb),5d4s23
     $       isb=1,nsymb),nlzzg                                         5d4s23
        ehfbsofarl0=bc(ibcoff+nddy-1)                                    5d7s21
        do isb=1,nsymb                                                  5d4s23
         ipack4=idoubg(isb)                                             5d4s23
         icanog(isb)=ipack2(2)                                          5d4s23
        end do
        read(1)((iapairg(j,i),j=1,3),i=1,natomg)                        10d27s17
        isstorg=ibcoff                                                   10d23s17
        ibstorg=isstorg+nallg                                           10d23s17
        ibcoff=ibstorg+nallg                                            10d23s17
        call enough('parahf.  9',bc,ibc)                                3d13s23
        nbdat=ngausg*9+nwcontg                                              9d12s17
        read(1)(bc(ibdatg+i),i=0,nbdat-1)                               1d28s19
        read(1)(ibc(ibstorg+i),ibc(isstorg+i),i=0,nallg-1)              10d23s17
        do isb=1,nsymbg                                                 10d23s17
         ibcode(isb)=ibcoff                                             10d23s17
         ibcoff=ibcoff+ibc(inbasg+isb-1)                                10d23s17
         call enough('parahf. 10',bc,ibc)                               3d13s23
         read(1)(ibc(ibcode(isb)+i),i=0,ibc(inbasg+isb-1)-1)            10d23s17
        end do                                                          10d23s17
        nn=0                                                            2d18s16
        do isb=0,nsymbg-1                                               2d18s16
         nn=nn+ibc(inbasg+isb)*ibc(inbasgp+isb)                         1d28s19
        end do                                                          2d18s16
        nn=nn*ncompg                                                    3d14s20
        ivguess=ibcoff                                                  2d18s16
        ibcoff=ivguess+nn                                               2d18s16
        ivdum=ibcoff                                                    2d18s16
        ibcoff=ivdum+nn                                                 2d18s16
        call enough('parahf. 11',bc,ibc)                                3d13s23
        jvguess=ivguess                                                 2d18s16
        do isb=0,nsymbg-1                                               2d18s16
         nh=ncompg*ibc(inbasgp+isb)*ibc(inbasg+isb)                     3d14s20
         write(6,*)('myguessi '),myguessi
         if(iabs(myguessi).eq.1.or.ngausg.ne.ngaus)then                 4d28s21
          read(1)(bc(jvguess+i),i=0,nh-1)                               2d18s16
          jvguess=jvguess+nh                                            2d15s19
         else                                                           2d18s16
          read(1)(bc(ivdum+i),i=0,nh-1)                                 2d18s16
         end if                                                         2d18s16
         nh=ibc(inbasg+isb)*ibc(inbasgc+isb)                            2d1s23
         if(myguessi.eq.2.and.ngausg.eq.ngaus)then                      4d28s21
          read(1)(bc(jvguess+i),i=0,nh-1)                               2d18s16
          jvguess=jvguess+nh                                            2d15s19
         else                                                           2d18s16
          read(1)(bc(ivdum+i),i=0,nh-1)                                 2d18s16
         end if                                                         2d18s16
         nb=ibc(inbasg+isb)                                             10d23s17
        end do                                                          2d18s16
        if(iwfromfile.ne.0)then                                         5d24s19
         write(6,*)('now grab initial weights from file: ')
         read(1)nstateg
         write(6,*)('nstateg vs nstate: '),nstateg,nstate               5d24s19
         if(nstateg.ne.nstate)then                                      5d24s19
          stop 'nstateg'                                                5d24s19
         end if                                                         5d24s19
         do i=1,nstate                                                  5d24s19
          read(1)iwh                                                    5d24s19
          write(6,*)('iwh vs isinfo: '),i,iwh,isinfo(4,i)
          if(iwh.ne.isinfo(4,i))stop 'isinfo'
          read(1)(eweight(j,i),j=1,isinfo(4,i))                         5d24s19
         end do                                                         5d24s19
        end if                                                          5d24s19
        ibcoff=ivdum                                                    2d18s16
        nsymbgx=nsymbg                                                  10d24s17
        call makeguess(idorel,idorelg,nsymb,nsymbg,ngaus,ngausg,        2d18s16
     $       bc(ibdat),bc(ibdatg),nbasdws,morb,bc(ivecr),               11d18s19
     $       bc(ivguess),bc(idwss),myguessi,ierr,idwsdeb,ibstor,        4d27s21
     $       isstor,                                                    2d15s19
     $       ibc(ibstorg),ibc(isstorg),ibcode,iptno,ipts,iapair,isym,   10d27s17
     $       iapairg,isymg,nhsz,ipao,morbp,nbasisp,0,idorel4c,ascale,   3d27s20
     $       ibc(iblstor),ibc(ibsstor),nbaslarge,nbassmall,ilsz,        3d27s20
     $       isnorm,nbasisc,nlzz,iorbsym,iorbsymz,morbc,morb,           5d7s21
     $       bc(iovr(1)),nrydb,nvguess,ircode,bc,ibc,makegbas,icanog)   5d3s23
        close(unit=1)                                                    3d19s07
       end if                                                            3d19s07
       xmyguessi=dfloat(myguessi)                                       6d1s21
       bc(ibcoff)=dfloat(myguessi)                                      10d11s22
       bc(ibcoff+1)=dfloat(nvguess)                                     10d11s22
       bc(ibcoff+2)=dfloat(ircode)                                      10d13s22
       call dws_bcast(bc(ibcoff),3)                                     10d13s22
       myguessi=nint(bc(ibcoff))                                        10d11s22
       nvguess=nint(bc(ibcoff+1))                                       10d11s22
       ircode=nint(bc(ibcoff+2))                                        10d13s22
       iarg1=nworb                                                       3d19s07
       if(ierr.ne.0)then                                                2d18s16
        bc(ivecr)=132d0                                                 2d18s16
       end if                                                           2d18s16
       if(idorel4c.ne.0)then                                            3d27s20
        iarg1=nbtot*nbaslarge*4                                         3d28s20
       end if                                                           3d27s20
       call dws_bcast(bc(ivecr),iarg1)                                  2d24s10
       if(abs(bc(ivecr)-132d0).lt.1d-8)then                             2d18s16
        write(6,*)('error from orbital guess routine ')                 2d18s16
        call dws_sync                                                   2d18s16
        call dws_finalize                                               2d18s16
        stop                                                            2d18s16
       end if                                                           2d18s16
       ivecold=ivecr                                                    1d11s16
       if(iwfromfile.ne.0)then                                          5d26s19
        nsend=maxst2*maxst1                                             5d26s19
        call dws_bcast(eweight,nsend)                                   5d26s19
       end if                                                           5d26s19
      end if                                                            3d19s07
      jorb=iorb                                                         6d28s04
      nsy1(1)=0                                                         6d28s04
      if(nlzz.ne.0)then                                                 6d7s21
       nxsb=0                                                           6d7s21
       do isb=1,nsymb                                                   6d7s21
        nxsb=nxsb+nbasisc(isb)                                          6d7s21
       end do                                                           6d7s21
       mxsb=7                                                           7d21s21
       ixsb=ibcoff                                                      6d7s21
       ibcoff=ixsb+nxsb*mxsb                                            6d7s21
       call enough('parahf. 12',bc,ibc)                                 3d13s23
       nxsb=0                                                           6d7s21
       do isb=1,nsymb                                                   6d7s21
        do i=0,nbasisc(isb)-1                                           6d7s21
         if(nlzz.eq.2)then                                              6d7s21
          itest=ibc(iorbsym(isb)+i)                                     6d7s21
         else                                                           6d7s21
          itest=ibc(iorbsymz(isb)+i)                                    6d7s21
         end if
         do j=0,nxsb-1                                                  6d7s21
          if(itest.eq.ibc(ixsb+j*mxsb))then                             6d7s21
           ibc(ixsb+j*mxsb+1)=ibc(ixsb+j*mxsb+1)+1                      6d7s21
           go to 3313                                                   6d7s21
          end if                                                        6d7s21
         end do                                                         6d7s21
         ibc(ixsb+nxsb*mxsb)=itest                                      6d7s21
         ibc(ixsb+nxsb*mxsb+1)=1                                        6d7s21
         ibc(ixsb+nxsb*mxsb+6)=0                                        7d21s21
         nxsb=nxsb+1                                                    6d7s21
 3313    continue                                                       6d7s21
        end do                                                          6d7s21
        do i=0,nbasisp(isb)*ncomp-1                                     7d21s21
         if(nlzz.eq.2)then                                              8d18s21
          itest=ibc(iorbsymao(isb)+i)                                   8d18s21
         else                                                           8d18s21
          itest=ibc(iorbsymao(isb)+i)/100                               8d18s21
         end if                                                         8d18s21
         do j=0,nxsb-1                                                  7d21s21
          if(itest.eq.ibc(ixsb+j*mxsb))then                             8d18s21
           ibc(ixsb+j*mxsb+6)=ibc(ixsb+j*mxsb+6)+1                      7d21s21
           go to 3413                                                   7d21s21
          end if                                                        7d21s21
         end do                                                         7d21s21
         write(6,*)('failed to match iorbsymao: '),i,isb,
     $        ibc(iorbsymao(isb)+i)
         write(6,*)('got: '),(ibc(ixsb+j*mxsb),j=0,nxsb-1)
         call dws_synca                                                 7d21s21
         call dws_finalize                                              7d21s21
         stop 'parahforbsymao'                                          7d21s21
 3413    continue                                                       7d21s21
        end do                                                          7d21s21
       end do                                                           6d7s21
       if(lprint)write(6,*)('number of unique extra symmetries: '),nxsb 8d23s21
       ibcoff=ixsb+mxsb*nxsb                                            6d7s21
       do j=0,nxsb-1                                                    6d7s21
        jp=j+1                                                          6d7s21
        if(nlzz.eq.6)then                                               6d7s21
         ndeg=2*ibc(ixsb+j*mxsb)+1                                      6d7s21
        else                                                            6d7s21
         if(ibc(ixsb+j*mxsb).eq.0.or.ibc(ixsb+j*mxsb).eq.50)then        10d1s21
          ndeg=1                                                        6d7s21
         else                                                           6d7s21
          ndeg=2                                                        6d7s21
         end if                                                         6d7s21
        end if                                                          6d7s21
        nall=ibc(ixsb+j*mxsb+1)                                         6d7s21
        nleft=mod(nall,ndeg)                                            6d7s21
        nhere=nall/ndeg                                                 6d7s21
        nhereao=ibc(ixsb+j*mxsb+6)/ndeg                                 7d21s21
        if(lprint)then
         if(j.eq.0)write(6,3933)                                        7d21s21
 3933    format(5x,'#',4x,'qn',2x,'nall',x,'nhere',x,'nleft',3x,'deg',  7d21s21
     $        x,'nhereao')                                                7d21s21
         if(nlzz.ne.2.or.nsymb.ne.8)then                                10d1s21
          write(6,3934)jp,ibc(ixsb+j*mxsb),ibc(ixsb+j*mxsb+1),           7d21s21
     $       nhere,nleft,ndeg,nhereao                                   7d21s21
         else                                                           10d1s21
          lsprim=ibc(ixsb+j*mxsb)                                       10d1s21
          if(lsprim.ge.50)then                                          10d1s21
           cdum='u'                                                     10d1s21
           lsprim=lsprim-50                                             10d1s21
          else                                                          10d1s21
           cdum='g'                                                     10d1s21
          end if                                                        10d1s21
          write(6,3935)jp,lsprim,cdum,ibc(ixsb+j*mxsb+1),               10d1s21
     $       nhere,nleft,ndeg,nhereao                                   7d21s21
         end if                                                         10d1s21
 3934    format(7i6)                                                    7d21s21
 3935    format(i6,i5,a1,5i6)                                           10d1s21
        end if                                                          7d21s21
        ibc(ixsb+mxsb*j+1)=nhere                                        6d7s21
        ibc(ixsb+mxsb*j+2)=ibcoff                                       6d7s21
        ibc(ixsb+mxsb*j+3)=0                                            6d7s21
        ibcoff=ibcoff+nhere*nhere                                       6d7s21
        ibc(ixsb+mxsb*j+4)=ibcoff                                       6d7s21
        ibcoff=ibcoff+nhere*nhereao                                     7d21s21
        ibc(ixsb+mxsb*j+5)=ibcoff                                       6d7s21
        ibcoff=ibcoff+nhere*nhere                                       6d7s21
        ibc(ixsb+mxsb*j+6)=nhereao                                      7d21s21
        call enough('parahf. 13',bc,ibc)                                3d13s23
       end do                                                           6d7s21
      end if                                                            6d7s21
      do isb=1,nsymb
       icanon(isb)=0                                                    11d28s22
       if(isb.gt.1)then                                                 6d28s04
        nsy1(isb)=nsy1(isb-1)+(nbasisp(isb-1)*ncomp)**2                 1d28s19
       end if                                                           6d28s04
       if(nbasisc(isb).gt.0)then                                        2d15s19
        nwds=(nbasdws(isb)*ncomp*(nbasdws(isb)*ncomp+1))/2              8d20s15
        iarg1=nwds
        if(idwsdeb.gt.5)then
         write(6,*)('overlap matrix for block '),isb,ioff
         call mpprnt2(bc(ioff),nbasisc(isb))                            2d15s19
        end if
        if(idorel4c.eq.0)then                                           2d21s20
         call square(bc(ioff),nbasisc(isb))                              2d15s19
        else                                                            2d21s20
        end if                                                          2d21s20
        nsymt=nsymt+(nbasisp(isb)*ncomp)**2                             2d14s19
        ieigs(isb)=ibcoff                                                6d28s04
        ivecs(isb)=ieigs(isb)+nbasisc(isb)                              2d15s19
c
c     first store orthogonalizing transformation, then overlap in       2d15s19
c     primitive basis.                                                  2d15s19
c     and if rel, orginal orthogonalizing transformation (as opposed to 2d9s20
c     the orthogonalizing transformation transformed to h0 electronic   2d9s20
c     states, which we will overwrite the start of ivecs)               2d9s20
c     if we have positron states, save h0 vectors in ob basis as well.  2d1s23
c
        if(idorel.eq.0)then                                             2d9s20
         ibcoff=ivecs(isb)                                              4d1s22
     $       +nbasisp(isb)*ncomp*(nbasisc(isb)+nbasisp(isb)*ncomp)      2d15s19
        else                                                            2d9s20
         ibcoff=ivecs(isb)                                              4d1s22
     $       +nbasisp(isb)*ncomp*(nbasisc(isb)*2+nbasisp(isb)*ncomp)    2d9s20
     $       +nbasisc(isb)*(nbasisc(isb)+nbasisp(isb)*ncomp)            2d1s23
        end if                                                          2d9s20
        idwst1=ibcoff                                                   4d1s22
        ibcoff=idwst1+nbasisc(isb)**2                                   2d15s19
        call enough('parahf. 14',bc,ibc)                                3d13s23
        nbas3=nbasisc(isb)                                              2d15s19
        nbas4=nbasisp(isb)*ncomp                                        5d2s18
        do 199 i=1,nbas3**2                                             2d15s19
         bc(idwst1+i-1)=bc(ioff+i-1)                                     6d28s04
  199   continue                                                         6d28s04
        isymdws(1)=ibcoff                                                9d10s07
        ibcoff=isymdws(1)+nbas3                                         2d15s19
        call enough('parahf. 15',bc,ibc)                                3d13s23
        if(nlzz.ne.0)then                                               4d19s21
         if(nlzz.eq.2)then                                              4d19s21
          do i=0,nbas3-1                                                4d19s21
           ibc(isymdws(1)+i)=ibc(iorbsym(isb)+i)                        4d19s21
          end do                                                        4d19s21
         else                                                           4d19s21
          do i=0,nbas3-1                                                4d19s21
           ibc(isymdws(1)+i)=ibc(iorbsym(isb)+i)                        4d19s21
     $          +100*ibc(iorbsymz(isb)+i)                               4d19s21
          end do                                                        4d19s21
         end if                                                         4d19s21
         igu=1                                                          8d31s23
         if(nsymb.eq.8)then                                             8d31s23
          if(isb.eq.2.or.isb.eq.3.or.isb.eq.5.or.isb.eq.8)igu=2         8d31s23
         end if                                                         8d31s23
         call diagy(nbas3,bc(idwst1),bc(ieigs(isb)),bc(ivecs(isb)),     4d19s21
     $       ibc(isymdws(1)),bc,ibc,nlzz,nlambda(1,igu),drangcut)       9d1s23
         inuq=ibcoff                                                    6d7s21
         nuq=0                                                          6d7s21
         do i=0,nbas3-1
          ibc(iobsym(isb)+i)=ibc(isymdws(1)+i)                          5d4s23
          do j=0,nuq-1                                                  6d7s21
           if(ibc(isymdws(1)+i).eq.ibc(inuq+j))go to 3136               6d7s21
          end do                                                        6d7s21
          ibcoff=inuq+nuq                                               6d7s21
          call enough('parahf. 16',bc,ibc)                              3d13s23
          ibc(inuq+nuq)=ibc(isymdws(1)+i)                               6d7s21
          nuq=nuq+1                                                     6d7s21
 3136     continue                                                      6d7s21
         end do
         ibcoff=inuq+nuq                                                6d7s21
         iherei=ibcoff                                                  6d7s21
         ibcoff=iherei+nuq*2                                            6d7s21
         call enough('parahf. 17',bc,ibc)                               3d13s23
         iherei0=ibcoff                                                 6d7s21
         do i=0,nuq-1                                                   6d7s21
          juse=-1                                                       6d7s21
          do j=0,nxsb-1                                                 6d7s21
           if(ibc(inuq+i).eq.ibc(ixsb+mxsb*j).or.                       6d7s21
     $        ibc(inuq+i)/100.eq.ibc(ixsb+mxsb*j))juse=j                6d7s21
          end do                                                        6d7s21
          if(juse.ge.0)then                                             6d7s21
          else                                                          6d7s21
           write(6,*)('could not match symmetries ')
           stop 'parahf'                                                6d7s21
          end if                                                        6d7s21
          ibc(iherei+i*2)=ibcoff                                        6d7s21
          ibc(iherei+i*2+1)=ibcoff                                      6d7s21
          ibcoff=ibcoff+ibc(ixsb+mxsb*juse+1)**2                        6d7s21
         end do                                                         6d7s21
         call enough('parahf. 18',bc,ibc)                               3d13s23
         do i=iherei0,ibcoff-1                                          6d7s21
          bc(i)=0d0                                                     6d7s21
         end do                                                         6d7s21
         do i=0,nbas3-1                                                 6d7s21
          if(nlzz.eq.2)then                                             6d7s21
           lsprim=ibc(iorbsym(isb)+i)                                   6d7s21
          else                                                          6d7s21
           lsprim=ibc(iorbsym(isb)+i)                                   6d7s21
     $          +100*ibc(iorbsymz(isb)+i)                               4d19s21
          end if                                                        6d7s21
          jxsbh=-1
          do j=0,nuq-1                                                  6d7s21
           if(lsprim.eq.ibc(inuq+j))then                                6d7s21
            jxsbh=iherei+j*2                                            6d7s21
            go to 3133                                                  6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
 3133     continue                                                      6d7s21
          do j=0,nbas3-1                                                6d7s21
           jjii=ivecs(isb)+i+nbas3*j                                    6d7s21b*
           if(ibc(isymdws(1)+j).eq.lsprim)then                          6d7s21
            bc(ibc(jxsbh))=bc(jjii)                                     6d7s21
            ibc(jxsbh)=ibc(jxsbh)+1                                     6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
         end do                                                         6d7s21
         inuqg=ibcoff                                                   6d7s21
         ibcoff=inuqg+nxsb                                              6d7s21
         call enough('parahf. 19',bc,ibc)                               3d13s23
         do i=0,nxsb-1                                                  6d7s21
          ibc(inuqg+i)=0                                                6d7s21
         end do                                                         6d7s21
         do i=0,nuq-1                                                   6d7s21
          if(ibc(iherei+2*i).ne.ibc(iherei+2*i+1))then                  6d7s21
           juse=-1                                                      6d7s21
           do j=0,nxsb-1                                                6d7s21
            if(ibc(inuq+i).eq.ibc(ixsb+j*mxsb).or.                      6d7s21
     $         ibc(inuq+i)/100.eq.ibc(ixsb+j*mxsb))juse=j               6d7s21
           end do                                                       6d7s21
           if(juse.ge.0)then                                            6d7s21
           else                                                         6d7s21
            write(6,*)('could not match symmetry ')
            stop 'parahf'                                               6d7s21
           end if                                                       6d7s21
           nsz=ibc(ixsb+mxsb*juse+1)                                    6d7s21
           if(ibc(ixsb+juse*mxsb+3).eq.0.and.ibc(inuqg+juse).eq.0)then  6d7s21
            do j=0,nsz*nsz-1                                            6d7s21
             bc(ibc(ixsb+mxsb*juse+2)+j)=bc(ibc(iherei+2*i+1)+j)        6d7s21
            end do                                                      6d7s21
            ibc(inuqg+juse)=1                                           6d7s21
           else                                                         6d7s21
            do j=0,nsz*nsz-1                                            6d7s21
             bc(ibc(iherei+2*i+1)+j)=bc(ibc(ixsb+mxsb*juse+2)+j)        6d7s21
            end do                                                      6d7s21
           end if                                                       6d7s21
           ibc(iherei+2*i)=ibc(iherei+2*i+1)                            6d7s21
          end if                                                        6d7s21
         end do                                                         6d7s21
         ibcoff=inuqg                                                   6d7s21
         do i=0,nbas3-1                                                 6d7s21
          if(nlzz.eq.2)then                                             6d7s21
           lsprim=ibc(iorbsym(isb)+i)                                   6d7s21
          else                                                          6d7s21
           lsprim=ibc(iorbsym(isb)+i)                                   6d7s21
     $          +100*ibc(iorbsymz(isb)+i)                               4d19s21
          end if                                                        6d7s21
          do j=0,nuq-1                                                  6d7s21
           if(lsprim.eq.ibc(inuq+j))then                                6d7s21
            jxsbh=iherei+j*2                                            6d7s21
            go to 3132                                                  6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
 3132     continue                                                      6d7s21
          do j=0,nbas3-1                                                6d7s21
           jjii=ivecs(isb)+i+nbas3*j                                    6d7s21b*
           if(ibc(isymdws(1)+j).eq.lsprim)then                          6d7s21
            bc(jjii)=bc(ibc(jxsbh))                                     6d7s21
            ibc(jxsbh)=ibc(jxsbh)+1                                     6d7s21
           else                                                         6d7s21
            bc(jjii)=0d0                                                6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
         end do                                                         6d7s21
         ibcoff=inuq                                                    6d7s21
        else                                                            4d19s21
         call diagx(nbas3,bc(idwst1),bc(ieigs(isb)),bc(ivecs(isb)),     4d19s21
     $       ibc(isymdws(1)),bc,ibc)                                    11d14s22
        end if                                                          4d19s21
        nwds=nbas3*(nbas3+1)                                            2d15s19
        call dws_bcast(bc(ieigs(isb)),nwds)                             4d22s14
        i18=1
         iarg2=nbas3                                                    2d15s19
  222    format(5es22.14)
         if(idwsdeb.gt.5)then
          write(6,*)('eigenvectors of s '),isb
          call prntm2(bc(ivecs(isb)),nbas3,nbas3,nbas3)
         end if
c
c     copy of vectors for gradient calculation
c
        ibcoff=idwst1                                                    6d28s04
        ivecso(isb)=ibcoff                                              4d1s22
        ibcoff=ivecso(isb)+nbas3*nbas3+nbas3*nbas4                      4d6s22
        ivecsop=ivecso(isb)+nbas3*nbas3                                 4d6s22
        ixinv(isb)=ibcoff                                               4d1s22
        ibcoff=ibcoff+2*nbasisc(isb)*nbasisc(isb)                       4d6s22
        call enough('parahf. 20',bc,ibc)                                3d13s23
        do i=0,nbas3*nbas3-1                                            4d4s22
         bc(ivecso(isb)+i)=bc(ivecs(isb)+i)                             4d1s22
         bc(ivecsop+i)=bc(ivecs(isb)+i)                                 4d6s22
        end do                                                          1d13s16
        if(idorel4c.eq.0)then                                           10d5s22
         call contractg(bc(ivecsop),ngaus,ibdat,ibstor,isstor,idorel,    4d6s22
     $       isb,nbasisp,nbasisc,iapair,0,idum,idum,idum,0,bc,ibc)      11d9s22
         jscopy=iscopy                                                   4d8s22
         do jsb=1,isb-1
          nn=(ncomp*nbasisp(jsb))**2                                     1d30s18
          jscopy=jscopy+nn
         end do
        end if                                                          10d5s22
        jeigs=ieigs(isb)-1                                               6d28s04
        if(lprint)write(6,445)bc(ieigs(isb)),bc(ieigs(isb)+nbas3-1)     12d2s22
  445   format('range of S eigenvalues: ',2es9.2)                       12d2s22
c     for debugging ...
        if(bc(ieigs(isb)).lt.sneglect)then                               1d18s07
         if(lprint)write(6,*)('for symmetry '),isb,                     5d4s23
     $       (', smallest eigenvalue of S is too small,'),              5d4s23
     $      ('so use cannonalical orthogonalization neglecting ')       5d2s23
        end if                                                           11d14s05
        do 26 i=1,nbasisc(isb)                                          2d15s19
         if(bc(jeigs+i).gt.sneglect)then                                 1d18s07
          bc(jeigs+i)=1d0/sqrt(bc(jeigs+i))                               2d23s04
         else
          if(lprint)then                                                9d12s24
           write(6,*)('eigenvalue no '),i,bc(jeigs+i)                   9d12s24
          end if                                                        9d12s24
          icanon(isb)=icanon(isb)+1                                     5d2s23
          bc(jeigs+i)=0d0                                                11d14s05
         end if
   26   continue
        ivect(isb)=ibcoff                                                6d28s04
        ibcoff=ivect(isb)+nbasisp(isb)*nbasisc(isb)*ncomp               2d15s19
        itmpdws=ibcoff
        ibcoff=itmpdws+nbasisc(isb)**2                                  2d15s19
        itmpdwsi=ibcoff                                                 4d1s22
        ibcoff=itmpdwsi+nbasisc(isb)**2                                 4d1s22
        call enough('parahf. 22',bc,ibc)                                3d13s23
        do i=1,nbasisc(isb)                                             2d15s19
         tmp=bc(jeigs+i)
         tmpi=1d0/tmp                                                   4d1s22
         do j=1,nbasisc(isb)                                            2d15s19
          iad1=itmpdws+j-1+nbas3*(i-1)                                  8d20s15
          iad2=ivecs(isb)+j-1+nbas3*(i-1)                               8d20s15
          iad3=itmpdwsi+j-1+nbas3*(i-1)                                 4d1s22
          bc(iad1)=bc(iad2)*tmp
          bc(iad3)=bc(iad2)*tmpi                                        4d1s22
         end do
        end do
        nbasdws(isb)=nbas3-icanon(isb)                                  5d2s23
        if(icanon(isb).eq.0)then                                        1d29s16
         call dgemm('n','t',nbas3,nbas3,nbas3,1d0,bc(itmpdws),nbas3,
     $       bc(ivecs(isb)),nbas3,0d0,bc(ivect(isb)),nbas3,
     d' parahf.  1')
         call dgemm('n','t',nbas3,nbas3,nbas3,1d0,bc(itmpdwsi),nbas3,   4d1s22
     $       bc(ivecs(isb)),nbas3,0d0,bc(ixinv(isb)),nbas3,             4d1s22
     d' parahf.  1')                                                    4d1s22
         ixinvp=ixinv(isb)+nbas3*nbas3                                  4d6s22
         do i=1,nbas3*nbas3
          bc(ivecs(isb)+i-1)=bc(ivect(isb)+i-1)
          bc(ixinvp+i-1)=bc(ivect(isb)+i-1)                             4d6s22
         end do
         if(idwsdeb.gt.5)then                                           5d31s22
          write(6,*)('x for '),ixinvp                                    4d6s22
          call prntm2(bc(ixinvp),nbas3,nbas3,nbas3)
         end if                                                         5d31s22
         if(nlzz.ne.0)then                                              5d4s23
          do i=0,nbas3-1                                                5d4s23
           ibc(iobsym(isb)+i)=ibc(iorbsym(isb)+i)                       5d4s23
          end do                                                        5d4s23
          if(nlzz.eq.6)then                                             5d4s23
           do i=0,nbas3-1                                                5d4s23
            ibc(iobsym(isb)+i)=ibc(iobsym(isb)+i)                       5d4s23
     $           +100*ibc(iorbsymz(isb)+i)                              5d4s23
           end do                                                        5d4s23
          end if                                                        5d4s23
         end if                                                         5d4s23
        else                                                             11d14s05
         itmpo=itmpdws+nbas3*icanon(isb)                                5d2s23
         if(nlzz.ne.0)then                                              9d12s24
          do i=0,nbasdws(isb)-1                                           5d4s23
           ip=i+icanon(isb)                                               5d4s23
           ibc(iobsym(isb)+i)=ibc(iobsym(isb)+ip)                        5d4s23
          end do                                                          5d4s23
         end if                                                         9d12s24
         do i=0,nbas3*nbasdws(isb)-1                                    5d2s23
          bc(ivecs(isb)+i)=bc(itmpo+i)                                  5d2s23
         end do                                                          11d14s05
         if(idwsdeb.gt.5)then
          write(6,*)('using Cannonical orthogalization ')
          write(6,*)('ivecs for isb = '),isb
          call prntm2(bc(ivecs(isb)),nbas3,nbasdws(isb),nbas3)
         end if
        end if                                                           11d14s05
        if(idorel4c.eq.0)then                                           2d21s20
         istmpx=ibcoff                                                   4d27s18
         ibcoff=istmpx+nbasisc(isb)*nbasisp(isb)*ncomp                   2d14s19
         call enough('parahf. 23',bc,ibc)                               3d13s23
         if(idwsdeb.gt.5)then                                           5d2s23
          jscopyp=iscopyp                                                 5d2s18
          do jsb=1,isb-1
           nn=(ncomp*nbasisp(jsb))**2                                     1d30s18
           jscopyp=jscopyp+nn
          end do
          if(nbasdws(isb).gt.0)then                                     7d3s23
           call dgemm('n','n',nbas3,nbasdws(isb),nbas3,1d0,               5d2s23
     $       bc(jscopyp),nbas3,bc(ivecs(isb)),nbas3,0d0,                1d30s18
     $       bc(istmpx),nbas3,                                          1d30s18
     d' parahf.  2')
          end if                                                        7d3s23
          write(6,*)('S*X ')
          call prntm2(bc(istmpx),nbas3,nbasdws(isb),nbas3)
          do i=0,nbasdws(isb)-1                                          5d2s23
           do j=0,nbas3-1                                                 1d30s18
            ji=istmpx+j+nbas3*i                                           1d30s18
            ij=jscopyp+i+nbasdws(isb)*j                                  5d2s23
            bc(ij)=bc(ji)
           end do
          end do
          write(6,*)('X*S ')
          call prntm2(bc(jscopyp),nbasdws(isb),nbas3,nbasdws(isb))
          if(nbasdws(isb).gt.0)then                                     7d3s23
           call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbas3,1d0,        5d2s23
     $       bc(jscopyp),nbasdws(isb),bc(ivecs(isb)),nbas3,0d0,         5d2s23
     $       bc(istmpx),nbasdws(isb),                                   5d2s23
     d' parahf.  3')
          end if                                                        7d3s23
          write(6,*)('what we have for istmpx=X*S*X ')
          call prntm2(bc(istmpx),nbasdws(isb),nbasdws(isb),nbasdws(isb))
          write(6,*)('vecs going into contractg ')
          call prntm2(bc(ivecs(isb)),nbas3,nbasdws(isb),nbas3)
         end if                                                         5d2s23
         call contractg(bc(ivecs(isb)),ngaus,ibdat,ibstor,isstor,idorel, 9d26s17
     $       isb,nbasisp,nbasisc,iapair,0,idum,idum,idum,0,bc,ibc)      11d9s22
         if(idwsdeb.gt.5)then                                           5d5s23
          write(6,*)('vecs coming out of contractg ')                   5d5s23
          call prntm2(bc(ivecs(isb)),nbasisp(isb)*ncomp,nbasisc(isb),   5d5s23
     $         nbasisp(isb)*ncomp)                                      5d5s23
         end if                                                         5d5s23
         jscopy=iscopy
         do jsb=1,isb-1
          nn=(ncomp*nbasisp(jsb))**2                                     1d30s18
          jscopy=jscopy+nn
         end do
         issave=ivecs(isb)+nbasisc(isb)*nbasisp(isb)*ncomp               2d15s19
         do i=0,nbas4*nbas4-1                                            2d15s19
          bc(issave+i)=bc(jscopy+i)                                      2d15s19
         end do                                                          2d15s19
         if(idorel.ne.0)then                                             2d9s20
          ixsave=ivecs(isb)+nbas4*(nbasisc(isb)+nbas4)                   2d9s20
          do i=0,nbasisc(isb)*nbas4-1                                    2d9s20
           bc(ixsave+i)=bc(ivecs(isb)+i)                                 2d9s20
          end do                                                         2d9s20
         end if                                                          2d9s20
         ibcoff=istmpx
        end if                                                          2d21s20
        do 27 i=1,nbasdws(isb)                                          5d2s23
         jvecs=ivecs(isb)-1+(i-1)*nbas4                                 1d30s18
         jvect=ivect(isb)+i-1-nbasdws(isb)                              5d2s23
         do 28 j=1,nbas4                                                1d30s18
          bc(jvect+j*nbasdws(isb))=bc(jvecs+j)                          5d2s23
   28    continue
   27   continue
        if(idwsdeb.gt.5)then                                            3d17s20
         writE(6,*)('orthogonalization matrix')
         iarg1=nbasisp(isb)*ncomp
         iarg2=nbasdws(isb)                                             5d2s23
         call prntm2(bc(ivecs(isb)),iarg1,iarg2,iarg1)                    6d28s04
         jprt=ivecs(isb)-1
        end if
        if(idorel4c.ne.0)then                                           2d21s20
         if(myguess)then                                                3d27s20
          jvecr=ivecr+nbaslarge*2*nbtot                                 3d27s20
         end if                                                         3d27s20
         call dhf(nbtot,ih0,nbb,ivecs,ilsz,idbly,idoub,iacto,isinfo,    3d1s20
     $       icas,idorel4c,eposcut,nnbb,idorel,ascale,nbaslarge,        3d1s20
     $       nbassmall,isnorm,iblstor,ibsstor,nsqbas,natom,ngaus,ibdat, 3d1s20
     $       nbasis,bc(iscopy),potdws,myguess,ivecr,lprint,bc,ibc)      11d9s22
         return                                                         11d17s22
        else                                                            2d21s20
        jdwsk=ioff+ih0-idwss                                            6d28s04
        iarg1=nwds
        call square(bc(jdwsk),nbas4)                                    1d30s18
        if(idwsdeb.gt.5)then
         write(6,*)('one electron ham '),jdwsk
         iarg1=nbasisp(isb)*ncomp                                       1d28s19
         call prntm2(bc(jdwsk),iarg1,iarg1,iarg1)                         6d28s04
         write(6,*)('ivecs ')
         call prntm2(bc(ivecs(isb)),nbas4,nbasdws(isb),nbas4)
        end if
        idwst1=ibcoff                                                    6d28s04
        ibcoff=idwst1+nbasisp(isb)*nbasisc(isb)*ncomp                   2d14s19
        call enough('parahf. 24',bc,ibc)                                3d13s23
        if(nbasdws(isb).gt.0)then                                       7d3s23
         call dgemm('n','n',nbas4,nbasdws(isb),nbas4,1d0,bc(jdwsk),     7d3s23
     $       nbas4,bc(ivecs(isb)),nbas4,0d0,bc(idwst1),nbas4,           7d3s23
     d' parahf.  4')
        end if                                                          7d3s23
        idwst3=ibcoff
        ibcoff=idwst3+nbas4*nbas4
        call enough('parahf. 25',bc,ibc)                                3d13s23
        if(nbasdws(isb).gt.0)then                                       7d3s23
         call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbas4,1d0,         5d2s23
     $       bc(ivect(isb)),nbasdws(isb),bc(idwst1),nbas4,0d0,          5d2s23
     $       bc(idwst3),nbasdws(isb),                                   5d2s23
     d' parahf.  5')
        end if                                                          7d3s23
        end if                                                          2d21s20
        if(idwsdeb.gt.5)then
         writE(6,*)('h0 in ob')
         iarg1=nbasdws(isb)                                             5d2s23
         call prntm2(bc(idwst3),iarg1,iarg1,iarg1)                      6d28s04
        end if
        isymdws(1)=ibcoff                                                9d10s07
        ibcoff=isymdws(1)+nbas4                                          9d10s07
        if(nlzz.ne.0)then                                               4d19s21
         do i=0,nbas3-1
          ibc(isymdws(1)+i)=ibc(iobsym(isb)+i)                          5d4s23
         end do
         call diagy(nbasdws(isb),bc(idwst3),bc(ieighs(isb)),bc(jorb),   5d2s23
     $       ibc(isymdws(1)),bc,ibc,0,idum,dum)                         9d1s23
         if(nlzz.eq.2.and.(isb.eq.1.or.isb.eq.5))then                   8d25s22
          nn00=0                                                        8d25s22
          do i=0,noc(isb)-1                                             8d25s22
           lcode=ibc(isymdws(1)+i)                                      8d25s22
           if(mod(lcode,50).ne.0)nn00=nn00+1                            8d25s22
          end do                                                        8d25s22
          isbp=isb+3                                                    8d25s22
          if(noc(isbp).eq.0.and.nn00.gt.0)then                          8d25s22
           if(lprint)then                                               8d25s22
            write(6,*)('we have a single component of Lambda gt 0')       8d25s22
            write(6,*)('swap for next Lambda=0 orbital ')                 8d25s22
           end if                                                       8d25s22
10101      continue                                                      8d25s22
           do i=0,noc(isb)-1                                             8d25s22
            lcodei=ibc(isymdws(1)+i)                                     8d25s22
            if(mod(lcodei,50).ne.0)then                                  8d25s22
             if(lprint)write(6,*)('swapping out orbital '),i+1          8d25s22
             do j=i+1,nbasdws(isb)-1                                    5d2s23
              lcodej=ibc(isymdws(1)+j)                                   8d25s22
              if(mod(lcodej,50).eq.0)then                                8d25s22
               if(lprint)write(6,*)('with orbital '),j+1                8d25s22
               ibc(isymdws(1)+j)=lcodei                                  8d25s22
               ibc(isymdws(1)+i)=lcodej                                  8d25s22
               do k=0,nbasdws(isb)-1                                    5d4s23
                iad=jorb+k+nbasdws(isb)*i                               5d4s23
                jad=jorb+k+nbasdws(isb)*j                               5d4s23
                cpy=bc(iad)                                              8d25s22
                bc(iad)=bc(jad)                                          8d25s22
                bc(jad)=cpy                                              8d25s22
               end do                                                    8d25s22
               go to 10101                                               8d25s22
              end if                                                     8d25s22
             end do                                                      8d25s22
            end if
           end do                                                        8d25s22
          end if                                                        8d25s22
         end if                                                         8d25s22
         inuq=ibcoff                                                    6d7s21
         nuq=0                                                          6d7s21
         do i=0,nbasdws(isb)-1                                          5d2s23
          do j=0,nuq-1                                                  6d7s21
           if(ibc(isymdws(1)+i).eq.ibc(inuq+j))go to 3137               6d7s21
          end do                                                        6d7s21
          ibcoff=inuq+nuq                                               6d7s21
          call enough('parahf. 26',bc,ibc)                              3d13s23
          ibc(inuq+nuq)=ibc(isymdws(1)+i)                               6d7s21
          nuq=nuq+1                                                     6d7s21
 3137     continue                                                      6d7s21
         end do
         ibcoff=inuq+nuq                                                6d7s21
         iherei=ibcoff                                                  6d7s21
         ibcoff=iherei+nuq*2                                            6d7s21
         call enough('parahf. 27',bc,ibc)                               3d13s23
         iherei0=ibcoff                                                 6d7s21
         do i=0,nuq-1                                                   6d7s21
          juse=-1                                                       6d7s21
          do j=0,nxsb-1                                                 6d7s21
           if(ibc(inuq+i).eq.ibc(ixsb+mxsb*j).or.                       6d7s21
     $        ibc(inuq+i)/100.eq.ibc(ixsb+mxsb*j))juse=j                6d7s21
          end do                                                        6d7s21
          if(juse.ge.0)then                                             6d7s21
          else                                                          6d7s21
           write(6,*)('could not match symmetries ')
           stop 'parahf'                                                6d7s21
          end if                                                        6d7s21
          ibc(iherei+i*2)=ibcoff                                        6d7s21
          ibc(iherei+i*2+1)=ibcoff                                      6d7s21
          ibcoff=ibcoff+ibc(ixsb+mxsb*juse+1)**2                        6d7s21
         end do                                                         6d7s21
         call enough('parahf. 28',bc,ibc)                               3d13s23
         do i=iherei0,ibcoff-1                                          6d7s21
          bc(i)=0d0                                                     6d7s21
         end do                                                         6d7s21
         do i=0,nbasdws(isb)-1                                          5d4s23
          lsprim=ibc(iobsym(isb)+i)                                     5d4s23
          do j=0,nuq-1                                                  6d7s21
           if(lsprim.eq.ibc(inuq+j))then                                6d7s21
            jxsbh=iherei+j*2                                            6d7s21
            go to 3134                                                  6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
 3134     continue                                                      6d7s21
          do j=0,nbasdws(isb)-1                                         5d2s23
           jjii=jorb+i+nbasdws(isb)*j                                   5d4s23
           if(ibc(isymdws(1)+j).eq.lsprim)then                          6d7s21
            bc(ibc(jxsbh))=bc(jjii)                                     6d7s21
            ibc(jxsbh)=ibc(jxsbh)+1                                     6d7s21
           else                                                         6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
         end do                                                         6d7s21
         inuqg=ibcoff                                                   6d7s21
         ibcoff=inuqg+nxsb                                              6d7s21
         call enough('parahf. 29',bc,ibc)                               3d13s23
         do i=0,nxsb-1                                                  6d7s21
          ibc(inuqg+i)=0                                                6d7s21
         end do                                                         6d7s21
         do i=0,nuq-1                                                   6d7s21
          if(ibc(iherei+2*i).ne.ibc(iherei+2*i+1))then                  6d7s21
           juse=-1                                                      6d7s21
           do j=0,nxsb-1                                                6d7s21
            if(ibc(inuq+i).eq.ibc(ixsb+j*mxsb).or.                      6d7s21
     $         ibc(inuq+i)/100.eq.ibc(ixsb+j*mxsb))juse=j               6d7s21
           end do                                                       6d7s21
           if(juse.ge.0)then                                            6d7s21
           else                                                         6d7s21
            write(6,*)('could not match symmetry ')
            stop 'parahf'                                               6d7s21
           end if                                                       6d7s21
           nsz=ibc(ixsb+mxsb*juse+1)                                       6d7s21
           if(ibc(ixsb+juse*mxsb+3).eq.0.and.ibc(inuqg+juse).eq.0)then  6d7s21
            do j=0,nsz*nsz-1                                            6d7s21
             bc(ibc(ixsb+mxsb*juse+4)+j)=bc(ibc(iherei+2*i+1)+j)        6d7s21
            end do                                                      6d7s21
            ibc(inuqg+juse)=1                                           6d7s21
            if(.not.(myguess.or.(icas.ne.0.and.icasguess.ne.1)))then    6d7s21
             ibc(ixsb+juse*mxsb+3)=1                                     6d7s21
            end if
           else                                                         6d7s21
            do j=0,nsz*nsz-1                                            6d7s21
             bc(ibc(iherei+2*i+1)+j)=bc(ibc(ixsb+mxsb*juse+4)+j)        6d7s21
            end do                                                      6d7s21
           end if                                                       6d7s21
           ibc(iherei+2*i)=ibc(iherei+2*i+1)                            6d7s21
          end if                                                        6d7s21
         end do                                                         6d7s21
         ibcoff=inuqg                                                   6d7s21
         do i=0,nbasdws(isb)-1                                          5d4s23
          lsprim=ibc(iobsym(isb)+i)                                     5d4s23
          do j=0,nuq-1                                                  6d7s21
           if(lsprim.eq.ibc(inuq+j))then                                6d7s21
            jxsbh=iherei+j*2                                            6d7s21
            go to 3135                                                  6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
 3135     continue                                                      6d7s21
          do j=0,nbasdws(isb)-1                                         5d2s23
           jjii=jorb+i+nbasdws(isb)*j                                   5d4s23
           if(ibc(isymdws(1)+j).eq.lsprim)then                          6d7s21
            bc(jjii)=bc(ibc(jxsbh))                                     6d7s21
            ibc(jxsbh)=ibc(jxsbh)+1                                     6d7s21
           else                                                         6d7s21
            bc(jjii)=0d0                                                      6d7s21
           end if                                                       6d7s21
          end do                                                        6d7s21
         end do                                                         6d7s21
         ibcoff=inuq                                                    6d7s21
         if(nlzz.eq.2)then                                              5d7s21
          do i=0,nbasdws(isb)-1                                         5d2s23
           ibc(iorbsym(isb)+i)=ibc(isymdws(1)+i)                        5d7s21
          end do                                                        5d7s21
         else                                                           5d7s21
          do i=0,nbasdws(isb)-1                                         5d2s23
           ibc(iorbsymz(isb)+i)=ibc(isymdws(1)+i)/100                   5d7s21
           ibc(iorbsym(isb)+i)=ibc(isymdws(1)+i)-100*ibc(iorbsymz(isb)) 5d7s21
          end do                                                        5d7s21
         end if                                                         5d7s21
        else                                                            4d19s21
         call diagx(nbasdws(isb),bc(idwst3),bc(ieighs(isb)),bc(jorb),   5d2s23
     $       ibc(isymdws(1)),bc,ibc)                                    11d14s22
        end if                                                          4d19s21
        iassgn=0                                                        12d15s22
        if(idorel.ne.0)then                                             8d21s15
         if(lprt)write(6,*)('for symmetry no. '),isb
         nposs(isb)=0                                                   5d2s23
         do i=0,nbasdws(isb)-1                                          5d2s23
          if(bc(ieighs(isb)+i).lt.eposcut)nposs(isb)=nposs(isb)+1       5d2s23
         end do                                                         2d14s19
         if(lprint)write(6,*)('number of positron states = '),nposs(isb)5d2s23
         npost=npost+nposs(isb)                                         5d2s23
         if(nposs(isb).eq.0)then                                        5d2s23
          if(lprint)then
           write(6,*)('we are good - there are no positron states! :)') 10d27s20
          end if                                                        10d27s20
          ismile(isb)=1                                                 10d27s20
          icodeao=0                                                     10d12s22
         else                                                           10d27s20
          ismile(isb)=0                                                 10d27s20
          if(lprint)write(6,*)('transform to electronic state basis')   10d27s20
          itmpe=ibcoff                                                   2d14s19
          itmpv=itmpe+nbasisc(isb)                                       2d15s19
          ibcoff=itmpv+nbasisc(isb)*nbasisc(isb)                         2d15s19
          if(nlzz.ne.0)then                                             12d15s22
           iltmp=ibcoff                                                 12d15s22
           ibcoff=iltmp+nbasisc(isb)                                    12d15s22
           jltmp=iltmp                                                  12d15s22
           iztmp=ibcoff                                                 12d15s22
           if(nlzz.ne.2)then                                            12d15s22
            ibcoff=iztmp+nbasisc(isb)                                   12d15s22
            jztmp=iztmp                                                 12d15s22
           end if                                                       12d15s22
          end if                                                        12d15s22
          call enough('parahf. 31',bc,ibc)                              3d13s23
          jtmpe=itmpe                                                    2d14s19
          jtmpv=itmpv                                                    2d14s19
          do i=0,nbasdws(isb)-1                                         5d2s23
           if(bc(ieighs(isb)+i).gt.eposcut)then                          2d15s19
            bc(jtmpe)=bc(ieighs(isb)+i)                                  2d15s19
            if(nlzz.ne.0)then                                           12d15s22
             ibc(jltmp)=ibc(iorbsym(isb)+i)                             12d15s22
             jltmp=jltmp+1                                              12d15s22
             if(nlzz.ne.2)then                                          12d15s22
              ibc(jztmp)=ibc(iorbsymz(isb)+i)                           12d15s22
              jztmp=jztmp+1                                             12d15s22
             end if                                                     12d15s22
            end if                                                      12d15s22
            jtmpe=jtmpe+1                                                2d14s19
            korb=jorb+nbas3*i                                            2d14s19
            do j=0,nbasisc(isb)-1                                        2d15s19
             bc(jtmpv+j)=bc(korb+j)                                      2d14s19
            end do                                                       2d14s19
            jtmpv=jtmpv+nbas3                                            2d14s19
           end if                                                        2d14s19
          end do                                                         2d14s19
          neles=nbas3-nposs(isb)                                        5d2s23
          icpy2=ivecs(isb)                                              2d1s23
     $         +nbasisp(isb)*ncomp*(nbasisc(isb)*2+nbasisp(isb)*ncomp)  2d1s23
          do i=0,nbas3*neles-1                                          2d1s23
           bc(icpy2+i)=bc(itmpv+i)                                      2d1s23
          end do                                                        2d1s23
          jtmpe=itmpe                                                    2d14s19
          jtmpv=itmpv                                                    2d14s19
          if(nlzz.ne.0)then                                             12d15s22
           jltmp=iltmp                                                  12d15s22
           if(nlzz.ne.2)jztmp=iztmp                                     12d15s22
          end if                                                        12d15s22
          do i=0,neles-1                                                 2d14s19
           bc(ieighs(isb)+i)=bc(jtmpe)                                   2d15s19
           if(nlzz.ne.0)then                                            12d15s22
            bc(iorbsym(isb)+i)=bc(jltmp)                                12d15s22
            jltmp=jltmp+1                                               12d15s22
            if(nlzz.ne.2)then                                           12d15s22
             bc(iorbsymz(isb)+i)=bc(jztmp)                              12d15s22
             jztmp=jztmp+1                                              12d15s22
            end if                                                      12d15s22
           end if                                                       12d15s22
           jtmpe=jtmpe+1                                                 2d14s19
           korb=jorb+nbas3*i                                             2d14s19
           do j=0,nbas3-1                                                2d15s19
            bc(korb+j)=bc(jtmpv+j)                                       2d14s19
           end do                                                        2d14s19
           jtmpv=jtmpv+nbas3                                             2d14s19
          end do                                                         2d14s19
          ibcoff=itmpe                                                   2d14s19
          nbasdws(isb)=neles                                             2d14s19
          if(lprint)write(6,*)('thus number of electronic states is '),  2d15s19
     $        nbasdws(isb)                                              2d15s19
          itmpx=ibcoff                                                   1d30s18
          ibcoff=itmpx+nbas4*nbasisc(isb)                                1d30s18
          call enough('parahf. 33',bc,ibc)                              3d13s23
          icpy2=ivecs(isb)                                              2d1s23
     $           +nbasisp(isb)*ncomp*(nbasisc(isb)*2+nbasisp(isb)*ncomp)2d1s23
     $       +nbasisc(isb)*nbasisc(isb)                                 2d1s23
          do i=0,nbas4*nbas3-1                                          2d1s23
           bc(icpy2+i)=bc(ivecs(isb)+i)                                 2d1s23
          end do                                                        2d1s23
          call dgemm('n','n',nbas4,nbasdws(isb),nbas3,1d0,              10d27s20
     $         bc(ivecs(isb)),nbas4,bc(jorb),nbas3,0d0,bc(itmpx),nbas4, 10d27s20
     d' parahf.  6')
          do i=0,nbas4*nbasdws(isb)-1                                    2d15s19
           bc(ivecs(isb)+i)=bc(itmpx+i)                                  1d30s18
          end do                                                         1d30s18
          if(myguess)then                                               12d15s22
           if(idwsdeb.gt.5)then                                         5d5s23
            write(6,*)('and mo guess vectors in ob basis? '),
     $          nbasisc(isb)
            call prntm2(bc(ivecold),nbasisc(isb),nbasdws(isb),
     $         nbasisc(isb))                                            12d15s22
           end if
           call dgemm('t','n',nbasdws(isb),nbasdws(isb),nbasisc(isb),   2d1s23
     $          1d0,bc(ivecold),nbasisc(isb),bc(ivecold),nbasisc(isb),  2d1s23
     $          0d0,bc(ibcoff),nbasdws(isb))                            12d15s22
           if(idwsdeb.gt.0)then                                         5d5s23
            write(6,*)('ortho test ')
            call prntm2(bc(ibcoff),nbasdws(isb),nbasdws(isb),
     $          nbasdws(isb))
           end if                                                       5d5s23
           sz=0d0                                                        12d15s22
           do i1=0,nbasdws(isb)-1                                        12d15s22
            do i2=0,i1-1                                                 12d15s22
             i12=ibcoff+i1+nbasdws(isb)*i2                               12d15s22
             i21=ibcoff+i2+nbasdws(isb)*i1                               12d15s22
             sz=sz+bc(i12)**2+bc(i21)**2                                 12d15s22
            end do                                                       12d15s22
            i11=ibcoff+i1*(nbasdws(isb)+1)                               12d15s22
            sz=sz+(bc(i11)-1d0)**2                                       12d15s22
           end do                                                        12d15s22
           sz=sqrt(sz/dfloat(nbasdws(isb)*nbasdws(isb)))                 12d15s22
           if(idwsdeb.gt.5)write(6,*)('rms size = '),sz
           if(sz.lt.1d-8)then                                            12d15s22
            igrabfrom=ivecs(isb)                                        2d1s23
     $           +nbasisp(isb)*ncomp*(nbasisc(isb)*2+nbasisp(isb)*ncomp)2d1s23
     $       +nbasisc(isb)*nbasisc(isb)                                 2d1s23
            call dgemm('n','n',nbas4,nbasdws(isb),nbasisc(isb),1d0,     2d1s23
     $         bc(igrabfrom),nbas4,bc(ivecold),nbasisc(isb),0d0,        2d1s23
     $         bc(itmpx),nbas4)
            if(idwsdeb.gt.5)then
             write(6,*)('mo guess vectors in ao basis')
             call prntm2(bc(itmpx),nbas4,nbasdws(isb),nbas4)
            end if
            ixh=ibcoff
            ibcoff=ixh+nbas4*nbasdws(isb)
            call dgemm('n','n',nbas4,nbasdws(isb),nbas4,1d0,
     $         bc(iovr(isb)),nbas4,bc(itmpx),nbas4,0d0,
     $         bc(ixh),nbas4)
            call dgemm('t','n',nbasdws(isb),nbasdws(isb),nbas4,1d0,
     $         bc(itmpx),nbas4,bc(ixh),nbas4,0d0,bc(ibcoff),
     $         nbasdws(isb))
            if(idwsdeb.gt.5)then
             write(6,*)('mo guess vectors in ao basis'),itmpx
             call prntm2(bc(itmpx),nbas4,nbasdws(isb),nbas4)
            end if
            ixr=ibcoff
            ibcoff=ixr+nbas4*nbasdws(isb)
            call enough('parahf.33a',bc,ibc)                            3d13s23
            call dgemm('n','n',nbas4,nbasdws(isb),nbas4,1d0,
     $         bc(iovr(isb)),nbas4,bc(ivecs(isb)),nbas4,0d0,
     $         bc(ixr),nbas4,'parahf.ixr')                              2d3s23
            ixrt=ibcoff
            ibcoff=ixrt+nbas4*nbasdws(isb)
            call enough('parahf.33b',bc,ibc)                            3d13s23
            do i1=0,nbasdws(isb)-1
             do i2=0,nbas4-1
              i12=ixrt+i1+nbasdws(isb)*i2
              i21=ixr+i2+nbas4*i1
              bc(i12)=bc(i21)
             end do
            end do
            call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbas4,1d0,
     $         bc(ixrt),nbasdws(isb),bc(itmpx),nbas4,0d0,bc(ibcoff),
     $         nbasdws(isb))
            if(idwsdeb.gt.5)then
             write(6,*)('re-generated guess vectors in mo basis ')
             call prntm2(bc(ibcoff),nbasdws(isb),nbasdws(isb),
     $           nbasdws(isb))
             write(6,*)('ircode '),ircode
            end if
            if(ircode.eq.4)then                                         2d3s23
             if(idwsdeb.gt.5)                                           5d5s23
     $           write(6,*)('new code parahf: overwrite vecold'),ivecold
             do ixi=0,nbasdws(isb)*nbasdws(isb)-1                        2d1s23
              bc(ivecold+ixi)=bc(ibcoff+ixi)                             2d1s23
             end do                                                      2d1s23
            end if                                                      2d3s23
            if(idwsdeb.gt.5)then                                        5d5s23
             write(6,*)('what is now in vecold: ')
             call prntm2(bc(ivecold),nbasdws(isb),nbasdws(isb),
     $           nbasdws(isb))
            end if                                                      5d5s23
            ihtmp=ibcoff                                                  12d15s22
            ibcoff=ihtmp+nbasdws(isb)                                     12d15s22
            itxcopy=ibcoff                                                12d15s22
            ibcoff=itxcopy+nbas4*nbasdws(isb)                             12d15s22
            call enough('parahf.29.9',bc,ibc)                           3d13s23
            do ixi=0,nbas4*nbasdws(isb)-1                                 12d15s22
             bc(itxcopy+ixi)=bc(itmpx+ixi)                                12d15s22
            end do                                                        12d15s22
            call assqn(bc(itmpx),nbasisp,nbasdws,isb,nlzz,                12d15s22
     $             ibc(ihtmp),natom,ngaus,ibdat,isym,iapair,ibstor,     12d8s22
     $             isstor,idorel,ascale,multh,nsymb,iorbsym,iorbsymz,bc,12d8s22
     $             ibc,epssym)                                          9d19s23
            itprod=ibcoff
            ibcoff=itprod+nbas4*nbasdws(isb)
            call enough('parahf.29.8',bc,ibc)                           3d13s23
            call dgemm('n','n',nbas4,nbasdws(isb),nbas4,1d0,              12d15s22
     $         bc(iovr(isb)),nbas4,bc(itmpx),nbas4,0d0,                 12d15s22
     $         bc(itprod),nbas4)                                        12d15s22
            itprodt=ibcoff                                                12d15s22
            ibcoff=itprodt+nbas4*nbasdws(isb)                             12d15s22
            call enough('parahf.29.7',bc,ibc)                           3d13s23
            do i1=0,nbasdws(isb)-1                                        12d15s22
             do i2=0,nbas4-1                                              12d15s22
              i12=itprod+i2+nbas4*i1                                      12d15s22
              i21=itprodt+i1+nbasdws(isb)*i2                              12d15s22
              bc(i21)=bc(i12)                                             12d15s22
             end do                                                       12d15s22
            end do                                                        12d15s22
            call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbas4,1d0,       12d15s22
     $         bc(itprodt),nbasdws(isb),bc(itxcopy),nbas4,0d0,          12d15s22
     $         bc(itprod),nbasdws(isb))                                 12d15s22
            do inew=0,nbasdws(isb)-1
             do iold=0,nbasdws(isb)-1
              iad=itprod+inew+nbasdws(isb)*iold                           12d15s22
              if(abs(bc(iad)-1d0).lt.1d-8)then                            12d15s22
               iadn=itprodt+inew*nbasdws(isb)                             12d15s22
               iado=ivecold+iold*nbasdws(isb)                             12d15s22
               do ix=0,nbasdws(isb)-1                                     12d15s22
                bc(iadn+ix)=bc(iado+ix)                                   12d15s22
               end do                                                     12d15s22
               go to 466                                                  12d15s22
              end if                                                      12d15s22
             end do                                                       12d15s22
             write(6,*)('could not match '),inew
             call dws_synca
             call dws_finalize
             stop
  466        continue
            end do                                                        12d15s22
            do ixi=0,nbasdws(isb)*nbasdws(isb)-1                        2d1s23
             bc(ivecold+ixi)=bc(itprodt+ixi)                              12d15s22
            end do                                                        12d15s22
            iassgn=1                                                      12d15s22
            ibcoff=ihtmp                                                  12d15s22
           end if                                                       12d15s22
          else                                                          12d15s22
           if(idwsdeb.gt.5)write(6,*)('filling jorb with unit matrix ')
           do i=0,nbasdws(isb)-1                                          2d15s19
            iad=jorb+nbasdws(isb)*i                                       2d15s19
            do j=0,nbasdws(isb)-1                                         2d15s19
             bc(iad+j)=0d0                                                1d30s18
            end do                                                        1d30s18
            bc(iad+i)=1d0                                                 1d30s18
           end do                                                         1d30s18
          end if                                                        12d15s22
         end if                                                         10d27s20
          if(idorel.ne.0.and.                                           10d11s22
     $         (myguess.or.(icas.ne.0.and.icasguess.ne.1)))then         10d11s22
           if(ircode.ne.4)then                                          10d12s22
            if(nvguess.gt.0)then                                        12d15s22
             icodeao=0                                                    10d12s22
             if(idwsdeb.gt.5)then
              write(6,*)('vecold is in ob basis')
              call prntm2(bc(ivecold),nbasisc(isb),nbasdws(isb),          12d15s22
     $          nbasisc(isb))                                           12d15s22
             end if                                                     5d5s23
             itmpt=ibcoff                                                 10d11s22
             itmx=itmpt+nbasisc(isb)*nbasdws(isb)                                10d11s22
             ibcoff=itmx+nbasdws(isb)*nbasdws(isb)                                      10d11s22
             call enough('parahf. 32',bc,ibc)                           3d13s23
             if(idwsdeb.gt.5)write(6,*)('filling jorb a')               5d5s23
             do i=0,nbasdws(isb)-1                                               10d11s22
              do j=0,nbasisc(isb)-1                                       10d11s22
               ji=jorb+j+nbasisc(isb)*i                                   10d11s22
               ij=itmpt+i+nbasdws(isb)*j                                12d15s22
               bc(ij)=bc(ji)                                              10d11s22
              end do                                                      10d11s22
             end do                                                       10d11s22
             call dgemm('n','n',nbasdws(isb),nbasdws(isb),              12d15s22
     $            nbasisc(isb),1d0,bc(itmpt),nbasdws(isb),bc(ivecold),  12d15s22
     $            nbasisc(isb),0d0,bc(itmx),nbasdws(isb),               12d15s22
     d'parahf.  4')
             do i=0,nbasdws(isb)-1                                               10d11s22
              iv=itmx+nbasdws(isb)*i                                             10d11s22
              do j=0,i-1                                                  10d11s22
               dot=0d0                                                    10d11s22
               jv=itmx+nbasdws(isb)*j                                            10d11s22
               do k=0,nbasdws(isb)-1                                             10d11s22
                dot=dot+bc(iv+k)*bc(jv+k)                                 10d11s22
               end do                                                     10d11s22
               do k=0,nbasdws(isb)-1                                             10d11s22
                bc(iv+k)=bc(iv+k)-dot*bc(jv+k)                            10d11s22
               end do                                                     10d11s22
              end do                                                      10d11s22
              dot=0d0                                                     10d11s22
              do k=0,nbasdws(isb)-1                                              10d11s22
               dot=dot+bc(iv+k)**2                                        10d11s22
              end do                                                      10d11s22
              dot=1d0/sqrt(dot)                                           10d11s22
              do k=0,nbasdws(isb)-1                                              10d11s22
               bc(iv+k)=bc(iv+k)*dot                                      10d11s22
              end do                                                      10d11s22
             end do                                                       10d11s22
             do i=0,nbasdws(isb)*nbasdws(isb)-1                                         10d11s22
              bc(ivecold+i)=bc(itmx+i)                                    10d11s22
             end do                                                       10d11s22
             if(idwsdeb.gt.5)then
              write(6,*)('re-write vecold with itmx')
              write(6,*)('guess vectors ')
              call prntm2(bc(ivecold),nbasdws(isb),nbasdws(isb),
     $           nbasdws(isb))
             end if
             ibcoff=itmpt                                                 10d11s22
            else if(nvguess.lt.0)then                                   12d15s22
             icodeao=1                                                    10d12s22
             if(idwsdeb.gt.5)then
              write(6,*)('nvguess lt 0 block ')
              write(6,*)('vecold is in ao basis'),nbas3                         10d11s22
              call prntm2(bc(ivecold),nbasisp(isb)*ncomp,nbasdws(isb),   2d2s23
     $          nbasisp(isb)*ncomp)                                     10d12s22
             end if
             if(nlzz.ne.0)then                                          12d7s22
              ihtmp=ibcoff                                              12d8s22
              ibcoff=ihtmp+nbasisc(isb)                                 12d8s22
              call enough('parahf.29.1',bc,ibc)                         3d13s23
              call assqn(bc(ivecold),nbasisp,nbasdws,isb,nlzz,          2d2s23
     $             ibc(ihtmp),natom,ngaus,ibdat,isym,iapair,ibstor,     12d8s22
     $             isstor,idorel,ascale,multh,nsymb,iorbsym,iorbsymz,bc,12d8s22
     $             ibc,epssym)                                          9d19s23
              iassgn=1                                                  12d15s22
              ibcoff=ihtmp                                              12d8s22
             end if                                                     12d7s22
             call contractg(bc(ivecold),ngaus,ibdat,ibstor,isstor,      10d12s22
     $            idorel,isb,nbasisp,nbasisc,iapair,nlzz,iorbsym,       10d12s22
     $            iorbsymz,nbasdws,1,bc,ibc)                            2d2s23
             itmx=ibcoff                                                 10d11s22
             ibcoff=itmx+nbas3*nbasisc(isb)                              10d11s22
             call enough('parahf. 30',bc,ibc)                           3d13s23
             call dgemm('n','n',nbas3,nbasdws(isb),nbas3,1d0,           2d2s23
     $           bc(ixinv(isb)),nbas3,bc(ivecold),nbas3,0d0,            10d12s22
     $           bc(itmx),nbas3,                                        10d11s22
     d'parahf.  3')
             do i=0,nbas3*nbasdws(isb)-1                                2d2s23
              bc(ivecold+i)=bc(itmx+i)                                   10d11s22
             end do                                                      10d11s22
             if(idwsdeb.gt.5)then
              write(6,*)('another overwrite vecold with tmpx')
              write(6,*)('in ob basis ')
              call prntm2(bc(ivecold),nbas3,nbasdws(isb),nbas3)          2d2s23
             end if                                                     5d5s23
             if(nposs(isb).ne.0)then                                    5d2s23
              igetfrom=ivecs(isb)                                        2d3s23
     $         +nbasisp(isb)*ncomp*(nbasisc(isb)*2+nbasisp(isb)*ncomp)  2d1s23
              itrnsx=ibcoff                                              2d3s23
              ibcoff=itrnsx+nbas3*nbasdws(isb)                           2d3s23
              call enough('parahf.trnsx',bc,ibc)                        3d13s23
              do i1=0,nbas3-1                                            2d3s23
               do i2=0,nbasdws(isb)-1                                    2d3s23
                i12=igetfrom+i1+nbas3*i2                                 2d3s23
                i21=itrnsx+i2+nbasdws(isb)*i1                            2d3s23
                bc(i21)=bc(i12)                                          2d3s23
               end do                                                    2d3s23
              end do                                                     2d3s23
              imxmx=ibcoff                                               2d3s23
              ibcoff=imxmx+nbasdws(isb)*nbasdws(isb)                     2d3s23
              call enough('parahf.mxmx',bc,ibc)                         3d13s23
              call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbas3,1d0,    2d3s23
     $            bc(itrnsx),nbasdws(isb),bc(ivecold),nbas3,0d0,        2d3s23
     $            bc(imxmx),nbasdws(isb),'parahf.mxmx')                 2d3s23
              if(idwsdeb.gt.5)then
               write(6,*)('vectors in h0 electronic basis: ')             2d3s23
               call prntm2(bc(imxmx),nbasdws(isb),nbasdws(isb),           2d3s23
     $            nbasdws(isb))                                         2d3s23
               write(6,*)('copy this to vecold ')
              end if
              do i12=0,nbasdws(isb)*nbasdws(isb)-1                       2d3s23
               bc(ivecold+i12)=bc(imxmx+i12)                             2d3s23
              end do                                                     2d3s23
c
c     modified graham-schmidt orthogonalization ...
c
              if(idwsdeb.gt.5)
     $             write(6,*)('orthogonalizing vecold in mo basis ')
              do i1=0,nbasdws(isb)-1                                     2d3s23
               iv1=ivecold+nbasdws(isb)*i1                               2d3s23
               do i2=0,i1-1                                              2d3s23
                iv2=ivecold+nbasdws(isb)*i2                              2d3s23
                dot=0d0                                                  2d3s23
                do i3=0,nbasdws(isb)-1                                   2d3s23
                 dot=dot+bc(iv1+i3)*bc(iv2+i3)                           2d3s23
                end do                                                   2d3s23
                do i3=0,nbasdws(isb)-1                                   2d3s23
                 bc(iv1+i3)=bc(iv1+i3)-dot*bc(iv2+i3)                    2d3s23
                end do                                                   2d3s23
               end do                                                    2d3s23
               dot=0d0                                                   2d3s23
               do i3=0,nbasdws(isb)-1                                    2d3s23
                dot=dot+bc(iv1+i3)**2                                    2d3s23
               end do                                                    2d3s23
               doti=1d0/sqrt(dot)                                        2d3s23
               do i3=0,nbasdws(isb)-1                                    2d3s23
                bc(iv1+i3)=bc(iv1+i3)*doti                               2d3s23
               end do                                                    2d3s23
              end do                                                     2d3s23
              if(idwsdeb.gt.5)then                                      5d5s23
               call dgemm('t','n',nbasdws(isb),nbasdws(isb),
     $             nbasdws(isb),
     $            1d0,bc(ivecold),nbasdws(isb),bc(ivecold),nbasdws(isb),
     $            0d0,bc(imxmx),nbasdws(isb),'parahf.otest')
               write(6,*)('ortho test ')
               call prntm2(bc(imxmx),nbasdws(isb),nbasdws(isb),
     $            nbasdws(isb))
              end if                                                    5d5s23
              ibcoff=itrnsx                                              2d3s23
             end if                                                     2d6s23
            end if                                                      10d12s22
           else                                                         4d27s23
            if(iabs(myguessi).eq.1.and..not.(nlzz.eq.6.or.idorel.eq.0)) 4d27s23
     $           then                                                   4d27s23
             if(idwsdeb.gt.5)                                           5d5s23
     $       write(6,*)('need to turn vecold from ao to ob vectors ...')4d27s23
             nrow=nbasisp(isb)*ncomp                                    4d27s23
             jvecs=ivecs(isb)+nrow*nbasisc(isb)                         4d27s23
             itmp1=ibcoff                                               4d27s23
             itmp2=itmp1+nrow*nbasdws(isb)                              4d27s23
             itmp3=itmp2+nrow*nbasdws(isb)                              4d27s23
             ibcoff=itmp3+nbasdws(isb)*nbasdws(isb)                     4d27s23
             call enough('parahf.tmp1',bc,ibc)                          4d27s23
             call dgemm('n','n',nrow,nbasdws(isb),nrow,1d0,bc(jvecs),   4d27s23
     $            nrow,bc(ivecs(isb)),nrow,0d0,bc(itmp1),nrow,          4d27s23
     $            'parahf.tmp1')                                        4d27s23
             do i=0,nrow-1                                              4d27s23
              do j=0,nbasdws(isb)-1                                     4d27s23
               ji=itmp2+j+nbasdws(isb)*i                                4d27s23
               ij=itmp1+i+nrow*j                                        4d27s23
               bc(ji)=bc(ij)                                            4d27s23
              end do                                                    4d27s23
             end do                                                     4d27s23
             call dgemm('n','n',nbasdws(isb),nbasdws(isb),nrow,1d0,     4d27s23
     $            bc(itmp2),nbasdws(isb),bc(ivecold),nrow,0d0,          4d27s23
     $            bc(itmp3),nbasdws(isb),'parahf.tmp3')                 4d27s23
             if(idwsdeb.gt.5)then
              write(6,*)('in ob')
              call prntm2(bc(itmp3),nbasdws(isb),nbasdws(isb),
     $            nbasdws(isb))
             end if
             do i=0,nbasdws(isb)*nbasdws(isb)-1                         4d27s23
              bc(ivecold+i)=bc(itmp3+i)                                 4d27s23
             end do                                                     4d27s23
             ibcoff=itmp1                                               4d27s23
            end if                                                      4d27s23
           end if                                                       10d12s22
          end if                                                        10d11s22
        end if                                                          8d21s15
        call dws_bcast(bc(ieighs(isb)),nbasdws(isb))                    2d15s19
        if(idorel.ne.0)then                                             1d30s18
         nwds=nbas4*nbasdws(isb)                                        2d15s19
         call dws_bcast(bc(ivecs(isb)),nwds)                            1d30s18
         do i=0,nbasdws(isb)-1                                          1d30s18
          do j=0,nbas4-1                                                1d30s18
           iad1=ivecs(isb)+j+nbas4*i                                    1d30s18
           iad2=ivect(isb)+i+nbasdws(isb)*j                             1d30s18
           bc(iad2)=bc(iad1)                                            1d30s18
          end do                                                        1d30s18
         end do                                                         1d30s18
        else                                                            1d30s18
         nwds=nbas3**2                                                   8d20s15
         call dws_bcast(bc(jorb),nwds)                                   4d22s14
        end if                                                          1d30s18
        if(myguess.or.(icas.ne.0.and.icasguess.ne.1))then               8d8s14
         if(lprint)then                                                  10d13s04
          write(6,*)('we''ve got guess from previous orbitals ')          9d20s04
          writE(6,*)('over write h0 eigenvectors with these ')            9d20s04
         end if                                                          10d13s04
         call enough('parahf. 34',bc,ibc)                               3d13s23
c
c     splush is the transformation for ao mos to ob mos
c
         if((myguessi.eq.1.and..not.(idorel.ne.0.and.myguess.or.         10d13s22
     $        (icas.ne.0.and.icasguess.ne.1))).or.nrydb.ne.0)then       11d2s22
          nrow=nbasisp(isb)*ncomp                                       2d15s19
          if(nlzz.eq.6.or.idorel.eq.0)then                              12d8s22
           ism=ivecs(isb)+nrow*nbasisc(isb)                              2d15s19
           itmp2=ibcoff                                                 5d5s23
            itmp1=ibcoff
            itmp2=itmp1+nrow*nbasdws(isb)
            ibcoff=itmp2+nrow*nbasdws(isb)
            call enough('parahf. 35',bc,ibc)                             3d13s23
            if(idwsdeb.gt.5)then                                         5d5s23
             write(6,*)('to get transformation from mo to ao basis ')
             write(6,*)('multiply ism ')
             call prntm2(bc(ism),nrow,nrow,nrow)
             write(6,*)('times vecs ')
             call prntm2(bc(ivecs(isb)),nrow,nbasdws(isb),nrow)
            end if                                                       5d5s23
            call dgemm('n','n',nrow,nbasdws(isb),nrow,1d0,bc(ism),nrow,   2d15s19
     $         bc(ivecs(isb)),nrow,0d0,bc(itmp1),nrow,                  2d15s19
     d' parahf.  7')
            if(idwsdeb.gt.5)then                                         5d5s23
             write(6,*)('to yield ')
             call prntm2(bc(itmp1),nrow,nbasdws(isb),nrow)
            end if                                                       5d5s23
            do i=0,nrow-1                                                 2d15s19
             do j=0,nbasdws(isb)-1                                        2d15s19
              ji=itmp2+j+nbasdws(isb)*i                                   2d15s19
              ij=itmp1+i+nrow*j                                           2d15s19
              bc(ji)=bc(ij)                                               2d15s19
             end do                                                       2d15s19
            end do                                                        2d15s19
            if(idwsdeb.gt.5)then
             write(6,*)
     $         ('just for kicks, multiply transpose of transformation ')
             write(6,*)('by vecs. this should be the unit matrix ')
             call dgemm('n','n',nbasdws(isb),nbasdws(isb),nrow,1d0,
     $           bc(itmp2),nbasdws(isb),bc(ivecs(isb)),nrow,0d0,
     $           bc(ibcoff),nbasdws(isb),'parahf.otest')
             call prntm2(bc(ibcoff),nbasdws(isb),nbasdws(isb),
     $           nbasdws(isb))
            end if                                                      5d5s23
           if(idwsdeb.gt.5)then                                         5d5s23
            write(6,*)('sort out sym the first '),mxsb,loc(mxsb)
            write(6,*)('vecold going into sortoutsym ')
            call prntm2(bc(ivecold),nrow,nbasdws(isb),nrow)             5d5s23
            itmpvo=ibcoff                                               5d5s23
            itmpov=itmpvo+nrow*nbasdws(isb)
            ibcoff=itmpov+nrow*nbasdws(isb)
            call enough('parahf.tmpovvo',bc,ibc)                        5d5s23
            call dgemm('n','n',nrow,nbasdws(isb),nrow,1d0,bc(ism),nrow,
     $           bc(ivecold),nrow,0d0,bc(itmpvo),nrow,'parahf.tmpvo')   5d5s23
            do i=0,nrow-1                                               5d5s23
             do j=0,nbasdws(isb)-1                                      5d5s23
              ji=itmpov+j+nbasdws(isb)*i                                5d5s23
              ij=itmpvo+i+nrow*j                                        5d5s23
              bc(ji)=bc(ij)                                             5d5s23
             end do                                                     5d5s23
            end do                                                      5d5s23
            call dgemm('n','n',nbasdws(isb),nbasdws(isb),nrow,1d0,      5d5s23
     $           bc(itmpov),nbasdws(isb),                               5d5s23
     $           bc(ivecold),nrow,0d0,bc(itmpvo),nbasdws(isb),          5d5s23
     $           'parahf.tmpov')                                        5d5s23
            write(6,*)('ortho test')
            call prntm2(bc(itmpvo),nbasdws(isb),nbasdws(isb),
     $           nbasdws(isb))
            ibcoff=itmpvo                                               5d5s23
           end if
           if(nlzz.ne.0)then                                            7d5s23
            call sortoutsym(nlzz,iorbsymao(isb),                          5d5s21
     $         bc(ivecs(isb)),bc(ivecold),nrow,nbasdws(isb),1,          6d7s21
     $         ibc(ixsb),nxsb,mxsb,natom,ngaus,ibdat,isym,iapair,       12d4s22
     $           ibstor,isstor,idorel,ascale,multh,nbasisp,isb,nsymb,bc,12d4s22
     $         ibc)                                                     12d4s22
           end if                                                       7d5s23
           if(idwsdeb.gt.5)then                                         5d5s23
            write(6,*)('back ')
            write(6,*)('vecold coming out of sortoutsym ')
            call prntm2(bc(ivecold),nrow,nbasdws(isb),nrow)             5d5s23
            write(6,*)('transformation from mo to ao basis: ')
            call prntm2(bc(itmp2),nbasdws(isb),nrow,nbasdws(isb))
           end if
            call dgemm('n','n',nbasdws(isb),nbasdws(isb),nrow,1d0,        2d15s19
     $         bc(itmp2),nbasdws(isb),bc(ivecold),nrow,0d0,bc(jorb),    2d15s19
     $         nbasdws(isb),                                            2d15s19
     d' parahf.  8')
           if(idwsdeb.gt.5)then                                         5d5s23
c     I=VtSV, Q=S*x, xT*Q=QT*x=xT*S*x=I, u=QTV=xT*S*V
c     uTu=VT*Q*xT*S*V
c     if uT*xT*S*x*u=uT*u=I and x*u=V, xT*S*x*u=xT*S*V=u
            write(6,*)('thus jorb is ')
            call prntm2(bc(jorb),nbasdws(isb),nbasdws(isb),nbasdws(isb))
            call dgemm('t','n',nbasdws(isb),nbasdws(isb),nbasdws(isb),
     $           1d0,bc(jorb),nbasdws(isb),bc(jorb),nbasdws(isb),0d0,
     $           bc(ibcoff),nbasdws(isb),'parahf.jorb.orth')
            write(6,*)('ortho test: ')
            call prntm2(bc(ibcoff),nbasdws(isb),nbasdws(isb),
     $           nbasdws(isb))
            write(6,*)('vecs*jorb should be vecold? ')
            call dgemm('n','n',nrow,nbasdws(isb),nbasdws(isb),1d0,
     $           bc(ivecs(isb)),nrow,bc(jorb),nbasdws(isb),0d0,
     $           bc(ibcoff),nrow,'parahf.vecold?')                      5d5s23
            call prntm2(bc(ibcoff),nrow,nbasdws(isb),nrow)
           end if                                                       5d5s23
           ibcoff=itmp1                                                 2d6s23
          else                                                          12d8s22
           do i=0,nbasdws(isb)*nbasdws(isb)-1                           12d8s22
            bc(jorb+i)=bc(ivecold+i)                                    12d8s22
           end do                                                       12d8s22
          end if                                                        12d8s22
          if(idorel.ne.0)then                                           6d1s21
c
c     if guess ao vectors do not span the same electronic space
c     as h0 electronic space, jorb will not be orthogonal.
c     so orthogonalize.
c
           do i=0,nbasdws(isb)-1                                          6d1s21
            ii=jorb+i*nbasdws(isb)                                        6d1s21
            do j=0,i-1                                                    6d1s21
             jj=jorb+j*nbasdws(isb)                                       6d1s21
             dot=0d0                                                      6d1s21
             do k=0,nbasdws(isb)-1                                        6d1s21
              dot=dot+bc(jj+k)*bc(ii+k)                                   6d1s21
             end do                                                       6d1s21
             do k=0,nbasdws(isb)-1                                        6d1s21
              bc(ii+k)=bc(ii+k)-bc(jj+k)*dot                              6d1s21
             end do                                                       6d1s21
            end do                                                        6d1s21
            dot=0d0                                                       6d1s21
            do k=0,nbasdws(isb)-1                                         6d1s21
             dot=dot+bc(ii+k)**2                                          6d1s21
            end do                                                        6d1s21
            dot=1d0/sqrt(dot)                                             6d1s21
            do k=0,nbasdws(isb)-1                                         6d1s21
             bc(ii+k)=bc(ii+k)*dot                                        6d1s21
            end do                                                        6d1s21
           end do
          end if                                                        6d1s21
          if(ircode.eq.4)then                                           2d22s23
           ivecold=ivecold+nrow*nbasdws(isb)                            2d22s23
          else                                                          2d22s23
           ivecold=ivecold+nrow*nbasisc(isb)                             2d2s23
          end if                                                        2d22s23
         else
          if(nlzz.ne.0)then                                             4d22s21
           if(idorel.eq.0.or.ismile(isb).ne.0)then                      4d23s21
            if(nlzz.eq.2)then                                           12d8s22
             ivtmp=ibcoff                                               12d8s22
             ibcoff=ivtmp+nbasisp(isb)*ncomp*nbasisc(isb)               12d8s22
             call enough('parahf.35.1',bc,ibc)                          3d13s23
             call dgemm('n','n',nbasisp(isb)*ncomp,nbasdws(isb),        5d4s23
     $            nbasdws(isb),1d0,bc(ivecs(isb)),nbasisp(isb)*ncomp,   5d4s23
     $            bc(ivecold),nbasdws(isb),0d0,bc(ivtmp),               5d4s23
     $            nbasisp(isb)*ncomp,                                   12d8s22
     d'parahf.8.1')                                                     12d8s22
             ihtmp=ibcoff                                               12d8s22
             ibcoff=ihtmp+nbasisc(isb)                                  12d8s22
             call enough('parahf.35.2',bc,ibc)                          3d13s23
             call assqn(bc(ivtmp),nbasisp,nbasdws,isb,nlzz,ibc(ihtmp),  5d4s23
     $            natom,ngaus,ibdat,isym,iapair,ibstor,isstor,idorel,   12d8s22
     $            ascale,multh,nsymb,iorbsym,iorbsymz,bc,ibc,epssym)    9d19s23
             jvtmp=ivtmp                                                12d8s22
             do i=0,nbasdws(isb)-1                                      5d4s23
              iad=ivecold+nbasdws(isb)*(ibc(ihtmp+i)-1)                 5d4s23
              do j=0,nbasdws(isb)-1                                     5d4s23
               bc(jvtmp+j)=bc(iad+j)                                    12d8s22
              end do                                                    12d8s22
              jvtmp=jvtmp+nbasdws(isb)                                  5d4s23
             end do                                                     12d8s22
             do i=0,nbasdws(isb)*nbasdws(isb)-1                         5d4s23
              bc(ivecold+i)=bc(ivtmp+i)                                 12d8s22
             end do                                                     12d8s22
             ibcoff=ivtmp                                               12d8s22
            else                                                        12d8s22
             call sortoutsym(nlzz,iorbsymao(isb),                        4d23s21
     $          bc(ivecs(isb)),bc(ivecold),nbas4,nbasdws(isb),0,        6d7s21
     $          ibc(ixsb),nxsb,mxsb,natom,ngaus,ibdat,isym,iapair,      12d4s22
     $           ibstor,isstor,idorel,ascale,multh,nbasisp,isb,nsymb,bc,12d4s22
     $           ibc)                                                   12d4s22
            end if                                                      12d8s22
           else                                                         4d23s21
            if(iassgn.eq.0)then                                         12d15s22
             call sortoutsym(nlzz,iorbsymao(isb),                        4d23s21
     $          bc(ivecs(isb)),bc(ivecold),nbasisc(isb),nbasdws(isb),0, 6d7s21
     $           ibc(ixsb),nxsb,mxsb,natom,ngaus,ibdat,isym,iapair,     12d4s22
     $           ibstor,isstor,idorel,ascale,multh,nbasisp,isb,nsymb,bc,12d4s22
     $           ibc)                                                   12d4s22
            end if                                                      12d15s22
           end if                                                       4d23s21
          end if                                                        4d22s21
          do i=0,nbasdws(isb)*nbasdws(isb)-1                             2d15s19
           bc(jorb+i)=bc(ivecold+i)
          end do
          if(nvguess.lt.0)then                                          10d12s22
           ivecold=ivecold+iabs(nvguess)*nbasisp(isb)*nbasdws(isb)      10d12s22
          else                                                          10d12s22
           if(nbasisc(isb).eq.nbasisp(isb)*2)then                       2d2s23
           ivecold=ivecold+iabs(nvguess)*(nbasdws(isb)**2)               10d11s22
           else                                                         2d2s23
            ivecold=ivecold+iabs(nvguess)*nbasdws(isb)*nbasisc(isb)      2d1s23
           end if                                                       2d2s23
          end if                                                        10d12s22
         end if
        end if                                                           9d20s04
        ibcoff=idwst3                                                    9d9s04
        jorb=jorb+nbasdws(isb)*nbas4                                    2d15s19
        if(idwsdeb.gt.10)then
         write(6,*)('eigenvalues of h0 '),nbas3
         iarg1=nbasdws(isb)                                             2d15s19
         call prntm2(bc(ieighs(isb)),i18,iarg1,i18)                     1d13s16
        end if
        if(noc(1).lt.0)then
         do idws=1,min(nbasdws(isb),not)                                2d15s19
          bc(jeigh)=bc(ieighs(isb)+idws-1)                              1d13s16
          ibc(jsymh)=isb
          jeigh=jeigh+1                                                  6d28s04
          jsymh=jsymh+1                                                  6d28s04
         end do                                                          6d28s04
        end if                                                           6d28s04
        ioff=ioff+nbas4*nbas4                                           8d20s15
       end if
      end do
      if(lprt)write(6,*)('stamp1')
      if(noc(1).lt.0)then                                               11d25s08
       if(lprint)write(6,*)('sorting h0 eigenvalues across blocks '),
     $     potdws
       iarg1=jeigh-ieigh
       isort=ibcoff                                                     6d28s04
       ibcoff=isort+iarg1                                               6d28s04
       call dsortdws(bc(ieigh),ibc(isort),iarg1)                        1d18s23
       do idws=1,nsymb                                                  6d28s04
        noc(idws)=0                                                     6d28s04
       end do                                                           6d28s04
       do idws=1,not                                                    6d28s04
        isb=ibc(isymh+ibc(isort+idws-1)-1)                              6d28s04
        if(lprint)write(6,6243)bc(ieigh+idws-1),isb                     10d13s04
 6243   format(1pe15.7,i5)                                              6d28s04
        noc(isb)=noc(isb)+1                                             6d28s04
       end do                                                           6d28s04
       ibcoff=isort                                                     8d1s12
c
c     now make this consistent across processors ...
c
       ixnoc=ibcoff                                                     4d16s14
       ibcoff=ixnoc+nsymb                                               4d16s14
       call enough('parahf. 36',bc,ibc)                                 3d13s23
       if(mynowprog.eq.0)then                                           4d16s14
        do idws=1,nsymb                                                 4d16s14
         bc(ixnoc+idws-1)=dfloat(noc(idws))+1d-3                        4d16s14
        end do                                                          4d16s14
       end if                                                           4d16s14
       call dws_bcast(bc(ixnoc),nsymb)                                  4d16s14
       if(mynowprog.ne.0)then                                           4d16s14
        do idws=1,nsymb                                                 4d16s14
         nocxy=bc(ixnoc+idws-1)                                         4d16s14
         if(nocxy.ne.noc(idws))then                                     4d16s14
          write(6,*)('my noc differs from proc 0 '),noc(idws),nocxy,    4d16s14
     $        (' resetting ...')                                        4d16s14
          noc(idws)=nocxy                                               4d16s14
         end if                                                         4d16s14
        end do                                                          4d16s14
       end if                                                           4d16s14
       ibcoff=ixnoc                                                     4d16s14
       do idws=1,nsymb                                                  8d20s14
        idbly(idws)=noc(idws)                                           8d20s14
       end do                                                           8d20s14
      end if                                                            6d28s04
      if(lprt)write(6,*)('stamp2')
      jorb=iorb                                                         6d28s04
      nsy2(1)=0                                                         6d28s04
      nsy3(1)=0
      do isb=1,nsymb                                                    6d28s04
       nvirt(isb)=nbasdws(isb)-noc(isb)                                 1d5s12
       if(isb.gt.1)then                                                 6d28s04
        nsy2(isb)=nsy2(isb-1)+nbasdws(isb-1)*noc(isb-1)                 6d28s04
        nsy3(isb)=nsy3(isb-1)+nbasdws(isb-1)                            6d28s04
       end if                                                           6d28s04
       if(noc(isb).gt.0.and.(idwsdeb.gt.5.or.icas.lt.0))then            8d8s14
        write(6,*)('initial orbitals in orth. basis for symmetry '),
     $      isb,jorb
        iarg1=noc(isb)                                                  6d28s04
        call prntm2(bc(jorb),nbasdws(isb),iarg1,nbasdws(isb))           1d30s18
       end if                                                           6d28s04
       iorbx(isb)=jorb                                                  6d28s04
       jorb=jorb+nbasdws(isb)*ncomp*nbasisp(isb)                        2d15s19
      end do                                                            6d28s04
      call second(time2)                                                10d19s04
      telap=time2-time1-tovr
      times(2,1)=times(2,1)+telap
      times(2,2)=times(2,2)+telap**2
      times(2,3)=times(2,3)+1d0
      ispaaaa=ibcoff
      ispaabb=ispaaaa+nsymb
      ispaabb=ispaabb+nsymb*nsymb
      ispabcd=ispaabb+nsymb*nsymb
      ibcoff=ispabcd+nsymb*nsymb
      do i=1,nsymb*nsymb*3+nsymb
       ibc(ispaaaa+i-1)=0
      end do
      call make1dm(nsymb,noc,nbasdws,lprint,maxddi,jmats,kmats,jmats,   1d30s18
     $     jdenpt,isblk,isblkk,isblkh,nsdlk,nsdlkk,nsdlkh,multh,idbk,   3d9s12
     $     isblk1,nsdlk1,isblkd,nsdlkd,0)                               7d5s21
      if(icas.ne.0)then                                                 8d20s14
       do idws=1,nsdlk                                                  8d20s14
        j2den(idws)=0                                                   11d6s17
        n1=isblk(1,idws)                                                8d20s14
        n2=isblk(2,idws)                                                8d20s14
        n3=isblk(3,idws)                                                8d20s14
        n4=isblk(4,idws)                                                8d20s14
        if(n1.eq.n2)then                                                8d20s14
         nnn=(iacto(n1)*(iacto(n1)+1))/2                                8d20s14
        else                                                            8d20s14
         nnn=iacto(n1)*iacto(n2)                                        8d20s14
        end if                                                          8d20s14
        if(n3.eq.n4)then                                                8d20s14
         mmm=(iacto(n3)*(iacto(n3)+1))/2                                8d20s14
        else                                                            8d20s14
         mmm=iacto(n3)*iacto(n4)                                        8d20s14
        end if                                                          8d20s14
        if(nnn*mmm.gt.0)then                                            8d20s14
         j2den(idws)=ibcoff                                             8d20s14
         ibcoff=j2den(idws)+nnn*mmm                                     8d20s14
         call enough('parahf. 37',bc,ibc)                               3d13s23
        end if                                                          8d20s14
       end do                                                           8d20s14
      end if                                                            8d20s14
      npos=0                                                            5d2s23
      do i=1,nsymb
       do j=1,nsymb
        iij=multh(i,j)
        do k=1,nsymb
         do l=1,nsymb
          ikl=multh(k,l)
          if(multh(iij,ikl).eq.1)then
           if(i.ge.j)then                                               1d31s06
            npos=npos+1                                                 5d2s23
            isblk4x(1,npos)=i                                           5d2s23
            isblk4x(2,npos)=j                                           5d2s23
            isblk4x(3,npos)=k                                           5d2s23
            isblk4x(4,npos)=l                                           5d2s23
           end if                                                       1d30s06
          end if
         end do
        end do
       end do
      end do
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
      n4x=npos                                                          5d2s23
      xlamb=1d0                                                         2d25s19
      iter1=0
      micitl=0                                                          1d18s20
      micit=0                                                           1d11s23
      change=0d0                                                        1d11s23
      ehfbsofarl=0d0                                                    1d11s23
      ehf=0d0                                                           1d11s23
      ehfo=1d0
      do isb=1,nsymb                                                    3d23s05
       iden1(isb)=ibcoff                                                6d28s04
       ibcoff=iden1(isb)+noc(isb)*noc(isb)                                 6d28s04
       do 58 i1=1,noc(isb)
        jden1=iden1(isb)-1+(i1-1)*noc(isb)
        do 59 i2=1,noc(isb)
         bc(jden1+i2)=0d0
   59   continue
        bc(jden1+i1)=2d0
   58  continue
       if(idwsdeb.gt.5)then
       writE(6,*)('density ')
       iarg1=noc(isb)
       call prntm2(bc(iden1(isb)),iarg1,iarg1,iarg1)
       end if
      end do                                                            3d23s05
   29 continue
      do isb=1,nsymb                                                    2d28s21
       ifockeig(isb)=ibcoff                                             2d28s21
       ibcoff=ifockeig(isb)+idbly(isb)                                  2d28s21
      end do                                                            2d28s21
      nbaspass=0
      ih0h=ibcoff
      ibcoff=ih0h+nsymt                                                 6d28s04
      ih0mo=ibcoff                                                      2d23s04
      ibcoff=ih0mo+nsymt                                                6d28s04
      call enough('parahf. 38',bc,ibc)                                  3d13s23
      call second(time1)                                                2d20s14
      telap=time1-time2-tovr
      times(14,1)=times(14,1)+telap
      times(14,2)=times(14,2)+telap**2
      times(14,3)=times(14,3)+1d0
      jackitup=0                                                        3d12s21
      xx20l=1d10                                                        3d12s21
      xx20o=1d10                                                        2d17s23
      xx1b=0d0                                                          2d1s23
      icol=ibcoff                                                       10d11s24
      ibcoff=icol+mynprocg                                              10d11s24
      imsg=ibcoff                                                       10d11s24
      ibcoff=imsg+mynnode                                               10d11s24
      memneed(1:4)='imsg'                                               10d11s24
      call enough('parahf. 40',bc,ibc)                                  10d11s24
      do isb=1,nsymb                                                    10d11s24
       nbas4=nbasdws(isb)                                               1d30s18
       nbas3=nbasisp(isb)*ncomp                                         4d30s18
       iorbn(isb)=ibcoff                                                10d11s24
       ibcoff=ibcoff+nbas3*nbas4                                        10d11s24
      end do                                                            10d11s24
      call enough('parahf.orbn',bc,ibc)                                 10d11s24
      norbparam=0                                                       10d22s24
      do isb=1,nsymb                                                    10d22s24
       norbparam=norbparam+noc(isb)*nvirt(isb)+idoub(isb)*iacto(isb)    10d22s24
      end do                                                            10d22s24
      if(lprint)write(6,*)('total no. orbital rotation parameters = '), 10d22s24
     $     norbparam                                                    10d22s24
      ibc404=ibcoff                                                     10d11s24
  404 continue
      if(idwsdeb.gt.10)write(6,*)('at 404 continue '),iorbx(1),ibcoff
      call second(time2)                                                2d20s14
      iter1=iter1+1
      if(iter1.eq.2.and.icas.ne.0)then                                  10d22s24
       do i=0,ncasvecsize-1                                             10d11s24
        bc(ibc404+i)=bc(icasvec+i)                                      10d11s24
       end do                                                           10d11s24
       icasvec=ibc404                                                   10d11s24
       ibc404=ibc404+ncasvecsize                                        10d11s24
      end if                                                            10d11s24
      ibcoff=ibc404                                                     10d11s24
      if(lprt)write(6,*)('starting iteration no.'),iter1
      do isb=1,nsymb                                                    6d28s04
       if(idwsdeb.gt.10.or.icas.lt.0)then
        write(6,*)('transforming orbs from orthogonal to ao basis')
        write(6,*)('nsdlk '),nsdlk
        write(6,*)('for block '),isb
       end if
       nbas4=nbasdws(isb)                                               1d30s18
       nbas3=nbasisp(isb)*ncomp                                         4d30s18
       if(nbas4.gt.0)then
        if(iter1.eq.1)then                                              9d9s04
         idwst1=iorbn(isb)                                              10d11s24
         if(idwsdeb.gt.5.or.icas.lt.0)then
          write(6,*)('S-1/2 '),ivecs(isb)
          call prntm2(bc(ivecs(isb)),nbas3,nbas3,nbas3)                 2d3s23
          write(6,*)('orbs in orthogonal basis '),iorbx(isb),ibcoff
          call prntm2(bc(iorbx(isb)),nbas4,nbas4,nbas4)                 1d30s18
         end if
         call dgemm('n','n',nbas3,nbas4,nbas4,1d0,bc(ivecs(isb)),nbas3,   6d28s04
     $           bc(iorbx(isb)),nbas4,0d0,bc(idwst1),nbas3,             6d28s04
     d' parahf.  9')
         if(idwsdeb.gt.5.or.icas.lt.0)then
          write(6,*)('orbs in ao basis '),idwst1,ibcoff,iorbx(isb)
          call prntm2(bc(idwst1),nbas3,nbas4,nbas3)                     1d30s18
         end if
         jorb=iorbx(isb)-1
         iorbn(isb)=idwst1                                               6d28s04
         jdwst1=idwst1-1                                                 6d28s04
         do i=1,nbas3
          ii=(i-1)*nbas4
          do j=1,nbas4
           jj=(j-1)*nbas3
           bc(jorb+j+ii)=bc(jdwst1+i+jj)
          end do                                                         6d28s04
         end do                                                          6d28s04
        else                                                            9d9s04
         if(nbasdws(isb).gt.0)then                                      3d22s12
          do i=1,nbas4                                                   9d9s04
           j1=iorbn(isb)-1+(i-1)*nbas3
           j2=iorbx(isb)+i-1-nbas4                                       9d9s04
           do j=1,nbas3                                                  9d9s04
            bc(j2+j*nbas4)=bc(j1+j)                                      9d9s04
           end do                                                        9d9s04
          end do                                                         9d9s04
         end if                                                         9d17s04
        end if                                                          9d9s04
       end if
      end do                                                            6d28s04
      if(idwsdeb.gt.10.or.icas.lt.0)then
       write(6,*)('transform h0 to mo basis ')
       write(6,*)('nsdlk '),nsdlk
      end if
      jh0mo=ih0mo                                                       5d8s18
      do isb=1,nsymb                                                    6d28s04
       nbas4=nbasdws(isb)                                               1d30s18
       nbas3=nbasisp(isb)*ncomp                                         5d8s18
       jdwsk=ih0+nsy1(isb)                                              6d28s04
       jh0h=ih0h+nsy1(isb)
       if(nbas4.gt.0)then
        if(idwsdeb.gt.5.or.icas.lt.0)then
         itryh0=jdwsk+7-1+nbas3*(15-1)
         write(6,*)('h0 for block '),isb,jdwsk,itryh0,bc(itryh0)
         call prntm2(bc(jdwsk),nbas3,nbas3,nbas3)
         write(6,*)('orbitals '),iorbn(isb),ibcoff,nsy1(isb)
         call prntm2(bc(iorbn(isb)),nbas3,nbasdws(isb),                 1d30s18
     $      nbas3)
        end if
        call dgemm('n','n',nbas3,nbas4,nbas3,1d0,bc(jdwsk),nbas3,        6d28s04
     $           bc(iorbn(isb)),nbas3,0d0,bc(jh0h),nbas3,               6d28s04
     d' parahf. 10')
        if(idwsdeb.gt.5.or.icas.lt.0)then
         write(6,*)('half transformed '),jh0h
         call prntm2(bc(jh0h),nbas3,nbas4,nbas3)
         write(6,*)('transpose of orbitals '),iorbx(isb)
         call prntm2(bc(iorbx(isb)),nbasdws(isb),nbas3,                 1d30s18
     $      nbasdws(isb))                                               1d30s18
        end if
        call dgemm('n','n',nbas4,nbas4,nbas3,1d0,bc(iorbx(isb)),nbas4,   6d28s04
     $           bc(jh0h),nbas3,0d0,bc(jh0mo),nbas4,                    2d23s04
     d' parahf. 11')
        if(idwsdeb.gt.5.or.icas.lt.0)then
         write(6,*)('in mo basis '),jh0mo,jh0mo-ih0mo
         call prntm2(bc(jh0mo),nbasdws(isb),nbasdws(isb),               1d30s18
     $        nbasdws(isb))                                             1d30s18
        end if
       end if
       jh0mo=jh0mo+nbas4*nbas4                                          5d8s18
      end do                                                            6d28s04
      nsdlp=0                                                           8d27s04
      do isa=1,nsymb
       nocx=noc(isa)
       ipkla(isa)=ibcoff
       nsdlp=nsdlp+1                                                    8d27s04
       ipk(nsdlp)=ibcoff                                                8d27s04
       isblp(1,nsdlp)=isa                                               8d27s04
       isblp(2,nsdlp)=isa                                               8d27s04
       isblp(3,nsdlp)=isa                                               8d27s04
       isblp(4,nsdlp)=isa                                               8d27s04
      end do
      do isa=1,nsymb
       do isb=1,nsymb
        if(isa.ne.isb)then
         noca=noc(isa)
         nocb=noc(isb)
         nsdlp=nsdlp+1                                                  8d27s04
         ipk(nsdlp)=ibcoff                                              8d27s04
         isblp(1,nsdlp)=isb                                             8d27s04
         isblp(2,nsdlp)=isb                                             8d27s04
         isblp(3,nsdlp)=isa                                             8d27s04
         isblp(4,nsdlp)=isa                                             8d27s04
        end if
       end do
      end do
      do isa=1,nsymb
       do isb=1,nsymb
        if(isa.ne.isb)then
         noca=noc(isa)
         nocb=noc(isb)
         nsdlp=nsdlp+1                                                  8d27s04
         ipk(nsdlp)=ibcoff                                              8d27s04
         isblp(1,nsdlp)=isb                                             8d27s04
         isblp(2,nsdlp)=isa                                             8d27s04
         isblp(3,nsdlp)=isb                                             8d27s04
         isblp(4,nsdlp)=isa                                             8d27s04
        end if
       end do
      end do
      do isa=1,nsymb
       do isb=1,nsymb
        if(isa.ne.isb)then
         noca=noc(isa)
         nocb=noc(isb)
         nsdlp=nsdlp+1                                                  8d27s04
         ipk(nsdlp)=ibcoff                                              9d15s04
         isblp(1,nsdlp)=isa                                             9d15s04
         isblp(2,nsdlp)=isb                                             8d27s04
         isblp(3,nsdlp)=isb                                             8d27s04
         isblp(4,nsdlp)=isa                                             8d27s04
        end if
       end do
      end do
      call second(time1)
      telap=time1-time2-tovr
      times(15,1)=times(15,1)+telap
      times(15,2)=times(15,2)+telap**2
      times(15,3)=times(15,3)+1d0
      do i=0,mynnode-1                                                  3d15s12
       ibc(imsg+i)=0                                                    3d15s12
      end do                                                            3d15s12
      idwsdeb=0
      if(lprt)write(6,*)('paraeri next')
      ibcb4=ibcoff                                                      7d7s21
      call paraeri(natom,ngaus,ibdat,nbasis,ihmat,iorbn,noc,            5d12s10
     $             ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,  11d29s12
     $     iter1,idwsdeb,idorel,ascale,nbasisp,2,0,idum,0,idum,idum,bc, 5d20s24
     $     ibc)                                                         5d20s24
      if(lprt)write(6,*)('back from paraeri')
      call second(time2)                                                10d19s04
      telap=time2-time1-tovr                                            10d19s04
      times(3,1)=times(3,1)+telap                                       10d19s04
      times(3,2)=times(3,2)+telap**2                                    10d19s04
      times(3,3)=times(3,3)+1d0                                         10d19s04
      time1=time2                                                       10d19s04
      idwsdeb=00
      n4x=0                                                             8d18s23
      call parajkfromh(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,         11d29s12
     $     iorbn,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,     1d5s12
     $     ioooo,ionex,noc,iter1,idwsdeb,ncomp,nvirt,i3x,0,             2d15s19
     $     idum,idum,dum,nbasisp,idorel,2,idum,idum,ibcb4,bc,ibc,       8d18s23
     $     ibc(isend),ibc(irecv),ibc(iptx),ibc(ipf),0,iff)              7d11s23
      if(lprt)write(6,*)('back from parajkfromh')
      idwsdeb=0
      if(idwsdeb.ne.0)then
       call dws_sync
       call dws_finalize
       stop
      end if
      if(myguessi.lt.0)then
       write(6,*)('aborting since mo vectors were not orthogonalized ') 3d14s16
       call dws_sync
       call dws_finalize
       stop
      end if
      call second(time2)                                                10d19s04
      telap=time2-time1-tovr                                            10d19s04
      times(4,1)=times(4,1)+telap                                       10d19s04
      times(4,2)=times(4,2)+telap**2                                    10d19s04
      times(4,3)=times(4,3)+1d0                                         10d19s04
      time1=time2                                                       10d19s04
      call second(timeq1)                                               3d7s07
      if(icas.eq.0)then                                                 8d8s14
      else                                                              8d8s14
       call second(time1)                                               4d5s18
       if(nlzz.ne.0)then                                                12d28s19
        call propcas(natom,ngaus,ibdat,isym,iapair,ibstor,isstor,idorel, 12d20s19
     $       ascale,multh,nbasisp,iorbn,ixlzz,noc,islz,nlzz,0,iorbsym,1,4d19s21
     $      iorbsymz,bc,ibc,epssym)                                     9d19s23
        call orbsymspace(nsymb,iorbsym,idoub,iacto,noc,                 4d19s21
     $       nvirt,mlx,inbr,ipbr,bc,ibc)                                11d9s22
       end if                                                           12d28s19
       call cas0(ih0mo,ioooo,noc,numpro,nowpro,potdws,                  8d8s14
     $      icallcas,iden1,j2den,inewl,nlzz,multh,ehf,ncomp,enobreit,   4d6s18
     $      dynw,iconv,icasvec,lprint,eavg2,jmats,kmats,nvirt,          2d22s19
     $      ibc(isend),ibc(irecv),ibc(iptx),ibc(ipf),0,ixlzz,islz,      2d20s20
     $      iptr,dum,dum,idum,idum,ibasisp,                             8d2s22
     $      icsfp,nfcnp,nctp,ipoint2p,idum,ibc(ismx),ibc(irelx),        8d2s22
     $      iorbn,morb,nbasisp,nryd,ixw1p,ixw2p,iptrbitp,iorbsym,       8d2s22
     $      iorbsymz,nextradata,extradata,1,idum,idum,idum,idum,dum,1,  6d9s22
     $      idum,idum,idum,idum,idum,idum,idum,dum,dum,bc,ibc,dum,      4d18s23
     $      ibc(nveccnt),icanon,ncasvecsize)                            10d11s24
       if(iprtr(5).ne.0)write(6,*)('starting energy .... '),ehf,eavg2
       call second(time2)                                               4d5s18
       telap=time2-time1-tovr
       times(5,1)=times(5,1)+telap
       times(5,2)=times(5,2)+telap**2
       times(5,3)=times(5,3)+1d0
      end if                                                            8d8s14
       ibcreset=ibcoff                                                   2d19s14
       call genergy(nsymb,ih0mo,idbly,iact,nbasdws,iden1,isblk,nsdlk,   2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo,ehf,ncomp,idwsdeb,enobreit,lprint,       3d3s17
     $     iconv,bc,ibc)                                                11d11s22
       ehfbsofar=ehf                                                    2d16s18
      if(iconv.eq.1)then                                                2d19s14
       go to 1404                                                       2d19s14
      end if                                                            2d19s14
      ehf2o=ehf                                                         4d19s13
      ehfrdb=ehf                                                        1d11s20
      defrombl=1d0                                                      1d11s20
      if(ehfo.gt.0d0.and.iter1.lt.itmax)then                            3d17s04
       if(lprint)write(6,447)ehf                                        12d2s22
  447  format('energy from starting orbitals: ',f16.6)                  12d2s22
       if(lprint.and.iter1.eq.1)then                                     4d20s18
        if(nrydb.eq.0)then
         if(icas.ne.0)then                                              8d4s22
          write(6,1066)                                                 8d4s22
 1066     format('iter.',7x,'  energy',2x,                                4d20s18
     $       '    change  approx-vari micro diag  xlam      grad orb')  1d13s20
         else                                                           8d4s22
          write(6,1166)                                                 8d4s22
 1166     format('iter.',7x,'  energy',2x,                              8d4s22
     $       '    change     micro diag  xlam')                         8d4s22
         end if                                                         8d4s22
        end if                                                          8d4s22
       end if                                                            4d20s18
      else
       if(ehfo.eq.1d0)then                                              1d13s17
        change=ehf+potdws                                               1d13s17
       else                                                             1d13s17
        change=ehf-ehfo
       end if                                                           1d13s17
       micitl=micit                                                     1d18s20
       if(lprint)then                                                   4d20s18
        write(ostng,1551)ehf                                            4d20s18
        if(change.ne.0d0)then                                           5d1s18
         ndig=nint(-log10(abs(change)))+1                                4d20s18
         ist=7+ndig                                                      4d20s18
          do k=ist,19                                                     4d20s18
          ostng(k:k)=' '                                                 4d20s18
         end do                                                          4d20s18
        end if                                                          5d1s18
        nrydb=nryd(1)                                                   3d3s20
        do isb=2,nsymb                                                  3d3s20
         nrydb=nrydb+nryd(isb)                                          3d3s20
        end do                                                          3d3s20
        if(nrydb.eq.0)then                                              8d4s22
         if(icas.ne.0)then                                              8d4s22
          write(6,1552)iter1,ostng,change,ehfbsofarl-ehf,micit,ncas0,   8d4s22
     $       xlamb,xx20o,xx1b                                           3d16s21
 1552     format(i4,1x,a19,2es9.1,2x,2i5,f8.2,es9.1,' to',2es9.1)           1d14s20
         else                                                           8d4s22
          write(6,2552)iter1,ostng,change,micit,ncas0,xlamb             8d4s22
 2552     format(i4,1x,a19,es9.1,2x,2i5,f8.2)                           8d4s22
         end if                                                         8d4s22
        end if                                                          8d4s22
 1551   format(f19.12)                                                  4d20s18
        if(abs(ehfbsofarl-ehf).gt.1d-4)then                             3d12s21
         jackitup=1                                                     3d12s21
        else                                                            3d12s21
         jackitup=0                                                     3d12s21
        end if                                                          3d12s21
       end if                                                           4d20s18
       change=sign(min(abs(change),abs(ehfbsofarl-ehf)),change)         4d20s18
       xx20l=xx20o                                                      3d12s21
       if(change.gt.1d4)then                                            2d1s19
        if(lprint)then
         write(6,*)('calculations diverging')                           3d21s07
         write(6,*)('reload orbitals for last iteration '),             3d21s07
     $        ('and decrease trust radius ')                            3d21s07
        end if                                                          3d21s07
        write(6,*)('getting orb from vecr '),ivecr,bc(ivecr)
        xlam=xlam*2d0                                                   3d21s07
        jvecr=ivecr-1                                                   3d21s07
        do isb=1,nsymb                                                  3d21s07
         nn=nbasdws(isb)*nbasisp(isb)*ncomp                             2d1s19
         iorb=iorbx(isb)-1                                              3d21s07
         write(6,*)isb,iorb,jvecr,nn
         do i=1,nn                                                      3d21s07
          bc(iorb+i)=bc(jvecr+i)                                        3d21s07
         end do                                                         3d21s07
         jvecr=jvecr+nn                                                 3d21s07
        end do                                                          3d21s07
        go to 404                                                       2d25s19
        ehf=ebestsofar                                                  6d19s17
       else                                                             3d21s07
        jvecr=ivecr-1                                                   3d21s07
        do isb=1,nsymb                                                  3d21s07
         nn=nbasisp(isb)*nbasdws (isb)                                  2d1s19
         iorb=iorbx(isb)-1                                              3d21s07
         do i=1,nn
          bc(jvecr+i)=bc(iorb+i)                                        3d21s07
         end do                                                         3d21s07
         jvecr=jvecr+nn                                                 3d21s07
        end do                                                          3d21s07
        ebestsofar=ehf                                                  6d19s17
       end if
       if(icas.ne.0)then                                                5d1s18
        epscnv=casdat(2)                                                5d1s18
        epswfn=casdat(3)                                                5d1s18
       else                                                             5d1s18
        epscnv=1d-12
        epswfn=1d-7                                                      4d5s18
       end if                                                           5d1s18
       if(iprtr(5).ne.0)write(6,*)('rmsa vs epswfn '),rmsa,epswfn
       if((abs(change).lt.epscnv.and.min(xx20o,rmsa).lt.epswfn.and.     4d21s23
     $      xlamb.lt.2d0)                                               4d21s23
     $      .or.iter1.ge.itmax)then                                     9d21s09
        if(abs(change).lt.epscnv.and.min(xx20o,rmsa).lt.epswfn)then     9d22s22
         if(lprint)write(6,*)('calculations converged')
         if(lprint)then                                                 9d22s22
          write(6,446)rmsa                                              12d2s22
  446     format('last change in orbitals: ',es9.2)                     12d2s22
          if(rmsa.gt.casdat(6))then                                     12d19s22
           write(6,*)                                                   6d20s24
           write(6,*)('this change is rather large!'),                  11d2s22
     $         (' perhaps an active orbital is doubly occupied?')       11d2s22
           write(6,*)('fix this by moving orbital from active space to')11d2s22
     $          ,(' doubly occupied space,')                            11d2s22
           write(6,*)('or by adding excited states to the calculation.')11d2s22
           write(6,*)                                                   6d20s24
          end if                                                        11d2s22
         end if                                                         9d22s22
         if(lprint)then                                                 4d20s18
          if(xx1b.gt.1d-4)write(6,*)('but grado is still tooooo large!') 1s1s20
          if(inocan.eq.0)then                                            12d6s21
           do isb=1,nsymb                                                 4d5s18
            if(max(idbly(isb),iacto(isb)).gt.0)                           2d28s21
     $         write(6,*)('for symmetry '),isb                          2d28s21
            if(idbly(isb).gt.0)then                                       2d28s21
             write(6,*)                                                   2d28s21
     $         ('fock matrix eigenvalues for doubly occupied orbs')     2d28s21
             call prntm2(bc(ifockeig(isb)),1,idbly(isb),1)                2d28s21
            end if                                                        2d28s21
            if(iacto(isb).gt.0)then                                       4d5s18
             write(6,*)('n.o. occupation numbers')                        2d28s21
             call prntm2(bc(noocc(isb)),1,iacto(isb),1)                   4d5s18
            end if                                                        4d5s18
           end do                                                         4d5s18
          end if                                                        12d6s21
         end if                                                         4d20s18
         ehf=return(8)                                                  3d13s23
         iconv=1                                                        2d19s14
         ibcoff=ibcreset                                                2d19s14
         go to 404                                                      2d19s14
        end if
        if(iter1.ge.itmax.and.nrydb.eq.0)then                           1d22s20
         if(lprint)write(6,*)('maximum number of iterations reached ')            3d17s04
         if(lprint)write(6,*)('last change in orbitals: '),rmsa         4d5s18
        end if                                                          3d17s04
 1404   continue                                                        2d19s14
        call second(time2)                                              2d20s14
        call dws_gsumf(telapo,15)                                       4d26s18
        times(11,1)=timb(1,1)                                           11d26s04
        times(11,2)=timb(1,2)                                           11d26s04
        times(11,3)=timb(1,3)                                           11d26s04
        telap=time2-tstart-tovr                                         10d19s04
        times(12,1)=times(12,1)+telap
        times(12,2)=times(12,2)+telap**2
        times(12,3)=times(12,3)+1d0                                     10d19s04
        iarg1=(loc(times(17,3))-loc(times(1,1)))/8                      2d20s14
        call dws_gsumf(times,iarg1)                                     3d15s12
        if(lprint)write(6,*)('total telapo = '),telapo                    3d2s17
        times(3,1)=times(3,1)-times(11,1)                               3d7s07
        times(3,2)=times(3,2)-times(11,2)                               3d7s07
        times(8,1)=times(8,1)-times(5,1)-times(6,1)-times(7,1)          4d23s18
        times(8,2)=times(8,2)-times(5,2)-times(6,2)-times(7,2)          4d23s18
        do 1401 i=1,17
         if(times(i,3).gt.0)then
          times(i,1)=times(i,1)/times(i,3)
          times(i,2)=times(i,2)/times(i,3)
          dev=sqrt(abs(times(i,2)-times(i,1)**2))
          if(lprint)then
          write(6,1402)i,timlab(i),times(i,1),dev,
     $         times(i,1)*times(i,3),times(i,3)
          end if
 1402     format('timer ',i4,1x,a8,' ave time ',1pe10.3,
     $         ' dev ',1pe10.3,' total ',1pe10.3,
     $         ' no. calls ',0pf10.0)
         end if
 1401   continue
        if(iconv.ne.0.and.xlamb.gt.2d0)then                             8d27s19
         iconv=0                                                        8d27s19
         if(lprint)                                                     8d27s19
     $        write(6,*)('xlamb is too large - force unconverged flag') 8d27s19
        end if                                                          8d27s19
        if(lprint)then                                                  3d2s17
         if((iconv.ne.0.and.xx1b.lt.1d-4.and.min(xx20o,rmsa).le.epswfn) 12d23s22
     $       .or.nrydb.ne.0)then                                        11d2s22
          write(6,*)('opening file forbs ')                               2d18s16
         else                                                           8d26s19
          write(6,*)('best orbitals so far ')                           8d26s19
          write(6,*)('opening file borbs ')                             8d26s19
         end if                                                         8d26s19
        end if                                                          3d2s17
        if(mynowprog.eq.0)then                                          3d2s17
         if((iconv.ne.0.and.xx1b.lt.1d-4.and.min(xx20o,rmsa).le.epswfn) 3d27s23
     $       .or.nrydb.ne.0)then                                        11d7s22
          open(unit=1,file='forbs',form='unformatted')                    2d18s16
          write(6,*)('final converged orbitals:')                       4d12s23
          if(molden.ne.0)then                                           11d22s19
           call moldenf(ibdat,iorbn,nsymb,ngaus,ibstor,isstor,noocc,    11d22s19
     $         idoub,iacto,nbasisp,bc,ibc)                              11d9s22
          end if                                                        11d22s19
         else                                                           8d26s19
          write(6,*)('opening borbs ')
          open(unit=1,file='borbs',form='unformatted')                  8d26s19
          write(6,*)('best orbitals so far:')                           4d12s23
         end if                                                         8d26s19
         ntag=0                                                         5d26s21
         if(iagrp(1).ne.0)then                                          5d26s21
          ntag=6+natom                                                  5d26s21
         end if                                                         5d26s21
         nextradatatag=nextradata+ntag                                  5d26s21
         write(1)nsymb,idorel,ngaus,natom,nwcont,numeminus,lmax,        11d20s19
     $        nbasallp,multh,ascale,potdws,ipropsym,(isinfo(j,1),j=1,11)1d2s20
     $        ,nextradatatag                                            5d26s21
         write(1)(nbasdws(isb),isb=1,nsymb),                            6d3s21
     $        (nlorcont(isb),isb=1,nsymb)                               6d3s21
         write(1)(nbasisc(isb),isb=1,nsymb)                             2d15s19
         write(1)(nbasisp(isb),isb=1,nsymb)                             1d28s19
         if(nrydb.ne.0)then                                             2d4s23
          ehfbsofarl=ehfbsofarl0                                        2d4s23
         end if                                                         2d4s23
         do isb=1,nsymb                                                 5d4s23
          ipack2(2)=icanon(isb)                                         5d4s23
          ipack2(1)=idoub(isb)                                          5d4s23
          idoub(isb)=ipack4                                             5d4s23
         end do                                                         5d4s23
         if(ntag.eq.0)then                                              5d26s21
          write(1)isym,(xcartx(i),atnumx(i),i=1,natom*3),                5d26s21
     $        (extradata(i),i=1,nextradata),ehfbsofarl,                 1d10s19
     $        (idoub(isb),iacto(isb),isb=1,nsymb),nlzz                  1d2s20
         else                                                           5d26s21
          itag=ibcoff                                                   5d26s21
          ibcoff=itag+ntag                                              5d26s21
          jtag=itag                                                     5d26s21
          do i=1,2                                                      5d26s21
           do j=1,3                                                     5d26s21
            bc(jtag)=cmx(j,i)                                           5d26s21
            jtag=jtag+1                                                 5d26s21
           end do                                                       5d26s21
          end do                                                        5d26s21
          do i=1,natom                                                  5d26s21
           ibc(jtag)=iagrp(i)                                           5d26s21
           jtag=jtag+1                                                  5d26s21
          end do                                                        5d26s21
          write(1)isym,(xcartx(i),atnumx(i),i=1,natom*3),                5d26s21
     $        (extradata(i),i=1,nextradata),
     $        (bc(itag+mtag),mtag=0,ntag-1),ehfbsofarl,                 5d26s21
     $        (idoub(isb),iacto(isb),isb=1,nsymb),nlzz                  1d2s20
          ibcoff=itag                                                   5d26s21
         end if                                                         5d26s21
         do isb=1,nsymb                                                 5d4s23
          ipack4=idoub(isb)
          idoub(isb)=ipack2(1)
         end do                                                         5d4s23
         write(1)((iapair(j,i),j=1,3),i=1,natom)                        10d27s17
         nbdat=ngaus*9+nwcont
         write(1)(bc(ibdat+i),i=0,nbdat-1)                              1d28s19
         write(1)(ibstor(i),isstor(i),i=1,nbasallp*3)                   11d21s19
         do isb=1,nsymb                                                 10d23s17
          ibcode(isb)=ibcoff                                            10d23s17
          ibcoff=ibcode(isb)+nbasdws(isb)                               1d30s18
          call enough('parahf. 41',bc,ibc)                              3d13s23
c
c     ibcode=0 for virtual orbs, 1 for doubly occupied orbs,            10d23s17
c     2 for active orbs, 3 for rydberg orbs...                                            10d23s17
c
          do i=0,nbasdws(isb)-1                                         1d30s18
           ibc(ibcode(isb)+i)=0                                         10d23s17
          end do                                                        10d23s17
          if(icas.eq.0)then                                             3d27s20
           do i=0,noc(isb)-1                                            10d23s17
            ibc(ibcode(isb)+i)=1                                        10d23s17
           end do                                                       10d23s17
          else                                                          10d23s17
           do i=0,idoub(isb)-1                                          3d27s20
            ibc(ibcode(isb)+i)=1                                        3d27s20
           end do                                                       3d27s20
           jbcode=ibcode(isb)+idoub(isb)                                3d27s20
           do i=0,iacto(isb)-1                                          3d27s20
            ibc(jbcode+i)=2                                             3d27s20
           end do                                                       3d27s20
           jbcode=jbcode+iacto(isb)                                     3d27s20
           do i=0,nryd(isb)-1                                           3d27s20
            ibc(jbcode+i)=3                                             3d27s20
           end do                                                       3d27s20
          end if                                                        10d23s17
          write(1)(ibc(ibcode(isb)+i),i=0,nbasdws(isb)-1)               1d30s18
          ibcoff=ibcode(isb)                                            10d23s17
         end do                                                         10d23s17
        end if                                                          3d2s17
        if(iacto(1).gt.0.and.nsymb.eq.1)then                            3d19s24
         if(lprint)then                                                 8d20s24
          write(6,*)('mo occupation nos. ')                              3d19s24
          call prntm2(bc(noocc(1)),1,iacto(1),1)                         3d19s24
         end if                                                         8d20s24
         is=0                                                           3d19s24
 1837    continue                                                       3d19s24
          ref=bc(noocc(1)+is)                                           3d19s24
          do i=is,iacto(1)-1                                            3d19s24
           if(abs(ref-bc(noocc(1)+i)).gt.1d-8)then                      3d19s24
            ie=i-1                                                      3d19s24
            go to 1387                                                  3d19s24
           end if                                                       3d19s24
          end do                                                        3d19s24
          ie=iacto(1)-1                                                 3d19s24
 1387     continue                                                      3d19s24
          ndeg=ie+1-is                                                  3d19s24
          if(ndeg.gt.1)then                                             3d19s24
           i0=is+idbly(1)                                               3d19s24
           jj=iorbn(1)+nbasisp(1)*ncomp*i0
           if(lprint)then                                               8d20s24
            write(6,*)('mo occ same for vectors '),is+1,('through '),    3d19s24
     $         ie+1                                                     3d19s24
            write(6,*)('starting orbs ')
            call prntm2(bc(jj),nbasisp(1)*ncomp,ndeg,nbasisp(1)*ncomp)
           end if                                                       8d20s24
           do i=i0,ie+idbly(1)-1                                        3d19s24
            do j=i+1,ie+idbly(1)                                        3d19s24
             if(lprint)write(6,*)('rotating '),j,(' and '),i            8d20s24
             ii=iorbn(1)+nbasisp(1)*ncomp*i
             jj=iorbn(1)+nbasisp(1)*ncomp*j
             do k=0,nbasisp(1)*ncomp-1
              if(abs(bc(ii+k)).gt.1d-8.and.abs(bc(jj+k)).gt.1d-8)then
               if(lprint)write(6,*)('looking at row '),k+1
               aa=bc(ii+k)
               bb=bc(jj+k)
               if(lprint)write(6,*)('a,b: '),aa,bb
               try1=-bb/aa
               try2=aa/bb
               if(lprint)write(6,*)('try1,try2 '),try1,try2
               if(abs(try1).lt.abs(try2))then
                if(lprint)write(6,*)('go with try1 '),try1
                try=try1
               else
                if(lprint)write(6,*)('go with try2 '),try2
                try=try2
               end if
               cc=1d0/sqrt(1d0+try*try)
               ss=cc*try
               if(lprint)write(6,*)('cc,ss '),cc,ss
               go to 1873
              end if
             end do
             cc=1d0                                                       3d19s24
             ss=0d0                                                       3d19s24
 1873        continue
             do k=0,nbasisp(1)*ncomp-1
              anew=bc(ii+k)*cc-bc(jj+k)*ss
              bnew=bc(ii+k)*ss+bc(jj+k)*cc
              bc(ii+k)=anew
              bc(jj+k)=bnew
             end do
            end do
           end do
           jj=iorbn(1)+nbasisp(1)*ncomp*i0
           if(lprint)then                                               8d20s24
            write(6,*)('ending orbs ')
            call prntm2(bc(jj),nbasisp(1)*ncomp,3,nbasisp(1)*ncomp)
           end if                                                       8d20s24
          end if                                                        3d19s24
          is=ie+1                                                       3d19s24
         if(is.le.iacto(1)-1)go to 1837                                 3d19s24
c                         567 890         705869
        end if                                                          3d19s24
        call mofromao(ivecs,ieigs,iorbn,morb,idorel,ixinv,nbasisp,      2d15s19
     $       nbasisc,ismile,bc,ibc,nposs)                               5d2s23
        jvecr=ivecr                                                     8d8s14
        jdwss=idwss
        do isb=1,nsymb                                                  8d8s14
         nrow=nbasisp(isb)*ncomp                                        2d3s19
         if(lprint.and.nbasdws(isb).gt.0)then                           1d30s18
          write(6,*)('for symmetry block '),isb                         11d17s22
          write(6,*)('ao basis ')
          if(nlzz.eq.2)then                                             4d20s21
           if(nsymb.eq.8)then                                           10d1s21
            if(isb.eq.2.or.isb.eq.3.or.isb.eq.5.or.isb.eq.8)then        10d1s21
             cdum='u'                                                   10d1s21
            else                                                        10d1s21
             cdum='g'                                                   10d1s21
            end if                                                      10d1s21
           else                                                         10d1s21
            cdum=' '                                                    10d1s21
           end if                                                       10d1s21
           write(6,*)('lz qns '),cdum                                   10d1s21
           write(6,88)(mod(ibc(iorbsym(isb)+i),50),i=0,nbasdws(isb)-1)  10d1s21
          else if(nlzz.ne.0)then                                        4d20s21
           write(6,*)('l over lz qns')                                  4d20s21
           write(6,88)(ibc(iorbsym(isb)+i),i=0,nbasdws(isb)-1)          4d20s21
           write(6,88)(ibc(iorbsymz(isb)+i),i=0,nbasdws(isb)-1)         4d20s21
          end if                                                        4d20s21
          if(idorel.eq.0)then                                           1d2s23
           call prntm2(bc(iorbn(isb)),nbasisp(isb),nbasdws(isb)*ncomp,    5d10s19
     $        nbasisp(isb))                                             5d10s19
          else                                                          1d2s23
           call prntm2r(bc(iorbn(isb)),nbasisp(isb),nbasdws(isb)*ncomp, 1d2s23
     $        nbasisp(isb))                                             5d10s19
          end if                                                        1d2s23
          write(6,*)('ob basis ')
          if(nposs(isb).ne.0)then                                       5d2s23
           call prntm2(bc(morb(isb)),nbasisc(isb),nbasdws(isb),          2d1s23
     $        nbasisc(isb))                                             2d1s23
          else                                                          5d2s23
           call prntm2(bc(morb(isb)),nbasdws(isb),nbasdws(isb),         5d2s23
     $        nbasdws(isb))                                             5d2s23
          end if                                                        5d2s23
         end if                                                          3d2s17
         if(mynowprog.eq.0)then                                          3d2s17
          if(nrow.gt.0)then                                             12d18s19
           write(1)(bc(iorbn(isb)+i),i=0,nbasdws(isb)*nrow-1)             2d15s19
          else                                                          12d18s19
           write(1)                                                     12d18s19
          end if                                                        12d18s19
         end if
         if(mynowprog.eq.0)then                                         3d2s17
          if(nbasdws(isb).gt.0)then                                     12d18s19
           write(1)(bc(morb(isb)+i),i=0,nbasisc(isb)*nbasdws(isb)-1)    5d4s23
          else                                                          12d18s19
           write(1)                                                     12d18s19
          end if                                                        12d18s19
         end if                                                         3d2s17
         jdwss=jdwss+nbasdws (isb)*nbasdws (isb)
         do i=0,nbasdws (isb)*nbasisp(isb)-1                            2d1s19
          bc(jvecr+i)=bc(iorbn(isb)+i)                                  8d8s14
         end do                                                         8d8s14
         jvecr=jvecr+nbasdws (isb)*nbasdws (isb)                        8d31s15
        end do                                                          8d8s14
        if(mynowprog.eq.0)then                                          8d27s19
         write(1)nstate
         do i=1,nstate                                                   5d24s19
          write(1)isinfo(4,i)                                            5d24s19
          write(1)(eweight(j,i),j=1,isinfo(4,i))                         5d24s19
          write(6,1351)(eweight(j,i),j=1,isinfo(4,i))
 1351     format(20es14.5)
         end do                                                          5d24s19
         if(idogrado.eq.0)then                                          12d2s22
          close(unit=1)                                                 12d2s22
         end if                                                         12d2s22
        end if                                                          8d27s19
        if(idogrado.ne.0)then                                           1d11s16
         if(iconv.eq.0.or.xx1b.gt.1d-4)then                             5d20s22
          if(lprint)write(6,*)                                          8d3s22
     $        ('orbitals not converged, so don''t compute gradients')
         else                                                           5d20s22
          if(icas.eq.0)then                                              5d4s22
           call grado(natom,ngaus,ibdat,nbasis,isym,iapair,               1d13s16
     $     ibstor,isstor,iso,nbb,idwsdeb,idorel,ascale,multh,ipropsym,  1d13s16
     $     ivecs,icanon,ieigs,iorbn,ih0,morb,idbly,ixinv,iooood,ionexd, 3d17s21
     $     jmatd,kmatd,ioooo,ionex,jmats,kmats,nbasdws ,nvirt,          2d15s19
     $     ih0mo,ipropmat,i3x,lprint,iact,nbasisp,nbasisc,ivecso,bc,ibc)11d10s22
          else                                                           5d4s22
           nart=0                                                            7d19s22
           do i=1,nstate                                                     7d19s22
            nart=nart+isinfo(4,i)                                            7d19s22
           end do                                                            7d19s22
           call gradcas(natom,ngaus,ibdat,nbasis,isym,iapair,            5d4s22
     $     ibstor,isstor,idorel,ascale,multh,ipropsym,ivecs,ieigs,iorbn,5d4s22
     $     ih0mo,morb,ixinv,ioooo,ionex,jmats,kmats,i3x,idoub,          5d19s22
     $     iact,noc,nbasisp,nbasisc,ivecso,ihessa,ihessb,ihessc,        5d4s22
     $     ihessd,iden1,j2den,lprint,inewl,nlzz,dynw,icasvec,isend,     5d16s22
     $         irecv,iptx,ipf,ixlzz,islz,iptr,bc(idorb),bc(isorb),      5d16s22
     $         mdon,mdoo,ibasis,ibc(icsf),nfcn,nct,nec,ibc(icsfpd),     5d16s22
     $         ismx,irelx,ixw1,ixw2,iptrbit,iorbsym,iorbsymz,           1d17s23
     $         ehfo,potdws,inocan,idogrado,nart,denergy,bc,ibc,nstate,  3d15s23
     $          isinfo,ibc(nveccnt))                                    4d18s23
          end if                                                         5d4s22
          close(unit=1)                                                 7d28s22
         end if                                                         5d20s22
         if(icas.eq.0)call paracis(-1,ioooo,jmats,kmats,bc(ih0mo),noc,  5d4s22
     $        nvirt,                                                    5d4s22
     $        nbasdws ,icore,multh,ieigv,ehf,ibc(ipropmat),natom,       4d28s22
     $        iapair,bc,ibc)                                            11d9s22
        end if                                                          1d11s16
        ibcoff=itopbc                                                   8d8s14
        return
       end if
      end if
      ehfo=ehf
      call second(time1)                                                10d19s04
      idwsdeb=00
      if(idwsdeb.gt.10)write(6,*)('on to amat! ')
      ibmat10=1                                                         4d23s18
      iamat10=1                                                         4d23s18
       if(icas.ne.0)then                                                10d23s17
        call second(time1)                                              4d5s18
        call buildhesscas(ioooo,ionex,jmats,kmats,bc(ih0mo),noc,ihessa,  4d25s17
     $      ihessb,ihessc,ihessd,                                       12d4s17
     $      nvirt,1,multh,iden1,j2den,idoub,iact,000,bc,ibc)            11d9s22
        call second(time2)                                              4d5s18
        telap=time2-time1-tovr                                          4d5s18
        times(7,1)=times(7,1)+telap                                     4d5s18
        times(7,2)=times(7,2)+telap**2                                  4d5s18
        times(7,3)=times(7,3)+1d0                                       4d5s18
        time1=time2                                                     4d5s18
        call buildcasgrad(ioooo,ionex,iamatx,idoub,iacto,noc,nvirt,     2d15s19
     $      1,ih0mo,multh,000,iden1,j2den,nsdlk,isblk,nsdlk1,isblk1,      2d16s18
     $      iamatb,namatt,rmssz0,1,1,idum,idum,bc,ibc)                  11d9s22
        xx20=rmssz0(1)                                                  1d2s20
        do isb=2,nsymb                                                  1d2s20
         xx20=max(xx20,rmssz0(isb))                                     1d2s20
        end do                                                          1d2s20
        if(iprtr(5).ne.0)write(6,*)('xx20 = '),xx20
        if(xx20.gt.xx20l*1.2d0)xlam=xlam*1d2                            3d12s21
      call dws_bcast(xlam,1)                                            3d12s21
        xx20o=xx20                                                      1d14s20
        if(c3.eq.0d0.and.xx20.lt.thirdthres.and.micitl.lt.50)then       8d12s22
         if(lprint)write(6,*)('switch to 3rd order ')                   1d17s20
         c0=-1d0                                                          11d29s17
         c1=1d0                                                           11d29s17
         cjk=-1d0                                                         11d29s17
         c3=1d0                                                           11d29s17
        end if                                                          1d17s20
        call second(time2)                                              4d5s18
        telap=time2-time1-tovr                                          4d5s18
        times(6,1)=times(6,1)+telap                                     4d5s18
        times(6,2)=times(6,2)+telap**2                                  4d5s18
        times(6,3)=times(6,3)+1d0                                       4d5s18
       else                                                             4d23s18
c
c     orbs fall into 2 groups: occupied and virtual,                    8d20s14
c     with occupied subdivided into doubly occupied and active.
c     we will build 2 sets of amat, bmat, and gmats, with
c     the first set for rotations between occupied and virtual orbs,
c     and the second set for rotations between double occupied and active
c     orbs. In the HF case, only the first set is required.
c
c     amat set 1:
c     amat: split into 2 parts, with part 1 being occ*occ part,
c     and part 2 being nvirt*occ.                                   8d26s04
c     A=HD+sum kl JklPlk
c     amat set 2: this is just the noc*idbly part of the first part      8d20s14
c     of amat set 1.                                                    8d20s14
c
      jh0mo=ih0mo                                                       8d26s04
      iamat10=ibcoff                                                    8d26s04
      do isb=1,nsymb                                                    8d26s04
       if(idwsdeb.gt.10)then
       write(6,*)('for symmetry block '),isb
       end if
       if(noc(isb).gt.0)then                                            8d27s04
        iamat1=ibcoff
        ibcoff=iamat1+noc(isb)*nbasdws(isb)                             2d1s19
        call enough('parahf. 42',bc,ibc)                                3d13s23
        do i=1,noc(isb)*nbasdws(isb)                                    2d1s19
         bc(iamat1+i-1)=0d0                                             8d30s04
        end do                                                          8d30s04
        nbas4=nbasdws(isb)                                              2d1s19
        do isa=1,nsymb                                                  2d14s12
         if(noc(isa).gt.0)then                                          2d14s12
          do isk=1,nsdlk                                                2d14s12
           if(isblk(1,isk).eq.isblk(2,isk))then
            nrow=(noc(isblk(1,isk))*(noc(isblk(1,isk))+1))/2
           else
            nrow=noc(isblk(1,isk))*noc(isblk(2,isk))
           end if
           ncol=noc(isblk(3,isk))*noc(isblk(4,isk))
           if(isb.ne.isa)then                                           9d18s14
            if(iact(isa).gt.0.and.iact(isb).gt.0)then                   9d18s14
             isab=multh(isa,isb)                                         9d18s14
             do isc=1,nsymb                                              9d18s14
              isd=multh(isab,isc)                                        9d18s14
              if(iact(isc).gt.0.and.iact(isd).gt.0.and.isc.ne.isa.and.
     $             isc.ne.isb)then                                      9d18s14
               if(isblk(1,isk).eq.isa.and.isblk(2,isk).eq.isb.and.      9d18s14
     $            isblk(3,isk).eq.isc.and.isblk(4,isk).eq.isd)then      9d18s14
                call ilimts(noc(isc),noc(isd),mynprocg,mynowprog,il,    9d18s14
     $                 ih,i1s,i1e,i2s,i2e)                              9d18s14
                i10=i1s                                                 9d18s14
                i1n=noc(isc)                                            9d18s14
                nrow=noc(isa)*noc(isb)                                  9d18s14
                ii=ioooo(isk)                                           9d18s14
                do i2=i2s,i2e                                           9d18s14
                 if(i2.eq.i2e)i1n=i1e                                   9d18s14
                 do i1=i10,i1n                                          9d18s14
                  if(i2.gt.idbly(isd).and.i1.gt.idbly(isc))then         9d18s14
                   do l=0,iact(isa)-1                                   9d18s14
                    do k=0,noc(isb)-1                                   9d18s14
                     ixn=ii+l+idbly(isa)+noc(isa)*k                     9d18s14
                     do i=0,iact(isb)-1                                 9d18s14
                      iada=iamat1+k+nbasdws(isb)*(i+idbly(isb))         2d1s19
                      iad=j2den(isk)+l+iact(isa)*(i+iact(isb)           9d18s14
     $                   *(i1-idbly(isc)-1+iact(isc)*(i2-idbly(isd)-1)))9d18s14
                      bc(iada)=bc(iada)+bc(ixn)*bc(iad)                 9d18s14
                     end do                                             9d18s14
                    end do                                              9d18s14
                   end do                                               9d18s14
                  end if
                  ii=ii+nrow                                            9d18s14
                 end do                                                 9d18s14
                 i10=1                                                  9d18s14
                end do                                                  9d18s14
               end if                                                   9d18s14
               if(isblk(1,isk).eq.isa.and.isblk(2,isk).eq.isc.and.      9d18s14
     $              isblk(3,isk).eq.isb.and.isblk(4,isk).eq.isd)then    9d18s14
                igotit=1                                                9d18s14
                call ilimts(noc(isb),noc(isd),mynprocg,mynowprog,il,    9d18s14
     $                 ih,i1s,i1e,i2s,i2e)                              9d18s14
                i10=i1s                                                 9d18s14
                i1n=noc(isb)                                            9d18s14
                nrow=noc(isa)*noc(isc)                                  9d18s14
                ii=ioooo(isk)                                           9d18s14
                do i2=i2s,i2e                                           9d18s14
                 if(i2.eq.i2e)i1n=i1e                                   9d18s14
                 do i1=i10,i1n                                          9d18s14
                  if(i2.gt.idbly(isd))then                              9d18s14
                   do l=0,iact(isa)-1                                   9d18s14
                    do k=0,iact(isc)-1                                  9d18s14
                     ixn=ii+l+idbly(isa)+noc(isa)*(k+idbly(isc))        9d18s14
                     do i=0,iact(isb)-1                                 9d18s14
                      iada=iamat1+i1-1+nbasdws(isb)*(i+idbly(isb))      2d1s19
                      iad=j2den(isk)+l+iact(isa)*(k+iact(isc)           9d18s14
     $                   *(i+iact(isb)*(i2-idbly(isd)-1)))              9d18s14
                      bc(iada)=bc(iada)+bc(ixn)*bc(iad)                 9d18s14
                     end do                                             9d18s14
                    end do                                              9d18s14
                   end do                                               9d18s14
                  end if
                  ii=ii+nrow                                            9d18s14
                 end do                                                 9d18s14
                 i10=1                                                  9d18s14
                end do                                                  9d18s14
               end if                                                   9d18s14
              end if                                                    9d18s14
             end do                                                      9d18s14
            end if                                                      9d18s14
           end if                                                       9d18s14
           if(isblk(1,isk).eq.isa.and.isblk(2,isk).eq.isa.and.          2d14s12
     $        isblk(3,isk).eq.isb.and.isblk(4,isk).eq.isb)then
            call ilimts(noc(isb),noc(isb),mynprocg,mynowprog,il,ih,i1s, 2d14s12
     $          i1e,i2s,i2e)                                            2d14s12
            i10=i1s                                                     2d14s12
            i1n=noc(isb)                                                2d14s12
            nrow=(noc(isa)*(noc(isa)+1))/2                              2d14s12
            nrow2=(iact(isa)*(iact(isa)+1))/2                           9d15s14
            ii=ioooo(isk)-1                                             2d14s12
            do i2=i2s,i2e                                               2d14s12
             if(i2.eq.i2e)i1n=i1e                                       2d14s12
             do i1=i10,i1n                                              2d14s12
              iad=iamat1+i1-1+nbasdws(isb)*(i2-1)                       2d1s19
              if(i2.le.idbly(isb))then                                  8d22s14
               do i3=1,idbly(isa)                                         8d20s14
                ixn=(i3*(i3+1))/2                                        2d14s12
                bc(iad)=bc(iad)+4d0*bc(ixn+ii)                           2d14s12
               end do                                                    2d14s12
               do k=1,iact(isa)                                         9d4s14
                kp=k+idbly(isa)                                         9d4s14
                do l=1,iact(isa)                                        9d4s14
                 lp=l+idbly(isa)                                        9d4s14
                 in=min(kp,lp)                                            9d4s14
                 ix=max(kp,lp)                                            9d4s14
                 ixn=((ix*(ix-1))/2)+in                                 9d4s14
                 iadd=iden1(isa)+l-1+iact(isa)*(k-1)                    9d4s14
                 bc(iad)=bc(iad)+bc(ixn+ii)*bc(iadd)*2d0                9d5s14
                end do                                                  9d4s14
               end do                                                   9d4s14
              else                                                      8d22s14
               do l=1,idbly(isa)                                        9d5s14
                do i=1,iact(isb)                                        9d5s14
                 ixn=(l*(l+1))/2                                        9d5s14
                 iad=iamat1+i1-1+nbasdws(isb)*(i+idbly(isb)-1)          2d1s19
                 iadd=iden1(isb)+i2-idbly(isb)-1+iact(isb)*(i-1)        9d5s14
                 bc(iad)=bc(iad)+bc(ixn+ii)*bc(iadd)*2d0                9d5s14
                end do                                                  9d5s14
               end do                                                   9d5s14
               do l=1,iact(isa)                                         9d5s14
                lp=l+idbly(isa)                                         9d5s14
                do k=1,l
                 kp=k+idbly(isa)                                        9d5s14
                 in=min(lp,kp)                                          9d5s14
                 ix=max(lp,kp)                                          9d5s14
                 ixn=((ix*(ix-1))/2)+in
                 in=min(l,k)                                            9d15s14
                 ix=max(l,k)                                            9d15s14
                 ixn2=((ix*(ix-1))/2)+in                                9d15s14
                 do i=1,iact(isb)                                       9d5s14
                  iad=iamat1+i1-1+nbasdws(isb)*(i+idbly(isb)-1)         2d1s19
                  fmul=1d0                                               9d15s14
                  if(i.eq.i2-idbly(isb))fmul=2d0                        9d15s14
                  in=min(i2-idbly(isb),i)                               9d5s14
                  ix=max(i2-idbly(isb),i)                               9d5s14
                  iadd=j2den(isk)+ixn2-1+nrow2*(((ix*(ix-1))/2)+in-1)   9d15s14
                  bc(iad)=bc(iad)+fmul*bc(ixn+ii)*bc(iadd)              9d15s14
                 end do                                                  8d22s14
                end do                                                   8d22s14
               end do                                                   9d5s14
              end if                                                    8d22s14
              ii=ii+nrow                                                2d14s12
             end do                                                     2d14s12
             i10=1                                                      2d14s12
            end do                                                      2d14s12
           end if                                                       2d14s12
           if(isa.gt.isb)then                                            2d14s12
            if(isblk(1,isk).eq.isa.and.isblk(2,isk).eq.isb.and.          2d14s12
     $         isblk(3,isk).eq.isa.and.isblk(4,isk).eq.isb)then          2d14s12
             call ilimts(noc(isa),noc(isb),mynprocg,mynowprog,il,ih,i1s,2d14s12
     $           i1e,i2s,i2e)                                           2d14s12
             nrow=noc(isa)*noc(isb)                                     2d14s12
             i10=i1s                                                    2d14s12
             i1n=noc(isa)                                               2d14s12
             ii=ioooo(isk)-1-noc(isa)                                   2d14s12
             do i2=i2s,i2e                                              2d14s12
              if(i2.eq.i2e)i1n=i1e                                      2d14s12
              do i1=i10,i1n                                             2d14s12
               if(i1.gt.idbly(isa))then                                 9d22s14
                do l=0,iact(isa)-1                                       9d18s14
                 do i=0,idbly(isb)-1                                      9d18s14
                  ixn=ii+l+idbly(isa)+1+noc(isa)*(i+1)                    9d22s14
                  iada=iamat1+i2-1+nbasdws(isb)*i                       2d1s19
                  iadd=iden1(isa)+i1-idbly(isa)-1+iact(isa)*l             9d22s14
                  bc(iada)=bc(iada)-bc(ixn)*bc(iadd)                    9d29s14
                 end do                                                  9d18s14
                end do                                                   9d18s14
               else
                do l=0,iact(isb)-1                                       9d18s14
                 do i=0,iact(isb)-1                                     9d22s14
                  ixn=ii+i1+noc(isa)*(l+idbly(isb)+1)                   9d26s14
                  iada=iamat1+i2-1+nbasdws(isb)*(i+idbly(isb))          2d1s19
                  iadd=iden1(isb)+i+iact(isb)*l                         9d22s14
                  bc(iada)=bc(iada)-bc(ixn)*bc(iadd)                    9d22s14
                 end do                                                  9d18s14
                end do                                                   9d18s14
               end if                                                   9d22s14
               if(i1.le.idbly(isa))then                                 9d18s14
                if(i2.le.idbly(isb))then                                9d18s14
                do i3=1,noc(isb)                                        9d18s14
                 iad=iamat1+i3-1+nbasdws(isb)*(i2-1)                    2d1s19
                 iad2=ii+i1+i3*noc(isa)                                  2d14s12
                 bc(iad)=bc(iad)-2d0*bc(iad2)                            2d14s12
                end do                                                   2d14s12
                end if
               else                                                     9d18s14
                do k=0,iact(isa)-1                                      9d18s14
                 do l=0,iact(isb)-1                                     9d18s14
                  ixn=ii+1+k+idbly(isa)+noc(isa)*(l+1+idbly(isb))       9d18s14
                  do i=0,iact(isb)-1                                    9d18s14
                   iada=iamat1+i2-1+nbasdws(isb)*(i+idbly(isb))         2d1s19
                   iadd=j2den(isk)+k+iact(isa)*(l+iact(isb)*            9d18s14
     $                  (i1-idbly(isa)-1+iact(isa)*i))                  9d18s14
                   bc(iada)=bc(iada)+bc(iadd)*bc(ixn)                   9d18s14
                  end do                                                9d18s14
                 end do                                                 9d18s14
                end do                                                  9d18s14
               end if                                                   8d21s14
               ii=ii+nrow                                               2d14s12
              end do                                                    2d14s12
              i10=1                                                     2d14s12
             end do                                                     2d14s12
            end if                                                      2d14s12
           else if(isa.eq.isb)then                                      2d14s12
            if(isblk(1,isk).eq.isblk(2,isk).and.                        2d14s12
     $         isblk(2,isk).eq.isblk(3,isk).and.                        2d14s12
     $         isblk(3,isk).eq.isblk(4,isk).and.isblk(4,isk).eq.isb)then
             call ilimts(noc(isb),noc(isb),mynprocg,mynowprog,il,ih,i1s,2d14s12
     $           i1e,i2s,i2e)                                           2d14s12
             nrow=(noc(isb)*(noc(isb)+1))/2                             2d14s12
             i10=i1s                                                    2d14s12
             i1n=noc(isb)                                               2d14s12
             ii=ioooo(isk)-1                                            2d14s12
             do i2=i2s,i2e                                              2d14s12
              if(i2.eq.i2e)i1n=i1e                                      2d14s12
              do i1=i10,i1n                                             2d14s12
               if(i2.le.idbly(isb))then                                 9d5s14
                if(i1.le.idbly(isb))then                                9d9s14
                 do i3=1,noc(isb)                                        8d21s14
                  in=min(i2,i3)                                           2d14s12
                  ix=max(i2,i3)                                           2d14s12
                  ixn=((ix*(ix-1))/2)+in
                  iad=iamat1+i3-1+nbasdws(isb)*(i1-1)                   2d1s19
                  iad2=ii+ixn                                             2d14s12
                  bc(iad)=bc(iad)-2d0*bc(iad2)                            2d14s12
                 end do                                                   2d14s12
                end if                                                  9d9s14
                do k=1,iact(isb)                                         9d5s14
                 do i=1,iact(isb)                                        9d5s14
                  in=min(i2,k+idbly(isb))                               9d5s14
                  ix=max(i2,k+idbly(isb))                               9d5s14
                  ixn=((ix*(ix-1))/2)+in                                9d5s14
                  iadd=iden1(isb)+i-1+iact(isb)*(k-1)                   9d5s14
                  iada=iamat1+i1-1+nbasdws(isb)*(i+idbly(isb)-1)        2d1s19
                  bc(iada)=bc(iada)-0.5d0*bc(ixn+ii)*bc(iadd)             9d5s14
                 end do                                                  9d5s14
                end do                                                   9d5s14
                do i=1,iact(isb)                                        9d5s14
                 ip=i+idbly(isb)                                        9d5s14
                 do l=1,iact(isb)                                       9d5s14
                  lp=l+idbly(isb)                                       9d5s14
                  ixn=((lp*(lp-1))/2)+i2                                9d17s14
                  iadd=iden1(isb)+i-1+iact(isb)*(l-1)                   9d5s14
                  iada=iamat1+i1-1+nbasdws(isb)                         2d1s19
     $                 *(i+idbly(isb)-1)                                9d17s14
                  bc(iada)=bc(iada)-0.5d0*bc(ixn+ii)*bc(iadd)             9d5s14
                 end do                                                 9d5s14
                end do                                                  9d5s14
               else                                                     9d5s14
                do l=1,iact(isb)                                        9d5s14
                 lp=l+idbly(isb)                                        9d5s14
                 do i=1,idbly(isb)                                      9d5s14
                  ixn=((lp*(lp-1))/2)+i                                 9d5s14
                  iadd=iden1(isb)+i2-idbly(isb)-1+iact(isb)*(l-1)       9d5s14
                  iada=iamat1+i1-1+nbasdws(isb)*(i-1)                   2d1s19
                  bc(iada)=bc(iada)-0.5d0*bc(iadd)*bc(ixn+ii)             9d5s14
                 end do                                                 9d5s14
                end do                                                  9d5s14
                do i=1,idbly(isb)                                       9d5s14
                 do k=1,iact(isb)                                       9d5s14
                  kp=k+idbly(isb)                                       9d5s14
                  ixn=((kp*(kp-1))/2)+i                                 9d5s14
                  iadd=iden1(isb)+i2-idbly(isb)-1+iact(isb)*(k-1)       9d5s14
                  iada=iamat1+i1-1+nbasdws(isb)*(i-1)                   2d1s19
                  bc(iada)=bc(iada)-0.5d0*bc(iadd)*bc(ixn+ii)             9d5s14
                 end do                                                 9d5s14
                end do                                                  9d5s14
               end if                                                   8d21s14
               ii=ii+nrow                                               2d14s12
              end do                                                    2d14s12
              i10=1                                                     2d14s12
             end do                                                     2d14s12
            end if                                                      2d14s12
           else if(isa.lt.isb)then                                      2d14s12
            if(isblk(1,isk).eq.isb.and.isblk(2,isk).eq.isa.and.          2d14s12
     $         isblk(3,isk).eq.isb.and.isblk(4,isk).eq.isa)then          2d14s12
             call ilimts(noc(isb),noc(isa),mynprocg,mynowprog,il,ih,i1s,2d14s12
     $           i1e,i2s,i2e)                                           2d14s12
             nrow=noc(isa)*noc(isb)                                     2d14s12
             i10=i1s                                                    2d14s12
             i1n=noc(isb)                                               2d14s12
             ii=ioooo(isk)-1-noc(isb)                                   1d14s14
             do i2=i2s,i2e                                              2d14s12
              if(i2.eq.i2e)i1n=i1e                                      2d14s12
              do i1=i10,i1n                                             2d14s12
               if(i1.le.idbly(isb).and.i2.le.idbly(isa))then            9d18s14
                do i3=1,noc(isb)                                        8d21s14
                 iad=iamat1+i3-1+nbasdws(isb)*(i1-1)                    2d1s19
                 iad2=ii+i3+i2*noc(isb)                                  1d14s14
                 bc(iad)=bc(iad)-2d0*bc(iad2)                            2d14s12
                end do                                                   2d14s12
               end if                                                   8d21s14
               if(i2.le.idbly(isa))then                                 9d18s14
                do l=0,iact(isb)-1                                      9d18s14
                 do i=0,iact(isb)-1                                     9d18s14
                  iada=iamat1+i1-1+nbasdws(isb)*(i+idbly(isb))          2d1s19
                  ixn=l+idbly(isb)+1+noc(isb)*i2+ii                     9d18s14
                  iadd=iden1(isb)+i+iact(isb)*l                         9d18s14
                  bc(iada)=bc(iada)-bc(ixn)*bc(iadd)                    9d18s14
                 end do                                                 9d18s14
                end do                                                  9d18s14
               else                                                     9d26s14
                do i=0,idbly(isb)-1                                     9d29s14
                 do k=0,iact(isa)-1                                     9d29s14
                  iada=iamat1+i1-1+nbasdws(isb)*i                       2d1s19
                  ixn=ii+i+1+noc(isb)*(idbly(isa)+k+1)                  9d29s14
                  iadd=iden1(isa)+i2-idbly(isa)-1+iact(isa)*k           9d29s14
                  bc(iada)=bc(iada)-bc(ixn)*bc(iadd)                    9d29s14
                 end do                                                 9d26s14
                end do                                                  9d26s14
               end if                                                   9d18s14
               if(i2.gt.idbly(isa))then                                 9d18s14
                do k=0,iact(isb)-1                                      9d18s14
                 do l=0,iact(isa)-1                                     9d18s14
                  ixn=ii+k+idbly(isb)+1+noc(isb)*(l+idbly(isa)+1)       9d18s14
                  do i=0,iact(isb)-1                                     9d18s14
                   iada=iamat1+i1-1+nbasdws(isb)*(i+idbly(isb))         2d1s19
                   iadd=j2den(isk)+i+iact(isb)*(i2-idbly(isa)-1         9d18s14
     $                  +iact(isa)*(k+iact(isb)*l))                     9d18s14
                   bc(iada)=bc(iada)+bc(ixn)*bc(iadd)                   9d18s14
                  end do                                                9d18s14
                 end do                                                 9d18s14
                end do                                                  9d18s14
               end if                                                   9d18s14
               ii=ii+nrow                                               2d14s12
              end do                                                    2d14s12
              i10=1                                                     2d14s12
             end do                                                     2d14s12
            end if
           end if                                                        2d14s12
          end do                                                         2d14s12
          do isk=1,nsdlk1                                               3d9s12
           if(iact(isblk1(1,isk))*iact(isblk1(2,isk))                   9d30s14
     $          *iact(isblk1(3,isk)).gt.0.and.isblk1(1,isk).ne.         9d30s14
     $          isblk1(2,isk).and.isblk1(1,isk).ne.isblk1(3,isk).and.   9d30s14
     $          isblk1(1,isk).ne.isblk1(4,isk).and.                     9d30s14
     $          isblk1(4,isk).eq.isb.and.isa.eq.1)then                  9d30s14
            do jsk=1,nsdlk                                              9d30s14
             if(isblk1(1,isk).eq.isblk(1,jsk).and.isblk1(2,isk).eq.     9d30s14
     $          isblk(2,jsk).and.isblk1(3,isk).eq.isblk(3,jsk).and.     9d30s14
     $          isblk1(4,isk).eq.isblk(4,jsk))then                      9d30s14
              j2sk=jsk                                                  9d30s14
              irule=1                                                   9d30s14
              go to 1022                                                9d30s14
             end if                                                     9d30s14
             if(isblk1(1,isk).eq.isblk(1,jsk).and.isblk1(2,isk).eq.     9d30s14
     $          isblk(2,jsk).and.isblk1(3,isk).eq.isblk(4,jsk).and.     9d30s14
     $          isblk1(4,isk).eq.isblk(3,jsk))then                      9d30s14
              j2sk=jsk                                                  9d30s14
              irule=2                                                   9d30s14
              go to 1022                                                9d30s14
             end if                                                     9d30s14
            end do                                                      9d30s14
            write(6,*)('2 part density not matched for block abcd '),
     $           (isblk1(j,isk),j=1,4)
            call dws_sync
            call dws_finalize
            stop
 1022       continue                                                    9d30s14
            call ilimts(noc(isblk1(3,isk)),nvirt(isb),mynprocg,         2d15s19
     $           mynowprog,il,ih,i1s,i1e,i2s,i2e)                       9d30s14
            i10=i1s                                                     9d30s14
            i1n=noc(isblk1(3,isk))                                      9d30s14
            nrow=noc(isblk1(1,isk))*noc(isblk1(2,isk))                  9d30s14
            nrow2=iact(isblk1(1,isk))*iact(isblk1(2,isk))               10d3s14
            ii=ionex(isk)                                               9d30s14
            do i2=i2s,i2e                                               9d30s14
             if(i2.eq.i2e)i1n=i1e                                       9d30s14
             i2p=i2+noc(isb)                                            9d30s14
             do i1=i10,i1n                                              9d30s14
              if(i1.gt.idbly(isblk1(3,isk)))then                        9d30s14
               do i=0,iact(isb)-1                                       10d3s14
                if(irule.eq.1)then                                       9d30s14
                 j2plus=nrow2*(i1-idbly(isblk1(3,isk))-1                 9d30s14
     $              +iact(isblk1(3,isk))*i)                             10d3s14
                else                                                     9d30s14
                 j2plus=nrow2*(i+iact(isb)*(i1-idbly(isblk1(3,isk))-1)) 9d30s14
                end if                                                   9d30s14
                do k=0,iact(isblk1(1,isk))-1                            10d3s14
                 do l=0,iact(isblk1(2,isk))-1                           10d3s14
                  iada=iamat1+i2p-1+nbasdws(isb)*(i+idbly(isb))         2d1s19
                  i2e=ii+k+idbly(isblk1(1,isk))+noc(isblk1(1,isk))*      9d30s14
     $                 (l+idbly(isblk1(2,isk)))                          9d30s14
                  iadd=j2den(j2sk)+k+iact(isblk1(1,isk))*l+j2plus       10d3s14
                  bc(iada)=bc(iada)+bc(i2e)*bc(iadd)                     9d30s14
                 end do                                                  9d30s14
                end do                                                   9d30s14
               end do                                                    9d30s14
              end if                                                    9d30s14
              ii=ii+nrow                                                9d30s14
             end do                                                     9d30s14
             i10=1                                                      9d30s14
            end do                                                      9d30s14
           end if                                                       9d30s14
           ncol=noc(isblk1(3,isk))*nvirt(isblk1(4,isk))                 2d15s19
           if(isblk1(1,isk).eq.isblk1(2,isk))then
            nrow=(noc(isblk1(1,isk))*(noc(isblk1(1,isk))+1))/2
           else
            nrow=noc(isblk1(1,isk))*noc(isblk1(2,isk))
           end if
           if(isblk1(1,isk).eq.isa.and.isblk1(2,isk).eq.isa.and.        3d9s12
     $        isblk1(3,isk).eq.isb.and.isblk1(4,isk).eq.isb)then
            call ilimts(noc(isb),nvirt(isb),mynprocg,mynowprog,il,ih,   2d15s19
     $           i1s,i1e,i2s,i2e)                                        2d14s12
            do jsk=1,nsdlk                                              9d30s14
             if(isblk1(1,isk).eq.isblk(1,jsk).and.isblk1(2,isk).eq.     9d30s14
     $          isblk(2,jsk).and.isblk1(3,isk).eq.isblk(3,jsk).and.     9d30s14
     $          isblk1(4,isk).eq.isblk(4,jsk))then                      9d30s14
              j2sk=jsk                                                  9d30s14
              go to 1020                                                9d30s14
             end if                                                     9d30s14
            end do                                                      9d30s14
            write(6,*)('2 part density not matched for block 1 '),
     $           (isblk1(j,isk),j=1,4)
            call dws_sync
            call dws_finalize
            stop
 1020       continue                                                    9d30s14
            i10=i1s                                                     2d14s12
            i1n=noc(isb)                                                2d14s12
            nrow=(noc(isa)*(noc(isa)+1))/2                              8d4s14
            nrow2=(iact(isa)*(iact(isa)+1))/2                           9d29s14
            ii=ionex(isk)-1                                             8d4s14
            do i2=i2s,i2e                                               2d14s12
             i2p=i2+noc(isb)                                            2d14s12
             if(i2.eq.i2e)i1n=i1e                                       2d14s12
             do i1=i10,i1n                                              2d14s12
              if(i1.le.idbly(isb))then                                  8d22s14
               iad=iamat1+i2p-1+nbasdws(isb)*(i1-1)                     2d1s19
               do l=1,idbly(isa)                                        4d4s17
                ixn=(l*(l+1))/2                                         4d4s17
                bc(iad)=bc(iad)+4d0*bc(ixn+ii)                           2d14s12
               end do                                                    2d14s12
               do k=0,iact(isa)-1                                       10d3s14
                kp=k+idbly(isa)+1                                       10d3s14
                fmul=4d0                                                10d3s14
                do l=0,k                                                10d3s14
                 if(l.eq.k)fmul=2d0                                     10d3s14
                 lp=l+idbly(isa)+1                                      10d3s14
                 j2e=((kp*(kp-1))/2)+lp+ii                              4d4s17
                 iadd=iden1(isa)+l+iact(isa)*k                          10d3s14
                 bc(iad)=bc(iad)+bc(j2e)*bc(iadd)*fmul                  4d4s17
                end do                                                  10d3s14
               end do                                                   10d3s14
              else                                                      9d29s14
               do i=0,iact(isb)-1                                       10d3s14
                do l=0,idbly(isa)-1                                     10d3s14
                 lp=l+1                                                 10d3s14
                 j2e=((lp*(lp+1))/2)+ii                                 4d4s17
                 iada=iamat1+i2p-1+nbasdws(isb)*(i+idbly(isb))          2d1s19
                 iadd=iden1(isb)+i+iact(isb)*(i1-idbly(isb)-1)          10d3s14
                 bc(iada)=bc(iada)+2d0*bc(j2e)*bc(iadd)                 4d4s17
                end do                                                  10d3s14
               end do                                                   10d3s14
               do i=0,iact(isb)-1
                if(i1.eq.i+idbly(isb)+1)then                            9d30s14
                 fmul=2d0                                                9d29s14
                else                                                     9d29s14
                 fmul=1d0                                                9d29s14
                end if                                                   9d29s14
                do k=0,iact(isa)-1
                 do l=0,k                                                9d29s14
                  iad=iamat1+i2p-1+nbasdws(isb)*(i+idbly(isb))          2d1s19
                  ix=max(k,l)+idbly(isa)+1                               9d29s14
                  in=min(k,l)+idbly(isa)+1                               9d29s14
                  ixn=((ix*(ix-1))/2)+in+ii                              9d29s14
                  jx=max(k,l)+1                                         9d29s14
                  jn=min(k,l)+1                                         9d29s14
                  kx=max(i1-idbly(isb),i+1)                             9d29s14
                  kn=min(i1-idbly(isb),i+1)                             9d29s14
                  iadd=j2den(j2sk)+((jx*(jx-1))/2)+jn-1                 9d30s14
     $                 +nrow2*(((kx*(kx-1))/2)+kn-1)                    9d29s14
                  bc(iad)=bc(iad)+fmul*bc(iadd)*bc(ixn)                  9d29s14
                 end do                                                 9d29s14
                end do
               end do
              end if                                                    8d22s14
              ii=ii+nrow                                                2d14s12
             end do                                                     2d14s12
             i10=1                                                      2d14s12
            end do                                                      2d14s12
           end if                                                       2d14s12
           if(isa.ge.isb)then                                            2d14s12
            if(isblk1(1,isk).eq.isa.and.isblk1(2,isk).eq.isb.and.       3d9s12
     $         isblk1(3,isk).eq.isa.and.isblk1(4,isk).eq.isb)then       8d21s14
             call ilimts(noc(isa),nvirt(isb),mynprocg,mynowprog,il,ih,  2d15s19
     $           i1s,i1e,i2s,i2e)                                       2d14s12
            do jsk=1,nsdlk                                              9d30s14
             if(isblk1(1,isk).eq.isblk(1,jsk).and.isblk1(2,isk).eq.     9d30s14
     $          isblk(2,jsk).and.isblk1(3,isk).eq.isblk(3,jsk).and.     9d30s14
     $          isblk1(4,isk).eq.isblk(4,jsk))then                      9d30s14
              j2sk=jsk                                                  9d30s14
              go to 1021                                                9d30s14
             end if                                                     9d30s14
            end do                                                      9d30s14
            write(6,*)('2 part density not matched for block 2 '),
     $           (isblk1(j,isk),j=1,4)
            call dws_sync
            call dws_finalize
            stop
 1021       continue                                                    9d30s14
             if(isa.eq.isb)then                                         8d4s14
              nrow=(noc(isa)*(noc(isa)+1))/2                            8d4s14
              ii=ionex(isk)-1                                           8d4s14
             else                                                       8d4s14
              nrow=noc(isa)*noc(isb)                                     2d14s12
              ii=ionex(isk)-1-noc(isa)                                   2d14s12
             end if                                                     8d4s14
             i10=i1s                                                    2d14s12
             i1n=noc(isa)                                               2d14s12
             do i2=i2s,i2e                                              2d14s12
              i2p=i2+noc(isb)                                           2d14s12
              if(i2.eq.i2e)i1n=i1e                                      2d14s12
              if(isa.eq.isb)then                                        8d4s14
               do i1=i10,i1n                                             2d14s12
                if(i1.le.idbly(isa))then                                8d21s14
                 do i3=1,idbly(isb)                                         2d14s12
                  in=min(i1,i3)                                          8d4s14
                  ix=max(i1,i3)                                          8d4s14
                  iad2=((ix*(ix-1))/2)+in+ii                             8d4s14
                  iad=iamat1+i2p-1+nbasdws(isb)*(i3-1)                  2d1s19
                  bc(iad)=bc(iad)-2d0*bc(iad2)                            2d14s12
                 end do                                                   2d14s12
                 do i=0,iact(isb)-1                                       10d3s14
                  do k=0,iact(isb)-1                                      10d3s14
                   kp=k+idbly(isb)+1                                      10d3s14
                   j2e=((kp*(kp-1))/2)+i1+ii                            4d4s17
                   iada=iamat1+i2p-1+nbasdws(isb)*(i+idbly(isb))        2d1s19
                   iadd=iden1(isb)+i+iact(isb)*k                          10d3s14
                   bc(iada)=bc(iada)-bc(j2e)*bc(iadd)                   4d4s17
                  end do                                                  10d3s14
                 end do                                                   10d3s14
                else
                 do i=0,idbly(isb)-1                                      10d3s14
                  ip=i+1                                                  10d3s14
                  do l=0,iact(isa)-1                                      10d3s14
                   lp=l+idbly(isa)+1                                      10d3s14
                   j2e=((lp*(lp-1))/2)+ip+ii                            4d4s17
                   iada=iamat1+i2p-1+nbasdws(isb)*i                     2d1s19
                   iadd=iden1(isa)+i1-idbly(isa)-1+iact(isa)*l            10d3s14
                   bc(iada)=bc(iada)-bc(j2e)*bc(iadd)                   4d4s17
                  end do                                                  10d3s14
                 end do                                                   10d3s14
                end if                                                  8d21s14
                ii=ii+nrow                                               2d14s12
               end do                                                    2d14s12
              else                                                      8d4s14
               do i1=i10,i1n                                             2d14s12
                if(i1.le.idbly(isa))then                                8d21s14
                 do i3=1,idbly(isb)                                     8d21s14
                  iad=iamat1+i2p-1+nbasdws(isb)*(i3-1)                  2d1s19
                  iad2=ii+i1+i3*noc(isa)                                  2d14s12
                  bc(iad)=bc(iad)-2d0*bc(iad2)                            2d14s12
                 end do                                                   2d14s12
                 do i=0,iact(isb)-1                                       10d3s14
                  do k=0,iact(isb)-1                                      10d3s14
                   kp=k+idbly(isb)+1                                      10d3s14
                   j2e=i1+noc(isa)*kp+ii                                4d4s17
                   iada=iamat1+i2p-1+nbasdws(isb)*(i+idbly(isb))        2d1s19
                   iadd=iden1(isb)+i+iact(isb)*k                          10d3s14
                   bc(iada)=bc(iada)-bc(j2e)*bc(iadd)                   4d4s17
                  end do                                                  10d3s14
                 end do                                                   10d3s14
                else                                                    9d29s14
                 do i=0,idbly(isb)-1                                      10d3s14
                  ip=i+1                                                  10d3s14
                  do l=0,iact(isa)-1                                      10d3s14
                   lp=l+idbly(isa)+1                                      10d3s14
                   j2e=lp+noc(isa)*ip+ii                                4d4s17
                   iada=iamat1+i2p-1+nbasdws(isb)*i                     2d1s19
                   iadd=iden1(isa)+i1-idbly(isa)-1+iact(isa)*l            10d3s14
                   bc(iada)=bc(iada)-bc(j2e)*bc(iadd)                   4d4s17
                  end do                                                  10d3s14
                 end do                                                   10d3s14
                 do i=0,iact(isb)-1                                     10d3s14
                  iad=iamat1+i2p-1+nbasdws(isb)*(i+idoub(isb))          2d1s19
                  do k=0,iact(isb)-1                                    10d3s14
                   do l=0,iact(isa)-1                                   10d3s14
                    iadd=j2den(j2sk)+l+iact(isa)*(k+iact(isb)           10d3s14
     $                   *(i1-idoub(isa)-1+iact(isa)*i))                10d3s14
                    j2e=ii+l+idoub(isa)+1+noc(isa)*(k+idoub(isb)+1)     4d4s17
                    bc(iad)=bc(iad)+bc(j2e)*bc(iadd)                    4d4s17
                   end do                                               9d29s14
                  end do                                                9d29s14
                 end do                                                 9d29s14
                end if                                                  8d21s14
                ii=ii+nrow                                               2d14s12
               end do                                                    2d14s12
              end if                                                    8d4s14
              i10=1                                                     2d14s12
             end do                                                     2d14s12
            end if                                                      2d14s12
           else if(isa.lt.isb)then                                      2d14s12
            if(isblk1(1,isk).eq.isb.and.isblk1(2,isk).eq.isa.and.       3d9s12
     $         isblk1(3,isk).eq.isa.and.isblk1(4,isk).eq.isb)then       3d9s12
             call ilimts(noc(isa),nvirt(isb),mynprocg,mynowprog,il,ih,  2d15s19
     $           i1s,i1e,i2s,i2e)                                       2d14s12
            do jsk=1,nsdlk                                              9d30s14
             if(isblk1(1,isk).eq.isblk(1,jsk).and.isblk1(2,isk).eq.     9d30s14
     $          isblk(2,jsk).and.isblk1(3,isk).eq.isblk(4,jsk).and.     10d3s14
     $          isblk1(4,isk).eq.isblk(3,jsk))then                      10d3s14
              j2sk=jsk                                                  9d30s14
              go to 1023                                                9d30s14
             end if                                                     9d30s14
            end do                                                      9d30s14
            write(6,*)('2 part density not matched for block 2 '),
     $           (isblk1(j,isk),j=1,4)
            call dws_sync
            call dws_finalize
            stop
 1023       continue                                                    9d30s14
             nrow=noc(isa)*noc(isb)                                     2d14s12
             i10=i1s                                                    2d14s12
             i1n=noc(isa)                                               3d9s12
             ii=ionex(isk)-1-noc(isb)                                   1d14s14
             do i2=i2s,i2e                                              2d14s12
              i2p=i2+noc(isb)                                           2d14s12
              if(i2.eq.i2e)i1n=i1e                                      2d14s12
              do i1=i10,i1n                                             2d14s12
               if(i1.le.idbly(isa))then                                 8d21s14
                do i3=1,idbly(isb)                                         2d14s12
                 iad=iamat1+i2p-1+nbasdws(isb)*(i3-1)                   2d1s19
                 iad2=ii+i3+i1*noc(isb)                                  1d14s14
                 bc(iad)=bc(iad)-2d0*bc(iad2)                            2d14s12
                end do                                                   2d14s12
                do i=0,iact(isb)-1                                      10d3s14
                 do l=0,iact(isb)-1                                     10d3s14
                  lp=l+idbly(isb)+1                                     10d3s14
                  j2e=lp+noc(isb)*i1+ii                                 4d4s17
                  iada=iamat1+i2p-1+nbasdws(isb)*(i+idbly(isb))         2d1s19
                  iadd=iden1(isb)+l+iact(isb)*i                         10d3s14
                  bc(iada)=bc(iada)-bc(j2e)*bc(iadd)                    4d4s17
                 end do                                                 10d3s14
                end do                                                  10d3s14
               else                                                     10d3s14
                do i=0,idbly(isb)-1                                     10d3s14
                 do l=0,iact(isa)-1                                     10d3s14
                  j2e=i+1+noc(isb)*(l+idbly(isa)+1)+ii                  4d4s17
                  iada=iamat1+i2p-1+nbasdws(isb)*i                      2d1s19
                  iadd=iden1(isa)+i1-idbly(isa)-1+iact(isa)*l           10d3s14
                  bc(iada)=bc(iada)-bc(j2e)*bc(iadd)                    4d4s17
                 end do                                                 10d3s14
                end do                                                  10d3s14
                 do i=0,iact(isb)-1                                     10d3s14
                  iad=iamat1+i2p-1+nbasdws(isb)*(i+idoub(isb))          2d1s19
                  do k=0,iact(isb)-1                                    10d3s14
                   do l=0,iact(isa)-1                                   10d3s14
                    iadd=j2den(j2sk)+k+iact(isb)*(l+iact(isa)           10d3s14
     $                   *(i+iact(isb)*(i1-idoub(isa)-1)))              10d3s14
                    j2e=ii+k+idoub(isb)+1+noc(isb)*(l+idoub(isa)+1)     4d4s17
                    bc(iad)=bc(iad)+bc(j2e)*bc(iadd)                    4d4s17
                   end do                                               9d29s14
                  end do                                                9d29s14
                 end do                                                 9d29s14
               end if                                                   8d21s14
               ii=ii+nrow                                               2d14s12
              end do                                                    2d14s12
              i10=1                                                     2d14s12
             end do                                                     2d14s12
            end if                                                      2d14s12
           end if                                                       2d14s12
          end do                                                        2d14s12
         end if                                                         2d14s12
        end do                                                          2d14s12
        call dws_sync                                                   8d30s04
        nwds=noc(isb)*nbasdws(isb)                                      2d1s19
        iarg1=nwds
        call dws_gsumf(bc(iamat1),iarg1)                                 4d12s07
        call dws_sync                                                   1d5s12
        if(idwsdeb.gt.10)then
        iarg1=noc(isb)
        iarg2=nbasdws(isb)                                              2d1s19
        write(6,*)('J*P ')
        call prntm2(bc(iamat1),iarg2,iarg1,iarg2)
        if(nsymb.eq.1)then
         call printa(bc(iamat1),nbasdws,0,1,0,noc,0,1,0,bc(ibcoff))     2d1s19
         end if
        end if
       end if
        do i=1,idbly(isb)                                               8d20s14
         do j=1,nbas4                                                   8d20s14
          iad1=iamat1+j-1+nbas4*(i-1)                                    8d20s14
          iad2=jh0mo+j-1+nbas4*(i-1)                                    8d20s14
          bc(iad1)=bc(iad1)+2d0*bc(iad2)                                8d20s14
         end do                                                         8d20s14
        end do                                                          8d20s14
        if(iact(isb).gt.0)then                                          8d20s14
         jh0moa=jh0mo+nbas4*idbly(isb)                                  4d4s17
         iamat1a=iamat1+idbly(isb)*nbas4                                4d4s17
         write(6,*)('dgemmx')
         call dgemm('n','n',nbas4,iact(isb),iact(isb),1d0,              4d4s17
     $              bc(jh0moa),nbas4,bc(iden1(isb)),iact(isb),1d0,      8d20s14
     $              bc(iamat1a),nbas4,                                  8d20s14
     d' parahf. 12')
        end if                                                          8d20s14
        if(idwsdeb.gt.10)then
         iarg2=nbas4
         iarg1=idbly(isb)
        write(6,*)('full amat')
        call prntm2(bc(iamat1),iarg2,iarg1,iarg2)
        if(nsymb.eq.1)then
         call printa(bc(iamat1),nbasdws,0,1,0,noc,0,1,0,bc(ibcoff))     2d1s19
        end if
       end if                                                           8d27s04
       jh0mo=jh0mo+nbasdws(isb)*nbasdws(isb)                            2d1s19
       call msl('amat',bc(iamat1),noc(isb),nbasdws(isb))                2d1s19
      end do                                                            8d27s04
      call second(time2)
      telap=time2-time1-tovr
      times(6,1)=times(6,1)+telap
      times(6,2)=times(6,2)+telap**2
      times(6,3)=times(6,3)+1d0
      end if                                                            4d23s18
c
c     initialize umat
c
      do isb=1,nsymb
       if(icas.eq.0)then                                                1d3s18
       iumat(isb)=ibcoff
       itmat(isb)=iumat(isb)+nbasdws(isb)*nbasdws(isb)                  2d1s19
       ibcoff=itmat(isb)+nbasdws(isb)*nbasdws(isb)                      2d1s19
       call enough('parahf. 43',bc,ibc)                                 3d13s23
       do 187 i=1,nbasdws(isb)                                          2d1s19
        jumat=iumat(isb)-1+(i-1)*nbasdws(isb)                           2d1s19
        jtmat=itmat(isb)-1+(i-1)*nbasdws(isb)                           2d1s19
        do 188 j=1,nbasdws(isb)                                         2d1s19
         bc(jumat+j)=0d0
         bc(jtmat+j)=0d0
  188   continue
        bc(jumat+i)=1d0
  187  continue
       else                                                             1d3s18
c                                                                       1d3s18
c     for Werner updates, we need to store occ-virt transformation      1d3s18
c     separately from doub-active transformation. We will store the     1d3s18
c     former in umat and the latter in tmat.                            1d3s18
c     in itotm will be the product of the two.                          1d9s18
c                                                                       1d3s18
        iumat(isb)=ibcoff                                               1d3s18
        itmat(isb)=iumat(isb)+nbasdws(isb)*nbasdws(isb)                 2d1s19
        itotm(isb)=itmat(isb)+noc(isb)*noc(isb)                         1d9s18
        ibcoff=itotm(isb)+nbasdws(isb)*nbasdws(isb)                     2d1s19
        iumatb(isb)=ibcoff                                              2d16s18
        itmatb(isb)=iumatb(isb)+nbasdws(isb)*nbasdws(isb)               2d1s19
        itotmb(isb)=itmatb(isb)+noc(isb)*noc(isb)                       4d5s18
        ibcoff=itotmb(isb)+nbasdws(isb)*nbasdws(isb)                    2d1s19
        call enough('parahf. 44',bc,ibc)                                3d13s23
        jumat=iumat(isb)                                                1d3s18
        jtot=itotm(isb)                                                 4d5s18
        do i=0,nbasdws(isb)-1                                           2d1s19
         do j=0,nbasdws(isb)-1                                          2d1s19
          bc(jumat+j)=0d0                                               1d3s18
          bc(jtot+j)=0d0                                                4d5s18
         end do                                                         1d3s18
         bc(jumat+i)=1d0                                                1d3s18
         jumat=jumat+nbasdws(isb)                                       2d1s19
         bc(jtot+i)=1d0                                                 4d5s18
         jtot=jtot+nbasdws(isb)                                         2d1s19
        end do                                                          1d3s18
        jtmat=itmat(isb)                                                1d3s18
        do i=0,noc(isb)-1                                               1d3s18
         do j=0,noc(isb)-1                                              1d3s18
          bc(jtmat+j)=0d0                                               1d3s18
         end do                                                         1d3s18
         bc(jtmat+i)=1d0                                                1d3s18
         jtmat=jtmat+noc(isb)                                           1d3s18
        end do                                                          1d3s18
        do i=0,2*nbasdws(isb)*nbasdws(isb)+noc(isb)*noc(isb)-1          2d1s19
         bc(iumatb(isb)+i)=bc(iumat(isb)+i)                             2d16s18
        end do                                                          2d16s18
       end if                                                           1d3s18
      end do                                                            8d31s04
c
c     bmat. Same drill as amat, except store transpose of ov
c
      if(icas.eq.0)then                                                 4d23s18
      if(idwsdeb.gt.10)then
      write(6,*)('bmat!')
      end if
      jh0mo=ih0mo                                                       8d31s04
      ibmat10=ibcoff                                                    8d31s04
      iamat=iamat10                                                     8d31s04
      do isb=1,nsymb                                                    8d31s04
       if(idwsdeb.gt.10)write(6,*)('for symmetry block '),isb
       if(noc(isb).gt.0)then                                            8d31s04
        ibmat1=ibcoff                                                    8d31s04
        ibmat2=ibmat1+noc(isb)*noc(isb)                                 9d9s04
        ibcoff=ibmat2+noc(isb)*nvirt(isb)                               2d15s19
        call enough('parahf. 45',bc,ibc)                                3d13s23
        do i=1,noc(isb)                                                 9d9s04
         j1=ibmat1-1+(i-1)*noc(isb)                                     9d9s04
         j2=iamat-1+(i-1)*nbasdws(isb)                                  2d1s19
         do j=1,noc(isb)                                                9d9s04
          bc(j1+j)=bc(j2+j)                                             9d9s04
         end do                                                         9d9s04
        end do                                                          9d9s04
        do j=1,noc(isb)
         j1=ibmat2-1+(j-1)*nvirt(isb)                                   2d15s19
         j2=iamat+noc(isb)-1+(j-1)*nbasdws(isb)                         2d1s19
         do i=1,nvirt(isb)                                              2d15s19
          bc(j1+i)=bc(j2+i)
         end do
        end do
        if(idwsdeb.gt.10)then
         write(6,*)('bmat1 for symmetry block '),isb,ibmat1
         call prntm2(bc(ibmat1),noc(isb),noc(isb),noc(isb))
         if(nsymb.eq.1)then
          call printa(bc(ibmat1),noc,0,1,0,noc,0,1,0,bc(ibcoff))
         end if
         write(6,*)('bmat2 for symmetry block'),isb,ibmat2
         iarg1=nvirt (isb)                                              8d31s15
         iarg2=noc(isb)
         call prntm2(bc(ibmat2),iarg1,iarg2,iarg1)
         if(nsymb.eq.1)then
          call printa(bc(ibmat2),nvirt ,noc,1,0,noc,0,1,0,bc(ibcoff))
         end if
        end if
        iamat=iamat+noc(isb)*nbasdws(isb)                               2d1s19
        call msl('bmat',bc(ibmat1),noc(isb),nbasdws(isb))               2d1s19
       end if                                                           8d31s04
       nwds=nbasdws(isb)*nbasdws(isb)                                   2d1s19
       jh0mo=jh0mo+nwds                                                 8d31s04
      end do                                                            8d31s04
      if(idwsdeb.gt.10)write(6,*)('after bmat')
      idwsdeb=00
c
c     gmat. distributed across processors
c     dimension is not same as exchange integrals
c
c     note: for optimization step, we only need occ*virt part of g,
c     however for Werner method, we need to transform g, thus we need
c     full occ*nbas part of g. Thus we will distribute and store
c     full g, but then copy only occ*virt part into ready memory.
c
c     this strategy is
c     on shared memory machine, that is a good strategy. However on
c     numa machine, we will store occ*occ part separately from occ*virt
c     part.
c
c
      call second(time1)
      telap=time1-time2-tovr
      times(16,1)=times(16,1)+telap
      times(16,2)=times(16,2)+telap**2
      times(16,3)=times(16,3)+1d0
      if(idwsdeb.gt.10)write(6,*)('gmat!')
      jh0mo=ih0mo                                                       8d31s04
      do isb=1,nsymb                                                    8d31s04
       do isa=1,isb                                                     2d5s14
        igmat(isa,isb)=-1                                               2d5s14
        igmat(isb,isa)=-1                                               2d5s14
        if(noc(isb)*noc(isa).gt.0)then                                  9d15s04
         n1=nvirt (isb)                                                 8d31s15
         n2=nvirt (isa)                                                 8d31s15
         call ilimts(n1,n2,numpro,nowpro,ilj,ihj,i1s,i1e,i2s,i2e)        8d31s04
         nhere=ihj+1-ilj                                                3d9s12
         nk1=noc(isa)*noc(isb)                                          3d14s12
         if(idwsdeb.gt.10)then
         end if
         jsze=nk1*nhere                                                 3d9s12
         igmat(isa,isb)=ibcoff                                           9d15s04
         ibcoff=igmat(isa,isb)+jsze                                      9d15s04
         call enough('parahf. 46',bc,ibc)                               3d13s23
         jgmat=igmat(isa,isb)                                           3d9s12
         i10=i1s                                                        3d9s12
         i1n=n1                                                         3d9s12
         do i2=i2s,i2e                                                  3d9s12
          i2p=i2+noc(isa)                                               2d4s14
          if(i2.eq.i2e)i1n=i1e                                          3d9s12
          do i1=i10,i1n                                                 3d9s12
           i1p=i1+noc(isb)                                              2d4s14
           if(isa.eq.isb)then                                            9d15s04
            jdwsk=jh0mo+i1p-1+(i2p-1)*nbasdws(isb)                      2d1s19
            do i3=1,idbly(isa)                                          4d4s17
             do i4=1,idbly(isa)                                         4d4s17
              if(i4.eq.i3)then                                          4d4s17
               bc(jgmat)=bc(jdwsk)*2d0                                  4d4s17
              else                                                      4d4s17
               bc(jgmat)=0d0                                            4d4s17
              end if                                                    4d4s17
              jgmat=jgmat+1
             end do                                                     4d4s17
             do i4=1,iacto(isa)                                         4d4s17
              bc(jgmat)=0d0                                             4d4s17
              jgmat=jgmat+1                                             4d4s17
             end do                                                     4d4s17
            end do                                                      4d4s17
            jden1=iden1(isb)                                            4d4s17
            do i3=1,iacto(isa)                                          4d4s17
             do i4=1,idbly(isa)                                         4d4s17
              bc(jgmat)=0d0                                             4d4s17
              jgmat=jgmat+1                                             4d4s17
             end do                                                     4d4s17
             do i4=1,iacto(isa)                                         4d4s17
              bc(jgmat)=bc(jdwsk)*bc(jden1)                             4d4s17
              jgmat=jgmat+1                                             4d4s17
              jden1=jden1+1                                             4d4s17
             end do                                                     4d4s17
            end do                                                      4d4s17
           else                                                          9d15s04
            do i3=1,noc(isb)
             do i4=1,noc(isa)                                           9d15s04
              bc(jgmat)=0d0                                              9d15s04
              jgmat=jgmat+1                                              9d15s04
             end do                                                      9d15s04
            end do                                                       9d15s04
           end if                                                        9d15s04
          end do                                                        3d9s12
          i10=1                                                         3d9s12
         end do                                                         3d9s12
         fact=1d0                                                       9d15s04
         if(isa.gt.isb)fact=2d0                                         9d15s04
         if(idwsdeb.gt.10)then
         write(6,*)('gmat after h*Dij '),igmat(isa,isb)
         iarg1=nk1
         iarg2=(ihj+1-ilj)                                               8d31s04
         call prntm2(bc(igmat(isa,isb)),iarg1,iarg2,iarg1)
         end if
         ihit=0                                                         2d4s14
         do idws1=1,nsdlk                                                  8d31s04
          if(isblk(3,idws1).eq.isb.and.isblk(4,idws1).eq.isa)then       2d4s14
           if(isa.eq.isb)then                                           3d9s12
            ihit=1                                                      2d4s14
            if(noc(isblk(1,idws1)).gt.0)then                            3d9s12
             nrow=(noc(isblk(1,idws1))*(noc(isblk(1,idws1))+1))/2       3d9s12
             i10=i1s                                                    3d9s12
             i1n=nvirt (isa)                                            8d31s15
             nocp=noc(isa)+1                                            3d12s12
             jgmat=igmat(isa,isb)-nocp                                  3d9s12
             jj=jmats(idws1)-1                                          3d9s12
             do i2=i2s,i2e                                              3d9s12
              if(i2.eq.i2e)i1n=i1e                                      3d9s12
              do i1=i10,i1n                                             3d9s12
               sum=0d0                                                  3d12s12
               do k=1,idbly(isblk(1,idws1))                             4d4s17
                iad2=jj+((k*(k+1))/2)                                   3d9s12
                sum=sum+4d0*bc(iad2)                                    3d12s12
               end do                                                   3d9s12
               do k=1,idbly(isa)                                        4d4s17
                iad1=jgmat+k*nocp                                       3d9s12
                bc(iad1)=bc(iad1)+sum                                   3d12s12
               end do                                                   3d12s12
               jgmat=jgmat+nk1                                          3d9s12
               jj=jj+nrow                                               3d9s12
              end do                                                    3d9s12
              i10=1                                                     3d9s12
             end do                                                     3d9s12
            end if                                                      3d9s12
           end if                                                       3d9s12
           if(isblk(1,idws1).eq.isa.and.isblk(2,idws1).eq.isb)then      3d9s12
            if(isa.eq.isb)then                                          3d9s12
             ihit=2                                                     2d4s14
             nrow=(noc(isa)*(noc(isa)+1))/2                             3d9s12
             jgmat=igmat(isa,isb)-1-noc(isa)                            3d9s12
             jj=jmats(idws1)-1                                          3d9s12
             do i12=1,nhere                                             3d9s12
              do i3=1,noc(isa)                                          3d9s12
               do i4=1,noc(isa)                                         3d9s12
                in=min(i3,i4)                                           3d9s12
                ix=max(i3,i4)                                           3d9s12
                inx=((ix*(ix-1))/2)+in+jj                               3d9s12
                iad=jgmat+i4+i3*noc(isa)                                3d9s12
                bc(iad)=bc(iad)-bc(inx)*2d0                             3d9s12
               end do                                                   3d9s12
              end do                                                    3d9s12
              jgmat=jgmat+nk1                                           3d9s12
              jj=jj+nrow                                                3d9s12
             end do                                                     3d9s12
            else                                                        3d9s12
             ihit=3                                                     2d4s14
             do i1234=0,jsze-1                                          3d9s12
              bc(igmat(isa,isb)+i1234)=bc(igmat(isa,isb)+i1234)         3d9s12
     $             -bc(jmats(idws1)+i1234)*2d0                          3d9s12
             end do                                                     3d9s12
            end if                                                      3d9s12
           end if                                                       3d9s12
           if(isblk(1,idws1).eq.isb.and.isblk(2,idws1).eq.isa.and.      3d9s12
     $          isa.ne.isb)then                                         3d9s12
            jgmat=igmat(isa,isb)-1-noc(isa)                             3d9s12
            jj=jmats(idws1)-1-noc(isb)                                  3d9s12
            ihit=4                                                      2d4s14
            do i12=1,nhere                                              3d9s12
             do i3=1,noc(isb)                                           3d9s12
              do i4=1,noc(isa)                                          3d9s12
               iad1=jgmat+i4+i3*noc(isa)                                3d9s12
               iad2=jj+i3+i4*noc(isb)                                   3d9s12
               bc(iad1)=bc(iad1)-bc(iad2)*2d0                           3d9s12
              end do                                                    3d9s12
             end do                                                     3d9s12
             jgmat=jgmat+nk1                                            3d9s12
             jj=jj+nk1                                                  3d9s12
            end do                                                      3d9s12
           end if                                                       3d9s12
          else if(isblk(3,idws1).eq.isa.and.isblk(4,idws1).eq.isb)then  2d5s14
           if(isblk(1,idws1).eq.isa.and.isblk(2,idws1).eq.isb)then      3d9s12
            ihit=-4
           end if
          end if                                                        3d9s12
         end do                                                         3d9s12
         if(idwsdeb.gt.10)then
         write(6,*)('gmat after JP '),igmat(isa,isb)
         write(6,*)('ihit = '),ihit                                     2d4s14
         iarg1=nk1
         iarg2=(ihj+1-ilj)                                               8d31s04
         call prntm2(bc(igmat(isa,isb)),iarg1,iarg2,iarg1)
         end if
         ihit=0                                                         2d4s14
         do idws1=1,nsdlkk                                                  8d31s04
          if(isblkk(3,idws1).eq.isb.and.isblkk(4,idws1).eq.isa)then     2d4s14
           if(isa.eq.isb)then                                           3d9s12
            ihit=1
            if(noc(isblkk(1,idws1)).gt.0)then                            3d9s12
             nrow=noc(isblkk(1,idws1))**2                               3d12s12
             i10=i1s                                                    3d9s12
             i1n=nvirt (isb)                                            8d31s15
             nocp=noc(isa)+1                                            3d12s12
             nocpk=noc(isblkk(1,idws1))+1                               3d12s12
             jgmat=igmat(isa,isb)-nocp                                  3d9s12
             jj=kmats(idws1)-nocpk                                      3d12s12
             do i2=i2s,i2e                                              3d9s12
              if(i2.eq.i2e)i1n=i1e                                      3d9s12
              do i1=i10,i1n                                             3d9s12
               sum=0d0                                                  3d12s12
               do k=1,noc(isblkk(1,idws1))                               3d9s12
                iad2=jj+k*nocpk                                         3d12s12
                sum=sum+bc(iad2)                                        3d12s12
               end do                                                   3d9s12
               sum=-2d0*sum                                             3d12s12
               do k=1,noc(isa)                                          3d12s12
                iad1=jgmat+k*nocp                                       3d9s12
                bc(iad1)=bc(iad1)+sum                                   3d12s12
               end do                                                   3d12s12
               jgmat=jgmat+nk1                                          3d9s12
               jj=jj+nrow                                               3d9s12
              end do                                                    3d9s12
              i10=1                                                     3d9s12
             end do                                                     3d9s12
            end if                                                      3d9s12
           end if                                                       3d9s12
           if(isblkk(1,idws1).eq.isa.and.isblkk(2,idws1).eq.isb)then      3d9s12
           ihit=2                                                       2d4s14
             fuse=8d0                                                   3d12s12
            do i1234=0,jsze-1                                           3d12s12
             bc(igmat(isa,isb)+i1234)=bc(igmat(isa,isb)+i1234)          3d12s12
     $            +bc(kmats(idws1)+i1234)*fuse                          3d12s12
            end do                                                      3d12s12
           end if                                                       3d9s12
           if(isblkk(1,idws1).eq.isb.and.isblkk(2,idws1).eq.isa)then    3d12s12
           ihit=3                                                       2d4s14
             fuse=-2d0                                                  3d12s12
            jgmat=igmat(isa,isb)-1-noc(isa)                             3d9s12
            jj=kmats(idws1)-1-noc(isb)                                  3d9s12
            do i12=1,nhere                                              3d9s12
             do i3=1,noc(isb)                                           3d9s12
              do i4=1,noc(isa)                                          3d9s12
               iad1=jgmat+i4+i3*noc(isa)                                3d9s12
               iad2=jj+i3+i4*noc(isb)                                   3d9s12
               bc(iad1)=bc(iad1)+bc(iad2)*fuse                          3d12s12
              end do                                                    3d9s12
             end do                                                     3d9s12
             jgmat=jgmat+nk1                                            3d9s12
             jj=jj+nk1                                                  3d9s12
            end do                                                      3d9s12
           end if                                                       3d9s12
          else if(isblkk(3,idws1).eq.isa.and.isblkk(4,idws1).eq.isb)then2d5s14
           if(isblkk(1,idws1).eq.isa.and.isblkk(2,idws1).eq.isb)then    2d5s14
            ihit=-2
           end if
          end if                                                        3d9s12
         end do                                                         3d9s12
         if(idwsdeb.gt.10)then
          write(6,*)('ihit = '),ihit                                    2d4s14
         write(6,*)('gmat after 2KQ '),igmat(isa,isb)
         writE(6,*)('block '),isa,isb
         iarg1=nk1
         iarg2=(ihj+1-ilj)                                               8d31s04
         call prntm2(bc(igmat(isa,isb)),iarg1,iarg2,iarg1)
         end if
         call msl('gmat',bc(igmat(isa,isb)),nk1,ihj+1-ilj)
        end if                                                           8d31s04
       end do
       jh0mo=jh0mo+nbasdws(isb)*nbasdws(isb)                            2d1s19
      end do
      call second(time2)
      telap=time2-time1-tovr
      times(7,1)=times(7,1)+telap
      times(7,2)=times(7,2)+telap**2
      times(7,3)=times(7,3)+1d0
      timeah1=time2                                                     10d19s04
      end if                                                            4d23s18
      iter2=0
  402 continue
      maxd=50                                                           4d23s18
      tchange=0d0                                                       9d14s04
      if(icas.eq.0)then                                                 4d23s18
      do isb=1,nsymb                                                    9d8s04
       nv4(isb)=0                                                       4d1s22
       if(noc(isb).gt.0)then                                            9d8s04
        call ilimts(nvirt (isb),nvirt (isb),numpro,nowpro,ilj,ihj,i1s,  8d31s15
     $      i1e,i2s,i2e)                                                8d31s15
        ibmat2=ibmat1+noc(isb)*noc(isb)                                 9d9s04
c
c
c     diagonals of gmat in mo basis
c
        igdig(isb)=ibcoff                                               9d14s04
        ibcoff=igdig(isb)+nvirt (isb)*noc(isb)                          8d31s15
        call enough('parahf. 47',bc,ibc)                                3d13s23
        do idws=0,nvirt (isb)*noc(isb)-1                                8d31s15
         bc(igdig(isb)+idws)=0d0                                        3d12s12
        end do
        nrow=noc(isb)**2                                                3d12s12
        i10=i1s                                                         3d12s12
        i1n=nvirt (isb)                                                 8d31s15
        nocp=noc(isb)+1                                                 3d12s12
        ii=igmat(isb,isb)-nocp                                          3d12s12
        do i2=i2s,i2e                                                   3d12s12
         if(i2.eq.i2e)i1n=i1e                                           3d12s12
         do i1=i10,i1n                                                  3d12s12
          if(i1.eq.i2)then                                              3d12s12
           jgdig=igdig(isb)-1+noc(isb)*(i1-1)                           3d12s12
           do i4=1,noc(isb)                                             3d12s12
            bc(jgdig+i4)=bc(ii+i4*nocp)                                 3d12s12
           end do                                                       3d12s12
          end if                                                        3d12s12
          ii=ii+nrow                                                    3d12s12
         end do                                                         3d12s12
         i10=1                                                          3d12s12
        end do                                                          3d12s12
        if(idwsdeb.gt.10)then
         write(6,*)('this processors contribution to diag g: '),
     $        noc(isb),nvirt (isb)                                      8d31s15
         iarg1=noc(isb)
         iarg2=nvirt (isb)                                              8d31s15
         call prntm2(bc(igdig(isb)),iarg1,iarg2,iarg1)                  3d12s12
        end if
        call dws_sync                                                   3d12s12
        iarg1=noc(isb)*nvirt (isb)                                      8d31s15
        call dws_gsumf(bc(igdig(isb)),iarg1)                            3d12s12
        call dws_sync                                                   3d12s12
        if(idwsdeb.gt.10)then
         write(6,*)('after global sum ')
         iarg1=noc(isb)
         iarg2=nvirt (isb)                                              8d31s15
         call prntm2(bc(igdig(isb)),iarg1,iarg2,iarg1)                          9d9s04
         idwsdeb=0                                                      1d3s17
        end if
        nv4(isb)=noc(isb)*nvirt (isb)+1                                 8d31s15
        iyc(isb)=ibcoff
        ixyc(isb)=iyc(isb)+nv4(isb)
        ibig2(isb)=ixyc(isb)+nv4(isb)*maxd
        ibcoff=ibig2(isb)+nv4(isb)*maxd
        call enough('parahf. 48',bc,ibc)                                3d13s23
       end if                                                            9d14s04
      end do                                                            9d14s04
      iter2=iter2+1
      if(iter2.gt.10*nsymb)stop 'iter2'                                 9d14s04
      call second(timeq1)                                               3d7s07
      telap=timeq1-time2-tovr
      times(17,1)=times(17,1)+telap
      times(17,2)=times(17,2)+telap**2
      times(17,3)=times(17,3)+1d0
      end if                                                            4d23s18
      if(icas.ne.0)then                                                 4d23s18
       do isb=1,nsymb                                                   4d23s18
        nv4(isb)=idoub(isb)*iacto(isb)+noc(isb)*nvirt (isb)             4d23s18
        if(isb.eq.1)nv4(isb)=nv4(isb)+1                                 4d23s18
        iyc(isb)=ibcoff                                                 4d23s18
        ixyc(isb)=iyc(isb)+nv4(isb)                                     4d23s18
        ibig2(isb)=ixyc(isb)+nv4(isb)*maxd                              4d23s18
        ibcoff=ibig2(isb)+nv4(isb)*maxd                                 4d23s18
        call enough('parahf. 49',bc,ibc)                                3d13s23
       end do                                                           4d23s18
      end if                                                            4d23s18
      iaseward=ibcoff                                                   3d14s12
      iucpy=ibcoff                                                      8d5s14
      nucpy=0                                                           8d5s14
      do isb=1,nsymb                                                    8d5s14
       nucpy=nucpy+nbasdws(isb)*nbasdws(isb)                            2d1s19
      end do                                                            8d5s14
      ibcoff=iucpy+nucpy                                                8d5s14
      call enough('parahf. 50',bc,ibc)                                  3d13s23
      micit=0                                                           8d4s14
      ncas0=1                                                           4d5s18
      xlam0=xlam                                                        8d5s14
      xlam00=xlam0                                                      8d5s14
      ios=ibcoff                                                        2d26s19
      ior=ios+mynprocg*idbk                                             2d27s19
      ibcoff=ior+mynprocg*idbk                                          2d27s19
      call enough('parahf. 51',bc,ibc)                                  3d13s23
      ibctop=ibcoff                                                     3d30s18
      deltaeci=1d10                                                     4d3s18
      npullback=0                                                       4d23s18
      if(ifirsto.ne.0.and.iter1.le.ifirsto)then                         8d27s19
       iforcef=0                                                         8d26s19
      else                                                              8d27s19
       iforcef=10                                                       8d27s19
      end if                                                            8d27s19
      idwsdeb=0
      ngodown=0                                                         2d25s19
 2020 continue                                                          8d4s14
      ibcoff=ibctop                                                     3d30s18
      micit=micit+1                                                     8d4s14
      if(iprtr(5).ne.0)write(6,*)('for micro iteration '),micit,xlam,   5d7s21
     $     ibcoff,ngodown                                                       2d25s19
      xlaml=xlam                                                        2d25s19
      if(micit.gt.50)then
            do isb=1,nsymb                                                1d3s20
             if(iprtr(5).ne.0.and.nbasdws(isb).gt.0)then                5d7s21
              write(6,*)('for symmetry block '),isb                     5d7s21
              write(6,*)('itotmb is ')                                  5d7s21
              call prntm2(bc(itotmb(isb)),nbasdws(isb),nbasdws(isb),    5d7s21
     $             nbasdws(isb))                                        5d7s21
             end if                                                     5d7s21
             do i=0,nbasdws(isb)*nbasdws(isb)-1                          1d8s19
              bc(iumat(isb)+i)=bc(itotmb(isb)+i)                        1d11s20
             end do                                                       1d3s20
            end do                                                        1d3s20
       go to 2021
      end if                                                            1d11s20
      call opto(noc,iyc,ixyc,nv4,ibig2,ibmat10,numpro,nowpro,igmat,     3d12s12
     $     times,lprint,maxd,iorbn,iumat,iamat10,itmat,mgmat,           8d5s14
     $     iorbx,ih0mo,tovr,tchange,i404,ipk,jmats,iqk,kmats,iden1,
     $     igdig,iaseward,thrs,thrsa,xlam,idwsdeb,idavopt,nbasdws,      2d1s19
     $     nvirt ,eiglastgood,1,ihessa,ihessb,ihessc,ihessd,icas,       1d30s18
     $     idoub,iacto,iamatx,smallest,nlzz,iorbsym,mlx,inbr,ipbr,bc,   11d9s22
     $     ibc)                                                         11d9s22
      if(icas.ne.0)then                                                 1d3s18
       call second(time1)                                               4d5s18
       call updateg(c0,c1,cjk,c3,ih0mo,ioooo,ionex,jmats,kmats,i3x,ih02, 1d9s18
     $      ioooo2,idoub,iacto,noc,nvirt ,nbasdws,iamatx,iden1,j2den,   2d1s19
     $      itmat,iumat,itotm,rmssz,iter1,idwsdeb,ibc(ios),ibc(ior),    2d27s19
     $      mynprocg,bc,ibc)                                            11d9s22
       call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk,nsdlk,    2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo2,ehf,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
        ediff=ehf-ehfbsofar                                              4d9s18
        call dws_bcast(ediff,1)                                         2d18s20
       call second(time2)                                               4d5s18
       telap=time2-time1-tovr                                           4d5s18
       times(6,1)=times(6,1)+telap                                      4d5s18
       times(6,2)=times(6,2)+telap**2                                   4d5s18
       times(6,3)=times(6,3)+1d0                                        4d5s18
       xx1=rmssz(1)                                                     4d4s18
       xx2=rmssz0(1)                                                    4d4s18
       epswfn=casdat(3)                                                 3d3s20
       esq=sqrt(epswfn)                                                 3d3s20
       iok=0                                                            3d3s20
       if(xx1.lt.epswfn)iok=1                                           3d3s20
       do isb=2,nsymb                                                   4d4s18
        xx1=max(xx1,rmssz(isb))                                         1d2s20
        if(rmssz(isb).lt.epswfn)iok=iok+1                               3d3s20
        xx2=max(xx2,rmssz0(isb))                                        1d2s20
       end do                                                           4d4s18
       iok=0
       if(iprtr(5).gt.0)write(6,*)('pre test iok: '),iok                3d3s20
       if(iok.eq.nsymb.and.micit.gt.1)then                              3d3s20
        if(iprtr(5).gt.0)write(6,*)('jump to converged ')               3d3s20
        do isb=1,nsymb                                                  3d3s20
         do i=0,nbasdws(isb)*nbasdws(isb)-1                             3d3s20
          bc(iumat(isb)+i)=bc(itotmb(isb)+i)                            3d3s20
         end do                                                         3d3s20
        end do                                                          3d3s20
        go to 2021                                                      3d3s20
       end if                                                           3d3s20
       sz=0d0                                                           12d28s19
       isz=0                                                            12d28s19
       do isb=1,nsymb                                                   12d28s19
        isz=isz+nbasdws(isb)*nbasdws(isb)                               12d28s19
        do i=0,nbasdws(isb)-1                                           12d28s19
         do j=0,i-1                                                     12d28s19
          ji=itotm(isb)+j+nbasdws(isb)*i                                12d28s19
          ij=itotm(isb)+i+nbasdws(isb)*j                                12d28s19
          sz=sz+bc(ji)**2+bc(ij)**2                                     12d28s19
         end do                                                         12d28s19
        end do                                                          12d28s19
       end do                                                           12d28s19
       sz=sqrt(sz/dfloat(isz))                                          12d28s19
       if(iprtr(5).ne.0)write(6,*)micit,xlam,xx1,ehf,xx2/xx1,ediff,sz      1d3s20
       if(idwsdeb.ge.1)write(6,*)('max grads: old '),xx2,(' new '),xx1                  4d4s18
       if(xx1.gt.xx2.and.ediff.lt.1d-10.and.iforcef.eq.0)then           5d2s18
        if(iprtr(5).ne.0)write(6,*)('force moving ahead ')              2d13s20
        xx2=2d0*xx1                                                     5d2s18
        iforcef=1                                                       5d2s18
       end if
       call dws_bcast(xx1,1)                                            2d18s20
       call dws_bcast(xx2,1)                                            2d18s20
       if(xx1.lt.xx2*0.999d0.or.xlam.gt.5d1)then                        1d11s20
        ngodown=ngodown+1                                               2d25s19
        call second(timedown)
        if(iprtr(5).ne.0)                                               1d3s20
     $   write(6,*)('excellent! the size of the gradient has decreased')1d3s20
        npullback=0                                                     4d23s18
        if(iprtr(5).ne.0)write(6,*)('ediff,deltaeci '),ediff,deltaeci
          if(ediff.gt.1d-9)then                                         11d17s20
           if(iprtr(5).ne.0)write(6,*)
     $         ('energy went up - increase xlam with best previous')
           xlam=xlam*3d0                                                1d13s20
            do isb=1,nsymb                                                1d3s20
           do i=0,2*nbasdws(isb)*nbasdws(isb)+noc(isb)*noc(isb)-1       1d3s20
              bc(iumat(isb)+i)=bc(iumatb(isb)+i)                        1d11s20
             end do                                                       1d3s20
            end do                                                        1d3s20
            go to 2020                                                  1d11s20
          end if                                                        1d11s20
        xx1b=xx1                                                        2d25s19
        xlamb=xlam                                                      2d25s19
        do isb=1,nsymb                                                  4d2s18
         rmssz0(isb)=min(rmssz0(isb),rmssz(isb))                        12d28s19
        end do                                                          4d2s18
        ehfbsofar=ehf                                                   3d30s18
        ehfbsofarl=ehf                                                  4d5s18
        epswfn=casdat(3)                                                5d1s18
        esq=sqrt(epswfn)                                                4d3s18
        iok=0
        do isb=1,nsymb                                                  4d2s18
         if(rmssz0(isb).le.epswfn)iok=iok+1                             4d3s18
        end do                                                          4d2s18
        if(xx20o.le.epswfn)then                                         1d14s20
         iok=nsymb                                                      1d14s20
         xlamb=1d0                                                      1d14s20
        end if                                                          1d14s20
        if(iok.ne.nsymb.or.abs(ediff).gt.casdat(2).or.micit.le.2)then   2d10s20
         iok=0                                                           4d3s18
         do isb=1,nsymb                                                  4d3s18
          if(rmssz0(isb).lt.esq)iok=iok+1                                4d3s18
         end do                                                          4d3s18
         if(iprtr(5).ne.0)write(6,*)('iok,ediff,deltaeci: '),iok,ediff,
     $        deltaeci
         if(iok.eq.nsymb.and.abs(ediff).lt.abs(deltaeci))then           4d3s18
          ncas0=ncas0+1                                                 4d5s18
          call second(time1)                                            4d5s18
          call cas0(ih02,ioooo2,noc,numpro,nowpro,potdws,                  8d8s14
     $      icallcas,iden1,j2den,inewl,nlzz,multh,ehf,ncomp,enobreit,   4d6s18
     $      dynw,0,icasvec,lprint,eavg2,jmats,kmats,nvirt,ibc(isend),   2d25s19
     $      ibc(irecv),ibc(iptx),ibc(ipf),0,ixlzz,islz,iptr,            2d20s20
     $         idum,idum,idum,idum,ibasisp,icsfp,nfcnp,nctp,            8d2s22
     $         ipoint2p,idum,ibc(ismx),ibc(irelx),iorbn,morb,nbasisp,   8d2s22
     $         nryd,ixw1p,ixw2p,iptrbitp,iorbsym,iorbsymz,nextradata,      5d3s21
     $        extradata,1,idum,idum,idum,idum,dum,1,idum,idum,idum,idum,7d1s22
     $         idum,idum,idum,dum,dum,bc,ibc,dum,ibc(nveccnt),icanon,   10d11s24
     $     ncasvecsize)                                                 10d11s24
          if(iprtr(5).ne.0)write(6,*)('after cas0: '),ehf,eavg2
          call second(time2)                                            4d5s18
          telap=time2-time1-tovr
          times(5,1)=times(5,1)+telap
          times(5,2)=times(5,2)+telap**2
          times(5,3)=times(5,3)+1d0
          time1=time2                                                   4d5s18
          call updateg(c0,c1,cjk,c3,ih0mo,ioooo,ionex,jmats,kmats,i3x,
     $         ih02,ioooo2,idoub,iacto,noc,nvirt ,nbasdws,iamatx,iden1, 2d1s19
     $         j2den,itmat,iumat,itotm,rmssz,iter1,idwsdeb,ibc(ios),    2d27s19
     $         ibc(ior),mynprocg,bc,ibc)                                11d9s22
          call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk,nsdlk, 2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo2,ehf,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
          if(iprtr(5).ne.0)write(6,*)('ehf from genergy '),ehf
          call second(time2)                                               4d5s18
          telap=time2-time1-tovr                                           4d5s18
          times(6,1)=times(6,1)+telap                                      4d5s18
          times(6,2)=times(6,2)+telap**2                                   4d5s18
          times(6,3)=times(6,3)+1d0                                        4d5s18
          deltaeci=abs(ehf-ehfbsofar)                                   1d18s20
          call dws_bcast(deltaeci,1)                                    2d18s20
          if(iprtr(5).ne.0)write(6,*)('with new weights,'),             1d11s20
     $         (' change in e from new densities: '),deltaeci,ehf,      1d11s20
     $         ehf-ehf2o,ehfbsofar                                                1d11s20
          end if                                                        1d17s20
          defromb=ediff                                                 1d18s20
          if(iprtr(5).ne..0)write(6,*)('defromb: '),defromb             1d11s20
          if(defromb.gt.1d-5)then                                       1d18s20
           if(iprtr(5).ne.0)write(6,*)
     $         ('energy went up - increase xlam with best previous')
           xlam=xlam*3d0                                                1d13s20
            do isb=1,nsymb                                                1d3s20
           do i=0,2*nbasdws(isb)*nbasdws(isb)+noc(isb)*noc(isb)-1       1d3s20
              bc(iumat(isb)+i)=bc(iumatb(isb)+i)                        1d11s20
             end do                                                       1d3s20
            end do                                                        1d3s20
            go to 2020                                                  1d11s20
          end if                                                        1d11s20
          ehfrdb=ehf                                                    1d11s20
          defrombl=defromb
          ehfbsofar=ehf                                                    4d3s18
          ehfbsofarl=ehf                                                4d5s18
          xx1b=xx1                                                      2d25s19
          xlamb=xlam                                                    2d25s19
          xx2n=rmssz(1)                                                 1d2s20
          do isb=1,nsymb
           rmssz0(isb)=rmssz(isb)
           xx2n=max(xx2n,rmssz(isb))                                    1d2s20
          end do
          call dws_bcast(xx2n,1)                                        2d18s20
          if(xx2n.gt.xx20)then                                          1d3s20
           if(lprint.and.iprtr(5).ne.0)then                             1d9s18
            write(6,*)('gradient went up ! : '),xx2n,(' vs. '),xx20
            write(6,*)('and change in energy is '),deltaeci
            write(6,*)('what''s up with that? ')
            write(6,*)('do not pass go - go straight to next iteration')1d3s20
           end if                                                       1d3s20
           xlamb=xlam                                                   1d8s19
           xlam0=xlam                                                   1d8s19
           xlam00=xlam                                                  1d8s19
          end if                                                        1d2s20
          if(iprtr(5).ne.0)write(6,*)('new gradient '),xx2n             1d8s19
          xx20=xx2n                                                     1d8s19
          if(iprtr(5).ne.0)write(6,*)('copy umat to umatb')             1d3s20
          do isb=1,nsymb                                                1d3s20
           do i=0,2*nbasdws(isb)*nbasdws(isb)+noc(isb)*noc(isb)-1       1d3s20
            bc(iumatb(isb)+i)=bc(iumat(isb)+i)                          1d3s20
           end do                                                       1d3s20
          end do                                                        1d3s20
         go to 2020                                                     4d3s18
        end if
       else                                                             3d30s18
        npullback=npullback+1                                           4d23s18
        xlam=xlam*3d0                                                   4d5s18
        if(iprtr(5).ne.0)                                               1d9s18
     $       write(6,*)('where we have xlam size test: '),xlam,ngodown  1d9s18
        if(xlam.gt.5d2)then                                             1d9s18
         if(iprtr(5).ne.0)write(6,*)('xlam is really large, so provided'
     $       ),('gradient was computed correctly, energy will go down') 1d9s18
            do isb=1,nsymb                                                1d3s20
             do i=0,nbasdws(isb)*nbasdws(isb)-1                          1d8s19
              bc(iumat(isb)+i)=bc(itotm(isb)+i)                         1d8s19
             end do                                                       1d3s20
            end do                                                        1d3s20
           go to 2021                                                   1d3s20
        end if                                                          1d9s18
        if((xlam.gt.1d3.and.ngodown.gt.2).or.xlam.gt.1d6)then           2d25s19
         nbaspass=nbaspass+1
         if(iprtr(5).ne.0)write(6,*)('copy itotmb to iumat ')
         do isb=1,nsymb                                                 4d5s18
          do i=0,nbasdws(isb)*nbasdws(isb)-1                            2d1s19
           bc(iumat(isb)+i)=bc(itotmb(isb)+i)                           4d5s18
          end do                                                        4d5s18
         end do                                                         4d5s18
         go to 2021                                                     4d3s18
        end if                                                          3d30s18
        if(iprtr(5).ne.0)write(6,*)('copy umatb to umat')
        do isb=1,nsymb                                                  3d30s18
         do i=0,2*nbasdws(isb)*nbasdws(isb)+noc(isb)*noc(isb)-1         2d1s19
          bc(iumat(isb)+i)=bc(iumatb(isb)+i)                            3d30s18
         end do                                                         3d30s18
        end do                                                          3d30s18
        go to 2020                                                      3d30s18
       end if                                                           3d30s18
       if(iprtr(5).ne.0)write(6,*)('loading itotm into umat ')          1d3s20
       do isb=1,nsymb
        do i=0,nbasdws(isb)*nbasdws(isb)-1                              2d1s19
         bc(iumat(isb)+i)=bc(itotm(isb)+i)
        end do
       end do
       go to 2021
  212  continue
       ifcnp=ibcoff
       ifcnm=ifcnp+8
       ibcoff=ifcnm+8
       call enough('parahf. 52',bc,ibc)                                 3d13s23
       do isb=1,nsymb                                                   1d5s18
        write(6,*)('for symmetry block '),isb
        nn=noc(isb)*noc(isb)                                            1d5s18
        isav=ibcoff                                                     1d5s18
        ipert=isav+nn                                                   1d5s18
        ibcoff=ipert+nn
        call enough('parahf. 53',bc,ibc)                                3d13s23
        do i=0,nn-1                                                     1d5s18
         bc(isav+i)=bc(itmat(isb)+i)                                    1d5s18
        end do                                                          1d5s18
        nda=idoub(isb)*iacto(isb)
        idert=ibcoff                                                    1d5s18
        ibcoff=idert+nda*8                                             1d16s17
        call enough('parahf. 54',bc,ibc)                                3d13s23
        do i=0,nda*8-1
         bc(idert+i)=0d0
        end do
        do id=0,idoub(isb)-1                                            1d5s18
         do ia=0,iacto(isb)-1                                           1d5s18
          write(6,*)('numerical ders for ia = '),ia,(' id = '),id       1d5s18
          step=2d-4                                                     1d5s18
          wgt=-1d0/3d0                                                  1d5s18
          do ipass=1,2                                                  1d5s18
           write(6,*)('for ipass = '),ipass
           do i=0,nn-1                                                   1d5s18
            bc(ipert+i)=0d0                                              1d5s18
           end do                                                        1d5s18
           iad1=ipert+ia+idoub(isb)+noc(isb)*id                         1d5s18
           bc(iad1)=step                                                1d5s18
           iad1=ipert+id+noc(isb)*(idoub(isb)+ia)                       1d5s18
           bc(iad1)=-step                                               1d5s18
           call pertrb(bc(ipert),noc(isb))
           if(noc(isb).gt.0)then
           call dgemm('n','n',noc(isb),noc(isb),noc(isb),1d0,bc(isav),  1d5s18
     $          noc(isb),bc(ipert),noc(isb),0d0,bc(itmat(isb)),noc(isb),1d5s18
     d' parahf. 13')
           end if
           ibcup=ibcoff                                                 1d5s18
       call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,ionex,jmats,   1d3s18
     $     kmats,i3x,ih02,ioooo2,idoub,iacto,noc,nvirt ,multh)          1d3s18
       call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk,nsdlk,    2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo2,ehf,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
           ibcoff=ibcup
           do i=1,8
            bc(ifcnp+i-1)=return(i)
           end do
           do i=0,nn-1                                                   1d5s18
            bc(ipert+i)=0d0                                              1d5s18
           end do                                                        1d5s18
           iad1=ipert+ia+idoub(isb)+noc(isb)*id                         1d5s18
           bc(iad1)=-step                                                1d5s18
           iad1=ipert+id+noc(isb)*(idoub(isb)+ia)                       1d5s18
           bc(iad1)=step                                                1d5s18
           call pertrb(bc(ipert),noc(isb))
           if(noc(isb).gt.0)then
           call dgemm('n','n',noc(isb),noc(isb),noc(isb),1d0,bc(isav),  1d5s18
     $          noc(isb),bc(ipert),noc(isb),0d0,bc(itmat(isb)),noc(isb),1d5s18
     d' parahf. 14')
           end if
           ibcup=ibcoff                                                 1d5s18
       call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,ionex,jmats,   1d3s18
     $     kmats,i3x,ih02,ioooo2,idoub,iacto,noc,nvirt ,multh)          1d3s18
       call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk,nsdlk,    2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo2,ehf,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
           ibcoff=ibcup
           do i=1,8
            der=0.5d0*(bc(ifcnp+i-1)-return(i))/step                    1d5s18
            write(6,*)i,der
            iad=idert+ia+iacto(isb)*id+nda*(i-1)
            bc(iad)=bc(iad)+wgt*der
           end do
           step=step*0.5d0                                              1d5s18
           wgt=4d0/3d0                                                  1d5s18
          end do                                                        1d5s18
         end do                                                         1d5s18
        end do                                                          1d5s18
        if(nda.gt.0)then
         write(6,*)('for symmetry block '),isb
         do i=0,7
          ip=i+1
          write(6,*)('der matrix for term no. '),ip
          iad=idert+nda*i
          call prntm2(bc(iad),iacto(isb),idoub(isb),iacto(isb))
         end do
         do i=0,nn-1
          bc(itmat(isb)+i)=bc(isav+i)
         end do
        end if
        ibcoff=isav
        nn=nbasdws(isb)*nbasdws(isb)                                    2d1s19
        isav=ibcoff                                                     1d5s18
        ipert=isav+nn                                                   1d5s18
        ibcoff=ipert+nn
        call enough('parahf. 55',bc,ibc)                                3d13s23
        do i=0,nn-1                                                     1d5s18
         bc(isav+i)=bc(iumat(isb)+i)                                    1d5s18
        end do                                                          1d5s18
        nda=noc(isb)*nvirt (isb)
        idert=ibcoff                                                    1d5s18
        ibcoff=idert+nda*8
        call enough('parahf. 56',bc,ibc)                                3d13s23
        do i=0,nda*8-1
         bc(idert+i)=0d0
        end do
        do id=0,noc(isb)-1                                              1d5s18
         do ia=0,nvirt (isb)-1                                           1d5s18
          write(6,*)('numerical ders for iv = '),ia,(' io = '),id       1d5s18
          step=2d-4                                                     1d5s18
          wgt=-1d0/3d0                                                  1d5s18
          do ipass=1,2                                                  1d5s18
           write(6,*)('for ipass = '),ipass
           do i=0,nn-1                                                   1d5s18
            bc(ipert+i)=0d0                                              1d5s18
           end do                                                        1d5s18
           iad1=ipert+ia+noc(isb)+nbasdws(isb)*id                       2d1s19
           bc(iad1)=step                                                1d5s18
           iad1=ipert+id+nbasdws(isb)*(noc(isb)+ia)                     2d1s19
           bc(iad1)=-step                                               1d5s18
           call pertrb(bc(ipert),nbasdws(isb))                          2d1s19
           if(nbasdws(isb).gt.0)then                                    2d25s19
           call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasdws(isb)    2d1s19
     $          ,1d0,bc(isav),nbasdws(isb),bc(ipert),nbasdws(isb),0d0,  2d1s19
     $          bc(iumat(isb)),nbasdws(isb),                            2d1s19
     d' parahf. 15')
           end if                                                       2d25s19
           ibcup=ibcoff                                                 1d5s18
           write(6,*)('calling updatei')
       call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,ionex,jmats,   1d3s18
     $     kmats,i3x,ih02,ioooo2,idoub,iacto,noc,nvirt ,multh)          1d3s18
       write(6,*)('back')
       call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk,nsdlk,    2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo2,ehf,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
           ibcoff=ibcup
           do i=1,8
            bc(ifcnp+i-1)=return(i)
           end do
           do i=0,nn-1                                                   1d5s18
            bc(ipert+i)=0d0                                              1d5s18
           end do                                                        1d5s18
           iad1=ipert+ia+noc(isb)+nbasdws(isb)*id                       2d1s19
           bc(iad1)=-step                                                1d5s18
           iad1=ipert+id+nbasdws(isb)*(noc(isb)+ia)                     2d1s19
           bc(iad1)=step                                                1d5s18
           call pertrb(bc(ipert),nbasdws(isb))                          2d1s19
           if(nbasdws(isb).gt.0)then                                    2d25s19
           call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasdws(isb),   2d1s19
     $          1d0,bc(isav),nbasdws(isb),bc(ipert),nbasdws(isb),0d0,   2d1s19
     $          bc(iumat(isb)),nbasdws(isb),                            2d1s19
     d' parahf. 16')
           end if                                                       2d25s19
           ibcup=ibcoff                                                 1d5s18
       call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,ionex,jmats,   1d3s18
     $     kmats,i3x,ih02,ioooo2,idoub,iacto,noc,nvirt ,multh)          1d3s18
       call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk,nsdlk,    2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo2,ehf,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
           ibcoff=ibcup
           do i=1,8
            der=0.5d0*(bc(ifcnp+i-1)-return(i))/step                    1d5s18
            write(6,*)i,der
            iad=idert+ia+nvirt (isb)*id+nda*(i-1)
            bc(iad)=bc(iad)+wgt*der
           end do
           step=step*0.5d0                                              1d5s18
           wgt=4d0/3d0                                                  1d5s18
          end do                                                        1d5s18
         end do                                                         1d5s18
        end do                                                          1d5s18
        if(nda.gt.0)then
         write(6,*)('for symmetry block '),isb
         do i=0,7
          ip=i+1
          write(6,*)('der matrix for term no. '),ip
          iad=idert+nda*i
          call prntm2(bc(iad),nvirt (isb),noc(isb),nvirt (isb))
         end do
         do i=0,nn-1
          bc(iumat(isb)+i)=bc(isav+i)
         end do
        end if
        ibcoff=isav
       end do                                                           1d5s18
       write(6,*)('now for second derivatives ...')
       do isb=1,nsymb                                                   5d9s22
        isavu(isb)=ibcoff                                               5d9s22
        ibcoff=isavu(isb)+nbasdws(isb)**2                               5d9s22
        call enough('parahf. 57',bc,ibc)                                3d13s23
        do i=0,nbasdws(isb)**2-1                                        5d9s22
         bc(isavu(isb)+i)=bc(iumat(isb)+i)                              5d9s22
        end do                                                          5d9s22
       end do                                                           5d9s22
       ifcn0=ibcoff
       ibcoff=ifcn0+8
       ifcnpp=ibcoff                                                    5d9s22
       ifcnmp=ifcnpp+8                                                  5d9s22
       ifcnpm=ifcnmp+8                                                  5d9s22
       ifcnmm=ifcnpm+8                                                  5d9s22
       ibcoff=ifcnmm+8                                                  5d9s22
       call enough('parahf. 58',bc,ibc)                                 3d13s23
       ibcup=ibcoff                                                     5d9s22
       call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,ionex,jmats,   1d3s18
     $     kmats,i3x,ih02,ioooo2,idoub,iacto,noc,nvirt ,multh)          1d3s18
       call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk,nsdlk,    2d1s19
     $     j2den,potdws,                                                5d19s16
     $      numpro,nowpro,ioooo2,ehf,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
       ibcoff=ibcup
       do i=1,8
        bc(ifcn0+i-1)=return(i)*2d0                                     5d9s22
       end do
       do isk=1,nsymb
        ipertk=ibcoff                                                   5d9s22
        nnk=nbasdws(isk)*nbasdws(isk)                                   5d9s22
        ibcoff=ipertk+nnk
        call enough('parahf. 59',bc,ibc)                                3d13s23
        do isb=1,isk
         nnb=nbasdws(isb)**2                                            5d9s22
         ipertb=ibcoff                                                  5d9s22
         idera=ipertb+nnb                                               5d9s22
         ndera=idoub(isb)*idoub(isk)*iacto(isb)*iacto(isk)              5d9s22
         iderb=idera+ndera*8                                            5d9s22
         nderb=idoub(isb)*noc(isk)*iacto(isb)*nvirt(isk)                5d9s22
         iderc=iderb+nderb*8                                            5d9s22
         nderc=noc(isb)*noc(isk)*nvirt(isb)*nvirt(isk)                  5d9s22
         ibcoff=iderc+nderc*8                                           5d9s22
         call enough('parahf. 60',bc,ibc)                               3d13s23
         do iz=idera,ibcoff-1                                           5d9s22
          bc(iz)=0                                                      5d9s22
         end do                                                         5d9s22
         if(min(idoub(isb),idoub(isk),iacto(isb),iacto(isk)).gt.0)then  5d9s22
          write(6,*)('dera for syms '),isb,isk
          jdera=idera
          do i2=0,iacto(isk)-1
           i2p=i2+idoub(isk)
           do i1=0,iacto(isb)-1
            i1p=i1+idoub(isb)
            do i4=0,idoub(isb)-1
             do i3=0,idoub(isk)-1
              write(6,*)('for '),i3,i4,i1,i2
              if(isb.ne.isk.or.i1.ne.i2.or.i4.ne.i3)then
               step=2d-3                                                     1d5s18
               wgt=-1d0/3d0                                                  1d5s18
               do ipass=1,2
                do idptype=1,4
                 ifsav=ifcnpp+8*(idptype-1)-1                           5d9s22
                 do i=0,nnb-1
                  bc(ipertb+i)=0d0
                 end do
                 ij=ipertb+i4+nbasdws(isb)*i1p
                 ji=ipertb+i1p+nbasdws(isb)*i4
                 bc(ij)=-step*dptype(idptype,1)
                 bc(ji)=step*dptype(idptype,1)
                 if(isb.eq.isk)then
                  ij=ipertb+i3+nbasdws(isk)*i2p
                  ji=ipertb+i2p+nbasdws(isk)*i3
                  bc(ij)=-step*dptype(idptype,2)
                  bc(ji)=step*dptype(idptype,2)
                 end if
                 call pertrb(bc(ipertb),nbasdws(isb))
                 call dgemm('n','n',nbasdws(isb),nbasdws(isb),
     $                nbasdws(isb),1d0,bc(isavu(isb)),nbasdws(isb),
     $                bc(ipertb),nbasdws(isb),0d0,bc(iumat(isb)),
     $                nbasdws(isb),
     d' parahf. ppb')
                 if(isb.ne.isk)then
                  do i=0,nnk-1
                   bc(ipertk+i)=0d0
                  end do
                  ij=ipertk+i3+nbasdws(isk)*i2p
                  ji=ipertk+i2p+nbasdws(isk)*i3
                  bc(ij)=-step*dptype(idptype,2)
                  bc(ji)=step*dptype(idptype,2)
                  call pertrb(bc(ipertk),nbasdws(isk))
                  call dgemm('n','n',nbasdws(isk),nbasdws(isk),
     $               nbasdws(isk),1d0,bc(isavu(isk)),nbasdws(isk),
     $               bc(ipertk),nbasdws(isk),0d0,bc(iumat(isk)),
     $               nbasdws(isk),                                      5d9s22
     d' parahf. ppk')
                 end if
                 ibcup=ibcoff                                                     5d9s22
                 call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,     5d9s22
     $                ionex,jmats,kmats,i3x,ih02,ioooo2,idoub,iacto,noc,5d9s22
     $                nvirt ,multh)                                     5d9s22
                 call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk, 5d9s22
     $               nsdlk,j2den,potdws,numpro,nowpro,ioooo2,ehf,ncomp, 5d9s22
     $               idwsdeb,enobreit,lprint,iconv,bc,ibc)              11d11s22
                 ibcoff=ibcup
                 do i=1,8                                               5d9s22
                  bc(ifsav+i)=return(i)                                 5d9s22
                 end do                                                 5d9s22
                end do
                ff=0.25d0/(step*step)                                   5d9s22
                kdera=jdera                                             5d9s22
                do i=0,7                                                5d9s22
                 der=(bc(ifcnpp+i)+bc(ifcnmm+i)-bc(ifcnmp+i)
     $                -bc(ifcnpm+i))*ff
                 bc(ifcnpp+i)=der
                 bc(kdera)=bc(kdera)+der*wgt                            5d9s22
                 bc(ifcnpm+i)=bc(kdera)                                 5d9s22
                 kdera=kdera+ndera                                      5d9s22
                end do                                                  5d9s22
                write(6,1852)(bc(ifcnpp+i),i=0,7)
 1852           format(8es15.7,a1)
                if(ipass.eq.2)write(6,1852)(bc(ifcnpm+i),i=0,7)         5d9s22
                step=step*0.5d0                                              1d5s18
                wgt=4d0/3d0                                                  1d5s18
               end do
              else
               step=2d-3                                                     1d5s18
               wgt=-1d0/3d0                                                  1d5s18
               do ipass=1,2
                do idptype=1,2
                 ifsav=ifcnpp+8*(idptype-1)-1                           5d9s22
                 do i=0,nnb-1
                  bc(ipertb+i)=0d0
                 end do
                 ij=ipertb+i4+nbasdws(isb)*i1p
                 ji=ipertb+i1p+nbasdws(isb)*i4
                 bc(ij)=-step*dptype(idptype,1)
                 bc(ji)=step*dptype(idptype,1)
                 call pertrb(bc(ipertb),nbasdws(isb))
                 call dgemm('n','n',nbasdws(isb),nbasdws(isb),
     $                nbasdws(isb),1d0,bc(isavu(isb)),nbasdws(isb),
     $                bc(ipertb),nbasdws(isb),0d0,bc(iumat(isb)),
     $                nbasdws(isb),
     d' parahf. ppb')
                 ibcup=ibcoff                                                     5d9s22
                 call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,     5d9s22
     $                ionex,jmats,kmats,i3x,ih02,ioooo2,idoub,iacto,noc,5d9s22
     $                nvirt ,multh)                                     5d9s22
                 call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk, 5d9s22
     $               nsdlk,j2den,potdws,numpro,nowpro,ioooo2,ehf,ncomp, 5d9s22
     $               idwsdeb,enobreit,lprint,iconv,bc,ibc)              11d11s22
                 ibcoff=ibcup
                 do i=1,8                                               5d9s22
                  bc(ifsav+i)=return(i)                                 5d9s22
                 end do                                                 5d9s22
                end do
                ff=1d0/(step*step)
                kdera=jdera                                             5d9s22
                do i=0,7                                                5d9s22
                 der=(bc(ifcnpp+i)+bc(ifcnmp+i)-bc(ifcn0+i))*ff         5d9s22
                 bc(ifcnpp+i)=der
                 bc(kdera)=bc(kdera)+der*wgt                            5d9s22
                 bc(ifcnpm+i)=bc(kdera)                                 5d9s22
                 kdera=kdera+ndera                                      5d9s22
                end do                                                  5d9s22
                write(6,1852)(bc(ifcnpp+i),i=0,7)
                if(ipass.eq.2)write(6,1852)(bc(ifcnpm+i),i=0,7)         5d9s22
                step=step*0.5d0                                              1d5s18
                wgt=4d0/3d0                                                  1d5s18
               end do
              end if
              jdera=jdera+1
             end do
            end do
           end do
          end do
          nrowa=idoub(isb)*idoub(isk)
          ncola=iacto(isb)*iacto(isk)
          do i=1,8
           write(6,*)('dera for fcn '),i
           jdera=idera+ndera*(i-1)
           call prntm2(bc(jdera),nrowa,ncola,nrowa)
          end do
         end if                                                         5d9s22
         if(min(idoub(isb),noc(isk),iacto(isb),nvirt(isk)).gt.0)then    5d9s22
          write(6,*)('derb for syms '),isb,isk
          jderb=iderb
          do i2=0,nvirt(isk)-1
           i2p=i2+noc(isk)
           do i1=0,noc(isk)-1
            do i4=0,idoub(isb)-1
             do i3=0,iacto(isb)-1
              i3p=i3+idoub(isb)
              write(6,*)('for '),i3,i4,i1,i2
              step=2d-3                                                     1d5s18
              wgt=-1d0/3d0                                                  1d5s18
              do ipass=1,2
               do idptype=1,4
                ifsav=ifcnpp+8*(idptype-1)-1                            5d9s22
                do i=0,nnb-1
                 bc(ipertb+i)=0d0
                end do
                ij=ipertb+i4+nbasdws(isb)*i3p
                ji=ipertb+i3p+nbasdws(isb)*i4
                bc(ij)=-step*dptype(idptype,1)
                bc(ji)=step*dptype(idptype,1)
                if(isb.eq.isk)then
                 ij=ipertb+i1+nbasdws(isk)*i2p
                 ji=ipertb+i2p+nbasdws(isk)*i1
                 bc(ij)=-step*dptype(idptype,2)
                 bc(ji)=step*dptype(idptype,2)
                end if
                call pertrb(bc(ipertb),nbasdws(isb))
                call dgemm('n','n',nbasdws(isb),nbasdws(isb),
     $               nbasdws(isb),1d0,bc(isavu(isb)),nbasdws(isb),
     $               bc(ipertb),nbasdws(isb),0d0,bc(iumat(isb)),
     $               nbasdws(isb),
     d' parahf. ppb')
                if(isb.ne.isk)then
                 do i=0,nnk-1
                  bc(ipertk+i)=0d0
                 end do
                 ij=ipertk+i1+nbasdws(isk)*i2p
                 ji=ipertk+i2p+nbasdws(isk)*i1
                 bc(ij)=-step*dptype(idptype,2)
                 bc(ji)=step*dptype(idptype,2)
                 call pertrb(bc(ipertk),nbasdws(isk))
                 call dgemm('n','n',nbasdws(isk),nbasdws(isk),
     $               nbasdws(isk),1d0,bc(isavu(isk)),nbasdws(isk),
     $               bc(ipertk),nbasdws(isk),0d0,bc(iumat(isk)),
     $               nbasdws(isk),                                      5d9s22
     d' parahf. ppk')
                end if
                ibcup=ibcoff                                                     5d9s22
                call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,      5d9s22
     $                ionex,jmats,kmats,i3x,ih02,ioooo2,idoub,iacto,noc,5d9s22
     $                nvirt ,multh)                                     5d9s22
                call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk, 5d9s22
     $               nsdlk,j2den,potdws,numpro,nowpro,ioooo2,ehf,ncomp, 5d9s22
     $               idwsdeb,enobreit,lprint,iconv,bc,ibc)              11d11s22
                ibcoff=ibcup
                do i=1,8                                                5d9s22
                 bc(ifsav+i)=return(i)                                  5d9s22
                end do                                                  5d9s22
               end do
               ff=0.25d0/(step*step)                                    5d9s22
               kderb=jderb                                              5d9s22
               do i=0,7                                                 5d9s22
                der=(bc(ifcnpp+i)+bc(ifcnmm+i)-bc(ifcnmp+i)
     $                -bc(ifcnpm+i))*ff
                bc(ifcnpp+i)=der
                bc(kderb)=bc(kderb)+der*wgt                             5d9s22
                bc(ifcnpm+i)=bc(kderb)                                  5d9s22
                kderb=kderb+nderb                                       5d9s22
               end do                                                   5d9s22
               write(6,1852)(bc(ifcnpp+i),i=0,7)
               if(ipass.eq.2)write(6,1852)(bc(ifcnpm+i),i=0,7)          5d9s22
               step=step*0.5d0                                              1d5s18
               wgt=4d0/3d0                                                  1d5s18
              end do
              jderb=jderb+1
             end do
            end do
           end do
          end do
          nrowb=idoub(isb)*iacto(isb)
          ncolb=noc(isk)*nvirt(isk)
          do i=1,8
           write(6,*)('derb for fcn '),i
           jderb=iderb+nderb*(i-1)
           call prntm2(bc(jderb),nrowb,ncolb,nrowb)
          end do
         end if                                                         5d9s22
         if(min(noc(isb),noc(isk),nvirt(isb),nvirt(isk)).gt.0)then      5d9s22
          write(6,*)('derc for syms '),isb,isk
          jderc=iderc
          do i1=0,nvirt(isb)-1
           i1p=i1+noc(isb)
           do i2=0,nvirt(isk)-1
            i2p=i2+noc(isk)
            do i3=0,noc(isk)-1
             do i4=0,noc(isb)-1
              write(6,*)('for '),i3,i4,i1,i2
              lbug=.false.
              if(isb.ne.isk.or.i1.ne.i2.or.i4.ne.i3)then
               step=2d-3                                                     1d5s18
               wgt=-1d0/3d0                                                  1d5s18
               do ipass=1,2
                do idptype=1,4
                 if(lbug)write(6,*)('for idptype '),idptype
                 ifsav=ifcnpp+8*(idptype-1)-1                           5d9s22
                 do i=0,nnb-1
                  bc(ipertb+i)=0d0
                 end do
                 ij=ipertb+i4+nbasdws(isb)*i1p
                 ji=ipertb+i1p+nbasdws(isb)*i4
                 bc(ij)=-step*dptype(idptype,1)
                 bc(ji)=step*dptype(idptype,1)
                 if(isb.eq.isk)then                                     5d13s22
                  ij=ipertb+i3+nbasdws(isk)*i2p
                  ji=ipertb+i2p+nbasdws(isk)*i3
                  bc(ij)=-step*dptype(idptype,2)
                  bc(ji)=step*dptype(idptype,2)
                 end if                                                 5d13s22
                 call pertrb(bc(ipertb),nbasdws(isb))
                 if(lbug)then
                  write(6,*)('pertb ')
                  call prntm2(bc(ipertb),nbasdws(isb),nbasdws(isb),
     $                nbasdws(isb))
                 end if
                 call dgemm('n','n',nbasdws(isb),nbasdws(isb),
     $                nbasdws(isb),1d0,bc(isavu(isb)),nbasdws(isb),
     $                bc(ipertb),nbasdws(isb),0d0,bc(iumat(isb)),
     $                nbasdws(isb),
     d' parahf. ppb')
                 if(lbug)then
                  write(6,*)('savu*pertb ')
                  call prntm2(bc(iumat(isb)),nbasdws(isb),nbasdws(isb),
     $                 nbasdws(isb))
                 end if
                 if(isb.ne.isk)then
                  do i=0,nnk-1
                   bc(ipertk+i)=0d0
                  end do
                  ij=ipertk+i3+nbasdws(isk)*i2p
                  ji=ipertk+i2p+nbasdws(isk)*i3
                  bc(ij)=-step*dptype(idptype,2)
                  bc(ji)=step*dptype(idptype,2)
                  call pertrb(bc(ipertk),nbasdws(isk))
                  if(lbug)then
                   write(6,*)('pertk ')
                   call prntm2(bc(ipertk),nbasdws(isk),nbasdws(isk),
     $                nbasdws(isk))
                  end if
                  call dgemm('n','n',nbasdws(isk),nbasdws(isk),
     $               nbasdws(isk),1d0,bc(isavu(isk)),nbasdws(isk),
     $               bc(ipertk),nbasdws(isk),0d0,bc(iumat(isk)),
     $               nbasdws(isk),                                      5d9s22
     d' parahf. ppk')
                 end if
                 if(lbug)then
                  write(6,*)('orbitals to use ')
                  call prntm2(bc(iumat(isk)),nbasdws(isk),nbasdws(isk),
     $                 nbasdws(isk))
                 end if
                 ibcup=ibcoff                                                     5d9s22
                 call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,     5d9s22
     $                ionex,jmats,kmats,i3x,ih02,ioooo2,idoub,iacto,noc,  5d9s22
     $                nvirt ,multh)                                     5d9s22
                 call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk, 5d9s22
     $               nsdlk,j2den,potdws,numpro,nowpro,ioooo2,ehf,ncomp, 5d9s22
     $               idwsdeb,enobreit,lprint,iconv,bc,ibc)              11d11s22
                 if(lbug)then
                  write(6,1852)return
                 end if
                 ibcoff=ibcup
                 do i=1,8                                               5d9s22
                  bc(ifsav+i)=return(i)                                 5d9s22
                 end do                                                 5d9s22
                end do
                ff=0.25d0/(step*step)                                   5d9s22
                kderc=jderc                                             5d9s22
                do i=0,7                                                5d9s22
                 der=(bc(ifcnpp+i)+bc(ifcnmm+i)-bc(ifcnmp+i)
     $                -bc(ifcnpm+i))*ff
                 bc(ifcnpp+i)=der
                 bc(kderc)=bc(kderc)+der*wgt                            5d9s22
                 bc(ifcnpm+i)=bc(kderc)                                 5d9s22
                 kderc=kderc+nderc                                      5d9s22
                end do                                                  5d9s22
                write(6,1852)(bc(ifcnpp+i),i=0,7),'a'
                if(ipass.eq.2)write(6,1852)(bc(ifcnpm+i),i=0,7),'a'         5d9s22
                step=step*0.5d0                                              1d5s18
                wgt=4d0/3d0                                                  1d5s18
               end do
              else
               step=2d-3                                                     1d5s18
               wgt=-1d0/3d0                                                  1d5s18
               do ipass=1,2
                do idptype=1,2
                 ifsav=ifcnpp+8*(idptype-1)-1                           5d9s22
                 do i=0,nnb-1
                  bc(ipertb+i)=0d0
                 end do
                 ij=ipertb+i4+nbasdws(isb)*i1p
                 ji=ipertb+i1p+nbasdws(isb)*i4
                 bc(ij)=-step*dptype(idptype,1)
                 bc(ji)=step*dptype(idptype,1)
                 call pertrb(bc(ipertb),nbasdws(isb))
                 call dgemm('n','n',nbasdws(isb),nbasdws(isb),
     $                nbasdws(isb),1d0,bc(isavu(isb)),nbasdws(isb),
     $                bc(ipertb),nbasdws(isb),0d0,bc(iumat(isb)),
     $                nbasdws(isb),
     d' parahf. ppb')
                 ibcup=ibcoff                                                     5d9s22
                 call updatei(c0,c1,cjk,c3,iumat,itmat,ih0mo,ioooo,     5d9s22
     $               ionex,jmats,kmats,i3x,ih02,ioooo2,idoub,iacto,noc,  5d9s22
     $                nvirt ,multh)                                     5d9s22
                 call genergy(nsymb,ih02,idbly,iact,nbasdws,iden1,isblk, 5d9s22
     $               nsdlk,j2den,potdws,numpro,nowpro,ioooo2,ehf,ncomp, 5d9s22
     $               idwsdeb,enobreit,lprint,iconv,bc,ibc)              11d11s22
                 ibcoff=ibcup
                 do i=1,8                                               5d9s22
                  bc(ifsav+i)=return(i)                                 5d9s22
                 end do                                                 5d9s22
                end do
                ff=1d0/(step*step)
                kderc=jderc                                             5d9s22
                do i=0,7                                                5d9s22
                 der=(bc(ifcnpp+i)+bc(ifcnmp+i)-bc(ifcn0+i))*ff         5d9s22
                 write(6,*)bc(ifcnpp+i),bc(ifcnmp+i),bc(ifcn0+i)*0.5d0
                 bc(ifcnpp+i)=der
                 bc(kderc)=bc(kderc)+der*wgt                            5d9s22
                 bc(ifcnpm+i)=bc(kderc)                                 5d9s22
                 kderc=kderc+nderc                                      5d9s22
                end do                                                  5d9s22
                write(6,1852)(bc(ifcnpp+i),i=0,7),'s'
                if(ipass.eq.2)write(6,1852)(bc(ifcnpm+i),i=0,7),'s'         5d9s22
                step=step*0.5d0                                              1d5s18
                wgt=4d0/3d0                                                  1d5s18
               end do
              end if
              jderc=jderc+1
             end do
            end do
           end do
          end do
          nrowc=noc(isb)*noc(isk)
          ncolc=nvirt(isb)*nvirt(isk)
          do i=1,8
           write(6,*)('derc for fcn '),i
           jderc=iderc+nderc*(i-1)
           call prntm2(bc(jderc),nrowc,ncolc,nrowc)
          end do
         end if                                                         5d9s22
         ibcoff=ipertb
        end do
        ibcoff=ipertk                                                   5d9s22
       end do
       if(ibcoff.ne.-132)then
        call dws_sync
        call dws_finalize
        stop
       end if
      else                                                              1d3s18
      call secoi(iptoh,noc,multh,iumat,ioooo,itmat,ih0mo,ih02,ioooo2,   8d31s15
     $     nbasdws,bc,ibc)                                              11d9s22
      call genergy(nsymb,ih02,idbly,iact,idbly,iden1,isblk,nsdlk,       3d3s21
     $     j2den,0d0,                                                   3d3s21
     $      numpro,nowpro,ioooo2,ehf2nd,ncomp,idwsdeb,enobreit,lprint,     1d5s18
     $      iconv,bc,ibc)                                               11d11s22
      ehfbsofarl=ehf2nd+potdws                                          4d24s18
      end if                                                            1d3s18
      if(micit.eq.1)then                                                8d5s14
       ehfmicl=ehf2nd                                                   8d5s14
       jucpy=iucpy                                                      8d5s14
       do isb=1,nsymb                                                   8d5s14
        do i=0,nbasdws(isb)*nbasdws(isb)-1                              2d1s19
         bc(jucpy+i)=bc(iumat(isb)+i)                                   8d5s14
        end do                                                          8d5s14
        jucpy=jucpy+nbasdws(isb)**2                                     2d1s19
       end do                                                           8d5s14
      end if                                                            8d5s14
      echan=ehf2nd-ehfmicl                                              8d5s14
      call dws_bcast(echan,1)                                           8d5s14
      if(lprint)then
      end if
      if(echan.gt.1d-13)then                                            8d5s14
       write(6,*)('diverging ... ')                                     8d5s14
       write(6,*)('re-load last ok umat ')
       jucpy=iucpy                                                      8d5s14
       do isb=1,nsymb                                                   8d5s14
        do i=0,nbasdws(isb)*nbasdws(isb)-1                              2d1s19
         bc(iumat(isb)+i)=bc(jucpy+i)                                   8d5s14
        end do                                                          8d5s14
        jucpy=jucpy+nbasdws(isb)*nbasdws(isb)                           2d1s19
       end do                                                           8d5s14
       xlam=xlam*2d0                                                    8d5s14
       write(6,*)('and increase xlam to '),xlam                         8d5s14
       if(abs(xlam0-xlam).lt.1d-6)then                                  8d5s14
        if(lprint)write(6,*)('we are cycling ... ')                               8d5s14
        xlam=xlam*2d0
        if(lprint)write(6,*)('and increase xlam to '),xlam                         8d5s14
       end if                                                           8d5s14
       xlam0=xlam                                                       8d5s14
       if(xlam.gt.1d3)then
        if(lprint)write(6,*)('this is a rediculous value of xlam ')               8d5s14
        go to 2021                                                      8d5s14
       end if                                                           8d5s14
       go to 2020                                                       8d5s14
      end if                                                            8d5s14
      call werup(iumat,itmat,ibmat10,iamat10,jmats,kmats,iqk,ipk,iden1, 4d19s13
     $     ih0mo,numpro,nowpro,igmat,noc,ionex,ioooo,nbasdws,nvirt,bc,  11d9s22
     $     ibc)                                                         11d9s22
      ehfmicl=ehf2nd                                                    8d5s14
      if((abs(echan).gt.1d-12.or.micit.eq.1).and.iwerner.eq.1)then      8d5s14
       go to 2020                                                       8d4s14
      end if                                                            8d4s14
 2021 continue                                                          8d5s14
      xlam=min(xlam,xlam0,xlam00)                                       8d5s14
       do isb=1,nsymb
        if(noc(isb).gt.0)then
         nbas4=nbasdws(isb)                                             2d1s19
         nbas3=nbasisp(isb)*ncomp                                       4d30s18
         iarg1=nbas4
         iorb=iorbn(isb)                                                9d16s04
         if(iprtr(5).ne.0)then
         write(6,*)('multiply old orbs for symmetry '),isb
         call prntm2(bc(iorb),nbas3,iarg1,nbas3)
         write(6,*)('times umat ')
         call prntm2(bc(iumat(isb)),iarg1,iarg1,iarg1)
         end if
         idwst1=ibcoff                                                  8d5s14
         ibcoff=idwst1+nbas3*nbas4                                      8d5s14
         call enough('parahf. 61',bc,ibc)                               3d13s23
         if(nbas3.gt.0.and.nbas4.gt.0)then                              2d25s19
         call dgemm('n','n',nbas3,nbas4,nbas4,1d0,bc(iorb),nbas3,
     $        bc(iumat(isb)),nbas4,0d0,bc(idwst1),nbas3,
     d' parahf. 17')
         end if                                                         2d25s19
         if(iprtr(5).ne.0)then                                          5d7s21
         writE(6,*)('to get new orbs ')
         call prntm2(bc(idwst1),nbas3,iarg1,nbas3)
         end if
         iarg1=nbas3*nbas4                                              4d24s06
         call dws_bcast(bc(idwst1),iarg1)                               3d12s12
         do idws=1,nbas3*nbas4
          bc(iorbn(isb)+idws-1)=bc(idwst1+idws-1)                       3d15s12
         end do
         ibcoff=idwst1                                                  8d5s14
        end if
       end do
       call second(timec)                                                2d20s14
       if(inocan.eq.0)then                                              5d2s16
        if(iprtr(5).ne.0)write(6,*)('calling cannon')                   5d7s21
        call cannon(ih0mo,numpro,nowpro,kmats,jmats,ioooo,noc,iorbn,      3d19s12
     $       idwsdeb,ieigv,ionex,nbasdws,nvirt ,idorel,idbly,iacto,     2d1s19
     $       iden1,noocc,nbasisp,iovr,ifockeig,nlzz,iorbsym,iorbsymz,   6d7s21
     $       iorbsymc,bc,ibc)                                           11d9s22
        if(nlzz.ne.0)then                                               7d21s21
         do j=0,nxsb-1                                                  7d21s21
          ibc(ixsb+j*mxsb+3)=0                                          7d21s21
         end do                                                         7d21s21
         do isb=1,nsymb                                                  7d21s21
          jorb=iorbn(isb)                                                7d21s21
          nbas4=nbasisp(isb)*ncomp                                      7d21s21
          nbas3=nbasdws(isb)                                            3d2s22
          if(mynowprog.eq.0.and.nbas3.gt.0)then                         7d3s23
           isymdws(1)=ibcoff                                                9d10s07
           ibcoff=isymdws(1)+nbas3                                          9d10s07
           inuq=ibcoff                                                    6d7s21
           nuq=0                                                          6d7s21
           if(nlzz.eq.2)then                                            8d23s21
            do i=0,nbas3-1                                                4d19s21
             ibc(isymdws(1)+i)=ibc(iorbsymc(isb)+i)                       7d21s21
            end do                                                        4d19s21
           else                                                         8d23s21
            do i=0,nbas3-1                                                4d19s21
             ibc(isymdws(1)+i)=ibc(iorbsymc(isb)+i)/100                 8d23s21
            end do                                                        4d19s21
           end if                                                       8d23s21
           do i=0,nbas3-1
            do j=0,nuq-1                                                  6d7s21
             if(ibc(iorbsymc(isb)+i).eq.ibc(inuq+j))go to 3237          11d17s21
            end do                                                        6d7s21
            ibcoff=inuq+nuq                                               6d7s21
            call enough('parahf. 62',bc,ibc)                            3d13s23
            ibc(inuq+nuq)=ibc(iorbsymc(isb)+i)                          11d17s21
            nuq=nuq+1                                                     6d7s21
 3237       continue                                                      7d21s21
           end do
           ibcoff=inuq+nuq                                              11d17s21
           iherei=ibcoff                                                  6d7s21
           ibcoff=iherei+nuq*2                                            6d7s21
           call enough('parahf. 63',bc,ibc)                             3d13s23
           iherei0=ibcoff                                                 6d7s21
           do i=0,nuq-1                                                   6d7s21
            juse=-1                                                       6d7s21
            itest=ibc(inuq+i)                                           11d17s21
            if(nlzz.eq.6)itest=itest/100                                11d17s21
            do j=0,nxsb-1                                                 6d7s21
             if(itest.eq.ibc(ixsb+mxsb*j))juse=j                        11d17s21
            end do                                                        6d7s21
            if(juse.ge.0)then                                             6d7s21
            else                                                          6d7s21
             write(6,*)('could not match symmetries ')
             stop 'parahf'                                                6d7s21
            end if                                                        6d7s21
            ibc(iherei+i*2)=ibcoff                                        6d7s21
            ibc(iherei+i*2+1)=ibcoff                                      6d7s21
            ibcoff=ibcoff+ibc(ixsb+mxsb*juse+1)*ibc(ixsb+mxsb*juse+6)    7d21s21
           end do                                                         6d7s21
           call enough('parahf. 64',bc,ibc)                             3d13s23
           do i=iherei0,ibcoff-1                                          6d7s21
            bc(i)=0d0                                                     6d7s21
           end do                                                         6d7s21
           do i=0,nbas4-1                                                 6d7s21
            lsprim=ibc(iorbsymao(isb)+i)                                 7d21s21
            do j=0,nuq-1                                                  6d7s21
             if(lsprim.eq.ibc(inuq+j))then                              11d17s21
              jxsbh=iherei+j*2                                            6d7s21
              jxsbh0=ibc(iherei+j*2+1)
              go to 3234                                                  7d21s21
             end if                                                       6d7s21
            end do                                                        6d7s21
            write(6,*)('oh no! we had no match! '),lsprim
            write(6,*)('isb, nbas3: '),isb,nbas3                        7d3s23
            write(6,*)(ibc(inuq+j),j=0,nuq-1)
 3234       continue                                                      7d21s21
            do j=0,nbas3-1                                                6d7s21
             jjii=jorb+i+nbas4*j                                         7d21s21
             if(ibc(iorbsymc(isb)+j).eq.lsprim)then                     11d17s21
              bc(ibc(jxsbh))=bc(jjii)                                     6d7s21
              ibc(jxsbh)=ibc(jxsbh)+1                                     6d7s21
             else                                                         6d7s21
             end if                                                       6d7s21
            end do                                                        6d7s21
           end do                                                         6d7s21
           inuqg=ibcoff                                                   6d7s21
           ibcoff=inuqg+nxsb                                              6d7s21
           call enough('parahf. 65',bc,ibc)                             3d13s23
           do i=0,nxsb-1                                                  6d7s21
            ibc(inuqg+i)=0                                                6d7s21
           end do                                                         6d7s21
           do i=0,nuq-1                                                   6d7s21
            if(ibc(iherei+2*i).ne.ibc(iherei+2*i+1))then                  6d7s21
             itest=ibc(inuq+i)                                          11d17s21
             if(nlzz.eq.6)itest=itest/100                               11d17s21
             juse=-1                                                      6d7s21
             do j=0,nxsb-1                                                6d7s21
              if(itest.eq.ibc(ixsb+j*mxsb))juse=j                       11d17s21
             end do                                                       6d7s21
             if(juse.ge.0)then                                            6d7s21
             else                                                         6d7s21
              write(6,*)('could not match symmetry ')
              stop 'parahf'                                               6d7s21
             end if                                                       6d7s21
             nsz=ibc(ixsb+mxsb*juse+1)*ibc(ixsb+mxsb*juse+6)             7d21s21
             ndwsr=ibc(ixsb+mxsb*juse+6)
             ndwsc=ibc(ixsb+mxsb*juse+1)
             if(ibc(ixsb+juse*mxsb+3).eq.0.and.ibc(inuqg+juse).eq.0)then  6d7s21
              do j=0,nsz-1                                               7d21s21
               bc(ibc(ixsb+mxsb*juse+4)+j)=bc(ibc(iherei+2*i+1)+j)        6d7s21
              end do                                                      6d7s21
              ibc(inuqg+juse)=1                                           6d7s21
              ibc(ixsb+juse*mxsb+3)=1                                    7d21s21
             else                                                         6d7s21
              do j=0,nsz-1                                               7d21s21
               bc(ibc(iherei+2*i+1)+j)=bc(ibc(ixsb+mxsb*juse+4)+j)        6d7s21
              end do                                                      6d7s21
             end if                                                       6d7s21
             ibc(iherei+2*i)=ibc(iherei+2*i+1)                            6d7s21
            end if                                                        6d7s21
           end do                                                         6d7s21
           ibcoff=inuqg                                                   6d7s21
           do i=0,nbas4-1                                                7d21s21
            lsprim=ibc(iorbsymao(isb)+i)                                 7d21s21
            do j=0,nuq-1                                                  6d7s21
             if(lsprim.eq.ibc(inuq+j))then                                6d7s21
              jxsbh=iherei+j*2                                            6d7s21
              go to 3235                                                  6d7s21
             end if                                                       6d7s21
            end do                                                        6d7s21
 3235       continue                                                      6d7s21
            do j=0,nbas3-1                                                6d7s21
             jjii=jorb+i+nbas4*j                                         7d21s21
             if(ibc(iorbsymc(isb)+j).eq.lsprim)then                          6d7s21
              bc(jjii)=bc(ibc(jxsbh))                                     6d7s21
              ibc(jxsbh)=ibc(jxsbh)+1                                     6d7s21
             else                                                         6d7s21
              bc(jjii)=0d0                                                      6d7s21
             end if                                                       6d7s21
            end do                                                        6d7s21
           end do                                                         6d7s21
           ibcoff=inuq                                                    6d7s21
           ibcoff=isymdws(1)                                             7d21s21
          end if                                                        7d21s21
          if(min(nbas3,nbas4).gt.0)then                                 1d11s23
           nsnd=nbas3*nbas4                                              7d21s21
           call dws_bcast(bc(jorb),nsnd)                                 7d21s21
          end if                                                        1d11s23
         end do                                                          7d21s21
        end if                                                          7d21s21
   88    format(20(10i1,x))
       end if                                                           3d28s16
      call second(timec2)                                               2d20s14
      telap=timec2-timec-tovr                                           2d20s14
      times(13,1)=times(13,1)+telap
      times(13,2)=times(13,2)+telap**2
      times(13,3)=times(13,3)+1d0
      rmsa=0d0                                                          9d21s09
      nbasisps=0                                                         3d15s12
      do isb=1,nsymb                                                    9d21s09
       nbasisps=nbasisps+nbasdws(isb)                                   2d1s19
       rms=0d0                                                          9d21s09
       if(nbasdws(isb).gt.0)then                                        2d1s19
        do i=1,nbasdws(isb)                                             2d1s19
         iad=iumat(isb)+(i-1)*(nbasdws(isb)+1)                          2d1s19
         sav=bc(iad)                                                    9d21s09
         bc(iad)=sav-1d0                                                9d21s09
         do j=1,nbasdws(isb)                                            2d1s19
          iad=iumat(isb)+j-1+nbasdws(isb)*(i-1)                         2d1s19
          rms=rms+bc(iad)**2                                            9d21s09
         end do                                                         9d21s09
         bc(iad)=sav                                                    9d21s09
        end do                                                          9d21s09
        rmsa=rmsa+rms                                                   9d21s09
        rms=sqrt(rms/dfloat(nbasdws(isb)**2))                           2d1s19
       end if                                                           9d21s09
      end do                                                            9d21s09
      rmsa=sqrt(rmsa/dfloat(nbasisps**2))                                9d21s09
      call second(timeq2)                                               3d7s07
      telap=timeq2-timeq1-tovr                                          3d7s07
      times(8,1)=times(8,1)+telap                                       3d7s07
      times(8,2)=times(8,2)+telap**2                                    3d7s07
      times(8,3)=times(8,3)+1d0                                         3d7s07
      if(i404.eq.1)then                                                 2d19s14
       ibcoff=ibcreset                                                  2d19s14
       go to 404
      end if                                                            2d19s14
      ibcoff=iamat10                                                    11d24s04
      go to 402                                                         9d14s04
      return
      end
       subroutine readin(a,iunit,irow,icol,idone,id,itri,itrans)
       implicit real*8 (a-h,o-z)
       character*1 dcar
       character*80 line
c
c     read in matrix.
c
       dimension a(1)
       write(6,*)('in readin. iunit = '),iunit,itrans
       idone=1
       read(iunit,101,end=2)line
  101  format(a80)
       do i=1,80
        if(line(i:i).eq.'@')then
	 read(line(1:i-1),*)irow,icol
	 go to 102
	end if
       end do
  102  continue
        iread=5
       if(irow*iabs(icol).gt.id*id)stop 'id'
    1  format(1x,2i5)
       write(6,*)('irow,icol '),irow,icol
       if(icol.lt.0.or.irow.ne.icol)then
        itri=0
        icol=iabs(icol)
        npart=icol/iread
        if(iread*npart.ne.icol)npart=npart+1
        call flush(6)
        i0=1
        do 3 ipart=1,npart
         ie=min(i0+iread-1,icol)
         read(iunit,4)idum
    4    format(11x,i5)
         do 5 il=1,irow
          if(itrans.eq.0)then
           read(iunit,6,err=555)(a(il+(i-1)*irow),i=i0,ie)
          else
           read(iunit,6,err=555)(a(i+(il-1)*irow),i=i0,ie)
          end if
    6     format(6x,1p5e14.6)
    5    continue
         i0=ie+1
    3   continue
       else
        icol=iabs(icol)
        itri=1
        npart=icol/iread
        if(iread*npart.ne.icol)npart=npart+1
        i0=1
        do 13 ipart=1,npart
         read(iunit,4)idum
         do 15 il=i0,irow
          ie=i0+min(iread-1,il-i0)
          if(itrans.eq.0)then
           read(iunit,16,err=555)(a(il+(i-1)*irow),i=i0,ie)
          else
           read(iunit,16,err=555)(a(i+(il-1)*irow),i=i0,ie)
          end if
   16     format(6x,1p9e14.6)
   15    continue
         i0=ie+1
   13   continue
        do 17 ir=1,icol
         ii=(ir-1)*irow
         iii=ir-irow
         do 18 il=1,ir
          a(il+ii)=a(iii+il*irow)
   18    continue
   17   continue
       end if
       read(iunit,1532)dcar
 1532  format(1x,a1)
       return
    2  continue
       write(6,20)iunit
   20  format(/1x,'end of file encounter on unit ',i5)
       stop
  555  continue
       write(6,*)('i/o error on read')
       write(6,*)('perhaps you forgot to substitute "tiny" for 0.00?')
       stop
       end
