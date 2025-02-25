c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine basisz(ibdat,ibstor,isstor,pwd,ifname,ncore,nmn,bc,ibc,12d13s22
     $     cfile,ncfile)                                                12d13s22
c
c     input for Massively
c               Parallel
c               Electronic
c               Correlation
c
c     CAS and MRCI.
c
c     below, we give the input in "" to show spacing as required.
c     the actual input should not contain the "'s.
c     line 1:
c     "/u/schwenke/mpec/source/basis_library"
c      this is the name of the directory containing basis set libraries.
c      the above entry is appropriate for columbia.
c     the next lines specify atoms, positions and basis set
c     see below.
c     the following lines specify the calculation options.
c
c
c     what's under ibdat
c     for each gaussian...
c     1: l
c     2: exponential parameter
c     3: normalization coefficient
c     4: offset in total basis function list
c     5: x of center
c     6: y of center
c     7: z of center
c     8: nhere -> address in iapair of possible symmetry equivalent
c        center
c     9: contraction index
c
c     what's under iapair:
c     1: if nonzero, iabs is symmetry equivalent basis fcn
c     2: group operation relating symmetry equivalent basis fcns
c     3: phase?
c
      implicit real*8 (a-h,o-z)
      integer*2 idwstry                                                 5d18s06
c                                                                       1d10s19
c     iddx is the maximum number of orbital files for the diabatic      1d10s19
c     orbital option.                                                   1d10s19
c                                                                       1d10s19
      parameter (iddx=10)                                               1d10s19
      include "common.hf"
      include "common.cas"
      include "common.store"
      include "common.basis"
      include "common.input"                                            3d1s10
      include "common.spher"
      include "common.mrci"                                             5d24s18
      include "common.print"                                            1d3s20
      logical lcart,leof                                                8d3s21
      integer*1 idogrado1(4)                                            6d18s22
      integer*2 nlzzq2(2)                                               12d5s22
      integer*8 ipack8                                                  12d13s22
      dimension ipack4(2)                                               12d13s22
      equivalence (ipack8,ipack4)                                       12d13s22
      equivalence (nlzzq,nlzzq2)                                        12d5s22
      character*1 wedge                                                 3d23s23
      character*26 ends(ida,2)                                          4d12s19
      character*(*) cfile                                               12d13s22
      character*4 atomname                                              4d12s19
      character*2 sym,asym(2)                                           8d8s22
      character*80 line,bline,lineauto,opline                           1d11s16
      character*(*) pwd                                                 2d1s19
      character*80 bpath
      character*200 line200,bline200                                    5d15s19
      parameter(ide=101)                                                1d3s24
      character*3 stype(8,8),element(ide)                               12d24s19
      character*1 digit(11),btype(5)                                    2d26s21
      character*5 symatoms,itype                                        8d17s15
      parameter (idcg=40,idctr=40,idll=7)
      parameter (maxpack=100)                                           2d18s10
      dimension jsh(1),ctry(3),isop(7),noc(8),lkeep(11),ladd(11),       10d2s24
     $     inoc(11)                                                     10d2s24
      parameter (idpl=100)                                              11d15s05
      dimension iplist(2,idpl),jplist(idpl),lplist(idpl),ihit(2),       12d22s05
     $     nblock(8),nlptr(ida),ippack(maxpack),cartb(3),fact(3)        2d18s10
      dimension japair(ida),igrp(8),istinfo(11),nbaspre(8),             3d3s20
     $     mapd2h(8),ngausa(ida),ismul(3),irunc(8),irunc2(8)            1d28s19
      dimension isymt(3),isymsub(8),isymd2h(3,8),itmp(8),               3d17s21
     $     insym(8),negs(3),rjac(ida,3),iends(ida,2),xm(ida),mya(ida),  4d12s19
     $     mycart(maxpack,2),icoreg(4,ida),ivalg(4,ida),ncoreg(8),      6d11s19
     $     nvalg(8),jcoreg(4,ida),jvalg(4,ida),ibcode(8,iddx),          1d10s19
     $     iorb(8,iddx),iorbao(8,iddx),icanon(8),isymu(8)               12d7s23
      equivalence (idogrado,idogrado1)                                  6d18s22
      character*7 calc                                                  3d22s12
      data isym/0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1/ 4d5s10
      data btype/'s','p','d','f','g'/                                   11d22s19
      data stype(1,1)/'A  '/                                            12d24s19
      data stype(1,2),stype(2,2)/'A'' ','A'''''/                        12d24s19
      data (stype(i,4),i=1,4)/'A1 ','B2 ','B1 ','A2 '/                  8d16s21
      data (stype(i,8),i=1,8)                                           12d24s19
     $     /'Ag ','B3u','B2u','B1g','B1u','B2g','B3g','Au '/            12d24s19
      data digit/'0','1','2','3','4','5','6','7','8','9','-'/           2d26s21
      data element/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ',     4d6s18
     $     'F  ','Ne ','Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ', 4d6s18
     $     'K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ', 4d6s18
     $     'Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ', 4d6s18
     $     'Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ', 4d6s18
     $     'In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ', 4d6s18
     $     'Pr ','Nd ','Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ', 4d6s18
     $     'Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ', 4d6s18
     $     'Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ', 4d6s18
     $     'Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ', 4d6s18
     $     'Es ','Fm ','dmy'/                                           1d3s24
      data isymu/8*1/                                                   12d7s23
      do i=1,maxst1                                                     7d3s23
       nlambda(i,1)=-1                                                  8d31s23
       nlambda(i,2)=-1                                                  8d31s23
      end do                                                            7d3s23
      do ia=1,ida                                                       3d23s23
       xm(ia)=0d0                                                       3d23s23
       xmaxs(ia)=0d0                                                    3d23s23
       do ib=1,ida                                                      3d23s23
        coord(ib,ia)=0d0                                                3d23s23
       end do                                                           3d23s23
      end do                                                            3d23s23
      icdatsave=0                                                       1d11s23
      nmn=0                                                             8d12s22
      atom(1)='us  '                                                    3d2s22
      iusecartr=0                                                       2d21s20
      idiglow=ichar(digit(1))                                           12d24s19
      idighig=ichar(digit(1))                                           12d24s19
      do i=2,10                                                         12d24s19
       idiglow=min(idiglow,ichar(digit(i)))                             12d24s19
       idighig=max(idighig,ichar(digit(i)))                             12d24s19
      end do                                                            12d24s19
      try=0.5d0                                                         6d11s19
 1792 continue                                                          6d11s19
       try=try*0.5d0                                                    6d11s19
       test=1d0+try                                                     6d11s19
      if(test.ne.1d0)go to 1792                                         6d11s19
      smallest=try*2d0                                                  6d11s19
      write(6,1793)smallest                                             11d17s22
 1793 format('smallest number one can add to 1 and not get 1 is ',      11d17s22
     $     es10.3)                                                      11d17s22
      do i=1,8                                                          6d11s19
       ncoreg(i)=0                                                      6d11s19
       nvalg(i)=0                                                       6d11s19
      end do                                                            6d11s19
      do i=1,idprt                                                      1d3s20
       iprtr(i)=0                                                       1d3s20
      end do                                                            1d3s20
      nprtr=1
c     1
      prtrname(nprtr)(1:12)='stepwisecsf '                                  1d3s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='gencsf3     '                              1d3s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='gencsf2     '                              1d3s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='mrci        '                              1d3s20
      nprtr=nprtr+1
c     5
      prtrname(nprtr)(1:12)='casopt      '                              1d3s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='hccsfd      '                              1d3s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='parahf      '                              1d3s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='cas0        '                              1d3s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='cont2ecsf   '                              1d3s20
c     10
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='parajkfromh '                              1d6s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='updatex     '                              1d6s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='intcsf      '                              1d6s20
      nprtr=nprtr+1
      prtrname(nprtr)(1:12)='contractg   '                              1d22s20
      nprtr=nprtr+1                                                     2d13s20
      prtrname(nprtr)(1:12)='updateg     '                              2d13s20
c     15
      nprtr=nprtr+1                                                     2d14s20
      prtrname(nprtr)(1:12)='cupo2       '                              2d14s20
      nprtr=nprtr+1                                                     2d14s20
      prtrname(nprtr)(1:12)='genergy     '                              3d3s20
      nprtr=nprtr+1                                                     2d14s20
      prtrname(nprtr)(1:12)='opto        '                              3d3s20
      nprtr=nprtr+1                                                     2d14s20
      prtrname(nprtr)(1:12)='paraeri     '                              5d4s20
      nprtr=nprtr+1                                                     2d14s20
      prtrname(nprtr)(1:12)='genryd      '                              5d6s20
c     20
      nprtr=nprtr+1                                                     1d26s21
      prtrname(nprtr)(1:12)='genfcn      '                              1d26s21
      nprtr=nprtr+1                                                     1d26s21
      prtrname(nprtr)(1:12)='sortoutsym  '                              4d26s21
      nprtr=nprtr+1                                                     1d26s21
      prtrname(nprtr)(1:12)='mofromao    '                              4d28s21
      nprtr=nprtr+1                                                     5d3s21
      prtrname(nprtr)(1:12)='makeguess   '                              5d3s21
      nprtr=nprtr+1                                                     5d3s21
      prtrname(nprtr)(1:12)='dotdvecs    '                              5d7s21
c     25
      nprtr=nprtr+1                                                     5d3s21
      prtrname(nprtr)(1:12)='cannon      '                              7d21s21
      nprtr=nprtr+1                                                     5d3s21
      prtrname(nprtr)(1:12)='psisopsi    '                              3d2s22
      nprtr=nprtr+1                                                     5d3s21
      prtrname(nprtr)(1:12)='prntm2      '                              5d19s22
      nprtr=nprtr+1                                                     5d3s21
      prtrname(nprtr)(1:12)='grad        '                              8d4s22
      nprtr=nprtr+1                                                     5d3s21
      prtrname(nprtr)(1:12)='restart     '                              6d23s23
c     30
      nprtr=nprtr+1                                                     3d29s24
      prtrname(nprtr)(1:12)='solin       '                              3d29s24
      nprtr=nprtr+1                                                     3d29s24
      prtrname(nprtr)(1:12)='orbs        '                              4d15s24
      nprtr=nprtr+1                                                     3d29s24
      prtrname(nprtr)(1:12)='denmx12     '                              6d14s24
      nprtr=nprtr+1                                                     3d29s24
      prtrname(nprtr)(1:12)='mp2         '                              8d9s24
      epsym=1d-8                                                        2d3s12
      pi=acos(-1d0)
      rad=pi/1.8d2                                                      4d12s19
      srpi=sqrt(pi)
      irtrn=0                                                           5d13s05
      nidaxx=0                                                          10d23s08
      nlzz=0                                                            8d11s14
      nvec=0                                                            4d12s19
      do i=1,ida                                                        11d25s08
       nlptr(i)=0                                                       11d25s08
      end do                                                            11d25s08
      nextradatx=0                                                      1d10s19
      do i=1,3                                                          1d5s20
       xcart(i,1)=0d0                                                   1d5s20
      end do                                                            1d5s20
      itmax=10
      idomp2=0                                                          8d8s14
      idorel=0                                                          8d17s15
      ianuc=-1                                                          8d17s15
      idohf=0                                                           8d8s14
      inocan=0                                                          3d28s16
      idocas=0                                                          8d8s14
      idomrci=0                                                         5d24s18
      ipack=0
      lmax=0
      ngaus=0
      natom=0                                                           1d22s10
      ncore=1                                                           5d31s21
      maxbct=maxbc                                                      1d31s21
      iloopit=0                                                         5d24s18
 5577 continue                                                          11d4s05
      read(5,100,end=1008)line                                            11d4s05
      if(line(1:1).eq.'#')go to 5577                                    12d23s22
  254 continue                                                          8d8s14
      iloopit=iloopit+1                                                 5d24s18
      is=1                                                              8d8s14
      call delim(line,is,ie)                                            8d8s14
      if(ie.lt.is)go to 5577                                            5d19s16
      if(line(is:is).eq.'#')go to 5577                                  10d24s17
      if(line(is:is+2).eq.'&BG')go to 110                               8d8s14
      if(line(is:is+4).eq.'&SOHF')go to 1252                            8d8s14
      if(line(is:is+3).eq.'&CAS')go to 251                              8d8s14
      if(line(is:is+4).eq.'&MRCI')go to 2251                            5d24s18
      if(line(is:is+4).eq.'&PROP')go to 2252                            11d21s19
      if(line(is:is+3).eq.'&MP2')go to 253                              8d8s14
      if(line(is:is+3).eq.'&REL')go to 1253                              8d17s15
      if(line(is:is+4).eq.'&DIAB')go to 2253                            1d10s19
      if(line(is:is+3).eq.'&SYS')go to 3253                             1d31s21
      write(6,*)('whats under line: "'),line,('"')
      write(6,*)('what we have for is: '),is
      write(6,*)('what we have for string: '),line(is:is+3)
      write(6,*)('unknown input - jump to end ')                        8d7s14
      go to 1008                                                        8d7s14
c
c     for each atom basis type, input (atomic units only)               8d25s17
c     Z=nucz (or chemical symbol)                                       4d6s18
c     x1 y1 z1                                                          8d25s17
c        .                                                              8d25s17
c        .                                                              8d25s17
c        .                                                              8d25s17
c     xn yn zn                                                          8d25s17
c     l,ngaus,ncont                                                     1d28s19
c     zeta(1) zeta(2) ... zeta(ngaus)
c     c(1,1) c(1,2) ... c(1,nzetau)
c     c(2,1) c(2,2) ... c(1,nzetau)
c     .
c     .
c     .
c     c(ncont,1) c(ncont,2) ... c(ncont,nzetau)
c     with last l through c(ncont,1) lines repeated as necessary.
c
c     the number of atoms with the same basis set and Z, n in the input 8d25s17
c     above, is determined by the program from the input: the xyz input 8d25s17
c     will(must!) contain 3 decimal points, while the l,nzeta,ncont     1d28s19
c     input will have none.                                             1d28s19
c     if ncont is omitted, no basis set contraction will be used.       9d11s17
c     if this is a non-relativistic calculation, nzetau=nzeta, but if   9d11s17
c     this is a relativistic calculation, nzetau=2*nzeta.               9d11s17
c     optionally, one can specify a file name instead of the basis
c     set lines. In this case, the file name must not begin with a
c     numeral, and will contain exactly the above information. note
c     only one basis set and atom per file, and currently no contraction
c     is allowed.
c
c     end of this part of the input is signified by a line containing
c     & and the keyword for the next program                            8d8s14
c
c     all lists are blank delimited!!!
c
c     optionally...
c     one can include the line sym=off before the first Z=nucz line
c     to force the use of no symmetry.
c
  110 continue                                                          8d8s14
      alpha=1d0/137.0359895d0                                           1d2s23
      ialpha=0                                                          1d2s23
      ndfld=0                                                           12d6s23
      nvec=0                                                            4d10s19
      nblba=0                                                           9d10s19
      nwcont=0                                                          9d12s17
      xlength=1d0                                                       8d31s21
    1 continue
       read(5,200,end=2)line200
       line=line200(1:80)                                               3d23s23
       is=1                                                             12d7s23
       call delim(line,is,ie)                                           12d7s23
       if(line(is:is+2).eq.'sym')then                                   12d7s23
        if(line(ie-2:ie).eq.'off')then                                  12d7s23
         do i=1,8                                                       12d7s23
          isymu(i)=0                                                    12d7s23
         end do                                                         12d7s23
         go to 1                                                        12d7s23
        end if                                                          12d7s23
 5665   continue                                                        12d7s23
        is=ie+1                                                         12d7s23
        call delim(line,is,ie)                                          12d7s23
        if(ie.ge.is)then                                                12d7s23
         if(line(is:is+2).eq.'off')then                                 3d14s24
          write(6,*)('symmetry is turned off ...')                      3d14s24
          do i=1,8                                                      3d14s24
           isymu(i)=0                                                   3d14s24
          end do                                                        3d14s24
          go to 1                                                       3d14s24
         end if                                                         3d14s24
         if(line(is:is).eq.'x')then                                     12d7s23
          isymu(2)=2                                                    12d7s23
         else if(line(is:is).eq.'y')then                                12d7s23
          isymu(3)=2                                                    12d7s23
         else if(line(is:is).eq.'z')then                                12d7s23
          isymu(5)=2                                                    12d7s23
         else                                                           12d7s23
          write(6,*)('I don''t understand "'),line(is:ie),              12d7s23
     $         ('as an argument to the sym keyword')                    12d7s23
          write(6,*)('help!')                                           12d7s23
          irtrn=1                                                       12d7s23
          return                                                        12d7s23
         end if                                                         12d7s23
         go to 5665                                                     12d7s23
        end if                                                          12d7s23
        go to 1                                                         12d7s23
       end if                                                           12d7s23
       is=1                                                             8d8s14
       call delim(line,is,ie)                                           8d8s14
       if(line(is:is).eq.'&'.or.line(is:is+1).eq.'#&')then              5d24s18
        go to 2                                                         5d24s18
       end if                                                           5d24s18
       if(line(is:is).eq.'#')go to 1                                    8d29s22
       if(line(is:is+2).eq.'len')then                                   8d31s21
        is=ie+1                                                         8d31s21
        call delim(line,is,ie)                                          8d31s21
        if(ie.lt.is)then                                                8d31s21
         write(6,*)('you have included the length key word without'),   8d31s21
     $       (' an argument')                                           8d31s21
         write(6,*)('you must specify an argument')                     8d31s21
         irtrn=1                                                        8d31s21
         return                                                         8d31s21
        end if                                                          8d31s21
        write(6,*)('read xlength from '),line(is:ie)                    8d31s21
        do i=is,ie                                                      8d31s21
         if(line(i:i).eq.'/')then                                       8d31s21
          im=i-1                                                        8d31s21
          read(line(is:im),*)xnumerator                                 8d31s21
          ip=i+1                                                        8d31s21
          read(line(ip:ie),*)denominator                                8d31s21
          xlength=xnumerator/denominator                                8d31s21
          go to 1166                                                    8d31s21
         end if                                                         8d31s21
        end do                                                          8d31s21
        read(line(is:ie),*)xlength                                      8d31s21
 1166   continue                                                        8d31s21
        write(6,*)('input distance units are equal to '),xlength,       8d31s21
     $       ('bohr')                                                   8d31s21
        go to 1                                                         9d6s21
       end if                                                           8d31s21
       if(line(is:is+5).eq.'ffield')then                                8d29s22
        write(6,*)                                                      8d29s22
        write(6,*)('finite field ...')                                  8d29s22
        is=ie+1                                                         8d29s22
        call delim(line,is,ie)                                          8d29s22
 2743   continue                                                        12d6s23
        if(ie.ge.is)then                                                8d29s22
         write(6,*)('>field type: "'),line(is:ie),('"')                  8d29s22
         ndfld=ndfld+1                                                  12d6s23
         if(ndfld.gt.idfld)then                                         12d6s23
          write(6,*)('toooo many field types!!!')                       12d6s23
          write(6,*)('the maximum number of field types is set by')     12d6s23
          write(6,*)('the parameter idfld in common.input ')            12d6s23
          irtrn=1                                                       12d6s23
          return                                                        12d6s23
         end if                                                         12d6s23
         if(ie.eq.is)then                                               8d29s22
          iftype(ndfld)=0                                               12d6s23
          if(line(is:is).eq.'x')then                                    8d29s22
           iftype(ndfld)=1                                              12d6s23
          else if(line(is:is).eq.'y')then                               8d29s22
           iftype(ndfld)=2                                              12d6s23
          else if(line(is:is).eq.'z')then                               8d29s22
           iftype(ndfld)=3                                              12d6s23
          end if                                                        8d29s22
         else if(ie.eq.is+1)then                                        8d29s22
          if(line(is:is+1).eq.'xx')then                                   8d29s22
           iftype(ndfld)=7                                              12d6s23
          else if(line(is:is+1).eq.'yy')then                              8d29s22
           iftype(ndfld)=8                                              12d6s23
          else if(line(is:is+1).eq.'zz')then                              8d29s22
           iftype(ndfld)=9                                              12d6s23
          else if(line(is:is+1).eq.'xy')then                              8d29s22
           iftype(ndfld)=4                                              12d6s23
          else if(line(is:is+1).eq.'xz')then                              8d29s22
           iftype(ndfld)=5                                              12d6s23
          else if(line(is:is+1).eq.'yz')then                              8d29s22
           iftype(ndfld)=6                                              12d6s23
          end if                                                        8d29s22
         end if                                                         8d29s22
         if(iftype(ndfld).eq.0)then                                     12d6s23
          write(6,*)                                                    8d29s22
     $        ('I''m sorry, but the field type you have specified'),    8d29s22
     $        (' has not been coded yet. :(')                           8d29s22
          irtrn=1                                                       8d29s22
          return                                                        8d29s22
         end if                                                         8d29s22
         write(6,*)('>field type index is '),iftype(ndfld)              12d6s23
         is=ie+1                                                        8d29s22
         call delim(line,is,ie)                                         8d29s22
         if(ie.ge.is)then                                               8d29s22
          read(line(is:ie),*)fstgth(ndfld)                              12d6s23
          write(6,*)('>field strength = '),fstgth(ndfld)                12d6s23
          write(6,*)                                                      8d29s22
          is=ie+1                                                       12d6s23
          call delim(line,is,ie)                                        12d6s23
          if(is.le.ie)go to 2743                                        12d6s23
         else                                                           8d29s22
          write(6,*)('help! you need to specify a field strength!')     8d29s22
          irtrn=1                                                       8d29s22
          return                                                        8d29s22
         end if                                                         8d29s22
        else                                                            8d29s22
         write(6,*)('help! you need to specify a field type!')          8d29s22
         irtrn=1                                                        8d29s22
         return                                                         8d29s22
        end if                                                          8d29s22
        go to 1                                                         8d29s22
       end if                                                           8d29s22
       if(line(is:is+4).eq.'alpha')then                                 1d2s23
        is=ie+1                                                         1d2s23
        call delim(line,is,ie)                                          1d2s23
        if(ie.lt.is)then                                                1d2s23
         write(6,*)('there is no argument for the keyword alpha! ')     1d2s23
         write(6,*)('help!')                                            1d2s23
         irtrn=1                                                        1d2s23
         return                                                         1d2s23
        end if                                                          1d2s23
        write(6,*)('resetting the fine-structure constant alpha ')      1d2s23
        write(6,*)('from '),alpha                                       1d2s23
        read(line(is:ie),*)alpha                                        1d2s23
        write(6,*)('to '),alpha                                         1d2s23
        if(ialpha.eq.1)ascale=0.25d0*alpha*alpha                        1d2s23
        go to 1                                                         1d2s23
       end if                                                           1d2s23
       if(line(is:is+2).eq.'rel')then                                   5d25s21
        ascale=0.25d0*alpha*alpha
        ialpha=1                                                        1d2s23
        idorel4c=0                                                       2d20s20
c
c     idorel=1: no s 2e integrals
c     idorel=2: SS 2e integrals
c     idorel=3: SSSS 2e integrals
c     idorel minus those numbers: include breit
c
        is=ie+1                                                         5d25s21
        call delim(line,is,ie)                                          5d25s21
        if(ie.gt.is)then                                                10d4s22
         if(line(is:is+2).eq.'off')go to 1                              10d4s22
        end if                                                          10d4s22
        idorel=2
        write(6,*)
        write(6,*)('This is going to be a relativistic calculation ')
        write(6,*)
        if(ie.ge.is)then                                                5d25s21
         write(6,*)('we have itype in input ')                          5d25s21
         if(ie.eq.is+1)then                                             5d25s21
          itype='ss   '                                                 8d26s15
          idorel=2                                                      8d26s15
         else if(ie.eq.is+2)then                                        8d26s15
          if(line(is:is+1).eq.'ss')then                                 10d5s22
           if(line(ie:ie).eq.'b')then                                    2d18s20
            itype='ssb  '                                                 8d26s15
            idorel=-2                                                     8d26s15
           else if(line(ie:ie).eq.'g')then                               2d18s20
            itype='ssg  '                                                2d18s20
            idorel=-4                                                    2d18s20
           else                                                          2d18s20
            write(6,*)('unknown tail end on integral type: '),          10d5s22
     $           line(ie:ie)                                            10d5s22
            irtrn=1                                                      2d18s20
            return                                                       2d18s20
           end if                                                        2d18s20
          else if(line(is:ie).eq.'dhf')then                             10d5s22
           write(6,*)
     $ ('a full four component calculation will be carried out, thus')  5d25s21
           write(6,*)('symmetry is turned off, ')                         2d20s20
           do jj=1,8                                                    12d7s23
            isymu(jj)=0                                                 12d7s23
           end do                                                       12d7s23
           write(6,*)('basis set will be uncontracted, ')                 2d20s20
           do jj=1,8                                                    12d19s22
            inoc(jj)=-1                                                 12d19s22
           end do                                                       12d19s22
           idorel4c=-1                                                    2d26s20
           write(6,*)                                                     2d20s20
     $        ('cartesian rather than spherical harmonics are used')    2d20s20
           iusecartr=1                                                    2d21s20
           is=ie+1                                                      10d5s22
           call delim(line,is,ie)                                       10d5s22
           if(ie.lt.is.or.ie.eq.is+1)then                               10d5s22
            itype='ss   '                                               10d5s22
            idorel=2                                                    10d5s22
           else if(line(is+2:is+2).eq.'g')then                          10d5s22
            itype='ss   '                                               10d5s22
            idorel=-4                                                   10d5s22
           else if(line(is+2:is+2).eq.'b')then                          10d5s22
            itype='ss   '                                               10d5s22
            idorel=-2                                                   10d5s22
           else if(line(is+4:is+4).eq.'g')then                          10d5s22
            itype='ssss '                                               10d5s22
            idorel=-5                                                   10d5s22
           else if(line(is+4:is+4).eq.'b')then                          10d5s22
            itype='ssss '                                               10d5s22
            idorel=-3                                                   10d5s22
           else                                                         10d5s22
            itype='ssss '                                               10d5s22
            idorel=3                                                    10d5s22
           end if                                                       10d5s22
          end if                                                        10d5s22
         else if(ie.eq.is+3)then                                        8d26s15
          itype='ssss '                                                 8d26s15
          idorel=3                                                      8d26s15
         else if(ie.eq.is+4)then                                        10d5s22
          if(line(ie:ie).eq.'b')then                                    10d5s22
           itype='ssssb'                                                 8d26s15
           idorel=-3                                                     8d26s15
          else if(line(ie:ie).eq.'g')then                               2d18s20
           itype='ssssg'                                                2d18s20
           idorel=-5                                                    2d18s20
          else                                                          2d18s20
           write(6,*)('unknown tail end on integral type: '),           2d18s20
     $          line(ie:ie)                                             2d18s20
           irtrn=1                                                      2d18s20
           return                                                       2d18s20
          end if                                                        2d18s20
         else                                                           8d26s15
          write(6,*)('unknow itype code: '),line(is:ie)                 8d26s15
          go to 1008                                                    8d26s15
         end if                                                         8d26s15
        else
         itype='ss   '                                                  5d25s21
        end if                                                          5d25s21
        write(6,*)('2 electron integral type = '),itype                 5d25s21
        go to 1                                                         5d25s21
       end if                                                           5d25s21
  100  format(a80)
c
c     vibrational coordinate input?
c
 2287  continue                                                         4d10s19
       if(line(is:is+3).eq.'ends')then                                  4d10s19
        is=ie+1                                                         4d10s19
        call delim(line,is,ie)                                          4d10s19
        read(line(is:ie),*)ivec                                         4d10s19
        nvec=max(nvec,ivec)                                             4d10s19
        is=ie+1                                                         4d10s19
        call delim(line,is,ie)                                          4d10s19
        nn=ie+1-is                                                      4d10s19
        ends(ivec,1)(1:nn)=line(is:ie)                                  4d10s19
        iends(ivec,1)=nn                                                  4d10s19
        is=ie+1                                                         4d10s19
        call delim(line,is,ie)                                          4d10s19
        if(line(is:ie).ne.'to')then                                     4d10s19
         write(6,*)('vector input error: expecting ''to'', got '),      4d10s19
     $       line(is:ie)                                                4d10s19
         stop 'to'                                                      4d10s19
        end if                                                          4d10s19
        is=ie+1                                                         4d10s19
        call delim(line,is,ie)                                          4d10s19
        nn=ie+1-is                                                      4d10s19
        ends(ivec,2)(1:nn)=line(is:ie)                                  4d10s19
        iends(ivec,2)=nn                                                4d12s19
        is=ie+1                                                         4d10s19
        call delim(line,is,ie)                                          4d10s19
        write(6,*)('vector '),ivec,(' is '),ends(ivec,1)                4d12s19
     $       (1:iends(ivec,1))                                          4d12s19
     $       ,(' to '),ends(ivec,2)(1:iends(ivec,2))                    4d10s19
        if(iends(ivec,1).eq.1.and.iends(ivec,2).eq.1)nblba=nblba+1      9d10s19
        if(ie.ge.is)go to 2287                                          4d10s19
        go to 1                                                         4d10s19
       end if                                                           4d10s19
       if(nvec.gt.0)then                                                4d12s19
c
c     mass number followed by element name?
c
        iok=0                                                           4d12s19
        do i=is,ie                                                       4d12s19
         do j=1,10                                                      4d12s19
          if(line(i:i).eq.digit(j))then                                 4d12s19
           iok=iok+1                                                     4d12s19
           go to 212                                                     4d12s19
          end if                                                         4d12s19
         end do                                                         4d12s19
c
c     not a number. look for atomic name                                4d12s19
c
         npos=ie+1-i                                                    4d12s19
         if(npos.gt.0.and.npos.le.3)then                                8d9s22
          atomname='    '                                               4d12s19
          atomname(1:npos)=line(i:ie)                                   4d12s19
          do j=1,ide                                                     4d12s19
           if(atomname(1:3).eq.element(j))then                          4d12s19
            nnew=0                                                       4d12s19
            igot=0                                                       4d12s19
            nucz=j                                                       4d12s19
            nuczu=nucz                                                  7d5s24
            go to 122                                                    4d12s19
           end if                                                        4d12s19
          end do                                                         4d12s19
         end if                                                         4d12s19
c
c     not an element. jump out of loop
c
         go to 221                                                      4d12s19
  122    continue                                                       4d12s19
         if(is.ne.i)then                                                10d6s22
          nx=i-is                                                       10d6s22
          atomname(npos+1:npos+nx)=line(is:i-1)                         10d6s22
         end if                                                         10d6s22
         write(6,*)(' ')                                                12d20s19
         call getm(atomname,1,xmass,1,pwd//'/atwgts')                   4d12s19
         write(6,*)('xmass from getma '),xmass,natom,nuczu
         if(xmass.eq.0d0)xmass=-1d0                                     6d20s24
c
c     now search for "is a" etc.                                        4d12s19
c
         is=ie+1                                                        4d12s19
         call delim(line,is,ie)                                         4d12s19
         nahere=0                                                       4d12s19
         if(line(is:ie).eq.'is')then                                    4d12s19
  222     continue                                                      4d12s19
          is=ie+1                                                       4d12s19
          call delim(line,is,ie)                                        4d12s19
          if(ie.ge.is)then                                              4d12s19
           itry1=ichar(line(is:is))-ichar('a')                          4d12s19
           itry2=ichar(line(is:is))-ichar('A')                          4d12s19
           if(itry1.ge.0.and.itry1.le.25)then                           4d12s19
            iatom=itry1+1                                               4d12s19
           else if(itry2.ge.0.and.itry2.le.25)then                      1s1s20
            iatom=itry2+1                                               4d12s19
           else                                                         4d12s19
            write(6,*)('don''t know how to parse '),line(is:ie),        4d12s19
     $           (' as a letter of the alphabet'),itry1,itry2                       4d12s19
            stop 'alphabet'
           end if                                                       4d12s19
           xm(iatom)=xmass                                              4d12s19
           atom(iatom)=atomname                                         5d26s21
           write(6,*)('set atom '),iatom,('to '),atomname,atom(iatom)
           nahere=nahere+1                                              4d12s19
           ista=natom+nahere                                            4d12s19
           write(6,*)('is atoms '),nahere,line(is:is),iatom
           mya(ista)=iatom                                              4d12s19
           go to 222                                                    4d12s19
          end if                                                        4d12s19
          go to 121                                                     4d12s19
         else                                                           4d12s19
          write(6,*)('input error: expecting ''is'' but got '),
     $         line(is:ie)
          write(6,*)('from ')
          write(6,*)line
          stop 'no is'                                                  4d12s19
         end if                                                         4d12s19
  212    continue                                                       4d12s19
        end do                                                           4d12s19
       else                                                             4d12s19
        do i=1,ide                                                       4d6s18
         if(line(1:3).eq.element(i))then                                 4d6s18
          nnew=0                                                         4d6s18
          igot=0                                                         4d6s18
          xm(natom+1)=0d0                                               1d11s24
          nucz=i                                                         4d6s18
          go to 121                                                      4d6s18
         end if                                                          4d6s18
        end do                                                           4d6s18
       end if                                                           4d12s19
  221  continue                                                         4d12s19
c
c     values of vibrational coordinates
c
       if(nvec.gt.0)then                                                4d10s19
        igrab=0                                                         1d8s19
        if(line(is:is+3).eq.'grab')then                                  4d10s19
         igrab=1                                                        1d8s19
         is=ie+1                                                         4d10s19
         call delim(line,is,ie)                                          4d10s19
         write(6,*)('file to grab from is "'),line(is:ie),('"')
         open(unit=43,file=line(is:ie))                                 1d8s19
 1101    continue                                                       1d8s19
         read(43,200,end=1102)line200                                   3d23s23
         is=1                                                           1d8s19
         call delim(line200,is,ie)                                      3d23s23
         if(line200(is:is).eq.'r')then                                  3d23s23
          do ivec=1,nvec                                                 4d12s19
           is=ie+1                                                        4d10s19
           call delim(line200,is,ie)                                    3d23s23
           if(ie.lt.is)then
            write(6,*)('not enough vector lengths read in '),ivec
            stop 'length '                                               4d12s19
           end if                                                        4d12s19
           read(line200(is:ie),*)rjac(ivec,1)                           3d23s23
           rjac(ivec,1)=rjac(ivec,1)*xlength                            8d31s21
          end do                                                         4d12s19
          go to 1101                                                    1d8s19
         end if
         if(line200(is:is+2).eq.'cth'.or.line200(is:is+2).eq.'the')then 3d23s23
          ikind=0                                                        4d12s19
          if(line200(ie:ie).eq.'d')then                                 3d23s23
           ikind=1                                                       4d12s19
          else if(line200(is:is).eq.'c')then                            3d23s23
           ikind=2                                                       4d12s19
          else if(line200(is+4:is+4).ne.'a')then                        3d23s23
           write(6,*)('don''t know how to parse '),line200(is:ie)       3d23s23
           write(6,*)('testing on "'),line200(is+4:is+4),('"')          3d23s23
           stop 'theta'                                                  4d12s19
          end if                                                         4d12s19
          do ivec=1,nvec-1                                               4d12s19
           is=ie+1                                                        4d10s19
           call delim(line200,is,ie)                                    3d23s23
           if(ie.lt.is)then
            write(6,*)('not enough vector thetas read in '),ivec         4d12s19
            stop 'theta'                                                 4d12s19
           end if                                                        4d12s19
           read(line200(is:ie),*)rjac(ivec,2)                           3d23s23
           if(ikind.eq.2)then                                            4d12s19
            rjac(ivec,2)=acos(rjac(ivec,2))                              4d12s19
           else if(ikind.eq.1)then                                       4d12s19
            rjac(ivec,2)=rjac(ivec,2)*rad                                4d12s19
           end if                                                        4d12s19
          end do                                                         4d12s19
          go to 1101                                                    1d8s19
         end if
         if(line200(is:is+2).eq.'phi')then                              3d23s23
          ikind=0                                                        4d12s19
          if(line200(ie:ie).eq.'d')then                                 3d23s23
           ikind=1                                                       4d12s19
          else if(ie+1-is.gt.3)then                                      4d15s19
           write(6,*)('don''t know how to parse '),line200(is:ie)       3d23s23
           stop 'phi'                                                    4d15s19
          end if                                                         4d12s19
          do ivec=1,nvec-2                                               4d12s19
           is=ie+1                                                        4d10s19
           call delim(line200,is,ie)                                    3d23s23
           if(ie.lt.is)then
            write(6,*)('not enough vector phis read in '),ivec           4d12s19
            stop 'phi'                                                   4d12s19
           end if                                                        4d12s19
           read(line200(is:ie),*)rjac(ivec,3)                           3d23s23
           if(ikind.eq.1)then                                            4d12s19
            rjac(ivec,3)=rjac(ivec,3)*rad                                4d12s19
           end if                                                        4d12s19
          end do                                                         4d12s19
          go to 1101                                                    1d8s19
         end if                                                         1d8s19
 1102    continue                                                       1d8s19
         close(unit=43)                                                 1d8s19
         go to 1                                                        1d8s19
        end if                                                           4d10s19
        if(line200(is:is).eq.'r')then                                   3d23s23
         do ivec=1,nvec                                                 4d12s19
          is=ie+1                                                        4d10s19
          call delim(line200,is,ie)                                     3d23s23
          if(ie.lt.is)then
           write(6,*)('not enough vector lengths read in '),ivec
           stop 'length '                                               4d12s19
          end if                                                        4d12s19
          write(6,*)('vector length: '),line200(is:ie)                  3d23s23
          read(line200(is:ie),*)rjac(ivec,1)                            3d23s23
          rjac(ivec,1)=rjac(ivec,1)*xlength                             8d31s21
         end do                                                         4d12s19
         go to 1                                                        4d12s19
        end if
        if(line200(is:is+2).eq.'cth'.or.line200(is:is+2).eq.'the')then  3d23s23
         ikind=0                                                        4d12s19
         if(line200(ie:ie).eq.'d')then                                  3d23s23
          ikind=1                                                       4d12s19
         else if(line200(is:is).eq.'c')then                             3d23s23
          ikind=2                                                       4d12s19
         else if(line200(is+4:is+4).ne.'a')then                         3d23s23
          write(6,*)('don''t know how to parse '),line(is:ie)           4d12s19
          stop 'theta'                                                  4d12s19
         end if                                                         4d12s19
         do ivec=1,nvec-1                                               4d12s19
          is=ie+1                                                        4d10s19
          call delim(line200,is,ie)                                     3d23s23
          if(ie.lt.is)then
           write(6,*)('not enough vector thetas read in '),ivec         4d12s19
           stop 'theta'                                                 4d12s19
          end if                                                        4d12s19
          read(line200(is:ie),*)rjac(ivec,2)                            3d23s23
          if(ikind.eq.2)then                                            4d12s19
           rjac(ivec,2)=acos(rjac(ivec,2))                              4d12s19
          else if(ikind.eq.1)then                                       4d12s19
           rjac(ivec,2)=rjac(ivec,2)*rad                                4d12s19
          end if                                                        4d12s19
         end do                                                         4d12s19
         go to 1                                                        4d12s19
        end if
        if(line200(is:is+2).eq.'phi')then                               3d23s23
         ikind=0                                                        4d12s19
         if(line200(ie:ie).eq.'d')then                                  3d23s23
          ikind=1                                                       4d12s19
         else if(ie+1-is.gt.3)then                                      4d15s19
          write(6,*)('don''t know how to parse '),line200(is:ie)        3d23s23
          stop 'phi'                                                    4d15s19
         end if                                                         4d12s19
         do ivec=1,nvec-2                                               4d12s19
          is=ie+1                                                        4d10s19
          call delim(line200,is,ie)                                     3d23s23
          if(ie.lt.is)then
           write(6,*)('not enough vector phis read in '),ivec           4d12s19
           stop 'phi'                                                   4d12s19
          end if                                                        4d12s19
          read(line200(is:ie),*)rjac(ivec,3)                            3d23s23
          if(ikind.eq.1)then                                            4d12s19
           rjac(ivec,3)=rjac(ivec,3)*rad                                4d12s19
          end if                                                        4d12s19
         end do                                                         4d12s19
          go to 1                                                       4d12s19
        end if
       end if                                                           4d10s19
       if(line(1:1).eq.'Z')then
        nnew=0
        go to 101
       end if
       igot=1
       go to 103
  101  continue
       igot=0
       is=3
       ie=nbr(line,is,80,0)
       if(ie.lt.0d0)stop 'error for nucz '
       read(line(is:ie),*)nucz
  121  continue                                                         4d6s18
       lcart=.true.                                                     4d12s19
       if(nvec.gt.0)then                                                4d12s19
        iok=0                                                           4d12s19
        do i=1,nvec+1                                                   4d12s19
         if(xm(i).ne.0d0)iok=iok+1                                      6d20s24
        end do                                                          4d12s19
        if(iok.ne.nvec+1)then                                           4d12s19
         lcart=.false.                                                  4d12s19
        else                                                            4d12s19
        end if                                                          4d12s19
       end if                                                           4d12s19
       ibunit=5                                                         2d18s10
       if(nvec.eq.0)then                                                4d12s19
        is=1
 1104   continue                                                        12d8s22
        read(5,100,end=2)line
        if(line(1:1).eq.'#')go to 1104                                  12d8s22
  104   continue
          iss=1                                                         1d5s20
          call delim(line,iss,iee)                                      1d5s20
          inok=0                                                        1d5s20
          do izz=iss,iee                                                4d11s22
           if(line(izz:izz).eq.'.')then                                 4d11s22
            inok=inok+1                                                 4d11s22
            go to 1216                                                  4d11s22
           end if                                                       4d11s22
           do i=1,11                                                     2d26s21
            if(line(izz:izz).eq.digit(i))then                           4d11s22
             inok=inok+1                                                4d11s22
             go to 1216                                                 4d11s22
            end if                                                      4d11s22
           end do                                                        1d5s20
           inok=0                                                       4d11s22
           go to 1217                                                   4d11s22
 1216      continue                                                     4d11s22
          end do                                                        4d11s22
 1217     continue                                                      4d11s22
          if(inok.eq.0)then                                             1d5s20
           jdorel4c=0                                                   2d20s20
          end if                                                        1d5s20
        do ixyz=1,3
         ie=nbr(line,is,80,0)
         if(ie.eq.-2)then                                                8d25s17
  443    continue                                                       12d23s22
          read(5,100)line
          if(line(1:1).eq.'&')go to 2                                   2d25s20
          if(line(1:1).eq.'#')go to 443                                 12d23s22
          ndecimal=0                                                       8d25s17
          do i=1,80                                                        8d25s17
           if(line(i:i).eq.'.')ndecimal=ndecimal+1                         8d25s17
          end do                                                           8d25s17
          if(ndecimal.eq.0)go to 1215                                   1d5s20
          is=1
          ie=nbr(line,is,80,0)
         end if
         if(ixyz.eq.1)then
          natom=natom+1
          nnew=nnew+1
          do l=1,4                                                      6d11s19
           icoreg(l,natom)=0                                            6d11s19
           ivalg(l,natom)=0                                             6d11s19
           jcoreg(l,natom)=0                                            6d11s19
           jvalg(l,natom)=0                                             6d11s19
          end do                                                        6d11s19
          if(xm(natom).lt.0d0)then                                      1d11s24
           nuczu=0                                                      1d5s24
          else                                                          1d5s24
           nuczu=nucz                                                   1d5s24
           if(nucz.le.2)then                                             6d11s19
            ivalg(1,natom)=1                                             6d11s19
           else if(nucz.le.10)then                                       6d11s19
            icoreg(1,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
           else if(nucz.le.18)then                                       6d11s19
            icoreg(1,natom)=2                                            6d11s19
            icoreg(2,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
           else if(nucz.le.36)then                                       6d11s19
            icoreg(1,natom)=3                                            6d11s19
            icoreg(2,natom)=2                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
           else if(nucz.le.54)then                                       6d11s19
            icoreg(1,natom)=4                                            6d11s19
            icoreg(2,natom)=3                                            6d11s19
            icoreg(3,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
           else if(nucz.le.86)then                                       6d11s19
            icoreg(1,natom)=5                                            6d11s19
            icoreg(2,natom)=4                                            6d11s19
            icoreg(3,natom)=2                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
            ivalg(4,natom)=1                                             6d11s19
           else                                                          6d11s19
            icoreg(1,natom)=6                                            6d11s19
            icoreg(2,natom)=5                                            6d11s19
            icoreg(3,natom)=3                                            6d11s19
            icoreg(4,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
            ivalg(4,natom)=1                                             6d11s19
           end if                                                        6d11s19
          end if                                                        1d5s24
          atnum(1,natom)=dfloat(nuczu)                                  1d5s24
          ngausa(natom)=0                                                4d8s10
          if(natom.gt.ida)stop 'ida'
         end if
         if(inok.eq.0)go to 1215                                        1d5s20
         read(line(is:ie),*)xcart(ixyz,natom)
         xcart(ixyz,natom)=xcart(ixyz,natom)*xlength                    8d31s21
         is=ie+1
        end do
        go to 104
  103   continue
       else                                                             4d12s19
        do i=1,nahere                                                   4d12s19
         natom=natom+1                                                  4d12s19
          do l=1,4                                                      6d11s19
           icoreg(l,natom)=0                                            6d11s19
           ivalg(l,natom)=0                                             6d11s19
           jcoreg(l,natom)=0                                            6d11s19
           jvalg(l,natom)=0                                             6d11s19
          end do                                                        6d11s19
          write(6,*)('xm of natom = '),natom,(' is '),xm(natom)
          if(xm(natom).lt.0d0)then                                      1d11s24
           nuczu=0                                                      1d5s24
          else                                                          1d5s24
           nuczu=nucz                                                   1d5s24
           if(nucz.le.2)then                                             6d11s19
            ivalg(1,natom)=1                                             6d11s19
           else if(nucz.le.10)then                                       6d11s19
            icoreg(1,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
           else if(nucz.le.18)then                                       6d11s19
            icoreg(1,natom)=2                                            6d11s19
            icoreg(2,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
           else if(nucz.le.36)then                                       6d11s19
            icoreg(1,natom)=3                                            6d11s19
            icoreg(2,natom)=2                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
           else if(nucz.le.54)then                                       6d11s19
            icoreg(1,natom)=4                                            6d11s19
            icoreg(2,natom)=3                                            6d11s19
            icoreg(3,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
           else if(nucz.le.86)then                                       6d11s19
            icoreg(1,natom)=5                                            6d11s19
            icoreg(2,natom)=4                                            6d11s19
            icoreg(3,natom)=2                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
            ivalg(4,natom)=1                                             6d11s19
           else                                                          6d11s19
            icoreg(1,natom)=6                                            6d11s19
            icoreg(2,natom)=5                                            6d11s19
            icoreg(3,natom)=3                                            6d11s19
            icoreg(4,natom)=1                                            6d11s19
            ivalg(1,natom)=1                                             6d11s19
            ivalg(2,natom)=1                                             6d11s19
            ivalg(3,natom)=1                                             6d11s19
            ivalg(4,natom)=1                                             6d11s19
           end if                                                        6d11s19
          end if                                                        1d5s24
         nnew=nnew+1                                                    4d12s19
          atnum(1,natom)=dfloat(nuczu)                                   1d2s24
         ngausa(natom)=0                                                4d12s19
        end do                                                          4d12s19
        igot=0                                                          4d12s19
        read(5,100)line                                                 4d12s19
        go to 1215                                                      4d12s19
       end if                                                           4d12s19
       if(igot.ne.1)then                                                2d18s10
        if(ibunit.eq.4)then                                             2d18s10
         read(ibunit,100,end=1)line                                     2d18s10
        else                                                            2d18s10
         read(ibunit,100,end=2)line                                     2d18s10
         if(line(1:1).eq.'Z')go to 1                                    2d18s10
        end if                                                          2d18s10
       end if
       if(ibunit.eq.5)then                                              2d18s10
        is=1
        ie=nbr(line,is,80,0)
        do i=1,10
         if(line(is:is).eq.digit(i))then                                 2d18s10
          go to 1066                                                     2d18s10
         end if                                                          2d18s10
        end do                                                           2d18s10
        ibunit=4                                                         2d18s10
        open(unit=4,file=bpath(1:ibpath)//line(is:ie),status='old',     2d18s10
     $       err=15)                                                    2d18s10
        read(4,100)line                                                  2d18s10
 1066   continue                                                         2d18s10
       end if                                                           2d18s10
 1215  continue                                                         8d25s17
       is=1                                                             9d11s17
       ie=nbr(line,is,80,1)                                               9d11s17
       if(line(is:is).eq.'#')go to 1                                    5d26s21
       do izz=is,ie                                                     4d11s22
        icheck=ichar(line(izz:izz))                                     4d11s22
        if(icheck.lt.ichar(digit(1)).or.icheck.gt.ichar(digit(10)))then 4d11s22
         go to 1218                                                     4d11s22
        end if                                                          4d11s22
       end do                                                           4d11s22
 1218  continue                                                         4d11s22
       noscont=0                                                        3d29s21
       if(icheck.lt.ichar(digit(1)).or.icheck.gt.ichar(digit(10)))then  2d27s19
        write(6,*)('we are getting basis set info from file '),
     $      line(is:ie)
        ixtr=ibcoff                                                     12d12s19
        mxtr=200                                                        12d12s19
        jxtr=ixtr+mxtr                                                  12d12s19
        icvxtr=jxtr+mxtr                                                12d13s22
        ibcoff=icvxtr+mxtr                                              12d13s22
        nxtr=0                                                          12d12s19
        call enough('basisz.  1',bc,ibc)
        open(unit=4,file=pwd//line(is:ie),status='old',iostat=ios)      2d27s19
        do jj=1,11                                                      10d2s24
         lkeep(jj)=1                                                    3d5s19
         ladd(jj)=0                                                     3d5s19
         inoc(jj)=0                                                     3d5s19
        end do                                                          3d5s19
        js=ie+1                                                         3d5s19
 6610   continue                                                        3d5s19
        if(js.lt.0)then                                                 2d20s20
         je=js-1                                                        2d20s20
        else                                                            2d20s20
         je=nbr(line,js,80,0)                                              3d5s19
        end if                                                          2d20s20
        if(je.ge.js.or.(idorel4c.ne.0.and.jdorel4c.eq.0))then           2d20s20
         if(je.ge.js)then                                               12d2s22
          write(6,*)('modifications to this basis: '),line(js:je)       12d2s22
         end if                                                         12d2s22
         if(line(js:js).eq.'&')then                                     12d12s19
          icv=0                                                         12d13s22
          do j=js+1,je-1                                                12d13s22
           if(line(j:j+1).eq.'cv')icv=1                                 12d13s22
          end do                                                        12d13s22
          if(line(js:js+5).eq.'&dwscv')then                             3d29s21
           write(6,*)('turning off contraction of s functions ')        3d29s21
           inoc(1)=-1                                                   12d12s22
           noscont=1                                                    12d12s22
          end if                                                        3d29s21
          write(6,*)('looking for additional functions in the file ')   12d12s19
          line(js:js)='.'                                               12d12s19
          write(6,*)('"'),line(is:ie),line(js:je),('"')                 12d12s19
          open(unit=44,file=pwd//line(is:ie)//line(js:je),status='old', 12d12s19
     $         iostat=ios)                                              12d12s19
          call getxtr(44,element(nucz),nxtr,bc(ixtr),ibc(jxtr),mxtr,    12d13s22
     $         ibc(icvxtr),icv)                                         12d13s22
          close(unit=44)                                                12d12s19
          js=je+1                                                        3d5s19
          go to 6610                                                     3d5s19
         end if                                                         12d12s19
         if(line(js:js).eq.'u'.or.(idorel4c.ne.0.and.jdorel4c.eq.0))then3d27s20
          if(je.eq.js.or.idorel4c.ne.0)then                             2d20s20
           write(6,*)('all of basis will be uncontracted ')             3d5s19
           jdorel4c=jdorel4c+1                                          2d20s20
           do jj=1,8                                                    3d5s19
            inoc(jj)=-1                                                 12d19s22
           end do                                                       3d5s19
          else                                                          3d5s19
           do jj=js+1,je                                                3d5s19
            write(6,*)('uncontract fcns of symmetry '),line(jj:jj)      3d5s19
            if(line(jj:jj).eq.'s'.or.line(jj:jj).eq.'S')then             3d5s19
             inoc(1)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'p'.or.line(jj:jj).eq.'P')then       3d5s19
             inoc(2)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'d'.or.line(jj:jj).eq.'D')then        3d5s19
             inoc(3)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'f'.or.line(jj:jj).eq.'F')then        3d5s19
             inoc(4)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'g'.or.line(jj:jj).eq.'G')then        3d5s19
             inoc(5)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'h'.or.line(jj:jj).eq.'H')then        3d5s19
             inoc(6)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'i'.or.line(jj:jj).eq.'I')then        3d5s19
             inoc(7)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'k'.or.line(jj:jj).eq.'K')then        3d5s19
             inoc(8)=-1                                                 3d5s19
            else if(line(jj:jj).eq.'l'.or.line(jj:jj).eq.'L')then       8d7s24
             inoc(9)=-1                                                 8d7s24
            else if(line(jj:jj).eq.'m'.or.line(jj:jj).eq.'M')then       8d7s24
             inoc(10)=-1                                                8d7s24
            else if(line(jj:jj).eq.'n'.or.line(jj:jj).eq.'N')then       10d2s24
             inoc(11)=-1                                                10d2s24
            else
             write(6,*)('unknown fcn type: '),line(jj:jj)                3d5s19
             stop                                                        3d5s19
            end if                                                       3d5s19
           end do                                                       3d5s19
          end if                                                        3d5s19
         end if                                                         3d27s20
         if(line(js:js).eq.'-')then                                     10d25s20
          write(6,*)('eliminate the functions '),line(js+1:je)          3d5s19
          do jj=js+1,je                                                 3d5s19
           if(line(jj:jj).eq.'p'.or.line(jj:jj).eq.'P')then             3d5s19
            lkeep(1)=-1                                                 3d5s19
           else if(line(jj:jj).eq.'d'.or.line(jj:jj).eq.'D')then        3d5s19
            lkeep(2)=-1                                                 3d5s19
           else if(line(jj:jj).eq.'f'.or.line(jj:jj).eq.'F')then        3d5s19
            lkeep(3)=-1                                                 3d5s19
           else if(line(jj:jj).eq.'g'.or.line(jj:jj).eq.'G')then        3d5s19
            lkeep(4)=-1                                                 3d5s19
           else if(line(jj:jj).eq.'h'.or.line(jj:jj).eq.'H')then        3d5s19
            lkeep(5)=-1                                                 3d5s19
           else if(line(jj:jj).eq.'i'.or.line(jj:jj).eq.'I')then        3d5s19
            lkeep(6)=-1                                                 3d5s19
           else if(line(jj:jj).eq.'k'.or.line(jj:jj).eq.'K')then        3d5s19
            lkeep(7)=-1                                                 3d5s19
           else if(line(jj:jj).eq.'l'.or.line(jj:jj).eq.'L')then        3d5s19
            lkeep(8)=-1                                                 8d7s24
           else if(line(jj:jj).eq.'m'.or.line(jj:jj).eq.'M')then        8d7s24
            lkeep(9)=-1                                                 8d7s24
           else if(line(jj:jj).eq.'n'.or.line(jj:jj).eq.'N')then        10d2s24
            lkeep(10)=-1                                                10d2s24
           else
            write(6,*)('unknown fcn type: '),line(jj:jj)                3d5s19
            stop                                                        3d5s19
           end if                                                       3d5s19
          end do                                                        3d5s19
         else if(line(js:js).eq.'+')then                                3d5s19
          write(6,*)('we are adding functions ')                        3d5s19
          nadd=1                                                        3d5s19
          if(je.eq.js)then                                              3d5s19
           write(6,*)('we are adding 1 diffuse fcn for each l')         3d5s19
           do jj=1,10                                                   8d7s24
            ladd(jj)=ladd(jj)+1                                         3d5s19
           end do                                                       3d5s19
          else                                                          3d5s19
           jsp=js+1                                                     3d5s19
           ival=ichar(line(jsp:jsp))                                    3d5s19
           if(ival.ge.ichar('0').and.ival.le.ichar('9'))then            3d5s19
            read(line(jsp:jsp),*)nadd                                   3d5s19
            write(6,*)('we are adding '),nadd,(' fcns ')                3d5s19
            jsp=jsp+1                                                   3d5s19
            if(jsp.gt.je)then                                           3d5s19
             write(6,*)('for each l')                                   3d5s19
             do jj=1,10                                                 8d7s24
              ladd(jj)=ladd(jj)+nadd                                    3d5s19
             end do                                                     3d5s19
            end if                                                      3d5s19
           end if                                                       3d5s19
           do jj=jsp,je                                                 3d5s19
            if(line(jj:jj).eq.'s'.or.line(jj:jj).eq.'S')then            3d5s19
             ladd(1)=ladd(1)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'p'.or.line(jj:jj).eq.'P')then            3d5s19
             ladd(2)=ladd(2)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'d'.or.line(jj:jj).eq.'D')then            3d5s19
             ladd(3)=ladd(3)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'f'.or.line(jj:jj).eq.'F')then            3d5s19
             ladd(4)=ladd(4)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'g'.or.line(jj:jj).eq.'G')then            3d5s19
             ladd(5)=ladd(5)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'h'.or.line(jj:jj).eq.'H')then            3d5s19
             ladd(6)=ladd(6)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'i'.or.line(jj:jj).eq.'I')then            3d5s19
             ladd(7)=ladd(7)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'k'.or.line(jj:jj).eq.'K')then            3d5s19
             ladd(8)=ladd(8)+nadd                                       3d5s19
            else if(line(jj:jj).eq.'l'.or.line(jj:jj).eq.'L')then       8d7s24
             ladd(9)=ladd(9)+nadd                                       8d7s24
            else if(line(jj:jj).eq.'n'.or.line(jj:jj).eq.'N')then       8d7s24
             ladd(10)=ladd(10)+nadd                                     8d7s24
            else if(line(jj:jj).eq.'m'.or.line(jj:jj).eq.'M')then       10d2s24
             ladd(11)=ladd(11)+nadd                                     10d2s24
            else                                                        3d5s19
             write(6,*)('unknown fcn type: '),line(jj:jj)               3d5s19
             stop                                                       3d5s19
            end if                                                      3d5s19
           end do                                                       3d5s19
          end if                                                        3d5s19
         end if                                                         3d5s19
         if(je.ge.js)then                                               12d2s22
          js=je+1                                                        3d5s19
          go to 6610                                                     3d5s19
         end if                                                         12d2s22
        end if                                                          3d5s19
        ise=1                                                            2d27s19
        iee=nbr(element(nucz),ise,3,0)                                    2d27s19
        ndige=iee+1-ise                                                 2d27s19
        if(ios.ne.0)stop 'basis iostat'
        last=-1                                                         2d27s19
11215   continue                                                        2d27s19
         read(4,200,end=1008)line200                                    5d15s19
31215    continue                                                       2d27s19
         is=1                                                           2d27s19
         ie=nbr(line200,is,200,1)                                         5d15s19
         ndig=ie+1-is                                                   2d27s19
         if(line200(is:is).eq.'#'.and.last.ge.0.or.                     3d1s21
     $        line200(is:is+2).eq.'END')then                            3d1s21
          close(unit=4)                                                 2d27s19
          kzeta=jzeta                                                   2d27s19
          do l=1,numl                                                   2d27s19
           lu=l-1                                                       2d27s19
           if(lu.eq.0)then                                              3d5s19
            llkeep=1                                                    3d5s19
           else                                                         3d5s19
            llkeep=lkeep(lu)                                            3d5s19
           end if                                                       3d5s19
           nrow=ibc(inzeta+lu)                                          2d27s19
           ncol=1+ibc(incont+lu)                                        2d27s19
           if(llkeep.gt.0)then                                          3d5s19
            nzeta=nrow                                                   2d27s19
            if(ncol.eq.nrow+1)then                                       2d27s19
             if(inoc(lu+1).eq.0.and.idorel.ne.0)then                    10d28s20
              ncont=1                                                   10d27s20
             else                                                       10d27s20
              ncont=0                                                     2d27s19
             end if                                                     10d27s20
            else                                                         2d27s19
             ncont=ncol-1                                                2d27s19
            end if                                                       2d27s19
            if(ncont.gt.0.and.idorel.ne.0)then                          1d15s19
c
c     the current value of ncont is basis on non-relativistic contraction
c     it may not be accurate in the relativistic case, for the uncontracted
c     functions need to be put in twice (one for L, one for S).         1d15s19
             do icol=1,nzeta*2                                          1d15s19
              iad=kzeta+nzeta*icol                                      1d15s19
              sz=0d0                                                    1d15s19
              do i=0,nzeta-1                                            1d15s19
               sz=sz+bc(iad+i)**2                                       1d15s19
              end do                                                    1d15s19
              sz=sqrt(sz/dfloat(nzeta))                                 1d15s19
              if(sz.ne.0d0)ncontgg=icol                                   1d15s19
             end do                                                     1d15s19
             if(ncontgg.lt.-10.or.ncontgg.gt.100)then
              write(6,*)('bad value!!!')
              irtrn=1
              return
             end if                                                     8d16s22
c
c     now consider any extra functions
c
              if(nxtr.gt.0.and.ladd(l).eq.0)then                        12d13s19
               nplus=0                                                  12d12s19
               do ixt=0,nxtr-1                                            12d12s19
                if(ibc(jxtr+ixt).eq.l-1)nplus=nplus+1                   12d12s19
               end do                                                   12d12s19
               if(nplus.gt.0)then                                       12d12s19
                ixpl=nzeta-nplus                                        1d16s20
                kzetap=kzeta+ixpl+nzeta*(ncontgg+1)                     1d16s20
                lzeta=kzeta+ixpl                                        1d16s20
                do ixt=0,nxtr-1                                           12d12s19
                 if(ibc(jxtr+ixt).eq.l-1)then                           12d12s19
                  bc(lzeta)=bc(ixtr+ixt)                                1d16s20
                  lzeta=lzeta+1                                         1d16s20
                  if(ibc(icvxtr+ixt).eq.0)then                          12d13s22
                   bc(kzetap)=1d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   bc(kzetap)=1d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   ncontgg=ncontgg+2                                     10d27s20
                   kzetap=kzetap+1                                      12d13s22
                  else                                                  12d13s22
                   bc(kzetap)=1d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   bc(kzetap)=0d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   ncontgg=ncontgg+2                                     10d27s20
                   bc(kzetap)=0d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   bc(kzetap)=1d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   ncontgg=ncontgg+2                                     10d27s20
                   kzetap=kzetap+1                                      12d13s22
                  end if                                                12d13s22
                 end if                                                 12d12s19
                end do                                                  12d12s19
               end if                                                   12d12s19
              else if(nxtr.gt.0)then                                    12d13s19
               nplus=0                                                  12d12s19
               do ixt=0,nxtr-1                                            12d12s19
                if(ibc(jxtr+ixt).eq.l-1)nplus=nplus+1                   12d12s19
               end do                                                   12d12s19
               ixpl=nzeta-nplus-ladd(l)                                 1d22s20
               kzetap=kzeta+ixpl+nzeta*(ncontgg+1)                      1d22s20
               lzeta=kzeta+ixpl                                         1d22s20
               if(nplus.gt.0)then                                       12d12s19
                do ixt=0,nxtr-1                                           12d12s19
                 if(ibc(jxtr+ixt).eq.l-1)then                           12d12s19
                  if(ibc(icvxtr+ixt).eq.0)then                          12d13s22
                   bc(lzeta)=bc(ixtr+ixt)                                1d22s20
                   lzeta=lzeta+1                                         1d22s20
                   bc(kzetap)=1d0                                        1d22s20
                   kzetap=kzetap+nzeta                                   10d27s20
                   bc(kzetap)=1d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   ncontgg=ncontgg+2                                     10d27s20
                   kzetap=kzetap+1                                       12d13s22
                  else                                                  12d13s22
                   bc(lzeta)=bc(ixtr+ixt)                                1d22s20
                   lzeta=lzeta+1                                         1d22s20
                   bc(kzetap)=1d0                                        1d22s20
                   kzetap=kzetap+nzeta                                   10d27s20
                   bc(kzetap)=0d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   ncontgg=ncontgg+2                                     10d27s20
                   bc(kzetap)=0d0                                        1d22s20
                   kzetap=kzetap+nzeta                                   10d27s20
                   bc(kzetap)=1d0                                        12d12s19
                   kzetap=kzetap+nzeta                                   10d27s20
                   ncontgg=ncontgg+2                                     10d27s20
                   kzetap=kzetap+1                                       12d13s22
                  end if                                                12d13s22
                 end if                                                 12d12s19
                end do                                                  12d12s19
               end if                                                   12d12s19
               ixpl=ixpl+nplus                                          12d13s19
               isorte=ibcoff                                            2d14s23
               isortf=isorte+ixpl                                       2d14s23
               ibcoff=isortf+ixpl                                       2d14s23
               do ixt=0,ixpl-1                                          2d14s23
                bc(isorte+ixt)=bc(kzeta+ixt)                            2d14s23
               end do                                                   2d14s23
               call dsortdws(bc(isorte),ibc(isortf),ixpl)               2d14s23
               xlast=bc(isorte)                                         2d14s23
               xpen=bc(isorte+1)                                        2d14s23
               ibcoff=isorte                                            2d14s23
               if(ixpl.gt.0)then                                        3d5s19
                ixpp=ixpl-1                                             3d5s19
                dfr=xpen/xlast                                          2d14s23
                write(6,*)('last two exponential parameters '),
     $               xpen,xlast                                         2d14s23
                write(6,*)('ratio = '),dfr                              3d5s19
               else                                                     3d5s19
                dfr=2.5d0
                write(6,*)('using default diffuse exponent ratio of '),
     $               dfr
               end if                                                   3d5s19
               kzetap=kzeta+ixpl+nzeta*(ncontgg+1)                      1d22s20
               do i=1,ladd(l)                                           3d5s19
                xlast=xlast/dfr                                         3d5s19
                bc(lzeta)=xlast                                         1d22s20
                lzeta=lzeta+1                                           1d22s20
                bc(kzetap)=1d0                                          1d22s20
                kzetap=kzetap+nzeta                                     10d27s20
                bc(kzetap)=1d0                                          1d22s20
                kzetap=kzetap+nzeta+1                                   6d7s21
                ncontgg=ncontgg+2                                       10d27s20
               end do                                                   3d5s19
              else if(ladd(l).eq.0)then                                 12d12s19
              else                                                      3d5s19
               ixpl=nzeta-ladd(l)                                       1d22s20
               lzeta=kzeta+ixpl                                         1d22s20
               if(ixpl.gt.0)then                                        3d5s19
                ixpp=ixpl-1                                             3d5s19
                dfr=bc(lzeta-2)/bc(lzeta-1)                             1d22s20
               else                                                     3d5s19
                dfr=2.5d0
                write(6,*)('using default diffuse exponent ratio of '),
     $               dfr
               end if                                                   3d5s19
               xlast=bc(lzeta-1)                                        1d22s20
               kzetap=kzeta+ixpl+nzeta*(ncontgg+1)                      1d22s20
               do i=1,ladd(l)                                           3d5s19
                xlast=xlast/dfr                                         3d5s19
                bc(lzeta)=xlast                                         1d22s20
                lzeta=lzeta+1                                           1d22s20
                bc(kzetap)=1d0                                          1d22s20
                kzetap=kzetap+nzeta                                     10d27s20
                bc(kzetap)=1d0                                          1d22s20
                kzetap=kzetap+nzeta                                     10d27s20
                ncontgg=ncontgg+2                                       10d27s20
                kzetap=kzetap+1                                         5d13s21
               end do                                                   3d5s19
              end if                                                    3d5s19
             ncont=ncontgg/2                                            1d15s19
            end if                                                      1d15s19
            if(idorel.eq.0)then                                         1d17s20
             nwcont=nwcont+ncont*nzeta                                    2d27s19
            else                                                        1d17s20
             nwcont=nwcont+ncont*nzeta*2                                1d17s20
            end if                                                      1d17s20
            igot=0
            noff=natom-nnew                                                  1d22s10
            do ia=1,nnew                                                     1d22s10
             ipack=ipack+1
             lmax=max(lmax,lu)
             ippack(ipack)=ibcoff
             ibcoff=ibcoff+2*(nzeta+1)+4
             nzetau=nzeta                                                    9d12s17
             if(idorel.ne.0)nzetau=nzeta*2                                   9d12s17
             ibcoff=ibcoff+1+ncont*nzetau                                    9d12s17
             call enough('basisz.  2',bc,ibc)
             isto=ippack(ipack)
             bc(isto)=dfloat(lu)
             isto=isto+1
             bc(isto)=dfloat(nzeta)
             isto=isto+1
             bc(isto)=dfloat(ncont)                                          9d11s17
             isto=isto+1                                                     9d11s17
             if(nvec.gt.0)then
              mycart(ipack,1)=isto                                       4d12s19
              mycart(ipack,2)=noff+ia                                   4d12s19
             end if
             do i=1,3
              bc(isto)=xcart(i,ia+noff)                                       1d22s10
              isto=isto+1
             end do                                                          1d22s10
             jsto=isto+nzeta
             if(idorel.ne.0.and.ncont.ne.0)then                         1d17s20
              isto1=isto                                                1d16s20
              do i=0,nzeta-1                                            1d16s20
               bc(isto)=bc(kzeta+i)                                     1d16s20
               isto=isto+1                                              1d16s20
              end do                                                    1d16s20
              isto=isto+nzeta                                           1d16s20
              lzeta=kzeta+nzeta                                         1d16s20
              do i=0,nzetau*ncont-1                                     1d16s20
               bc(isto)=bc(lzeta+i)                                     1d16s20
               isto=isto+1                                              1d16s20
              end do                                                    1d16s20
             else                                                       1d16s20
              if(ia.eq.1)then                                           1d16s20
               if(nxtr.gt.0.and.ladd(l).eq.0)then                        12d13s19
                nplus=0                                                  12d12s19
                do ixt=0,nxtr-1                                            12d12s19
                 if(ibc(jxtr+ixt).eq.l-1)nplus=nplus+1                   12d12s19
                end do                                                   12d12s19
                 do ixt=0,nzeta-nplus-1                                  12d13s19
                  bc(isto+ixt)=bc(kzeta+ixt)                                 12d12s19
                 end do                                                  12d12s19
                 ixpl=nzeta-1-nplus                                      12d13s19
                 istox=isto+ixpl                                         12d12s19
                 kzetap=kzeta+ixpl+nzeta*(ncont-nplus)                   12d12s19
                 do ixt=0,nxtr-1                                           12d12s19
                  if(ibc(jxtr+ixt).eq.l-1)then                           12d12s19
                   istox=istox+1                                         12d12s19
                   bc(istox)=bc(ixtr+ixt)                                  12d12s19
                   kzetap=kzetap+(nzeta+1)                               12d12s19
                   bc(kzetap)=1d0                                        12d12s19
                  end if                                                 12d12s19
                 end do                                                  12d12s19
               else if(nxtr.gt.0)then                                    12d13s19
                nplus=0                                                  12d12s19
                do ixt=0,nxtr-1                                            12d12s19
                 if(ibc(jxtr+ixt).eq.l-1)nplus=nplus+1                   12d12s19
                end do                                                   12d12s19
                ixpl=nzeta-1-nplus-ladd(l)                               12d13s19
                kzetap=kzeta+ixpl+nzeta*(ncont-nplus-ladd(l))            12d13s19
                if(nplus.gt.0)then                                       12d12s19
                 do ixt=0,nzeta-nplus-ladd(l)-1                            12d12s19
                  bc(isto+ixt)=bc(kzeta+ixt)                                 12d12s19
                 end do                                                  12d12s19
                 istox=isto+ixpl                                         12d12s19
                 do ixt=0,nxtr-1                                           12d12s19
                  if(ibc(jxtr+ixt).eq.l-1)then                           12d12s19
                   istox=istox+1                                         12d12s19
                   bc(istox)=bc(ixtr+ixt)                                  12d12s19
                   kzetap=kzetap+(nzeta+1)                               12d12s19
                   bc(kzetap)=1d0                                        12d12s19
                  end if                                                 12d12s19
                 end do                                                  12d12s19
                end if                                                   12d12s19
                ixpl=ixpl+nplus                                          12d13s19
                ixpland=ixpl+1                                          4d21s21
                isorte=ibcoff                                           4d21s21
                isortf=isorte+ixpland                                   4d21s21
                ibcoff=isortf+ixpland                                   4d21s21
                do ixt=0,ixpl                                           4d21s21
                 bc(isorte+ixt)=bc(isto+ixt)                            4d21s21
                end do                                                  4d20s21
                call dsortdws(bc(isorte),ibc(isortf),ixpland)           1d18s23
                xlast=bc(isorte)                                        4d21s21
                xpen=bc(isorte+1)                                       4d21s21
                ibcoff=isorte                                           4d21s21
                if(ixpl.gt.0)then                                        3d5s19
                 ixpp=ixpl-1                                             3d5s19
                 dfr=xpen/xlast                                         4d21s21
                else                                                     3d5s19
                 dfr=2.5d0
                 write(6,*)('using default diffuse exponent ratio of '),
     $               dfr
                end if                                                   3d5s19
                kzetap=kzeta+ixpl+nzeta*(ncont-ladd(l))                  3d5s19
                do i=1,ladd(l)                                           3d5s19
                 xlast=xlast/dfr                                         3d5s19
                 bc(isto+ixpl+i)=xlast                                   3d5s19
                 if(ncont.gt.0)then                                     8d16s22
                  iad=kzetap+i*(nzeta+1)                                  3d5s19
                  bc(iad)=1d0                                             3d5s19
                 end if                                                 8d16s22
                end do                                                   3d5s19
               else if(ladd(l).eq.0)then                                 12d12s19
                do i=0,nzeta-1                                             2d27s19
                 bc(isto+i)=bc(kzeta+i)                                    2d27s19
                end do                                                     2d27s19
               else                                                      3d5s19
                do i=0,nzeta-1-ladd(l)                                   3d5s19
                 bc(isto+i)=bc(kzeta+i)                                  3d5s19
                end do                                                   3d5s19
                ixpl=nzeta-1-ladd(l)                                     3d5s19
                ixppp=ixpl+1                                            4d22s21
                isorte=ibcoff                                           4d22s21
                isortf=isorte+ixppp                                     4d22s21
                ibcoff=isortf+ixppp                                     4d22s21
                call enough('basisz.  3',bc,ibc)
                do i=0,ixppp-1                                          4d22s21
                 bc(isorte+i)=bc(isto+i)                                4d22s21
                end do                                                  4d22s21
                call dsortdws(bc(isorte),bc(isortf),ixppp)              1d18s23
                xlast=bc(isorte)                                        4d22s21
                xpen=bc(isorte+1)                                       4d22s21
                ibcoff=isorte                                           4d22s21
                if(ixpl.gt.0)then                                        3d5s19
                 ixpp=ixpl-1                                             3d5s19
                 dfr=xpen/xlast                                         4d22s21
                else                                                     3d5s19
                 dfr=2.5d0
                 write(6,*)('using default diffuse exponent ratio of '),
     $               dfr
                end if                                                   3d5s19
                kzetap=kzeta+ixpl+nzeta*(ncont-ladd(l))                  3d5s19
                do i=1,ladd(l)                                           3d5s19
                 xlast=xlast/dfr                                         3d5s19
                 bc(isto+ixpl+i)=xlast                                   3d5s19
                 iad=kzetap+i*(nzeta+1)                                  3d5s19
                 bc(iad)=1d0                                             3d5s19
                end do                                                   3d5s19
               end if                                                    3d5s19
               lzeta=kzeta+nzeta                                          2d27s19
               isto1=isto                                                     1d22s10
               isto=isto+nzeta*2                                              9d12s17
               if(ncont.gt.0)then                                             9d11s17
                do icont=1,ncont                                              9d11s17
                 do i=0,nzetau-1                                           2d27s19
                  bc(isto+i)=bc(lzeta+i)                                   2d27s19
                 end do
                 isto=isto+nzetau                                             9d12s17
                 lzeta=lzeta+nzetau                                       2d27s19
                end do
               end if                                                         9d11s17
              else                                                            1d22s10
               ncopy=nzetau*ncont+nzeta*2                                     9d12s17
               do i=0,ncopy-1                                                 9d12s17
                bc(isto+i)=bc(isto1+i)                                        1d22s10
               end do                                                         1d22s10
              end if                                                          1d22s10
             end if                                                     1d16s20
             if(nvec.gt.0)then
              lm=l-1                                                     12d20s19
              if(l.eq.1)then                                            5d25s21
               write(6,*)                                               5d25s21
               write(6,1327)                                            6d12s23
 1327          format(2x,'l',2x,'Z',6x,'x',14x,'y',14x,'z',15x,         1d11s24
     $              'exponential parameters')                           6d12s23
              end if                                                    5d25s21
              write(6,1103)lm,nucz,(bc(isto1+i),i=0,nzeta-1)            12d20s19
             else                                                       4d12s19
              lm=lu-1                                                   12d20s19
              if(lu.eq.0)then
               write(6,*)                                               5d25s21
               write(6,1327)                                            6d12s23
              end if                                                    5d25s21
              write(6,3)lu,nucz,(xcart(i,ia+noff),i=1,3),
     $            (bc(isto1+i),i=0,nzeta-1)
             end if                                                     4d12s19
             if(ncont.gt.0)then                                              9d11s17
              ksto=isto1+2*nzeta                                             9d12s17
              if(idorel.eq.0)then                                            9d12s17
               do icont=1,ncont                                               9d11s17
                write(6,303)icont,(bc(ksto+i),i=0,nzeta-1)                    9d11s17
                ksto=ksto+nzeta                                               9d12s17
               end do                                                         9d11s17
              else                                                           9d12s17
               do icont=1,ncont                                              9d12s17
                write(6,304)icont,(bc(ksto+i),i=0,nzeta-1)                   9d12s17
                ksto=ksto+nzeta                                              9d12s17
                write(6,305)(bc(ksto+i),i=0,nzeta-1)                         9d12s17
                ksto=ksto+nzeta                                              9d12s17
               end do                                                        9d12s17
              end if                                                         9d12s17
             end if                                                          9d11s17
             ngaus=ngaus+nzeta
             ngausa(ia+noff)=ngausa(ia+noff)+nzeta                           4d8s10
             do i=1,nzeta
              im=i-1
              bc(jsto+im)=0.25d0*srpi/sqrt((2d0*bc(isto1+im))**3)
             end do
             do il=0,lu-1
              ff=0.5d0*(dfloat(il)+1.5d0)
              do i=1,nzeta
               im=i-1
               bc(jsto+im)=bc(jsto+im)*ff/bc(isto1+im)
              end do
             end do
             do i=1,nzeta
              im=i-1
              bc(jsto+im)=1d0/sqrt(bc(jsto+im))
             end do
            end do                                                           1d22s10
           end if                                                       3d5s19
           if(idorel.eq.0)then                                          1d15s19
            kzeta=kzeta+nrow*ncol                                        2d27s19
           else                                                         1d15s19
            kzeta=kzeta+nrow*(1+4*nrow)                                 1d22s20
           end if                                                       1d15s19
          end do                                                        2d27s19
          go to 1                                                       2d27s19
         end if                                                         2d27s19
         if(ndig.eq.ndige)then                                          2d27s19
          if(element(nucz)(ise:iee).eq.line200(is:ie))then              5d15s19
           is=ie+1                                                      2d27s19
           ie=nbr(line200,is,200,0)                                       5d15s19
           if(line200(is:is).eq.'s'.or.line200(is:is).eq.'S')then       5d15s19
            l=0                                                         2d27s19
           else if(line200(is:is).eq.'p'.or.line200(is:is).eq.'P')then        2d27s19
            l=1                                                         2d27s19
           else if(line200(is:is).eq.'d'.or.line200(is:is).eq.'D')then        2d27s19
            l=2                                                         2d27s19
           else if(line200(is:is).eq.'f'.or.line200(is:is).eq.'F')then        2d27s19
            l=3                                                         2d27s19
           else if(line200(is:is).eq.'g'.or.line200(is:is).eq.'G')then        2d27s19
            l=4                                                         2d27s19
           else if(line200(is:is).eq.'h'.or.line200(is:is).eq.'H')then        2d27s19
            l=5                                                         2d27s19
           else if(line200(is:is).eq.'i'.or.line200(is:is).eq.'I')then        2d27s19
            l=6                                                         2d27s19
           else if(line200(is:is).eq.'k'.or.line200(is:is).eq.'K')then        3d5s19
            l=7                                                         3d5s19
           else if(line200(is:is).eq.'l'.or.line200(is:is).eq.'L')then  8d7s24
            l=8                                                         8d7s24
           else if(line200(is:is).eq.'m'.or.line200(is:is).eq.'M')then  8d7s24
            l=9                                                         8d7s24
           else if(line200(is:is).eq.'n'.or.line200(is:is).eq.'N')then  8d7s24
            l=10                                                        10d2s24
           else                                                         2d27s19
            write(6,*)('don''t know how to decode '),line200(is:is)
            stop
           end if
           if(l.ne.last)then                                            2d27s19
            if(l.eq.0)then                                              2d27s19
             is=1                                                       2d27s19
             do i=1,5                                                   2d27s19
              ie=nbr(bline200,is,200,1)                                   5d15s19
              if(bline200(is:is).eq.'(')then                               2d27s19
               ncom=0                                                   2d27s19
               do j=is+1,ie-1                                            2d27s19
                if(bline200(j:j).eq.',')then                               2d27s19
                 bline200(j:j)=' '                                         2d27s19
                 ncom=ncom+1                                            2d27s19
                end if                                                  2d27s19
               end do                                                   2d27s19
               numl=ncom+1                                              2d27s19
               inzeta=ibcoff                                              2d27s19
               incont=inzeta+numl                                         2d27s19
               incread=incont+numl                                      5d15s19
               ibcoff=incread+numl                                      5d15s19
               call enough('basisz.  4',bc,ibc)
               js=is+1                                                  2d27s19
               jnzeta=inzeta                                            2d27s19
               do j=1,numl
                je=nbr(bline200,js,ie-1,0)                                   2d27s19
                jem=je-1                                                2d27s19
                read(bline200(js:jem),*)ibc(jnzeta)                        2d27s19
                ibc(jnzeta)=ibc(jnzeta)+ladd(j)                         3d5s19
                if(nxtr.gt.0)then                                       12d12s19
                 jm=j-1                                                 12d12s19
                 do ixt=0,nxtr-1                                          12d12s19
                  if(ibc(jxtr+ixt).eq.jm)ibc(jnzeta)=ibc(jnzeta)+1        12d12s19
                 end do                                                 12d12s19
                end if                                                  12d12s19
                jnzeta=jnzeta+1                                         2d27s19
                js=je+1                                                 2d27s19
               end do
              else if(bline200(is:is).eq.'[')then                               2d27s19
               js=is+1                                                  2d27s19
               jncont=incont                                            2d27s19
               jncread=incread                                          5d15s19
               do j=is,ie                                               2d27s19
                if(bline200(j:j).eq.',')bline200(j:j)=' '                     2d27s19
               end do                                                   2d27s19
               do j=1,numl
                je=nbr(bline200,js,ie-1,0)                                   2d27s19
                jem=je-1                                                2d27s19
                read(bline200(js:jem),*)ibc(jncont)                        2d27s19
                ibc(jncread)=ibc(jncont)                                5d15s19
                ibc(jncont)=ibc(jncont)+ladd(j)                         3d5s19
                if(nxtr.gt.0)then                                       12d12s19
                 jm=j-1                                                 12d12s19
                 do ixt=0,nxtr-1                                        12d12s19
                  if(ibc(jxtr+ixt).eq.jm)then                           12d13s22
                   istrtcv=ibc(jncont)
                   ibc(jncont)=ibc(jncont)+1                            12d13s22
c                                                                       12d13s22
c     if our extra functions are core-valence functions and we are      12d13s22
c     including relativity, we need both a L and S contraction.         12d13s22
c                                                                       12d13s22
                   if(idorel.ne.0.and.ibc(icvxtr+ixt).ne.0)             12d13s22
     $                  ibc(jncont)=ibc(jncont)+1                       12d13s22
                  end if                                                12d13s22
                 end do                                                 12d12s19
                end if                                                  12d12s19
                jncread=jncread+1                                       5d15s19
                jncont=jncont+1                                         2d27s19
                js=je+1                                                 2d27s19
               end do
               do j=0,numl-1                                            3d5s19
                if(inoc(j+1).ne.0)then                                  3d5s19
                 bc(incont+j)=bc(inzeta+j)                              3d5s19
                end if                                                   3d5s19
               end do                                                   3d5s19
              end if                                                    2d27s19
              is=ie+1
             end do                                                     2d27s19
             jzeta=ibcoff                                               2d27s19
             do i=1,numl                                                2d27s19
              im=i-1                                                    2d27s19
              if(idorel.eq.0)then                                       1d15s19
               nwds=ibc(inzeta+im)*(1+ibc(incont+im))                    2d27s19
              else                                                      1d15s19
c                                                                       1d15s19
c     at this point we do not know how many fcns will be uncontracted   1d15s19
c     these have to have 2 contraction coefficients: one for L and one  1d15s19
c     for S. So allocate space for worst case senario: all uncontracted.1d15s19
c
               nwds=ibc(inzeta+im)*(1+4*ibc(inzeta+im))                 1d22s20
              end if
              ibcoff=ibcoff+nwds                                        2d27s19
             end do                                                     2d27s19
             call enough('basisz.  5',bc,ibc)
             do i=jzeta,ibcoff                                          2d27s19
              bc(i)=0d0                                                 2d27s19
             end do                                                     2d27s19
             kzeta=jzeta                                                2d27s19
            else                                                        2d27s19
             if(idorel.eq.0)then                                        1d15s19
              kzeta=kzeta+ibc(inzeta+last)*(1+ibc(incont+last))          2d27s19
             else                                                       1d15s19
              kzeta=kzeta+ibc(inzeta+last)*(1+4*ibc(inzeta+last))       1d22s20
             end if                                                     1d15s19
            end if                                                      2d27s19
            nzzz=0                                                       2d27s19
            nczz=0                                                      2d27s19
            last=l                                                      2d27s19
           end if                                                       2d27s19
           ngroup=0                                                     2d27s19
21215      continue                                                     2d27s19
            read(4,200,end=1008)line200                                 5d15s19
            is=1                                                        2d27s19
            ie=nbr(line200,is,200,0)                                      5d15s19
            ndot=0                                                      2d27s19
            irelcc=0                                                    1d15s19
            do i=is,199                                                 1d11s23
             if(line200(i:i).eq.'.')ndot=ndot+1                         5d15s19
             if(line200(i:i+1).eq.'rc')then                             1d15s19
              if(inoc(l+1).eq.0)irelcc=i+2                              1d17s20
              go to 1891                                                1d15s19
             end if                                                     1d15s19
            end do                                                      2d27s19
 1891       continue                                                    5d15s19
            if(ndot.eq.0)then                                           2d27s19
             nczz=nczz+ngroup                                           2d27s19
             go to 31215                                                2d27s19
            end if                                                      2d27s19
            nzzz=nzzz+1                                                 2d27s19
            if(inoc(l+1).ne.0)ndot=1                                    3d5s19
            read(line200,*)zz,(bc(ibcoff+idot),idot=0,ndot-2)           5d15s19
            if(irelcc.gt.0.and.idorel.ne.0)then                         1d15s19
             nnonr=ndot-1                                               1d15s19
             inonr=ibcoff+nnonr                                         1d15s19
             read(line200(irelcc:200),*)(bc(inonr+idot),idot=0,         1d15s19
     $            nnonr*2-1)                                            1d15s19
            end if                                                      1d15s19
            iad1=kzeta+nzzz-1                                           2d27s19
            bc(iad1)=zz                                                 2d27s19
            if(idorel.ne.0)then                                         1d15s19
             iad1=iad1+nczz*ibc(inzeta+last)*2                          1d15s19
            else                                                        1d15s19
             iad1=iad1+nczz*ibc(inzeta+last)                             2d27s19
            end if                                                      1d15s19
            nhere=ndot-1                                                5d15s19
            if(mod(nhere,3).eq.0.and.nhere.gt.ibc(incread+l))then       5d15s19
             nhere3=nhere/3                                             5d15s19
            else                                                        5d15s19
             nhere3=nhere                                               5d15s19
            end if                                                      5d15s19
            ngroup=min(ndot-1,nhere3)                                   5d15s19
c
c     if relitivistic contraction, then ndot=3*ibc(incread+l)+1,       55d15s19
c     otherwise ndot=ibc(incread+l)+1.                                 55d15s19
c
           if(irelcc.gt.0.and.idorel.ne.0)then                          1d15s19
            do i=1,nnonr*2                                              1d15s19
             im=i-1                                                      2d27s19
             iad1=iad1+ibc(inzeta+last)                                  2d27s19
             bc(iad1)=bc(inonr+im)                                      1d15s19
            end do                                                      1d15s19
           else                                                         1d15s19
            if(idorel.ne.0)then                                         1d15s19
             delta=abs(bc(ibcoff)-1d0)                                  1d15s19
             if(delta.gt.1d-8.and.inoc(last+1).eq.0)then                1d17s20
              write(6,*)
     $            ('we need relativistic contraction coefficients')     1d15s19
              write(6,*)('I did not find any in the basis set file.')   1d15s19
              write(6,*)('I was looking at ')
              write(6,200)line200
              stop                                                      1d15s19
             end if                                                     1d15s19
             nnzeta=ibc(inzeta+last)                                    1d15s19
             nncont=ibc(incont+last)*2
             iad1=iad1+ibc(inzeta+last)                                 1d15s19
             bc(iad1)=1d0                                               1d15s19
             iad1=iad1+ibc(inzeta+last)                                 1d15s19
             bc(iad1)=1d0                                               6d1s21
             nncont=nnzeta                                              6d1s21
            else                                                        1d15s19
             do i=1,min(ndot-1,nhere3)                                    5d15s19
              im=i-1                                                      2d27s19
              iad1=iad1+ibc(inzeta+last)                                  2d27s19
              bc(iad1)=bc(ibcoff+im)                                      2d27s19
             end do                                                       2d27s19
            end if                                                      1d15s19
           end if                                                       1d15s19
           go to 21215                                                  2d27s19
          end if                                                        2d27s19
         else                                                           2d27s19
          bline200=line200                                              5d15s19
         end if                                                         2d27s19
         go to 11215                                                    2d27s19
       else                                                             2d27s19
        read(line(is:ie),*)l                                             9d11s17
        is=ie+1                                                          9d11s17
        ie=nbr(line,is,80,0)                                               9d11s17
        read(line(is:ie),*)nzeta                                         9d11s17
        is=ie+1                                                          9d11s17
        ie=nbr(line,is,80,0)
        if(ie.ge.is)then                                                 9d11s17
         read(line(is:ie),*)ncont                                        9d11s17
         inoc(l+1)=0                                                    12d19s22
        write(6,*)('number of contracted functions = '),ncont           9d11s17
         if(ncont.gt.nzeta)then                                          2d14s19
          write(6,*)('number of contracted functions read in '),ncont,   2d14s19
     $       (' exceeds number of primitive functions '),nzeta          2d14s19
          stop 'ncont gt nzeta'                                          2d14s19
         end if                                                          2d14s19
         if(idorel.eq.0)then                                             9d12s17
          nwcont=nwcont+ncont*nzeta                                      9d12s17
         else                                                            9d12s17
          nwcont=nwcont+ncont*nzeta*2                                    9d12s17
         end if                                                          9d12s17
        else                                                             9d11s17
         ncont=0                                                         9d11s17
         inoc(l+1)=-1                                                   12d19s22
        end if                                                           9d11s17
        igot=0
        noff=natom-nnew                                                  1d22s10
        do ia=1,nnew                                                     1d22s10
         ipack=ipack+1
         lmax=max(lmax,l)
         ippack(ipack)=ibcoff
         ibcoff=ibcoff+2*(nzeta+1)+4
         nzetau=nzeta                                                    9d12s17
         if(idorel.ne.0)nzetau=nzeta*2                                   9d12s17
         ibcoff=ibcoff+1+ncont*nzetau                                    9d12s17
         call enough('basisz.  6',bc,ibc)
         isto=ippack(ipack)
         bc(isto)=dfloat(l)
         isto=isto+1
         bc(isto)=dfloat(nzeta)
         isto=isto+1
         bc(isto)=dfloat(ncont)                                          9d11s17
         isto=isto+1                                                     9d11s17
         if(nvec.gt.0)then                                              12d6s21
          mycart(ipack,1)=isto                                          12d6s21
          mycart(ipack,2)=noff+ia                                       12d6s21
         end if                                                         12d6s21
         do i=1,3
          bc(isto)=xcart(i,ia+noff)                                       1d22s10
          isto=isto+1
         end do                                                          1d22s10
         jsto=isto+nzeta
         if(ia.eq.1)then                                                 1d22s10
          read(ibunit,*)(bc(isto+i),i=0,nzeta-1)                         2d18s10
          isto1=isto                                                     1d22s10
          isto=isto+nzeta*2                                              9d12s17
          if(ncont.gt.0)then                                             9d11s17
           do icont=1,ncont                                              9d11s17
            read(ibunit,*)(bc(isto+i),i=0,nzetau-1)                      9d12s17
            isto=isto+nzetau                                             9d12s17
           end do                                                        9d11s17
          end if                                                         9d11s17
         else                                                            1d22s10
          ncopy=nzetau*ncont+nzeta*2                                     9d12s17
          do i=0,ncopy-1                                                 9d12s17
           bc(isto+i)=bc(isto1+i)                                        1d22s10
          end do                                                         1d22s10
         end if                                                          1d22s10
         lm=l                                                           4d21s23
         lu=lm                                                          12d19s22
         if(l.eq.1)then
          write(6,*)(' ')                                               5d25s21
          write(6,*)(' l  Z       x              y              z'),
     $     ('               exponential parameters ')
         end if                                                         5d25s21
         if(nvec.gt.0)then                                              4d12s19
          write(6,1103)lm,nucz,(bc(isto1+i),i=0,nzeta-1)                12d20s19
 1103     format(2i3,3(6x,'tbd',6x),5x,40es9.2)                         4d12s19
         else                                                           4d12s19
          write(6,3)lm,nucz,(xcart(i,ia+noff),i=1,3),                   12d20s19
     $       (bc(isto1+i),i=0,nzeta-1)
    3    format(2i3,3f15.8,5x,1p80e9.2)                                 8d7s24
         end if                                                         4d12s19
         if(ncont.gt.0)then                                              9d11s17
          ksto=isto1+2*nzeta                                             9d12s17
          if(idorel.eq.0)then                                            9d12s17
           do icont=1,ncont                                               9d11s17
            write(6,303)icont,(bc(ksto+i),i=0,nzeta-1)                    9d11s17
            ksto=ksto+nzeta                                               9d12s17
  303       format(2x,'contraction no. ',i3,35x,40f9.5)                                9d11s17
           end do                                                         9d11s17
          else                                                           9d12s17
           do icont=1,ncont                                              9d12s17
            write(6,304)icont,(bc(ksto+i),i=0,nzeta-1)                   9d12s17
  304       format(2x,'contraction no. ',i3,32x,'k1',x,40f9.5)          1d24s23
            ksto=ksto+nzeta                                              9d12s17
            write(6,305)(bc(ksto+i),i=0,nzeta-1)                         9d12s17
            ksto=ksto+nzeta                                              9d12s17
  305       format(53x,'k2',x,40f9.5)                                   1d24s23
           end do                                                        9d12s17
          end if                                                         9d12s17
         end if                                                          9d11s17
         ngaus=ngaus+nzeta
         ngausa(ia+noff)=ngausa(ia+noff)+nzeta                           4d8s10
         do i=1,nzeta
          im=i-1
          bc(jsto+im)=0.25d0*srpi/sqrt((2d0*bc(isto1+im))**3)
         end do
         do il=0,lu-1                                                   4d11s22
          ff=0.5d0*(dfloat(il)+1.5d0)
          do i=1,nzeta
           im=i-1
           bc(jsto+im)=bc(jsto+im)*ff/bc(isto1+im)
          end do
         end do
         do i=1,nzeta
          im=i-1
          bc(jsto+im)=1d0/sqrt(bc(jsto+im))
         end do
        end do                                                           1d22s10
       end if                                                           2d27s19
       if(ibunit.eq.4)close(unit=4)                                     2d18s10
       go to 1
    2 continue
      if(nvec.gt.0)then                                                 4d12s19
       if(.not.lcart)then                                               4d12s19
        write(6,*)('I''m at end of atoms input, but not all masses ')   4d12s19
     $      ,('are defined. make sure all "ends" string match up '),    4d12s19
     $      ('with an atom')                                            4d12s19
        stop 'lcart'
       end if                                                           4d12s19
       iextradatad=ibcoff                                               1d10s19
       if(nvec.eq.1)then                                                1d10s19
        ibcoff=iextradatad+1                                            1d10s19
       else                                                             1d10s19
        ibcoff=iextradatad+3*nvec-3                                     1d10s19
       end if                                                           1d10s19
c
c     space in case of finite field ...
c
       if(ndfld.gt.0)then                                               12d6s23
        do i=1,ndfld                                                    12d6s23
         bc(ibcoff)=dfloat(iftype(i))                                   12d6s23
         bc(ibcoff+1)=fstgth(i)                                         12d6s23
         ibcoff=ibcoff+2                                                12d6s23
        end do                                                          12d6s23
       end if                                                           12d6s23
       call enough('basisz.  7',bc,ibc)
       nextradatx=ibcoff-iextradatad                                    1d10s19
       call cartsfromv(xm,ends,iends,rjac,ida,mya,nvec,nblba,           1d10s19
     $      bc(iextradatad))                                            5d26s21
       do i=1,ipack                                                     4d12s19
        isto=mycart(i,1)                                                4d12s19
        iato=mycart(i,2)                                                4d12s19
        do ixyz=1,3                                                     4d12s19
         bc(isto)=xcart(ixyz,iato)                                      4d12s19
         isto=isto+1                                                    4d12s19
        end do                                                          4d12s19
       end do                                                           4d12s19
      end if                                                            4d12s19
      numeminus=0                                                       11d16s19
      if(idorel.eq.0)then                                               8d17s15
       write(6,*)('atoms: ')                                             1d22s10
       write(6,311)                                                     3d23s23
  311  format(3x,'#',3x,'Z',13x,'x',14x,'y',14x,'z',2x,'(a.u.)')        3d23s23
      else
       write(6,*)('Nuclei will have finite radii ')
       write(6,*)('atoms: ')                                             1d22s10
       write(6,310)
  310  format(3x,'#',3x,'Z',13x,'x',14x,'y',14x,'z',7x,                 3d23s23
     $      'nuc weight(amu) nuc rad   nuc gaus (a.u.)')                9d12s17
      end if                                                            8d17s15
      potdws=0d0                                                        2d9s12
      ffdws=0d0                                                         8d3s23
      if(idorel.ne.0)then                                               8d17s15
       if(ianuc.lt.0)then                                               8d17s15
        ianuc=ibcoff                                                    8d17s15
        ibcoff=ianuc+natom                                              8d17s15
        call enough('basisz.  8',bc,ibc)
        call mostpop(bc(ianuc),atnum,natom,iretrn,pwd)                  2d1s19
        if(iretrn.ne.0)go to 1008
       end if                                                           8d17s15
      end if                                                            8d17s15
      thrd=1d0/3d0                                                      8d17s15
      sz=0d0                                                            3d23s23
      do ib=1,natom                                                     3d23s23
       do ia=1,natom                                                    3d23s23
        sz=sz+coord(ia,ib)**2                                           3d23s23
       end do                                                           3d23s23
      end do                                                            3d23s23
      sz=sqrt(sz)/dfloat(natom)                                         3d23s23
      if(sz.gt.1d-4)then                                                3d23s23
       wedge=' '                                                        3d23s23
      else                                                              3d23s23
       wedge='>'                                                        3d23s23
      end if                                                            3d23s23
      do i=1,natom                                                      1d22s10
       if(idorel.ne.0)then                                              8d17s15
        anuc=bc(ianuc+i-1)                                              8d17s15
        if(anuc.gt.1822d0)anuc=anuc*5.485 799 110 d-4                   8d17s15
        rrms=0.836d0*(anuc**thrd)+0.57d0
        rrms=rrms*1d-5
        rrms=rrms/0.529177249d0
        etanuc=1.5d0/(rrms**2)
        etanuc=etanuc*0.5d0
        xnormn=1d0/sqrt(0.25d0*srpi/sqrt((2d0*etanuc)**3))
        atnum(2,i)=etanuc                                               8d17s15
        atnum(3,i)=xnormn                                               8d17s15
        nuc=nint(atnum(1,i))                                            4d6s18
        nuce=nuc                                                        1d3s24
        if(nuc.eq.0)nuce=ide                                            1d3s24
        numeminus=numeminus+nuc
        write(6,10)wedge,i,atnum(1,i),element(nuce),                    1d3s24
     $       (xcart(j,i),j=1,3),anuc,rrms,etanuc                        4d6s18
        atom(i)=element(nuce)                                           1d3s24
       else                                                             8d17s15
        nuc=nint(atnum(1,i))                                            4d6s18
        numeminus=numeminus+nuc
        nuce=nuc                                                        1d3s24
        if(nuc.eq.0)nuce=ide                                            1d3s24
        write(6,10)wedge,i,atnum(1,i),element(nuce),(xcart(j,i),j=1,3)  1d3s24
        atom(i)=element(nuce)                                           1d3s24
        atnum(2,i)=0d0                                                  8d17s15
       end if
       do ifi=1,ndfld                                                   12d6s23
        if(iftype(ifi).eq.1)then                                        12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)                12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)                  12d6s23
        else if(iftype(ifi).eq.2)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(2,i)                12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(2,i)                  12d6s23
        else if(iftype(ifi).eq.3)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(3,i)                12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(3,i)                  12d6s23
        else if(iftype(ifi).eq.4)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)*xcart(2,i)     12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)*xcart(2,i)       12d6s23
        else if(iftype(ifi).eq.5)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)*xcart(3,i)     12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)*xcart(3,i)       12d6s23
        else if(iftype(ifi).eq.6)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(2,i)*xcart(3,i)     12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(2,i)*xcart(3,i)       12d6s23
        else if(iftype(ifi).eq.7)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)*xcart(1,i)     12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(1,i)*xcart(1,i)       12d6s23
        else if(iftype(ifi).eq.8)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(2,i)*xcart(2,i)     12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(2,i)*xcart(2,i)       12d6s23
        else if(iftype(ifi).eq.9)then                                   12d6s23
         potdws=potdws+atnum(1,i)*fstgth(ifi)*xcart(3,i)*xcart(3,i)     12d6s23
         ffdws=ffdws+atnum(1,i)*fstgth(ifi)*xcart(3,i)*xcart(3,i)       12d6s23
        end if
       end do                                                           12d6s23
   10  format(a1,i3,f5.0,x,a3,2x,3f15.8,2x,f10.4,2es11.3)               3d23s23
       do j=i+1,natom                                                   2d9s12
         dist=sqrt((xcart(1,i)-xcart(1,j))**2+(xcart(2,i)-xcart(2,j))**2 2d9s12
     $       +(xcart(3,i)-xcart(3,j))**2)                               2d9s12
         potdws=potdws+atnum(1,i)*atnum(1,j)/dist                       8d17s15
       end do                                                           2d9s12
      end do                                                            1d22s10
      write(6,*)(' ')
      if(natom.gt.1)then                                                3d1s21
       write(6,*)('atom-atom distances in a.u. ')                           2d28s21
       do i=1,natom                                                      2d28s21
        do j=i+1,natom                                                   2d28s21
         dist=sqrt((xcart(1,i)-xcart(1,j))**2+(xcart(2,i)-xcart(2,j))**2 2d9s12
     $       +(xcart(3,i)-xcart(3,j))**2)                               2d9s12
         bc(ibcoff+j)=dist                                               2d28s21
        end do                                                           2d28s21
        write(6,3343)(i,j,bc(ibcoff+j),j=i+1,natom)                      2d28s21
 3343   format(20(i2,'-',i2,1x,f14.8))                                  1d19s23
       end do                                                            2d28s21
       write(6,*)(' ')
       write(6,*)('atom-atom distances in Angstroms ')                           2d28s21
       do i=1,natom                                                      2d28s21
        do j=i+1,natom                                                   2d28s21
         dist=sqrt((xcart(1,i)-xcart(1,j))**2+(xcart(2,i)-xcart(2,j))**2 2d9s12
     $       +(xcart(3,i)-xcart(3,j))**2)                               2d9s12
         dist=dist*0.529177249d0                                         2d28s21
         bc(ibcoff+j)=dist                                               2d28s21
        end do                                                           2d28s21
        write(6,3343)(i,j,bc(ibcoff+j),j=i+1,natom)                      2d28s21
       end do                                                            2d28s21
       write(6,*)(' ')
      end if                                                            3d1s21
      write(6,*)('nuclear repulsion energy = '),potdws                  2d9s12
      if(potdws.ne.potdws.or.potdws.gt.1d10)then                        1d2s20
       write(6,*)('oh no! something is wrong with the atomic '),
     $     ('coordinates!')
       irtrn=1                                                          1d2s20
       return                                                           1d2s20
      end if                                                            1d2s20
      write(6,*)(' ')
      first(1)=potdws                                                   2d9s12
c
c     perform symmetry analysis
c                                                                       4d5s10
c     change from sphersub notation to multiplicative form              4d5s10
c                                                                       4d5s10
      do ipass=1,8                                                      4d5s10
       do ixyz=1,3                                                      4d5s10
        if(isym(ixyz,ipass).eq.0)then                                   4d5s10
         isym(ixyz,ipass)=1                                             4d5s10
        else                                                            4d5s10
         isym(ixyz,ipass)=-1                                            4d5s10
        end if                                                          4d5s10
       end do                                                           4d5s10
      end do                                                            4d5s10
      do i=1,8                                                          4d27s10
       do j=1,3                                                         4d27s10
        isymd2h(j,i)=isym(j,i)                                          4d27s10
       end do                                                           4d27s10
      end do                                                            4d27s10
c
c     look for symmetry elements                                        4d5s10
c
      ntwo=0                                                            12d7s23
      do i=1,8                                                          12d7s23
       if(isymu(i).eq.2)ntwo=ntwo+1                                     12d7s23
      end do                                                            12d7s23
      if(ntwo.gt.0)then                                                 12d7s23
       if(ntwo.eq.2)then                                                12d7s23
        if(isymu(2).eq.2.and.isymu(3).eq.2)isymu(4)=2                    12d7s23
        if(isymu(2).eq.2.and.isymu(5).eq.2)isymu(6)=2                    12d7s23
        if(isymu(3).eq.2.and.isymu(5).eq.2)isymu(7)=2                   12d7s23
       else if(ntwo.eq.3)then                                           12d7s23
        isymu(4)=2                                                      12d7s23
        isymu(6)=2                                                      12d7s23
        isymu(7)=2                                                      12d7s23
        isymu(8)=2                                                      12d7s23
       end if                                                           12d7s23
       do i=1,8                                                         12d7s23
        isymu(i)=isymu(i)-1                                             12d7s23
       end do                                                           12d7s23
      end if                                                            12d7s23
      do i=1,natom                                                      4d5s10
       iapair(1,i)=0                                                    4d5s10
       iapair(3,i)=4                                                    4d8s10
      end do                                                            4d5s10
      do i=1,3                                                          4d5s10
       isymt(i)=1                                                       4d5s10
      end do                                                            4d5s10
      ngrp=1                                                            4d5s10
      igrp(1)=1                                                         4d5s10
      do i=1,3
       negs(i)=0
      end do
      do ipass=2,8                                                      4d5s10
       nmatch=0                                                         4d5s10
       do i=1,natom                                                     4d5s10
        do ixyz=1,3                                                     4d5s10
         cartb(ixyz)=xcart(ixyz,i)*dfloat(isym(ixyz,ipass))             4d5s10
        end do                                                          4d5s10
        do j=1,natom                                                    4d5s10
         if(atnum(1,i).eq.atnum(1,j))then                               8d17s15
          diff=0d0                                                      4d5s10
          do ixyz=1,3                                                   4d5s10
           diff=diff+(cartb(ixyz)-xcart(ixyz,j))**2                     4d5s10
          end do                                                        4d5s10
          diff=sqrt(diff/3d0)                                           4d5s10
          if(diff.lt.epsym.and.isymu(ipass).eq.1)then                   12d7s23
           nmatch=nmatch+1                                              4d5s10
           japair(i)=j                                                  4d5s10
          end if                                                        4d5s10
         end if                                                         4d5s10
        end do                                                          4d5s10
       end do                                                           4d5s10
       if(nmatch.eq.natom)then                                          4d5s10
        nneg=0                                                          4d5s10
        do j=1,3                                                        4d5s10
         if(isym(j,ipass).eq.1)then                                     4d5s10
          isymt(j)=isymt(j)+1                                           4d5s10
         else                                                           4d5s10
          negs(j)=1
          nneg=nneg+1                                                   4d5s10
         end if                                                         4d5s10
        end do                                                          4d5s10
 3351   format(a3,x,3i2,5x,11i3)                                        12d24s19
        ngrp=ngrp+1                                                     4d5s10
        igrp(ngrp)=ipass                                                4d5s10
        do i=1,natom                                                    4d5s10
         if(japair(i).ne.i.and.nneg.lt.iapair(3,i))then                 4d8s10
          iapair(1,i)=japair(i)                                         4d5s10
          iapair(2,i)=ngrp                                              4d5s10
          iapair(3,i)=nneg                                              4d8s10
         end if                                                         4d5s10
        end do                                                          4d5s10
       end if                                                           4d5s10
      end do                                                            4d5s10
      write(6,*)('order of group: '),ngrp                               4d5s10
      do i=1,ngrp                                                       4d5s10
       multh(i,i)=1                                                     4d5s10
       do j=i+1,ngrp                                                    4d5s10
        do ixyz=1,3                                                     4d5s10
         mapd2h(ixyz)=isym(ixyz,igrp(i))*isym(ixyz,igrp(j))             4d5s10
        end do                                                          4d5s10
        do k=1,ngrp                                                     4d5s10
         if(mapd2h(1).eq.isym(1,igrp(k)).and.                           4d5s10
     $      mapd2h(2).eq.isym(2,igrp(k)).and.                           4d5s10
     $      mapd2h(3).eq.isym(3,igrp(k)))then                           4d5s10
          multh(i,j)=k                                                  4d5s10
          multh(j,i)=k                                                  4d5s10
          go to 3354                                                    4d5s10
         end if                                                         4d5s10
        end do                                                          4d5s10
        write(6,*)('group is not closed!!!')                            4d5s10
        stop                                                            4d5s10
 3354   continue                                                        4d5s10
       end do                                                           4d5s10
      end do                                                            4d5s10
      do i=1,3                                                          4d5s10
       if(isymt(i).eq.ngrp)then                                         4d5s10
        isymt(i)=0                                                      4d5s10
       else                                                             4d5s10
        isymt(i)=1                                                      4d5s10
       end if                                                           4d5s10
      end do                                                            4d5s10
      write(6,1244)                                                     12d28s19
 1244 format('IRR',3x,'gen.',2x,'multiplication table: ')               12d28s19
      do i=1,ngrp                                                       4d5s10
       do ixyz=1,3                                                      4d5s10
        isym(ixyz,i)=isym(ixyz,igrp(i))                                 4d5s10
       end do                                                           4d5s10
       is=11+ngrp*3+4                                                   12d24s19
       write(opline(1:is),3351)stype(i,ngrp),                           12d24s19
     $      (isym(j,i),j=1,3),(multh(i,j),j=1,ngrp)                     12d24s19
       is=is+1
       opline(is:is+5)='      '
       is=is+6
       if(isym(1,i)*negs(1).eq.-negs(1).and.negs(2)*isym(2,i).eq.negs(2)
     $      .and.negs(3)*isym(3,i).eq.negs(3))then
        opline(is:is+2)=' X'
        ipropsym(1)=i
        is=is+2
       end if
       if(negs(2)*isym(2,i).eq.-negs(2).and.negs(1)*isym(1,i).eq.negs(1)
     $      .and.negs(3)*isym(3,i).eq.negs(3))then
        opline(is:is+2)=' Y'
        ipropsym(2)=i
        is=is+2
       end if
       if(negs(3)*isym(3,i).eq.-negs(3).and.negs(2)*isym(2,i).eq.negs(2)
     $      .and.negs(1)*isym(1,i).eq.negs(1))then
        opline(is:is+2)=' Z'
        ipropsym(3)=i
        is=is+2
       end if
       if(negs(1)*isym(1,i).eq.-negs(1).and.negs(2)*isym(2,i).eq.
     $      -negs(2).and.negs(3)*isym(3,i).eq.negs(3))then
        opline(is:is+3)=' XY'
        ipropsym(4)=i
        is=is+3
       end if
       if(negs(1)*isym(1,i).eq.-negs(1).and.negs(2)*isym(2,i).eq.negs(2)
     $      .and.negs(3)*isym(3,i).eq.-negs(3))then
        opline(is:is+3)=' XZ'
        ipropsym(5)=i
        is=is+3
       end if
       if(negs(2)*isym(2,i).eq.-negs(2).and.negs(1)*isym(1,i).eq.negs(1)
     $      .and.negs(3)*isym(3,i).eq.-negs(3))then
        opline(is:is+3)=' YZ'
        ipropsym(6)=i
        is=is+3
       end if
       is=is-1
       write(6,2245)(opline(j:j),j=1,is)
 2245  format(80a1)
      end do                                                            4d5s10
      if(ngrp.eq.1)then                                                 4d8s10
       write(6,*)('C1 symmetry ')                                       4d8s10
       msym=1                                                           4d8s10
       do i=1,8                                                         4d27s10
        isymsub(i)=1                                                    4d27s10
       end do                                                           4d27s10
      else if(ngrp.eq.8)then                                            4d8s10
       write(6,*)('D2h symmetry ')                                      4d8s10
       do i=1,8                                                         4d27s10
        isymsub(i)=i                                                    4d27s10
       end do                                                           4d27s10
       msym=2                                                           4d8s10
      else if(ngrp.eq.2)then                                            4d8s10
       nneg=0                                                           4d8s10
       do i=1,3                                                         4d8s10
        if(isym(i,2).lt.0)nneg=nneg+1                                   4d8s10
       end do                                                           4d8s10
       if(nneg.eq.1)then                                                4d8s10
        write(6,*)('Cs symmetry ')                                      4d8s10
        do i=1,3                                                        4d27s10
         if(isym(i,2).lt.0)ngen=i                                       4d27s10
        end do                                                          4d27s10
        do i=1,8                                                        4d27s10
         if(isymd2h(ngen,i).gt.0)then                                   4d27s10
          isymsub(i)=1                                                  4d27s10
         else                                                           4d27s10
          isymsub(i)=2                                                  4d27s10
         end if                                                         4d30s10
        end do                                                          4d27s10
        msym=3                                                          4d8s10
       else if(nneg.eq.2)then                                           4d8s10
        write(6,*)('C2 symmetry ')                                      4d8s10
        ipd=1                                                           4d27s10
        do i=1,3                                                        4d27s10
         if(isym(i,2).lt.0)then                                         4d27s10
          itmp(ipd)=i                                                   4d27s10
          ipd=ipd+1                                                     4d27s10
         end if                                                         4d27s10
        end do                                                          4d27s10
        do i=1,8                                                        4d27s10
         if(isymd2h(itmp(1),i)*isymd2h(itmp(2),i).gt.0)then             4d27s10
          isymsub(i)=1                                                  4d27s10
         else                                                           4d27s10
          isymsub(i)=2                                                  4d27s10
         end if                                                         4d27s10
        end do                                                          4d27s10
        msym=4                                                          4d8s10
       else                                                             4d8s10
        write(6,*)('Ci symmetry ')                                      4d8s10
        msym=5                                                          4d8s10
        do i=1,8                                                        4d27s10
         if(isymd2h(1,i)*isymd2h(2,i)*isymd2h(3,i).gt.0)then            4d27s10
          isymsub(i)=1                                                  4d27s10
         else                                                           4d27s10
          isymsub(i)=2                                                  4d27s10
         end if                                                         4d27s10
        end do                                                          4d27s10
       end if                                                           4d8s10
      else if(ngrp.eq.4)then                                            4d8s10
       nnin=3                                                           4d8s10
       nnax=0                                                           4d8s10
       do i=2,4                                                         4d8s10
        nneg=0                                                          4d8s10
        do j=1,3                                                        4d8s10
         if(isym(j,i).lt.0)nneg=nneg+1                                  4d8s10
        end do                                                          4d8s10
        nnin=min(nnin,nneg)                                             4d8s10
        nnax=max(nnax,nneg)                                             4d8s10
       end do                                                           4d8s10
       if(nnin.eq.1.and.nnax.eq.2)then                                  4d8s10
        write(6,*)('C2v symmetry ')                                     4d8s10
        do i=1,3                                                        4d30s10
         isum=0                                                         4d30s10
         do j=2,4                                                       4d30s10
          isum=isum+isym(i,j)                                           4d30s10
         end do                                                         4d30s10
         if(isum.eq.3)then                                              4d30s10
          ismul(i)=0                                                    4d30s10
         else                                                           4d30s10
          ismul(i)=1                                                    4d30s10
         end if                                                         4d30s10
        end do                                                          4d30s10
        do i=1,8                                                        4d30s10
         do j=1,4                                                       4d30s10
          idff=iabs(isym(1,j)-isymd2h(1,i))*ismul(1)                    4d30s10
     $        +iabs(isym(2,j)-isymd2h(2,i))*ismul(2)                    4d30s10
     $        +iabs(isym(3,j)-isymd2h(3,i))*ismul(3)                    4d30s10
          if(idff.eq.0)then                                             4d30s10
           isymsub(i)=j                                                 4d30s10
          end if                                                        4d30s10
         end do                                                         4d30s10
        end do                                                          4d30s10
        msym=6                                                          4d8s10
       else if(nnax.eq.3)then                                           4d8s10
        write(6,*)('C2h symmetry ')                                     4d8s10
        write(6,*)('haven''t coded up isymsub rules yet... ')           4d30s10
        do i=1,4
         write(6,*)(isym(j,i),j=1,3)
        end do
        write(6,*)('we have an inversion and a reflection')
        write(6,*)('reflection ')
        do i=2,4
         isum=isym(1,i)+isym(2,i)+isym(3,i)
         if(isum.eq.1)ngen=i
        end do
        write(6,*)('generator: '),ngen
        do i=1,8                                                        4d30s18
         if(isymd2h(ngen,i).gt.0)then
          nref=1                                                        4d30s18
         else
          nref=2                                                        4d30s18
         end if                                                         4d30s18
         if(isymd2h(1,i)*isymd2h(2,i)*isymd2h(3,i).gt.0)then            4d30s18
          nver=1                                                        4d30s18
         else                                                           4d30s18
          nver=2                                                        4d30s18
         end if                                                         4d30s18
         isymsub(i)=nver+2*(nref-1)                                     4d30s18
        end do                                                          4d30s18
        msym=7                                                          4d8s10
       else if(nnin.eq.2.and.nnax.eq.2)then                             4d8s10
        msym=8                                                          4d8s10
        write(6,*)('D2 symmetry ')                                      4d8s10
        write(6,*)('I am sorry, but I never coded up D2 symmetry')      4d8s10
        stop                                                            4d8s10
       end if                                                           4d8s10
      end if                                                            4d8s10
      if(idorel4c.ne.0)lmax=lmax+1                                      2d20s20
      call spher(lmax,bc,ibc)                                           11d9s22
      ngaus=0                                                           4d8s10
      do i=1,natom                                                      4d8s10
       if(iapair(1,i).gt.0)then                                         10d27s17
        iapair(1,iapair(1,i))=-i                                        4d8s10
       end if                                                           4d8s10
       if(iapair(1,i).ge.0)ngaus=ngaus+ngausa(i)                        4d8s10
 1065  format(2i5,5x,3i5)
      end do                                                            4d8s10
      write(6,*)
      write(6,*)('total number of gaussians: '),ngaus
      if(nwcont.gt.0)then                                               2d3s19
       write(6,*)('number of words for contraction coefficients = '),    9d12s17
     $     nwcont                                                       2d22s19
      end if                                                            2d3s19
      ngaus2=ngaus*2
      ngaus3=ngaus*3                                                    1d22s10
      ngaus4=ngaus3+ngaus                                               1d22s10
      ngaus5=ngaus4+ngaus                                               1d22s10
      ngaus6=ngaus5+ngaus                                               1d22s10
      ngaus7=ngaus6+ngaus                                               4d5s10
      ngaus8=ngaus7+ngaus                                               9d12s17
      ibdat=ibcoff
      ibcoff=ibdat+9*ngaus+nwcont                                       1d28s19
      ibctop=ibcoff                                                     9d12s17
      icdat=9*ngaus                                                     9d12s17
c
c     put contraction coefficients after this, and add 10 pointer
c     to contraction coefficients relative to ibdat
c
      jbdat=ibdat
      ioff=0                                                            1d22s10
      czetal=1d10
      explast=0d0                                                       3d17s21
      do ip=1,ipack
       isto=ippack(ip)
       l=nint(bc(isto))
       isto=isto+1
       nzeta=nint(bc(isto))
       isto=isto+1
       ncont=max(1,nint(bc(isto)))                                      9d12s17
       isto=isto+1                                                      9d12s17
       isto0=isto                                                       4d12s19
       do i=1,3                                                         1d22s10
        cartb(i)=bc(isto)                                               1d22s10
        isto=isto+1                                                     1d22s10
       end do                                                           1d22s10
       do j=1,natom                                                     2d18s10
        rms=sqrt((xcart(1,j)-cartb(1))**2                               4d8s10
     $          +(xcart(2,j)-cartb(2))**2                               4d8s10
     $          +(xcart(3,j)-cartb(3))**2)                              4d8s10
        if(rms.lt.1d-14)then                                            2d18s10
         nhere=j                                                        2d18s10
         go to 1067                                                     2d18s10
        end if                                                          2d18s10
       end do                                                           2d18s10
       write(6,*)('cartesian error!!! '),isto0,cartb,ip,ippack(ip)                                2d
       stop                                                             2d18s10
 1067  continue                                                         2d18s10
       if(iapair(1,nhere).ge.0)then                                     4d8s10
        icstor=max(1,ncont)                                             9d12s17
        istorc=isto+nzeta*2                                             9d12s17
        if(explast.eq.bc(isto))then                                     3d17s21
         icdatu=icdatsave                                               9d26s17
        else                                                            9d12s17
         czetal=bc(istorc)                                              9d12s17
         icdatu=icdat                                                   9d12s17
         if(idorel.ne.0.and.abs(czetal-1d0).lt.1d-8)czetal=1d10         3d8s21
        end if                                                          9d12s17
        do i=1,nzeta
         ibc(jbdat)=l
         bc(jbdat+ngaus)=bc(isto)
         if(i.eq.1)explast=bc(isto)                                     3d17s21
         bc(jbdat+ngaus2)=bc(isto+nzeta)
         ibc(jbdat+ngaus3)=ioff                                          1d22s10
         bc(jbdat+ngaus4)=cartb(1)                                       1d22s10
         bc(jbdat+ngaus5)=cartb(2)                                       1d22s10
         bc(jbdat+ngaus6)=cartb(3)                                       1d22s10
         ibc(jbdat+ngaus7)=nhere                                        4d8s10
         if(i.eq.1.and.(ncont.gt.1.or.(inoc(l+1).eq.0.and.idorel.ne.0)))10d27s20
     $        then                                                      10d27s20
          ibc(jbdat+ngaus8)=icdatu                                      9d12s17
         else if(i.gt.ncont.and.ncont.gt.1)then                         1d28s19
          ibc(jbdat+ngaus8)=0                                           9d12s17
         else                                                           9d12s17
          ibc(jbdat+ngaus8)=-ncont                                      9d12s17
         end if                                                         9d12s17
         ioff=ioff+2*l+1                                                 1d22s10
         if(iapair(1,nhere).gt.0)then                                   5d3s10
          ioff=ioff+2*l+1                                               5d3s10
         end if                                                         5d3s10
         isto=isto+1
         jbdat=jbdat+1
        end do
        if((ncont.gt.1.or.(inoc(l+1).eq.0.and.idorel.ne.0)).            10d27s20
     $       and.icdatu.ne.icdatsave)then                               10d27s20
         if(idorel.eq.0)then                                             9d12s17
          nzetau=nzeta                                                   9d12s17
         else                                                            9d12s17
          nzetau=nzeta*2                                                 9d12s17
         end if                                                          9d12s17
         nmove=ncont*nzetau                                              9d12s17
         istosave=istorc
         icdatsave=icdat                                                9d26s17
         do imove=0,nmove-1                                              9d12s17
          if(ibdat+icdat.gt.ibctop)then
           write(6,*)('when storing contraction coefficients, ')        10d28s20
           write(6,*)('stopping, since ibctop is only '),ibctop
           write(6,*)('while ibdat+icdat = '),ibdat+icdat,icdat
           write(6,*)('this is for zeta = '),i
           stop 'ibctop'
          end if
          bc(ibdat+icdat)=bc(istorc+imove)                              9d12s17
          icdat=icdat+1                                                  9d12s17
         end do                                                          9d12s17
        end if                                                          9d12s17
       else                                                             4d8s10
        isto=isto+nzeta                                                 4d8s10
       end if                                                           4d8s10
      end do
      write(6,*)('gaussians: ')
      jbdat=ibdat-1
      nbasis=0                                                          2d18s10
      do i=1,ngrp                                                       4d30s10
       nbasb(i)=0                                                       4d30s10
       nbasc(i)=0                                                       4d25s18
       irunc2(i)=0                                                      4s23s18
       irunc(i)=0                                                       5d3s10
      end do                                                            4d30s10
c
c     contraction. gaussians i, i+1, ..., i+p will be contracted down to
c     new set of functions j,j+1,...,j+c, with c le p. for setting info
c     in ibdat, we need the gaussians, while for setting ibstor etc,
c     we need contracted functions. we will achieve this by adding an
c     extra index to the ngaus data: this index will be p*100+c, meaning
c     this corresponds to the first gaussian used in the contraction,
c     1, meaning this function is not contracted, or -k, 1 le l < p,
c     for fcn i-k in the contraction.
c     the first nbasisp entries of bstor and sstor are for the          4s23s18
c     uncontracted functions. the 2nd are for the contracted functions. 4s23s18
c
      do i=1,ngaus                                                      5d3s10
       ii=jbdat+i                                                       4d12s19
       lp=ibc(jbdat+i)+1                                                5d3s10
       nl=2*lp-1                                                        5d3s10
       nhere=ibc(jbdat+i+ngaus7)                                        5d3s10
       if(iapair(1,nhere).gt.0)nl=nl*2                                  5d3s10
       nbasis=nbasis+nl                                                 5d3s10
      end do                                                            5d3s10
c                                                                       2d14s19
c     idea was primitive pointers followed by contracted pointers       2d14s19
c     but if relativistic, the number of contracted functions can be    2d14s19
c     greater than the number of primitives - up to twice as many.      2d14s19
c     this is because of large and small component ...                  2d14s19
c                                                                       2d14s19
      ibstor=ibcoff                                                     5d3s10
      isstor=ibstor+nbasis*3                                            2d14s19
      ibcoff=isstor+nbasis*3                                            2d14s19
      call enough('basisz.  9',bc,ibc)
      jbstor=ibstor                                                     5d3s10
      jsstor=isstor                                                     5d3s10
       jbstor2=ibstor+nbasis                                            12d13s22
       jsstor2=isstor+nbasis                                            12d13s22
      do isb=1,ngrp                                                     3d3s20
       nbaspre(isb)=0                                                   3d3s20
      end do                                                            3d3s20
      write(6,*)('   #  l      exp.     atom(s) no. in symmetry blocks')5d3s10
      nleft0=0                                                          12d13s22
      ngothere=0                                                        1d18s23
      do i=1,ngaus
       nhere=ibc(jbdat+i+ngaus7)                                        4d8s10
       nccode=ibc(jbdat+i+ngaus8)                                       4s23s18
       if(nccode.gt.0)then                                              12d13s22
        do kk=1,nleft0                                                  12d13s22
         do id2hsym=1,8                                                   4d30s10
          if(ipt(id2hsym,lp).gt.0)then                                    4d30s10
           isto=ipt(id2hsym,lp)                                           4d30s10
           mhere=ibc(isto+1)                                              4d30s10
           isb=isymsub(id2hsym)                                           5d3s10
           do k=1,mhere                                                   12d13s22
            irunc2(isb)=irunc2(isb)+1                                     4s23s18
            ibc(jbstor2)=irunc2(isb)                                      4s23s18
            ibc(jsstor2)=isb                                              4s23s18
            jbstor2=jbstor2+1                                             4s23s18
            jsstor2=jsstor2+1                                             4s23s18
           end do                                                         12d13s22
          end if                                                        12d13s22
         end do                                                         12d13s22
         if(iapair(1,nhere).gt.0)then                                     12d15s22
          do id2hsym=1,8                                                  5d6s10
           if(ipt(id2hsym,lp).gt.0)then                                   5d6s10
            isto=ipt(id2hsym,lp)                                           4d30s10
            mhere=ibc(isto+1)                                              4d30s10
            isb=isymsub(id2hsym)                                           5d3s10
            isbo=multh(isb,iapair(2,nhere))                               5d3s10
            do k=1,mhere                                                  5d3s10
             irunc2(isbo)=irunc2(isbo)+1                                  5d2s18
             ibc(jbstor2)=irunc2(isbo)                                    5d2s18
             ibc(jsstor2)=isbo                                            4s23s18
             jbstor2=jbstor2+1                                            4s23s18
             jsstor2=jsstor2+1                                            4s23s18
            end do                                                        12d15s22
           end if                                                         12d15s22
          end do                                                        12d15s22
         end if                                                         12d15s22
        end do                                                          12d13s22
        nleft0=0                                                        12d13s22
        ngothere=0                                                      12d13s22
       end if                                                           12d13s22
       lp=ibc(jbdat+i)+1
       nl=2*lp-1                                                        4d8s10
       if(iapair(1,nhere).gt.0)nl=nl*2                                  4d8s10
       symatoms='     '                                                 4d8s10
       do j=1,ngrp                                                      4d30s10
        insym(j)=0                                                      4d30s10
       end do                                                           4d30s10
       ioff=0                                                           5d6s10
       do id2hsym=1,8                                                   4d30s10
        if(ipt(id2hsym,lp).gt.0)then                                    4d30s10
         mgothere=ngothere                                              12d13s22
         isto=ipt(id2hsym,lp)                                           4d30s10
         mhere=ibc(isto+1)                                              4d30s10
         isb=isymsub(id2hsym)                                           5d3s10
33561    format('for d2h sym ',i1,' becoming sym ',i1,' no. fcns ',i2)  4d30s10
         insym(isb)=insym(isb)+mhere                                    5d3s10
c
c     for primitives
c
         do k=1,mhere                                                   12d13s22
          irunc(isb)=irunc(isb)+1                                       12d13s22
          ibc(jbstor)=irunc(isb)                                        12d13s22
          jbstor=jbstor+1                                               5d3s10
          ibc(jsstor)=isb                                               12d13s22
          jsstor=jsstor+1                                               5d3s10
         end do                                                         12d13s22
         if(nccode.gt.0)then                                            12d13s22
c
c     first contracted fcn
c
          mgothere=mgothere+1                                            12d13s22
          do k=1,mhere                                                   12d13s22
           irunc2(isb)=irunc2(isb)+1                                     4s23s18
           ibc(jbstor2)=irunc2(isb)                                      4s23s18
           ibc(jsstor2)=isb                                              4s23s18
           jbstor2=jbstor2+1                                             4s23s18
           jsstor2=jsstor2+1                                             4s23s18
          end do                                                         12d13s22
         end if                                                          12d13s22
         if(nccode.lt.0)then                                             12d13s22
          if(nccode.eq.-1)then                                          12d13s22
           do k=1,mhere                                                 12d13s22
            irunc2(isb)=irunc2(isb)+1                                     4s23s18
            ibc(jbstor2)=irunc2(isb)                                      4s23s18
            ibc(jsstor2)=isb                                              4s23s18
            jbstor2=jbstor2+1                                             4s23s18
            jsstor2=jsstor2+1                                             4s23s18
           end do                                                         12d13s22
           if(idorel.ne.0)nbaspre(isb)=nbaspre(isb)+mhere               12d13s22
          else                                                          12d13s22
           nleft=-nccode-mgothere                                       12d13s22
           nleft0=nleft                                                 12d13s22
           if(nleft0.gt.1)nleft=1
           do kk=1,nleft                                                12d13s22
            do k=1,mhere                                                12d13s22
             irunc2(isb)=irunc2(isb)+1                                     4s23s18
             ibc(jbstor2)=irunc2(isb)                                      4s23s18
             ibc(jsstor2)=isb                                              4s23s18
             jbstor2=jbstor2+1                                             4s23s18
             jsstor2=jsstor2+1                                             4s23s18
            end do                                                         12d13s22
           end do                                                       12d13s22
           nleft0=nleft0-nleft                                          12d13s22
           mgothere=mgothere+nleft                                      12d13s22
          end if                                                        12d13s22
         end if                                                          12d13s22
        end if                                                          5d6s10
       enddo                                                            5d6s10
       if(iapair(1,nhere).gt.0)then                                     5d6s10
        do id2hsym=1,8                                                  5d6s10
         if(ipt(id2hsym,lp).gt.0)then                                   5d6s10
          mgothere=ngothere                                             12d15s22
          isto=ipt(id2hsym,lp)                                           4d30s10
          mhere=ibc(isto+1)                                              4d30s10
          isb=isymsub(id2hsym)                                           5d3s10
          isbo=multh(isb,iapair(2,nhere))                               5d3s10
          insym(isbo)=insym(isbo)+mhere                                 5d3s10
c
c     for primitives
c
          do k=1,mhere                                                  5d3s10
           irunc(isbo)=irunc(isbo)+1                                    5d3s10
           ibc(jbstor)=irunc(isbo)                                      5d3s10
           jbstor=jbstor+1                                              5d3s10
           ibc(jsstor)=isbo                                             5d3s10
           jsstor=jsstor+1                                              5d3s10
          end do                                                        5d3s10
          if(nccode.gt.0)then                                            12d13s22
c
c     first contracted fcn
c
           mgothere=mgothere+1                                            12d13s22
           do k=1,mhere                                                  5d3s10
            irunc2(isbo)=irunc2(isbo)+1                                 5d2s18
            ibc(jbstor2)=irunc2(isbo)                                   5d2s18
            ibc(jsstor2)=isbo                                           4s23s18
            jbstor2=jbstor2+1                                            4s23s18
            jsstor2=jsstor2+1                                            4s23s18
           end do                                                       12d15s22
          end if                                                        12d15s22
          if(nccode.lt.0)then                                             12d13s22
           if(nccode.eq.-1)then                                          12d13s22
            do k=1,mhere                                                 12d13s22
             irunc2(isbo)=irunc2(isbo)+1                                 5d2s18
             ibc(jbstor2)=irunc2(isbo)                                   5d2s18
             ibc(jsstor2)=isbo                                           4s23s18
             jbstor2=jbstor2+1                                            4s23s18
             jsstor2=jsstor2+1                                            4s23s18
            end do                                                      12d15s22
            if(idorel.ne.0)nbaspre(isbo)=nbaspre(isbo)+mhere            12d15s22
           else                                                         12d15s22
            nleft=-nccode-mgothere                                       12d13s22
            nleft0=nleft                                                 12d13s22
            if(nleft0.gt.1)nleft=1
            do kk=1,nleft                                                12d13s22
             do k=1,mhere                                                  5d3s10
              irunc2(isbo)=irunc2(isbo)+1                                 5d2s18
              ibc(jbstor2)=irunc2(isbo)                                   5d2s18
              ibc(jsstor2)=isbo                                           4s23s18
              jbstor2=jbstor2+1                                            4s23s18
              jsstor2=jsstor2+1                                            4s23s18
             end do                                                     12d15s22
            end do                                                      12d15s22
            nleft0=nleft0-nleft                                          12d13s22
            mgothere=mgothere+nleft                                      12d13s22
           end if                                                       12d15s22
          end if                                                        12d15s22
         end if                                                         5d3s10
        end do                                                          5d6s10
       end if                                                           5d6s10
       ngothere=mgothere                                                12d15s22
       if(iapair(1,nhere).gt.0)then                                     4d8s10
        write(symatoms,1063)iapair(1,nhere)                             4d8s10
 1063   format(' & ',i2)                                                4d8s10
       end if                                                           4d8s10
       write(6,4)i,ibc(jbdat+i),bc(jbdat+i+ngaus),                      1d22s10
     $      nhere,symatoms,(insym(j),j=1,ngrp)                          4d30s10
       ll=ibc(jbdat+i)+1                                                6d11s19
       if(ll.le.4)then                                                  6d11s19
        if(jcoreg(ll,nhere).lt.icoreg(ll,nhere))then                    6d11s19
         do k=1,ngrp                                                    6d11s19
          ncoreg(k)=ncoreg(k)+insym(k)                                  6d11s19
         end do                                                         6d11s19
         jcoreg(ll,nhere)=jcoreg(ll,nhere)+1                            6d11s19
        end if                                                          6d11s19
        if(jvalg(ll,nhere).lt.ivalg(ll,nhere))then                      6d11s19
         do k=1,ngrp                                                    6d11s19
          nvalg(k)=nvalg(k)+insym(k)                                    6d11s19
         end do                                                         6d11s19
         jvalg(ll,nhere)=jvalg(ll,nhere)+1                              6d11s19
        end if                                                          6d11s19
       end if                                                           6d11s19
       do j=1,ngrp                                                      4d30s10
        nbasb(j)=nbasb(j)+insym(j)                                      4d30s10
       end do                                                           4d30s10
    4  format(i5,i3,1p1e15.7,i3,a5,2x,8i3)                              4d8s10
      end do
      do kk=1,nleft0                                                    12d13s22
       do id2hsym=1,8                                                   4d30s10
        if(ipt(id2hsym,lp).gt.0)then                                    4d30s10
         isto=ipt(id2hsym,lp)                                           4d30s10
         mhere=ibc(isto+1)                                              4d30s10
         isb=isymsub(id2hsym)                                           5d3s10
         do k=1,mhere                                                   12d13s22
          irunc2(isb)=irunc2(isb)+1                                     4s23s18
          ibc(jbstor2)=irunc2(isb)                                      4s23s18
          ibc(jsstor2)=isb                                              4s23s18
          jbstor2=jbstor2+1                                             4s23s18
          jsstor2=jsstor2+1                                             4s23s18
         end do                                                         12d13s22
        end if                                                          12d13s22
       end do                                                           12d13s22
       if(iapair(1,nhere).gt.0)then                                     12d15s22
        do id2hsym=1,8                                                  5d6s10
         if(ipt(id2hsym,lp).gt.0)then                                   5d6s10
          isto=ipt(id2hsym,lp)                                           4d30s10
          mhere=ibc(isto+1)                                              4d30s10
          isb=isymsub(id2hsym)                                           5d3s10
          isbo=multh(isb,iapair(2,nhere))                               5d3s10
          do k=1,mhere                                                  5d3s10
           irunc2(isbo)=irunc2(isbo)+1                                  5d2s18
           ibc(jbstor2)=irunc2(isbo)                                    5d2s18
           ibc(jsstor2)=isbo                                            4s23s18
           jbstor2=jbstor2+1                                            4s23s18
           jsstor2=jsstor2+1                                            4s23s18
          end do                                                        12d15s22
         end if                                                         12d15s22
        end do                                                          12d15s22
       end if                                                           12d15s22
      end do                                                            12d13s22
      njbs2=jbstor2-ibstor-nbasis
      ncomp=1                                                           2d14s19
      if(idorel.ne.0)ncomp=2                                            2d14s19
      write(6,1328)(ncoreg(i),i=1,ngrp)                                 6d12s23
 1328 format('default closed orbs:',8i5)                                6d12s23
      write(6,1329)(nvalg(i),i=1,ngrp)                                  6d12s23
 1329 format('default active orbs:',8i5)                                6d12s23
      write(6,*)('total number of primitive basis functions = '),       2d14s19
     $     nbasis*ncomp                                                 2d14s19
      write(6,352)(nbasb(j)*ncomp,j=1,ngrp)                             2d14s19
      nbpt=0                                                            4d25s18
      do j=1,ngrp                                                       4d25s18
       nlorcont(j)=irunc2(j)                                            3d29s21
       irunc2(j)=irunc2(j)+nbaspre(j)                                   3d3s20
       nbpt=nbpt+irunc2(j)                                              3d3s20
       nbasc(j)=irunc2(j)                                               4d25s18
      end do                                                            4d25s18
      write(6,*)('total number of contracted basis functions = '),nbpt  4d25s18
      write(6,352)(irunc2(j),j=1,ngrp)
  352 format(' by symmetry block: ',8i4)
      nsymb=ngrp                                                        4d30s10
      nsymbb=nsymb                                                      3d3s10
      nbasall=nbasis                                                    2d24s10
      if(atom(1).eq.'us  ')then                                         3d2s22
       atom(1)=element(nint(atnum(1,1)))                                3d2s22
      end if                                                            3d2s22
      go to 254                                                         8d8s14
c                                                                       2d19s04
c     sohf input                                                        2d19s04
c                                                                       2d19s04
 1252 continue
      casdat(6)=1d-6                                                    12d19s22
      do i=1,8
       icore(i)=0
      end do                                                            9d23s04
      nusecsf=-1                                                        1d8s19
      molden=0                                                          1d8s19
      myguess=0                                                         6d26s08
      if(numeminus.eq.0)then                                            11d21s19
       open(unit=1,file='orbs',form='unformatted')                      5d25s21
       write(6,*)('getting input from orbs ')                           5d25s21
       call loadr(potdws,ibdat,ibstor,isstor,nsymb,istinfo,ncoreg,      1d2s20
     $      nvalg,nlzzq,iextradatd,nectradatx,0,nbasdws,element,idum,   7d26s22
     $      idum,idum,idum,1,bc,ibc,icanon)                             5d5s23
      close(unit=1)                                                     11d20s19
       ngrp=nsymb                                                       12d27s19
       first(1)=potdws                                                  11d21s19
      end if                                                            11d21s19
      iwerner=0                                                         7d9s13
      idavopt=0                                                         7d9s13
      idohf=1                                                           8d8s14
      idogrado=0                                                        1d11s16
      nocc(1)=0                                                         6d28s04
      icartg=0
      ichrg=0                                                           11d20s19
       nocc(1)=-numeminus/2                                             11d20s19
       if(idorel4c.ne.0)nocc(1)=numeminus                               10d5s22
  250 continue                                                          2d19s04
  252   continue
         read(5,100,end=1008)line                                       3d1s10
         is=1                                                           8d8s14
         call delim(line,is,ie)                                         8d8s14
        if(line(is:is+4).eq.'debug')then                                6d5s20
 6606   continue                                                        6d5s20
         is=ie+1                                                        6d5s20
         call delim(line,is,ie)                                         6d5s20
         if(ie.ge.is)then                                               6d5s20
          write(6,*)('debug print out turned on for '),line(is:ie)      6d5s20
          nchere=ie+1-is                                                6d5s20
          do i=1,nprtr                                                  6d5s20
           if(line(is:ie).eq.prtrname(i)(1:nchere))then                 6d5s20
            iprtr(i)=1                                                  6d5s20
            go to 6606                                                  6d5s20
           end if                                                       6d5s20
          end do                                                        6d5s20
          write(6,*)('no matching string in prtrname found ')           6d5s20
         end if                                                         6d5s20
         go to 252                                                      6d5s20
        end if                                                          6d5s20
         if(line(is:is+4).eq.'nocan')then                               3d28s16
          inocan=1                                                      3d28s16
          write(6,*)('orbitals will not be cannonicalized ')            3d28s16
          go to 252                                                     3d28s16
         end if                                                         3d28s16
         if(line(is:is+3).eq.'grad')then                                8d4s22
          write(6,*)('in grado block ')
c
c     options:
c     nothing: compute totally symmetric gradients of energy
c     bodc: compute Born-Oppenheimer diagonal correction and derivative
c           couplings if more than one state.
c     nonad: if HF wavefunction, compute BODC as well as second order
c            Born-Oppenheimer correction terms.
c     minimize: optimize geometry to minimize energy.
c           no further argument means minimize the state averaged energy
c           otherwise specify root number and state name.
c     stationary: optimize geometry to minimize gradients
c
          idogrado1(1)=1                                                6d18s22
          is=ie+1                                                       6d18s22
          call delim(line,is,ie)                                        6d18s22
          write(6,*)('is,ie '),is,ie
          if(ie.ge.is)then                                              6d18s22
           if(line(is:is).eq.'b')then                                   6d18s22
            idogrado1(2)=1                                              6d18s22
           else if(line(is:is).eq.'n')then                              6d18s22
            idogrado1(1)=3                                              6d18s22
           else if(line(is:is).eq.'m')then                              6d18s22
            idogrado1(3)=1                                              6d18s22
           else if(line(is:is).eq.'s')then                              6d18s22
            idogrado1(3)=2                                              6d18s22
           else                                                         6d18s22
            write(6,*)('I don''t understand grado argument '),
     $           line(is:ie),('!')                                      6d18s22
            irtrn=1                                                     6d18s22
            return                                                      6d18s22
           end if                                                       6d18s22
          end if                                                        6d18s22
        write(6,*)('what we have for idogrado '),idogrado,idogrado1
          go to 252                                                     1d11s16
         end if                                                         1d11s16
         if(line(is:is+3).eq.'gues')then                                11d20s19
c
c     the mygues card indicates we will read in a guess for orbitals.
c     parahf saves orbitals on a file which contains information about
c     the basis, the orbitals in the ao basis, the overlap, and the
c     orbitals in the orthogonal basis. One can use the orbitals in the
c     ao basis or orbitals in the orthogonal basis as initial guess.
c     one can also use non-relativistic orbitals as initial guess for
c     relativitistic orbitals, and use orbitals from another point
c     group or basis as an initial guess.
c     the default or myguess ao means use the orbitals in the ao basis.
c     using myguess ob will mean using the orbitals in the orthogonal
c     basis.
c     (the difference is in the ao case, the orbitals will be
c     orthogonalized using the current overlap first, unless noortho
c     is included on the mygues line. Warning: using ao and noortho
c     will yield bad energies !!!!!).
c
          is=ie+1                                                       2d18s16
          call delim(line,is,ie)                                        2d18s16
          write(6,*)('after mygues card: '),is,ie,line(is:ie)
          if(is.gt.ie)then
           myguess=2                                                    12d19s22
          else                                                          2d18s16
           if(line(is:is+1).eq.'h0')then                                12d19s22
           else if(line(is:is+3).eq.'read')then
            myguess=2                                                   12d19s22
            is=ie+1                                                     12d19s22
            call delim(line,is,ie)                                      12d19s22
            if(is.le.ie)then                                            12d19s22
             if(line(is:ie).eq.'ob')then                                  2d18s16
              myguess=2                                                   2d18s16
             else if(line(is:ie).eq.'ao')then                             2d18s16
              myguess=1                                                   2d18s16
              is=ie+1                                                     3d14s16
              call delim(line,is,ie)                                      3d14s16
              if(ie.gt.is)then                                            3d14s16
               if(line(is:is+3).eq.'noor')then                            3d14s16
                write(6,*)('ao orbitals not orthogonalized ')             3d14s16
                write(6,*)('*** W A R N I N G *** ')                      3d14s16
                write(6,*)('this will give B A D energies')               3d14s16
                myguess=-1                                                3d14s16
               end if                                                     3d14s16
              end if                                                      3d14s16
             else                                                         2d18s16
              write(6,*)('unknown mygues flag: '),line(is:ie)              2d18s16
              stop                                                        2d18s16
             end if                                                       2d18s16
            end if                                                      12d19s22
           end if                                                       12d19s22
          end if                                                        2d18s16
          write(6,*)('myguess = '),myguess
          go to 252                                                       6d26s08
         end if                                                           6d26s08
         if(line(is:is+3).eq.'wern')then                                8d8s14
          ix=ie+1                                                       4d2s18
          call delim(line,ix,ie)                                        11d29s17
          if(ie.gt.0)then                                               11d29s17
           read(line(ix:ie),*)iwerner                                   11d29s17
          else                                                          11d29s17
           iwerner=2                                                       7d9s13
          end if                                                        11d29s17
         end if                                                           7d9s13
         if(line(is:is+3).eq.'davi')then                                      7d9s13
          idavopt=1                                                       7d9s13
          call intin(line,5,len(line),idavopt)                            7d9s13
         end if                                                           7d9s13
         if(line(is:is+3).eq.'char')then                                11d20s19
          is=ie+1                                                       11d20s19
          call delim(line,is,ie)                                        11d20s19
          read(line(is:ie),*)ichrg                                      11d20s19
          write(6,*)('charge read in is '),ichrg                        11d20s19
          nhere=numeminus-ichrg                                         11d20s19
          if(mod(nhere,2).ne.0)then                                     11d20s19
           write(6,*)('number of electrons '),nhere,('is odd')          11d20s19
           write(6,*)('use cas rather than sohf in this case')          11d20s19
           irtrn=1                                                      11d20s19
           return                                                       11d20s19
          end if                                                        11d20s19
          nocc(1)=-nhere/2                                              11d20s19
          if(idorel4c.ne.0)nocc(1)=nhere                                10d5s22
         end if                                                         11d20s19
         if(line(is:is+3).eq.'nocc'.or.line(is:is+2).eq.'occ')then      8d12s22
c                                                                       6d28s04
c     to have program figure out distribution of occupied orbitals,     6d28s04
c     read in nocc(1) as the negitive of total number of occupied orbs. 6d28s04
c                                                                       6d28s04
         ncomma=0                                                       8d8s14
         do i=1,8                                                       8d8s14
          is=ie+1                                                        8d8s14
          call delim(line,is,ie)                                        8d8s14
          if(ie.ge.is)then                                              8d8s14
           ncomma=ncomma+1                                              8d8s14
           read(line(is:ie),*)nocc(ncomma)                              8d8s14
          else                                                          8d8s14
           go to 255                                                    8d8s14
          end if                                                        8d8s14
         end do                                                         8d8s14
  255    continue                                                       8d8s14
          if(ncomma.gt.nsymb)then                                       8d12s22
           write(6,*)('you specified '),ncomma,                         8d12s22
     $         ('symmetries in occ input')                              8d12s22
           write(6,*)('but point group only has '),nsymb,('symmetries') 8d12s22
           write(6,*)('perhaps you input the geometry incorrectly?')    8d12s22
           irtrn=1                                                      8d12s22
           return                                                       8d12s22
          end if                                                        8d12s22
          write(6,*)('nocc read in as '),(nocc(idws),idws=1,ncomma)
          go to 252
         end if
         if(line(is:is+4).eq.'itmax'.or.line(is:is+4).eq.'maxit')then   6d6s23
          is=ie+1                                                       8d8s14
          call delim(line,is,ie)                                        8d8s14
          read(line(is:ie),*)itmax                                                3d17s04
          go to 252                                                     3d17s04
         end if                                                         3d17s04
         if(line(is:is).eq.'&')go to 254                                8d8s14
        go to 252                                                        2d19s04
c
c     input for relativistic calculations.
c     to not include relativity either do not include &REL input
c     or set the speed of light to a very high value.
c
c     when &REL is present in the input, 2 component spin-free calculations
c     are done. The details can be modified by the following parameters:
c     alpha, the fine structure constant
c     anuc(i) the nuclear mass for atom i for computing finite nuclear radius
c       if more than 1837, this is assumed to be in au, otherwise it is
c       take to be in amu. Defaults to most probable nuclear isotope.
c       (remember, blank delimited)
c     itype, the kind of small component 2-e integrals to keep.
c       can be "", "ss", "ssb", "ssg", "ssss", "ssssb" or "ssssg".
c     If there is b in                                                  2d14s20
c     in the string, the Breit interaction integrals are included in the
c     property calculation.                                             10d24s20
c     If there is a g in string, then the Gaunt approximation to the    2d14s20
c     Breit interaction will be included in the property calculation.   10d24s20
c     The Breit interaction will not be included in the cas or mrci part10d24s20
c     of the calculation, only in the property calculation.             10d24s20
c     if the keyword fourc is present, a full four component calculation2d20s20
c     will be carried out. In this case symmetry is turned off,         2d20s20
c     contraction can not be used, and the spin can not be specified.
c     nonetheless, the cas code can be used MC calculations and multiple2d20s20
c     roots.                                                            2d20s20
c
 1253  continue
       alpha=1d0/137.0359895d0
       ascale=0.25d0*alpha*alpha
       idorel4c=0                                                       2d20s20
c
c     idorel=1: no s 2e integrals
c     idorel=2: SS 2e integrals
c     idorel=3: SSSS 2e integrals
c     idorel minus those numbers: include breit
c
       idorel=2
       write(6,*)
       write(6,*)('This is going to be a relativistic calculation ')
       write(6,*)
       itype='ss   '
       write(6,*)('2 electron integral type = '),itype
 1254  continue
        read(5,100,end=1008)line                                        8d17s15
        is=1                                                            8d17s15
        call delim(line,is,ie)                                          8d17s15
        if(line(is:is+4).eq.'itype')then                                8d26s15
         write(6,*)('we have itype in input '),ie
         is=ie+1                                                        8d26s15
         call delim(line,is,ie)                                         8d26s15
         if(ie.lt.is)then                                               8d26s15
          itype='     '                                                 8d26s15
          idorel=1                                                      8d26s15
         else if(ie.eq.is+1)then                                        8d26s15
          itype='ss   '                                                 8d26s15
          idorel=2                                                      8d26s15
         else if(ie.eq.is+2)then                                        8d26s15
          if(line(ie:ie).eq.'b')then                                    2d18s20
           itype='ssb  '                                                 8d26s15
           idorel=-2                                                     8d26s15
          else if(line(ie:ie).eq.'g')then                               2d18s20
           itype='ssg  '                                                2d18s20
           idorel=-4                                                    2d18s20
          else                                                          2d18s20
           write(6,*)('unknown tail end on integral type: '),line(ie:ie)2d18s20
           irtrn=1                                                      2d18s20
           return                                                       2d18s20
          end if                                                        2d18s20
         else if(ie.eq.is+3)then                                        8d26s15
          itype='ssss '                                                 8d26s15
          idorel=3                                                      8d26s15
         else if(ie.eq.is+4)then                                        8d26s15
          if(line(ie:ie).eq.'b')then                                    2d18s20
           itype='ssssb'                                                 8d26s15
           idorel=-3                                                     8d26s15
          else if(line(ie:ie).eq.'g')then                               2d18s20
           itype='ssssg'                                                2d18s20
           idorel=-5                                                    2d18s20
          else                                                          2d18s20
           write(6,*)('unknown tail end on integral type: '),           2d18s20
     $          line(ie:ie)                                             2d18s20
           irtrn=1                                                      2d18s20
           return                                                       2d18s20
          end if                                                        2d18s20
         else                                                           8d26s15
          write(6,*)('unknow itype code: '),line(is:ie)                 8d26s15
          go to 1008                                                    8d26s15
         end if                                                         8d26s15
         write(6,*)('2 electron integral type modified to '),itype      8d26s15
         go to 1254                                                     8d26s15
        end if                                                          8d26s15
        if(line(is:is+4).eq.'fourc')then                                2d20s20
         write(6,*)('a full four component calculation will be carried')2d20s20
     $       ,(' out, thus')                                            2d20s20
         write(6,*)('symmetry is turned off, ')                         2d20s20
         do i=1,8                                                       12d7s23
          isymu(i)=0                                                    12d7s23
         end do                                                         12d7s23
         write(6,*)('basis set will be uncontracted, ')                 2d20s20
         idorel4c=-1                                                    2d26s20
         write(6,*)                                                     2d20s20
     $        ('cartesian rather than spherical harmonics are used')    2d20s20
         iusecartr=1                                                    2d21s20
         is=ie+1                                                        2d26s20
         call delim(line,is,ie)                                         2d26s20
         if(ie.ge.is)then                                               2d26s20
          write(6,*)('we have something else on fourc line: '),         2d26s20
     $        line(is:ie)                                               2d26s20
          read(line(is:ie),*)idorel4c                                   2d26s20
         end if                                                         2d26s20
         go to 1254                                                     2d20s20
        end if                                                          2d20s20
        if(line(is:is+4).eq.'alpha')then                                8d17s15
         is=ie+1                                                        8d17s15
         write(6,*)('the fine structure constant is reset from its'),   8d17s15
     $        (' default value of '),alpha,(' ie c= '),1d0/alpha        8d17s15
         call delim(line,is,ie)                                         8d17s15
         read(line(is:ie),*)alpha                                       8d17s15
         write(6,*)('to '),alpha                                        8d17s15
         ascale=0.25d0*alpha*alpha                                      8d17s15
        end if                                                          8d17s15
        if(line(is:is+4).eq.'itype')then                                8d17s15
         is=ie+1                                                        8d17s15
         call delim(line,is,ie)                                         8d17s15
         write(6,*)('itype is reset from its default value of '),itype  8d17s15
         itype='     '                                                  8d17s15
         iup=min(4,ie-is)                                               8d17s15
         itype=line(is:is+iup)                                          8d17s15
         write(6,*)('to '),itype                                        8d17s15
        end if                                                          8d17s15
        if(line(is:is+3).eq.'anuc')then                                 8d17s15
         ianuc=ibcoff                                                   8d17s15
         nanuc=0                                                        8d17s15
         is=ie+1                                                        8d17s15
         call delim(line,is,ie)                                         8d17s15
 1251    continue                                                       8d17s15
         ivald=nbr(line,is,80,0)                                          8d17s15
         if(ivald.eq.-1)then                                            8d17s15
          read(5,100,end=1008)line                                      8d17s15
          is=1                                                          8d17s15
          call delim(line,is,ie)                                        8d17s15
          go to 1251                                                    8d17s15
         else if(ivald.eq.-2)then                                       8d17s15
          go to 1254                                                    8d17s15
         end if                                                         8d17s15
         read(line(is:ivald),*)bc(ianuc+nanuc)                          8d17s15
         nanuc=nanuc+1                                                  8d17s15
         write(6,*)('nuclear mass for nucleus '),nanuc,(' set to '),    8d17s15
     $        bc(ianuc+nanuc)                                           8d17s15
         ibcoff=ianuc+nanuc                                             8d17s15
         is=ivald+1                                                     8d17s15
         go to 1251                                                     8d17s15
        end if                                                          8d17s15
        if(line(is:is).eq.'&')then
         go to 254                                                      8d26s15
        end if                                                          8d26s15
       go to 1254                                                       8d17s15
c                                                                       8d8s14
c     mp2 input                                                         8d8s14
c                                                                       8d8s14
  253  continue                                                         8d8s14
       write(6,*)('do we have nsymb at this point? '),nsymb
       do i=1,nsymb                                                     8d9s24
        icore(i)=ncoreg(i)                                              8d9s24
       end do                                                           8d9s24
       idomp2=1                                                         8d8s14
       read(5,100,end=1008)line                                         8d8s14
       is=1                                                             8d8s14
       call delim(line,is,ie)                                           8d8s14
       if(line(is:is+3).eq.'core')then
        nc=0
        do i=1,8
         is=ie+1                                                        8d8s14
         call delim(line,is,ie)                                         8d8s14
         if(ie.ge.is)then                                               8d8s14
          nc=nc+1                                                       8d8s14
          read(line(is:ie),*)icore(nc)                                  8d8s14
         else                                                           8d8s14
          go to 256                                                     8d8s14
         end if                                                         8d8s14
        end do                                                          8d8s14
  256   continue                                                        8d8s14
        writE(6,*)('icore ')
        writE(6,*)(icore(i),i=1,nc)                                     8d8s14
        go to 5577                                                      8d8s14
       end if
c
c     cas input
c
c     grad
c     debug
c     csf
c     3rd
c     npdiag
c     molden
c     first
c     sum
c     maxps
c     ryd
c     dynw
c     nocan
c     wern
c     nowern
c     maxitmac
c     guess
c     maxitci
c     pthrs
c     econv
c     vconv
c     gconv
c     occ
c     act
c     doub
c     states
c     sneglect
c     nlambda                                                           7d3s23
c     drangcut
c     epssym
c
  251 continue                                                          11d4s05
       write(6,*)('numeminus = '),numeminus
      if(numeminus.eq.0)then                                            11d21s19
       open(unit=1,file='orbs',form='unformatted')                      2d24s21
       write(6,*)('getting input from orbs ')                           2d24s21
       call loadr(potdws,ibdat,ibstor,isstor,nsymb,istinfo,ncoreg,      1d2s20
     $      nvalg,nlzzq,iextradatd,nectradatx,1,nbasdws,element,        6d3s21
     $      nlorcont,idum,idum,idum,1,bc,ibc,icanon)                    5d5s23
       ngrp=nsymb                                                       2d24s21
       do isb=1,nsymb                                                   2d24s21
        idoub(isb)=ncoreg(isb)                                          2d24s21
        iacto(isb)=nvalg(isb)                                           2d24s21
       end do                                                           2d24s21
       close(unit=1)                                                     11d20s19
       first(1)=potdws                                                  11d21s19
      end if                                                            11d21s19
      do i=1,8
       nocc2(i)=0
      end do
      nsumcas=0                                                         5d3s21
      nusecsf=-1                                                        12d24s19
      molden=0                                                          11d21s19
      iwerner=2                                                         12d28s19
      ifirsto=0                                                         8d26s19
      idogrado=0                                                        3d17s21
      ifcio=1
      nedoub=0
      neall=1
      ntype=0                                                           8d8s14
      include "cas.par"
      ncact=0                                                           1d17s23
      idynw=0                                                           12d31s19
      iwfromfile=0                                                      5d24s19
      myguess=2                                                         12d31s19
      if(idorel4c.ne.0)myguess=1                                        2d27s20
      idefd=0                                                           1s1s20
      idefa=0                                                           1s1s20
      idefo=0                                                           1s1s20
      do isb=1,ngrp                                                     1s1s20
       idoub(isb)=ncoreg(isb)                                           1s1s20
       iacto(isb)=nvalg(isb)                                            1s1s20
       nocc2(isb)=ncoreg(isb)+nvalg(isb)                                1s1s20
      end do                                                            1d2s20
      icasguess=3                                                       12d28s19
      if(idorel4c.ne.0)icasguess=1                                      3d25s20
      write(6,*)('we will be performing a cas calculation')
      idocas=1
      nrydb=0                                                           5d7s18
      nstateset=0                                                       11d22s19
 1363 continue
       ieof=1                                                           5d19s16
       is=1                                                             5d24s18
       read(5,100,end=260)line
       call delim(line,is,ie)                                           6d5s20
       if(line(is:is).eq.'#')go to 1363                                 11d19s20
       if(line(is:is).eq.'&')then
        ieof=0                                                          5d19s16
        go to 260                                                       5d19s16
       end if
 1368  continue                                                         11d21s19
       if(line(is:is+5).eq.'epssym')then                                9d19s23
        is=ie+1                                                         9d19s23
        call delim(line,is,ie)                                          9d19s23
        if(ie.ge.is)then                                                9d19s23
         read(line(is:ie),*)epssym                                      9d19s23
         write(6,*)('epssym read in to be '),epssym                     9d19s23
        else                                                            9d19s23
         write(6,*)('there is no argument for epssym!!!')               9d19s23
         write(6,*)('help!')                                            9d19s23
         irtrn=1                                                        9d19s23
         return                                                         9d19s23
        end if                                                          9d19s23
        go to 1363                                                      9d19s23
       end if                                                           9d19s23
       if(line(is:is+3).eq.'grad')then                                  7d19s22
        write(6,*)('in grad block ')
c
c     options:
c     grade: compute energy gradient.
c     grado: compute orbital gradients as
c     nothing: compute totally symmetric gradients of energy
c     bodc: compute Born-Oppenheimer diagonal correction and derivative
c           couplings if more than one state.
c     nonad: if HF wavefunction, compute BODC as well as second order
c            Born-Oppenheimer correction terms.
c     minimize: optimize geometry to minimize energy.
c           no further argument means minimize the state averaged energy
c           otherwise specify root number and state name.
c     stationary: optimize geometry to minimize gradients
c     idogrado1(1) ne 0 if dE only
c     idogrado1(2) ne 0 if dorb
c     idogrado1(3) ne 0 if bodc
c     idogrado1(4) ne 0 if minimize or stationary.
c     notes:
c          minimize or stationary not compatible with bodc
c     grade turns on idogrado1(1), as do minimize or stationary
c     no arg turns on idogrado1(2)
c     bodc turns on idogrado1(2) and (3)
c     minimize or stationary turn on idogrado1(1) and (4).
c
        if(line(is+4:is+4).eq.'e')then                                  12d19s22
         idogrado1(1)=1                                                  6d18s22
        else                                                            7d25s22
         idogrado1(2)=1                                                 7d25s22
        end if                                                          7d25s22
        is=ie+1                                                         6d18s22
        call delim(line,is,ie)                                          6d18s22
        if(ie.ge.is)then                                                6d18s22
         if(line(is:is).eq.'e')then                                     12d19s22
          idogrado1(1)=1                                                12d19s22
          idogrado1(2)=0                                                12d19s22
          idogrado1(3)=0                                                12d19s22
          idogrado1(4)=0                                                12d19s22
         else if(line(is:is).eq.'b')then                                12d19s22
          idogrado1(1)=0                                                7d25s22
          idogrado1(2)=1                                                6d18s22
          idogrado1(3)=1                                                6d18s22
         else if(line(is:is).eq.'m')then                                6d18s22
          idogrado1(1)=1                                                7d25s22
          idogrado1(2)=0                                                7d25s22
          idogrado1(4)=1                                                6d18s22
         else if(line(is:is).eq.'s')then                                6d18s22
          idogrado1(1)=1                                                7d25s22
          idogrado1(2)=0                                                7d25s22
          idogrado1(4)=2                                                6d18s22
         else                                                           6d18s22
          write(6,*)('I don''t understand grado argument '),
     $           line(is:ie),('!')                                      6d18s22
          irtrn=1                                                       6d18s22
          return                                                        6d18s22
         end if                                                         6d18s22
        end if                                                          6d18s22
        write(6,*)('what we have for idogrado '),idogrado,idogrado1
        go to 1363                                                      7d14s22
       end if                                                           1d11s16
       if(line(is:is+4).eq.'drang')then                                 9d1s23
         is=ie+1                                                        9d1s23
         call delim(line,is,ie)                                         9d1s23
         if(ie.lt.is)then                                               9d1s23
          write(6,*)('missing argument to drangcut keyword!!')          9d1s23
          write(6,*)('help!')                                           9d1s23
          irtrn=1                                                       9d1s23
          return                                                        9d1s23
         else                                                           9d1s23
          read(line(is:ie),*)drangcut                                   9d1s23
          write(6,*)('drangcut read in to be '),drangcut                9d1s23
         end if                                                         9d1s23
         go to 1363                                                     9d1s23
        end if                                                          9d1s23
        if(line(is:is+3).eq.'nlam')then                                 7d3s23
         igu=1                                                          8d31s23
         if(line(ie:ie).eq.'u')igu=2                                    8d31s23
         lambda=0                                                       7d3s23
          write(6,*)                                                    7d3s23
     $      ('if this is a linear molecule, for each lambda the first') 7d3s23
          write(6,*)('lambda no. of functions')                         7d3s23
          if(line(ie:ie).eq.'g')then                                    9d5s23
           write(6,*)('for gerade: ')                                   9d5s23
          else if(line(ie:ie).eq.'u')then                               9d5s23
           write(6,*)('for ungerade: ')                                 9d5s23
          end if                                                        9d5s23
 6699    continue                                                       7d3s23
          is=ie+1                                                        7d3s23
          call delim(line,is,ie)                                        7d3s23
          if(ie.ge.is)then                                              7d3s23
           lambdap=lambda+1                                             7d3s23
           read(line(is:ie),*)nlambda(lambdap,igu)                      8d31s23
           write(6,*)lambda,nlambda(lambdap,igu)                        8d31s23
           lambda=lambdap                                               7d3s23
           go to 6699                                                   7d3s23
          end if
         go to 1363                                                     7d3s23
        end if                                                          7d3s23
        if(line(is:is+4).eq.'debug')then                                1d3s20
 6605   continue                                                        1d3s20
         is=ie+1                                                        1d3s20
         call delim(line,is,ie)                                         1d3s20
         if(ie.ge.is)then                                               1d3s20
          write(6,*)('debug print out turned on for '),line(is:ie)      1d3s20
          nchere=ie+1-is                                                 1d3s20
          do i=1,nprtr                                                   1d3s20
           if(line(is:ie).eq.prtrname(i)(1:nchere))then                 1d3s20
            iprtr(i)=1                                                   1d3s20
            go to 6605                                                  1d3s20
           end if                                                        1d3s20
          end do                                                         1d3s20
          write(6,*)('no matching string in prtrname found ')           1d3s20
         end if                                                         1d3s20
         go to 1363                                                     1d3s20
        end if                                                          1d3s20
       if(line(is:is+2).eq.'csf')then                                   12d24s19
        nmn=1                                                           8d24s22
        is=ie+1                                                         12d24s19
        call delim(line,is,ie)                                          12d24s19
        if(ie.ge.is)then                                                12d24s19
         read(line(is:ie),*)nusecsf                                      12d24s19
        else                                                            12d24s19
         nusecsf=64                                                     12d24s19
        end if                                                          12d24s19
        go to 1363                                                      12d24s19
       end if                                                           12d24s19
       if(line(is:is+4).eq.'snegl')then                                 4d17s23
        is=ie+1                                                         4d17s23
        call delim(line,is,ie)                                          4d17s23
        if(ie.ge.is)then                                                4d17s23
         read(line(is:ie),*)sneglect                                    4d17s23
         write(6,*)                                                     4d17s23
     $        ('threshold for neglecting overlap matrix eigenvalues '), 4d17s23
     $        ('read in to be '),sneglect                               4d17s23
        end if                                                          4d17s23
        go to 1363                                                      4d17s23
       end if                                                           4d17s23
       if(line(is:is+2).eq.'3rd')then                                   8d12s22
        is=ie+1                                                         8d12s22
        call delim(line,is,ie)                                          8d12s22
        if(ie.ge.is)then                                                8d12s22
         read(line(is:ie),*)thirdthres                                  8d12s22
         write(6,*)                                                     8d12s22
     $        ('gradient threshold for switching to 3rd order is now'), 8d12s22
     $        thirdthres                                                8d12s22
        else                                                            8d12s22
         write(6,*)('there is no argument for thirdthres! ')            8d12s22
         write(6,*)('please supply one!')                               8d12s22
         irtrn=1                                                        8d12s22
         return                                                         8d12s22
        end if                                                          8d12s22
        go to 1363                                                      12d24s19
       end if                                                           8d12s22
       if(line(is:is+5).eq.'npdiag')then                                8d12s22
        is=ie+1                                                         8d12s22
        call delim(line,is,ie)                                          8d12s22
        if(ie.ge.is)then                                                8d12s22
         read(line(is:ie),*)npdiag(1)                                   8d12s22
         write(6,*)                                                     8d12s22
     $        ('threshold for using parallel diagonalization is now'),  8d12s22
     $        npdiag(1)                                                 8d12s22
        else                                                            8d12s22
         write(6,*)('there is no argument for npdiag! ')                8d12s22
         write(6,*)('please supply one!')                               8d12s22
         irtrn=1                                                        8d12s22
         return                                                         8d12s22
        end if                                                          8d12s22
        go to 1363                                                      12d24s19
       end if                                                           8d12s22
       if(line(is:is+5).eq.'molden')then                                11d21s19
        molden=1                                                        11d21s19
        open(unit=3,file='molden')                                      11d22s19
        write(3,*)('[Molden Format]')                                   11d22s19
        write(3,*)('[Atoms] AU')                                        11d22s19
        do ia=1,natom                                                   11d22s19
         nuc=nint(atnum(1,ia))                                          11d22s19
         write(3,*)element(nuc),ia,nuc,(xcart(ixyz,ia),ixyz=1,3)        11d22s19
        end do                                                          11d22s19
c
c     since we dump mos in primitive basis, spit out uncontracted
c     info
c
        write(3,*)('[GTO]')                                             11d22s19
        do ia=1,natom                                                   11d22s19
         write(3,*)ia,0
         ig=0                                                           11d22s19
 3136    continue                                                       11d22s19
          iat=ibc(ibdat+ig+ngaus*7)                                     11d22s19
          if(iat.eq.ia.or.iat.eq.-iapair(1,ia))then                     11d22s19
           lhere=ibc(ibdat+ig)
           kbdat=ibc(ibdat+ig+ngaus*8)                                  11d22s19
           if(lhere.gt.4)then                                           11d22s19
            write(6,*)('this basis has l = '),lhere
            write(6,*)('but molden only goes up to g functions')        11d22s19
            irtrn=1                                                     11d22s19
            return                                                      11d22s19
           end if                                                       11d22s19
           if(kbdat.ne.-100000)then                                          11d22s19
            write(3,*)btype(lhere+1),1,('1.00')                         11d22s19
            write(3,*)bc(ibdat+ig+ngaus),1d0                            11d22s19
            ig=ig+1
           else                                                         11d22s19
            ncont=-ibc(ibdat+ig+1+ngaus*8)                              11d22s19
            do jg=ig+1,ngaus-1                                          11d22s19
             itry=ibc(ibdat+jg+ngaus*8)                                 11d22s19
             if(itry.gt.0.or.itry.eq.-1)then                            11d22s19
              npr=jg-ig                                                 11d22s19
              jgx=jg-1                                                  11d22s19
              go to 3236                                                11d22s19
             end if                                                     11d22s19
            end do                                                      11d22s19
            npr=ngaus-ig                                                11d22s19
            jgx=ngaus-1                                                 11d22s19
 3236       continue                                                    11d22s19
            icc=ibdat+kbdat                                             11d22s19
            do ipass=1,ncont                                            11d22s19
             write(3,*)btype(lhere+1),npr,('1.00')                       11d22s19
             do ixp=ig,jgx                                              11d22s19
              iad=icc+ixp-ig
              write(3,*)bc(ibdat+ixp+ngaus),bc(iad)                     11d22s19
             end do                                                     11d22s19
             icc=icc+npr                                                11d22s19
            end do
            ig=ig+npr                                                   11d22s19
           end if                                                       11d22s19
          else                                                          11d22s19
           ig=ig+1                                                      11d22s19
          end if                                                        11d22s19
         if(ig.lt.ngaus)go to 3136                                      11d22s19
         write(3,*)(' ')                                                11d22s19
        end do                                                          11d22s19
        go to 1363                                                      11d21s19
       end if                                                           11d21s19
       if(line(is:is+4).eq.'first')then                                 8d27s19
        is=ie+1                                                         8d27s19
        call delim(line,is,ie)                                          8d27s19
        read(line(is:ie),*)ifirsto                                      8d27s19
        write(6,*)('first iteration will follow energy down hill')      8d27s19
        write(6,*)('for '),ifirsto,(' iterations')                      8d27s19
        go to 1363                                                      8d27s19
       end if                                                           8d27s19
       if(line(is:is+2).eq.'sum')then                                   5d3s21
        write(6,*)('I will write to CAS summary file')                  5d3s21
        nsumcas=1                                                       5d3s21
        go to 1363                                                      8d27s19
       end if                                                           5d3s21
       if(line(is:is+4).eq.'maxps')then                                 4d15s19
        is=ie+1                                                         4d15s19
        call delim(line,is,ie)                                          4d15s19
        read(line(is:ie),*)maxps                                        4d15s19
        write(6,*)('maxps is now '),maxps                               4d15s19
        go to 1363                                                      4d15s19
       end if                                                           4d15s19
       if(line(is:is+2).eq.'ryd')then                                   5d6s20
c
c     this means generate rydberg orbitals. when this is turned on,     5d7s18
c     we automatically set nocan and maxitmac to 1 and guess to readin. 5d7s18
c     the "core states" (the cations) will be given by the states       5d7s18
c     input, and we will perform the weighted average implied by the    5d7s18
c     states input.                                                     5d7s18
c
        is=ie+1                                                         5d7s18
        write(6,*)('we are generating rydberg orbitals ')               5d7s18
        icasguess=3                                                     5d7s18
        itmax=1                                                         5d7s18
        myguess=1
        write(6,*)('setting myguess to 1 ')
        inocan=1                                                        5d7s18
        nrydb=1                                                         1d18s20
        do isb=1,nsymb                                                  9d10s21
         mxryd(isb)=-1                                                  9d10s21
        end do                                                          9d10s21
        do isb=1,nsymb                                                  9d10s21
         is=ie+1                                                        9d10s21
         call delim(line,is,ie)                                         9d10s21
         if(ie.ge.is)then                                               9d10s21
          read(line(is:ie),*)mxryd(isb)                                 9d10s21
          write(6,*)('For symmetry '),isb,                              9d10s21
     $         (' number of Rydberg orbitals to flag is '),mxryd(isb)   9d10s21
         else                                                           9d10s21
          go to 1363                                                    9d10s21
         end if                                                         9d10s21
        end do                                                          9d10s21
        go to 1363                                                      4d29s20
       end if                                                           5d7s18
       if(line(is:is+3).eq.'dynw')then                                  5d7s18
        is=ie+1                                                         4d6s18
        call delim(line,is,ie)                                          4d6s18
        if(ie.gt.0)then                                                 4d6s18
         read(line(is:ie),*)dynw(1)                                     9d20s21
         write(6,*)('dynamic weights turned on with parameter '),       7d27s23
     $        dynw(1),(' Eh')                                           7d27s23
         dynw(1)=1d0/dynw(1)                                            9d20s21
         is=ie+1                                                        12d31s19
         mdynw=2                                                        9d20s21
 2021    continue                                                       9d20s21
          call delim(line,is,ie)                                         12d31s19
          ndot=0                                                         9d20s21
          do ii=is,ie                                                    9d20s21
           if(line(ii:ii).eq.'.')ndot=ndot+1                             9d20s21
          end do                                                         9d20s21
          if(ndot.ne.0)then                                             9d20s21
           if(mdynw.eq.2)then                                           9d20s21
            write(6,*)                                                  9d20s21
     $     ('we will be using 4/(cc+exp(2*wde)) as weighting function') 9d20s21
            read(line(is:ie),*)dynw(2)                                  9d20s21
            write(6,*)('with absolute energy '),dynw(2)                 9d20s21
            mdynw=3                                                     9d20s21
            is=ie+1                                                     9d20s21
            go to 2021                                                  9d20s21
           else                                                         9d20s21
            read(line(is:ie),*)dynw(3)                                  9d20s21
            mdynw=4                                                     9d20s21
           end if                                                       9d20s21
          end if                                                        9d20s21
         if(dynw(2).ne.0d0)then                                         9d20s21
          if(mdynw.eq.3)then                                            9d20s21
           write(6,*)('with cc default value of '),dynw(3)               9d20s21
          else                                                          9d20s21
           write(6,*)('with cc read in to be '),dynw(3)                 9d20s21
          end if                                                        9d20s21
         end if                                                         9d20s21
         if(ie.lt.is)then                                               12d31s19
         write(6,*)('lowest energy will be taken to be the first root ')4d6s18
     $        ,('of the first symmetry specified under states.')         4d6s18
         idynw=0                                                        12d31s19
         else if(line(is:is).eq.'p')then                                12d31s19
          write(6,*)                                                    12d31s19
     $         ('the lowest energy for each symmetry will be used.')    12d31s19
          idynw=1                                                       12d31s19
         else if(line(is:is).eq.'l')then                                3d17s21
          write(6,*)('the root weights will be fixed at their initial') 3d17s21
          write(6,*)                                                    3d17s21
     $         ('values, but will be re-computed for the orbital file') 3d17s21
          idynw=2                                                       3d17s21
         end if                                                         12d31s19
        end if                                                          4d6s18
        go to 1363                                                      4d6s18
       end if                                                           4d6s18
       if(line(is:is+4).eq.'nocan')then                                 4d5s18
        inocan=1                                                        4d5s18
        write(6,*)('orbitals will not be cannonicalized ')              4d5s18
        go to 1363                                                      4d5s18
       end if                                                           4d5s18
       if(line(is:is+3).eq.'wern')then                                  5d19s16
          ix=ie+1                                                       11d29s17
          write(6,*)('werner test '),is,ie
          call delim(line,ix,ie)                                        11d29s17
          write(6,*)('back from delim: '),ix,ie
          if(ie.gt.0)then                                               11d29s17
           write(6,*)('"'),line(ix:ie),('"')
           read(line(ix:ie),*)iwerner                                   11d29s17
          else                                                          11d29s17
           iwerner=2                                                       7d9s13
          end if                                                        11d29s17
        go to 1363                                                      5d19s16
       end if                                                           7d9s13
       if(line(is:is+5).eq.'nowern')then                                5d19s16
        iwerner=0                                                       5d19s16
        write(6,*)('setting iwerner to zero because of nowern? ')
        go to 1363                                                      5d19s16
       end if                                                           7d9s13
         if(line(is:is+8).eq.'maxitmac')then                            4d5s18
          is=ie+1                                                       8d8s14
          call delim(line,is,ie)                                        8d8s14
          read(line(is:ie),*)itmax                                                3d17s04
          go to 1363                                                    4d3s18
         end if                                                         3d17s04
       if(line(is:is+4).eq.'guess')then                                 8d8s14
        is=ie+1                                                         8d8s14
        call delim(line,is,ie)                                          5d19s16
        if(line(is:is+1).eq.'h0'.or.ie.lt.is)then                       12d28s19
         icasguess=1                                                    8d7s06
        else if(line(is:is+1).eq.'hf')then                              8d8s14
         icasguess=2                                                    8d7s06
        else if(line(is:is+5).eq.'readin')then                          8d8s14
         icasguess=3                                                    8d7s06
         is=ie+1                                                        5d19s16
         call delim(line,is,ie)                                         5d19s16
         if(ie.ge.is)then                                               4d22s21
          if(line(is:ie).eq.'ob')then                                    5d19s16
           myguess=2                                                     5d19s16
          else                                                           5d19s16
           myguess=1                                                     5d19s16
          end if                                                         5d19s16
         end if                                                         4d22s21
         write(6,*)('icasguess, myguess '),icasguess,myguess
        end if                                                          8d7s06
        go to 1363                                                      8d7s06
       end if                                                           8d7s06
       if(line(is:is+6).eq.'maxitci')then                               4d5s18
        is=ie+1                                                         8d8s14
        call delim(line,is,ie)                                          8d8s14
        read(line(is:ie),*)maxit                                        8d8s14
        go to 1363
       end if
       if(line(is:is+4).eq.'pthrs')then
        is=ie+1                                                         8d8s14
        call delim(line,is,ie)                                          8d8s14
        read(line(is:ie),*)pthrs
        go to 1363
       end if
       if(line(is:is+4).eq.'econv')then
        is=ie+1
        call delim(line,is,ie)
        read(line(is:ie),*)econv
        go to 1363
       end if
       if(line(is:is+4).eq.'oconv')then                                 12d19s22
        is=ie+1                                                         12d19s22
        call delim(line,is,ie)                                          12d19s22
        read(line(is:ie),*)gconv                                        12d19s22
        go to 1363                                                      12d19s22
       end if                                                           12d19s22
       if(line(is:is+4).eq.'gconv')then
        is=ie+1
        call delim(line,is,ie)
        read(line(is:ie),*)gconv
        go to 1363
       end if
       if(line(is:is+4).eq.'vconv')then
        is=ie+1
        call delim(line,is,ie)
        read(line(is:ie),*)vconv
        go to 1363
       end if
       if(line(is:is+2).eq.'occ')then                                   11d18s19
        ntype=ntype+100                                                 8d8s14
        nc=0
        do i=1,8                                                        8d8s14
         is=ie+1                                                        8d8s14
         call delim(line,is,ie)                                         8d8s14
         if(ie.ge.is)then                                               8d8s14
          nc=nc+1                                                       8d8s14
          read(line(is:ie),*)nocc2(nc)                                  8d8s14
          idefo=1                                                       1s1s20
         else                                                           8d8s14
          if(i.eq.1)then                                                11d18s19
           do j=1,ngrp                                                  11d18s19
            nocc2(j)=ncoreg(j)+nvalg(j)                                 11d18s19
           end do                                                       11d18s19
          end if                                                        11d18s19
          go to 259                                                     8d8s14
         end if                                                         8d8s14
        end do
  259   continue                                                        8d8s14
        ncact=nc                                                        4d24s07
        write(6,*)('no. occupied orbitals'),(nocc2(i),i=1,nc)
        go to 1363                                                      8d8s14
       end if                                                           8d8s14
       if(line(is:is+2).eq.'act')then                                   5d19s16
        ntype=ntype+10                                                  8d8s14
        nc=0                                                            11d18s19
        do i=1,8                                                        8d8s14
         is=ie+1                                                        8d8s14
         call delim(line,is,ie)                                         8d8s14
         if(ie.ge.is)then                                               8d8s14
          nc=nc+1                                                       8d8s14
          read(line(is:ie),*)iacto(nc)                                  8d8s14
          idefa=1                                                       1s1s20
         else                                                           8d8s14
          if(i.eq.1)then                                                11d18s19
           do j=1,ngrp                                                  11d18s19
            iacto(j)=nvalg(j)                                           11d18s19
           end do                                                       11d18s19
          end if                                                        11d18s19
          go to 257                                                     8d8s14
         end if                                                         8d8s14
        end do
  257   continue                                                        8d8s14
        go to 1363                                                      8d8s14
       end if
       if(line(is:is+3).eq.'doub')then                                  5d19s16
        ntype=ntype+1                                                   8d8s14
        nc=0                                                            11d18s19
        do i=1,8
         is=ie+1                                                        8d8s14
         call delim(line,is,ie)                                         8d8s14
         if(ie.ge.is)then                                               8d8s14
          nc=nc+1                                                       8d8s14
          read(line(is:ie),*)idoub(nc)                                  8d8s14
          idefd=1                                                       1s1s20
         else                                                           8d8s14
          if(i.eq.1)then                                                11d18s19
           do j=1,ngrp                                                  11d18s19
            idoub(j)=ncoreg(j)                                          11d18s19
           end do                                                       11d18s19
          end if                                                        11d18s19
          go to 258                                                     8d8s14
         end if                                                         8d8s14
        end do
  258   continue                                                        8d8s14
        write(6,*)('no. of doubly occupied orbitals '),
     $        (idoub(i),i=1,nc)
        nc=max(nc,ncact)                                                4d24s07
        go to 1363                                                      8d8s14
       end if
       if(line(is:is+4).eq.'state')then                                 11d18s19
        nstate=1
        nstasub=1                                                       8d9s22
 1365   continue
         read(5,100,end=1367)line                                       11d21s19
         if(line(1:1).eq.'#')go to 1365                                 12d12s19
         is=1                                                           11d20s19
         call delim(line,is,ie)                                         11d20s19
         iok=0                                                          11d20s19
         do j=1,10                                                      11d20s19
          if(line(is:is).eq.digit(j))then                               11d20s19
           iok=1                                                        11d20s19
           go to 3165                                                   11d20s19
          end if                                                        11d20s19
         end do                                                         11d20s19
         if(line(is:is).eq.'+'.or.line(is:is).eq.'-')then               11d20s19
          iok=1                                                         11d20s19
         end if                                                         11d20s19
 3165    continue                                                       11d20s19
         if(iok.eq.0)then                                               11d21s19
         end if                                                         11d21s19
         if(iok.ne.0)then                                               11d20s19
          ipass=1
          is=1
          ntop=4
          if(nstate.gt.maxst1)then
           write(6,*)('trying to read in more state groups than')
           write(6,*)('dimensions = '),maxst1
           write(6,*)('increase maxst1 in common.cas')                  5d10s19
           irtrn=1                                                      8d7s06
           return                                                       8d7s06
          end if
          do j=1,maxst2                                                 6d18s20
           eweight(j,nstate)=1d0                                        6d18s20
          end do                                                        6d18s20
 1366     continue
          call delim(line,is,ie)                                        8d8s14
c
c     ultimately: charge,spin,sym,root,wgt, to mimic term symbols
c     if neutral, no charge needs to be given,                          11d18s19
c     if anion, give charge as -, or -1, or -2, etc                     11d18s19
c     if cation, give charge as +, or +1, or +2, etc                    11d18s19
c     one can also specify the charge as +0 or -0 to get neutral.       11d18s19
c     ipass = 1: - no. electrons                                        5d3s07
c     ipass = 2: term symbol
c     and in all cases                                                  5d3s07
c     ipass = 3: no. roots
c     ipass > 3 : weight for root
c     alternatively, instead of weights, one can put "file" which means 5d24s19
c     initial weights are read from orbital guess file.                 5d24s19
c     if no. roots is omitted, it defaults to one, and the weight to one
c     as well.                                                          12d24s19
c
          if(ipass.eq.1)then                                            11d18s19
           if(line(is:is).eq.'+')then                                   11d18s19
            if(is.eq.ie)then                                            11d18s19
             isinfo(1,nstate)=numeminus-1                               11d18s19
            else                                                        11d18s19
             read(line(is+1:ie),*)ichrg                                 11d18s19
             isinfo(1,nstate)=numeminus-ichrg                           11d18s19
            end if                                                      11d18s19
            is=ie+1                                                     11d18s19
            ipass=ipass+1                                               11d18s19
            go to 1366                                                  11d18s19
           else if(line(is:is).eq.'-')then                              11d18s19
            if(is.eq.ie)then                                            11d18s19
             isinfo(1,nstate)=numeminus+1                               11d18s19
            else                                                        11d18s19
             read(line(is+1:ie),*)ichrg                                 11d18s19
             isinfo(1,nstate)=numeminus+ichrg                           11d18s19
            end if                                                      11d18s19
            is=ie+1                                                     11d18s19
            ipass=ipass+1                                               11d18s19
            go to 1366                                                  11d18s19
           else                                                         11d18s19
            isinfo(1,nstate)=numeminus                                  11d18s19
            ipass=2                                                     11d18s19
           end if                                                       11d18s19
          end if                                                        11d18s19
          if(ipass.gt.ntop)stop
          if(ipass.le.ntop)then
           if(ipass.eq.2)then                                           12d24s19
c
c     if we are diatomic, we can also specify 2 atomic states ...
c
            if(natom.eq.2)then                                          8d8s22
c     either 2 or 3 characters for atomic symbol ... if 3 characters
c     and last two are not "Pi", then atomic symbol, unless S1.
c
             ncha=ie+1-is                                               8d8s22
             if(ncha.le.3)then                                          8d8s22
              if(ncha.eq.3)then                                         8d8s22
               do i=1,nsymb                                             4d21s23
                if(line(is+1:ie).eq.stype(i,nsymb))go to 1492           4d21s23
               end do                                                   4d21s23
               if(line(is+2:is+2).ne.'o')then                           5d2s24
                if(line(is+1:ie).eq.'Pi'.or.(line(is+1:is+1).eq.'S'.and. 1d2s24
     $              line(is+2:is+2).ne.'i'))then                        5d2s24
                 go to 1492                                              8d8s22
                end if                                                  5d2s24
               end if                                                   8d8s22
              end if                                                    8d8s22
              write(6,*)('yes! ')                                       8d8s22
              write(6,*)('first atomic asymtote: '),line(is:ie)
              read(line(is:is),*)ispina                                 8d8s22
              asym(1)='  '                                              8d8s22
              asym(1)(1:ncha-1)=line(is+1:ie)                              8d8s22
              is=ie+1                                                   8d8s22
              call delim(line,is,ie)                                    8d8s22
              if(ie.lt.is)then                                          8d8s22
               write(6,*)('sorry, I can not find second asymptote!')    8d8s22
               irtrn=1                                                  8d8s22
               return                                                   8d8s22
              end if                                                    8d8s22
              write(6,*)('second asymptote: '),line(is:ie)              8d8s22
              nchb=ie-is                                                8d8s22
              asym(2)='  '                                              8d8s22
              asym(2)(1:nchb)=line(is+1:ie)                             8d8s22
              read(line(is:is),*)ispinb                                 8d8s22
              is=ie+1                                                   8d8s22
              call delim(line,is,ie)                                    8d8s22
              if(ie.ge.is)then                                          8d8s22
               write(6,*)('weight will be read in from '),line(is:ie)   8d8s22
               read(line(is:ie),*)aww                                   8d8s22
              else                                                      8d8s22
               aww=1d0                                                  8d8s22
              end if                                                    8d8s22
              lspin=iabs(ispina-ispinb)
              ispin=ispina+ispinb                                       8d8s22
              nspin=(ispin-lspin)/2
              jspin=ibcoff                                              8d8s22
              ibcoff=jspin+nspin                                        8d8s22
              call enough('basisz. 10',bc,ibc)
              jspinm=jspin-1                                            8d8s22
              ibc(jspin)=lspin+1                                        8d8s22
              do i=2,nspin                                              8d8s22
               im=i-1                                                   8d8s22
               ibc(jspinm+i)=ibc(jspinm+im)+2                           8d8s22
              end do                                                    8d8s22
              write(6,*)('spins generated: '),(ibc(jspinm+i),           8d8s22
     $             i=1,nspin)                                           8d8s22
              isum=ispina+ispinb                                        8d8s22
              nec=isinfo(1,nstate)                                      8d8s22
              if(isinfo(2,nstate).eq.0)nstate=nstate-1                  8d9s22
              nlzz=2                                                    8d8s22
              call genr(asym,asym(2),nsymb.eq.8.,ibc(jspin),nspin,isum, 8d8s22
     $             aww,isinfo,nstate,maxst1,nec,eweight,maxst2,         8d8s22
     $             pspacex,pthrs,ispina)                                6d16s23
              nstasub=0                                                 8d9s22
              write(6,2022)                                             2d4s23
 2022         format('states so far:',/,                                2d4s23
     $             '     #    ne  Smult   sym nroot   Lz')                 2d4s23
              do ii=1,nstate
               do i=6,11                                                 8d8s22
                if(isinfo(i,ii).ne.0)then
                 line(i:i)=char(isinfo(i,ii))
                else
                 line(i:i)=' '
                end if
               end do
               write(6,2023)ii,(isinfo(j,ii),j=1,5),line(6:11)          2d4s23
 2023          format(6i6,1x,a5)                                        2d4s23
              end do
              go to 1365                                                8d8s22
             end if                                                     8d8s22
            end if                                                      8d8s22
 1492       continue                                                    8d8s22
            do ilook=is,ie                                              12d24s19
             itry=ichar(line(ilook:ilook))                              12d24s19
             if(itry.lt.idiglow.or.itry.gt.idighig)then                 12d24s19
              ilm=ilook-1                                               12d24s19
              read(line(is:ilm),*)issmult                               12d24s19
              is=ilook                                                  12d24s19
              go to 1047                                                12d24s19
             end if                                                     12d24s19
            end do                                                      12d24s19
            write(6,*)('could not find spin multiplicity in "'),        12d24s19
     $           line(is:ie),('"')                                      12d24s19
            stop                                                        12d24s19
 1047       continue                                                    12d24s19
            do j=6,11                                                    12d28s19
             isinfo(j,nstate)=0                                         12d28s19
            end do                                                      12d28s19
            nshere=ie+1-is                                              12d28s19
            do j=1,min(nshere,6)                                        12d28s19
             js=is+j-1                                                  12d28s19
             isinfo(5+j,nstate)=ichar(line(js:js))                      12d28s19
            end do                                                      12d28s19
            if(ie.eq.is.or.line(is+1:is+1).eq.'o')then                  12d31s19
             nlzz=6                                                     12d31s19
             if(line(is:is).eq.'S')then                                 12d31s19
              isinfo(5,nstate)=0                                        12d31s19
              if(is.eq.ie)then                                          12d31s19
               issym=1                                                  12d31s19
              else                                                      12d31s19
               issym=8                                                  12d31s19
              end if                                                    12d31s19
             else                                                       12d31s19
              if(is.eq.ie)then                                          12d31s19
               issym=6                                                  12d31s19
              else                                                      12d31s19
               issym=2                                                  12d31s19
              end if                                                    12d31s19
              if(line(is:is).eq.'P')then                                12d31s19
               isinfo(5,nstate)=1                                       12d31s19
               if(is.eq.ie)then                                         4d30s21
                issym=4                                                 4d30s21
               else                                                     4d30s21
                issym=5                                                 4d30s21
               end if                                                   4d30s21
              else if(line(is:is).eq.'D')then                           12d31s19
               isinfo(5,nstate)=2                                       12d31s19
               if(is.eq.ie)then                                         4d30s21
                issym=1                                                 4d30s21
               else                                                     4d30s21
                issym=8                                                 4d30s21
               end if                                                   4d30s21
              else if(line(is:is).eq.'F')then                           12d31s19
               isinfo(5,nstate)=3                                       12d31s19
               if(is.eq.ie)then                                         4d30s21
                issym=4                                                 4d30s21
               else                                                     4d30s21
                issym=5                                                 4d30s21
               end if                                                   4d30s21
              else if(line(is:is).eq.'G')then                           12d31s19
               isinfo(5,nstate)=4                                       12d31s19
               if(is.eq.ie)then                                         4d30s21
                issym=1                                                 4d30s21
               else                                                     4d30s21
                issym=8                                                 4d30s21
               end if                                                   4d30s21
              else                                                      12d31s19
               write(6,*)('unknown atomic symmetry: '),line(is:is)      12d31s19
               irtrn=1                                                  12d31s19
               return                                                   12d31s19
              end if                                                    12d31s19
             end if                                                     12d31s19
            else if(line(is:is+2).eq.'Sig')then                         12d31s19
             if(nsymb.eq.8.and.                                         10d3s21
     $            .not.(line(ie:ie).eq.'g'.or.line(ie:ie).eq.'u'))then  10d3s21
              write(6,*)
              write(6,*)('!!!!!looking for g or u at end of "'),             10d3s21
     $            line(is:ie),('" but can not find it. Help!')          10d3s21
              write(6,*)
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             iem=ie-1                                                   3d3s21
             nlzz=2                                                     12d28s19
             isinfo(5,nstate)=0                                         12d28s19
             if(.not.(line(ie:ie).eq.'-'.or.line(ie:ie).eq.'+'.or.      3d23s21
     $            line(iem:iem).eq.'-'.or.line(iem:iem).eq.'+'))then    3d23s21
              write(6,*)('I''m looking for + or - at the end of '),
     $            line(is:ie),('but I can''t find it!')                 3d19s21
              write(6,*)('this is an input error!!')                    3d19s21
              irtrn=1                                                   3d19s21
              return                                                    3d19s21
             end if                                                     3d19s21
             if(line(ie:ie).eq.'+'.or.line(iem:iem).eq.'+')then         3d3s21
              issym=1                                                   12d28s19
             else                                                       12d24s19
              issym=4                                                   12d28s19
             end if                                                     12d24s19
             if(line(ie:ie).eq.'u')issym=issym+4                        3d3s21
            else if(line(is:is+1).eq.'Pi')then                          12d28s19
             if(nsymb.eq.8.and.                                         10d3s21
     $            .not.(line(ie:ie).eq.'g'.or.line(ie:ie).eq.'u'))then  10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line(is:ie),('" but can not find it. Help!')          10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzz=2                                                     12d28s19
             isinfo(5,nstate)=1                                         12d28s19
             if(is+1.eq.ie)then                                         12d24s19
              issym=2                                                   12d28s19
             else if(line(is+2:is+2).eq.'u')then                        12d28s19
              issym=2                                                   12d28s19
             else                                                       12d24s19
              issym=6                                                   12d28s19
             end if                                                     12d24s19
            else if(line(is:is+2).eq.'Del')then                         12d28s19
             if(nsymb.eq.8.and.                                         10d3s21
     $            .not.(line(ie:ie).eq.'g'.or.line(ie:ie).eq.'u'))then  10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line(is:ie),('" but can not find it. Help!')          10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzz=2                                                     12d28s19
             isinfo(5,nstate)=2                                         12d28s19
             if(is+2.eq.ie.or.(line(ie:ie).ne.'g'.and.                  3d3s21
     $            line(ie:ie).ne.'u'))then                              3d3s21
              issym=1                                                   12d28s19
             else if(line(ie:ie).eq.'g')then                            3d3s21
              issym=1                                                   12d28s19
             else                                                       12d24s19
              issym=5                                                   12d28s19
             end if                                                     12d24s19
            else if(line(is:is+2).eq.'Phi')then                         12d28s19
             if(nsymb.eq.8.and.                                         10d3s21
     $            .not.(line(ie:ie).eq.'g'.or.line(ie:ie).eq.'u'))then  10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line(is:ie),('" but can not find it. Help!')          10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzz=2                                                     12d28s19
             isinfo(5,nstate)=3                                         12d28s19
             if(is+2.eq.ie)then                                         12d24s19
              issym=2                                                   12d28s19
             else if(line(is+3:is+3).eq.'u')then                        12d28s19
              issym=2                                                   12d28s19
             else                                                       12d24s19
              issym=6
             end if                                                     12d24s19
            else if(line(is:is+2).eq.'Gam')then                         12d28s19
             if(nsymb.eq.8.and.                                         10d3s21
     $            .not.(line(ie:ie).eq.'g'.or.line(ie:ie).eq.'u'))then  10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line(is:ie),('" but can not find it. Help!')          10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzz=2                                                     12d28s19
             isinfo(5,nstate)=4                                         12d28s19
             if(is+2.eq.ie.or.(line(ie:ie).ne.'g'.and.                  3d3s21
     $            line(ie:ie).ne.'u'))then                              3d3s21
              issym=1                                                   12d28s19
             else if(line(ie:ie).eq.'g')then                            3d3s21
              issym=1                                                   12d28s19
             else                                                       12d24s19
              issym=5                                                   12d28s19
             end if                                                     12d24s19
            else if(line(is:is).eq.'S')then                             12d28s19
             isp=is+1                                                   12d24s19
             read(line(isp:ie),*)issym                                  12d24s19
            else                                                        12d24s19
             if(ngrp.eq.1)then                                          12d24s19
              issym=1                                                   12d24s19
             else if(ngrp.eq.2)then                                     12d24s19
              if(line(is+2:is+2).eq.''''.or.line(is+1:is+1).eq.'"'.or.  12d24s19
     $             line(is+1:is+1).eq.'2')then                          12d24s19
               issym=2                                                  12d24s19
              else                                                      12d24s19
               issym=1                                                  12d24s19
              end if                                                    12d24s19
             else                                                       12d24s19
              do isb=1,ngrp                                              12d24s19
               if(line(is:is+2).eq.stype(isb,ngrp))then                 12d24s19
                issym=isb                                               12d24s19
                go to 2047                                              12d24s19
               end if                                                   12d24s19
              end do                                                     12d24s19
              write(6,*)('bad term symbol: '),line(is:ie)
              write(6,*)('we were searching from the list ')            1d14s20
              do isb=1,ngrp                                             1d14s20
               write(6,*)stype(isb,ngrp)                                1d14s20
              end do                                                    1d14s20
              stop                                                      12d24s19
             end if                                                     12d24s19
            end if                                                      12d24s19
 2047       continue                                                    12d24s19
c
c     if nlzz = 2, we must be along z axis. check this now.
c
            if(nlzz.eq.2)then                                           12d31s19
             rms=0d0                                                     12d31s19
             do ia=1,natom                                               12d31s19
              rms=rms+xcart(1,ia)**2+xcart(2,ia)**2                      12d31s19
             end do                                                      12d31s19
             rms=sqrt(rms/dfloat(2*natom))                               12d31s19
             if(rms.gt.1d-6)then                                        12d31s19
              write(6,*)('wait ... you have asked for linear molecule'),12d31s19
     $            (' symmetry, but either your geometry is non-linear,')12d31s19
              write(6,*)('or your molecule is not along the z axis.')   12d31s19
              irtrn=1                                                   12d31s19
              return                                                    12d31s19
             end if                                                     12d31s19
            end if                                                      12d31s19
c
c     if nlzz = 6, we must atom at origin. check this now.
c
            if(nlzz.eq.6)then                                           12d31s19
             rms=0d0                                                    12d31s19
             do ixyz=1,3                                                12d31s19
              rms=rms+xcart(ixyz,1)**2                                  12d31s19
             end do                                                      12d31s19
             rms=sqrt(rms/3d0)                                          12d31s19
             if(rms.gt.1d-6.or.natom.gt.1)then                          12d31s19
              write(6,*)('wait ... you have asked for spherical'),      12d31s19
     $            (' symmetry, but either you have more than one atom,')12d31s19
              write(6,*)('or your atom is not at the origin.')          12d31s19
              irtrn=1                                                   12d31s19
              return                                                    12d31s19
             end if                                                     12d31s19
            end if                                                      12d31s19
            isinfo(2,nstate)=issmult                                    12d24s19
            isinfo(3,nstate)=issym                                      12d24s19
            ipass=4                                                     12d24s19
            is=ie+1                                                     12d24s19
            go to 1366                                                  12d24s19
           end if                                                       12d24s19
           if(ipass.le.4)then                                           5d1s18
            if(ie.lt.is)then                                            12d24s19
             isinfo(4,nstate)=1                                         12d24s19
             eweight(1,nstate)=1d0                                      12d24s19
            else                                                        12d24s19
             if(line(is:is).eq.'f')then                                 1d18s20
              iwfromfile=1                                              1d18s20
              isinfo(4,nstate)=1                                        1d18s20
              eweight(1,nstate)=1d0                                     1d18s20
             else                                                       1d18s20
              write(6,*)('going after roots? '),line(is:ie)
              isrt=is
              isinfo(ipass,nstate)=0                                    4d24s21
              do iii=is+1,ie-1                                          4d24s21
               if(line(iii:iii).eq.'+')then                             4d24s21
                read(line(isrt:iii-1),*)irtp                            4d24s21
                isinfo(ipass,nstate)=isinfo(ipass,nstate)+irtp          4d24s21
                isrt=iii+1                                              4d24s21
               end if                                                   4d24s21
              end do                                                    4d24s21
              if(isrt.gt.is)then                                        4d24s21
               if(isrt.le.ie)then                                       4d24s21
                read(line(isrt:ie),*)irtp                               4d24s21
                isinfo(ipass,nstate)=isinfo(ipass,nstate)+irtp          4d24s21
                write(6,*)('total number of roots: '),
     $               isinfo(ipass,nstate)
               end if                                                   4d24s21
              else                                                      4d24s21
               read(line(is:ie),*)isinfo(ipass,nstate)
              end if                                                    4d24s21
             end if                                                     1d18s20
            end if                                                      12d24s19
           else                                                         5d1s18
            if(ie.ge.is)then                                            6d18s20
             if(line(is:is).eq.'f')then
              write(6,*)('eweights will be taken from orbital file')
              iwfromfile=1                                              5d24s19
             else                                                        5d24s19
              read(line(is:ie),*)eweight(ipass-4,nstate)                  5d1s18
             end if                                                      5d24s19
            end if                                                      5d24s19
           end if                                                       5d1s18
          end if
          if(ipass.eq.4)then
            noelec=iabs(isinfo(1,nstate))                               5d3s07
            write(6,*)('noelec for state '),noelec,nstate
            issmult=isinfo(2,nstate)                                    5d3s07
            issym=isinfo(3,nstate)                                      5d3s07
            write(6,*)('nstasub = '),nstasub                            10d12s22
            if(nstasub.eq.1)then                                        10d12s22
             isinfo(1,nstate)=issym                                      5d3s07
             isinfo(2,nstate)=issmult                                    5d3s07
             isinfo(3,nstate)=noelec                                     5d3s07
            end if                                                      10d12s22
            write(6,*)('isinfo '),(isinfo(j,nstate),j=1,3)
           ntop=ntop+isinfo(ipass,nstate)                               2d21s07
           if(isinfo(ipass,nstate).gt.maxst2)then
            write(6,*)('more roots requested '),isinfo(ipass,nstate)
            write(6,*)('than dimensions '),maxst2
            write(6,*)('increase maxst2 in common.cas')                 5d10s19
            irtrn=1                                                      8d7s06
            return                                                       8d7s06
           end if
          end if
          ipass=ipass+1
          is=ie+1
          if(ipass.le.ntop.and.ie.gt.0)go to 1366                       12d24s19
           if(nlzz.ne.0.and.isinfo(5,nstate).ne.0)then                  1d5s20
            write(6,*)('now we need to consider degeneracies')          1d5s20
            if(nlzz.eq.2)then                                           1d5s20
             nstatep=nstate+1                                           1d5s20
             do j=1,11                                                  1d5s20
              isinfo(j,nstatep)=isinfo(j,nstate)                        1d5s20
             end do                                                     1d5s20
             do j=1,isinfo(4,nstate)                                    1d5s20
              eweight(j,nstatep)=eweight(j,nstate)                      1d5s20
             end do                                                     1d5s20
             ntop=ntop+isinfo(4,nstate)                                 1d5s20
             if(mod(isinfo(5,nstate),2).eq.1)then                       1d5s20
              isinfo(1,nstatep)=isinfo(1,nstate)+1                      1d5s20
             else                                                       1d5s20
              isinfo(1,nstatep)=isinfo(1,nstate)+3                      1d5s20
             end if                                                     1d5s20
             pspacex(nstate)=pthrs                                      1d5s20
             nstate=nstatep                                             1d5s20
            else                                                        1d5s20
             if(isinfo(5,nstate).eq.1)then                              1d5s20
              ndegm=2                                                   1d5s20
             else                                                       1d5s20
              ndegm=3                                                   1d5s20
             end if                                                     1d5s20
             nstatem=nstate-1                                           1d5s20
             do j=1,ndegm                                               1d5s20
              pspacex(nstatem+j)=pthrs                                  1d5s20
              do l=1,11                                                 1d5s20
               isinfo(l,nstate+j)=isinfo(l,nstate)                      1d5s20
              end do                                                    1d5s20
             end do                                                     1d5s20
             do j=1,ndegm                                               1d5s20
              do l=1,isinfo(4,nstate)                                   1d5s20
               eweight(l,nstate+j)=eweight(l,nstate)                    1d5s20
              end do                                                    1d5s20
             end do                                                     1d5s20
             if(isinfo(5,nstate).eq.1)then                              1d5s20
              if(isinfo(1,nstate).eq.5)then                             4d30s21
               isinfo(1,nstate+1)=2                                     4d30s21
               isinfo(1,nstate+2)=3                                     4d30s21
              else                                                      1d5s20
               isinfo(1,nstate+1)=6                                     4d30s21
               isinfo(1,nstate+2)=7                                     4d30s21
              end if                                                    1d5s20
             else if(isinfo(5,nstate).eq.2)then                         1d5s20
              if(isinfo(1,nstate).eq.8)then                             4d30s21
               isinfo(1,nstate+1)=2                                     4d30s21
               isinfo(1,nstate+2)=3                                     4d30s21
               isinfo(1,nstate+3)=5                                     1d9s18
              else                                                      1d5s20
               isinfo(1,nstate+1)=6                                     4d30s21
               isinfo(1,nstate+2)=7                                     4d30s21
               isinfo(1,nstate+3)=4                                     4d30s21
              end if                                                    1d5s20
              ll=1
              do l=1,isinfo(4,nstate)                                   1d5s20
               do j=1,2                                                 1d5s20
                eweight(ll,nstate)=eweight(l,nstate+1)                  4d30s21
                ll=ll+1                                                 1d5s20
               end do                                                   1d5s20
              end do                                                    1d5s20
              isinfo(4,nstate)=2*isinfo(4,nstate)                       4d30s21
             else                                                       1d5s20
              write(6,*)('haven''t yet figured out degeneracies ...')   1d5s20
              irtrn=1                                                   1d5s20
              return                                                    1d5s20
             end if                                                     1d5s20
             nstate=nstate+ndegm                                        1d5s20
            end if                                                      1d5s20
           end if                                                       1d5s20
          write(6,*)('for state symmetry/spin group '),nstate
          write(6,*)('symmetry '),isinfo(1,nstate)
          write(6,*)('spin multiplicity '),isinfo(2,nstate)
          write(6,*)('no. active electrons '),isinfo(3,nstate)
          write(6,*)('no. roots '),isinfo(4,nstate)
          pspacex(nstate)=pthrs                                         5d1s18
          write(6,*)('weights '),(eweight(j,nstate),                    5d1s18
     $         j=1,isinfo(4,nstate))                                    5d1s18
          nstate=nstate+1
          go to 1365
         else
          nstate=nstate-nstasub                                         8d9s22
          if(nstasub.eq.0)then                                          8d9s22
           write(6,*)('nowa swap 1 and 3: ')
           do j=1,nstate                                                8d9s22
            icpy=isinfo(1,j)                                            8d9s22
            isinfo(1,j)=isinfo(3,j)                                     8d9s22
            isinfo(3,j)=icpy                                            8d9s22
           end do                                                       8d9s22
           if(isinfo(3,nstate).eq.0)nstate=nstate-1                     10d12s22
          end if                                                        8d9s22
          ieof=0                                                        11d21s19
         end if
        end if
        if(line(is:is).eq.'&')then                                      11d21s19
         ieof=0                                                         11d21s19
         go to 260                                                      11d21s19
        end if                                                          11d21s19
        if(iok.ne.0)go to 1363                                          11d21s19
        go to 260                                                       11d21s19
 1367   continue                                                        11d21s19
        nstate=nstate-nstasub                                           8d9s22
          if(nstasub.eq.0)then                                          8d9s22
           do j=1,nstate                                                8d9s22
            icpy=isinfo(1,j)                                            8d9s22
            isinfo(1,j)=isinfo(3,j)                                     8d9s22
            isinfo(3,j)=icpy                                            8d9s22
           end do                                                       8d9s22
           if(isinfo(3,nstate).eq.0)nstate=nstate-1                     10d12s22
          end if                                                        8d9s22
  260   continue                                                        8d8s14
        if(nstateset.eq.0)then                                          11d22s19
         nstateset=1                                                     11d22s19
         ngoto=ntype/100
         ntype=ntype-ngoto*100
         ngota=ntype/10
         ngotd=ntype-ngota*10
         if(ngoto.eq.0.and.idefo.eq.0)then                              1s1s20
          do i=1,nsymb
           nocc2(i)=iacto(i)+idoub(i)                                    8d8s14
          end do                                                         8d8s14
         end if                                                          8d8s14
         if(ngota.eq.0.and.idefa.eq.0)then                              1s1s20
          do i=1,nsymb                                                   8d8s14
           iacto(i)=nocc2(i)-idoub(i)                                    8d8s14
          end do                                                         8d8s14
         end if                                                          8d8s14
         if(ngotd.eq.0.and.idefd.eq.0)then                              1d2s20
          do i=1,nsymb                                                   8d8s14
           idoub(i)=nocc2(i)-iacto(i)                                    8d8s14
          end do                                                         8d8s14
         end if                                                          8d8s14
         nedoub=0                                                        8d8s14
         do i=1,nsymb                                                    8d8s14
          nedoub=nedoub+idoub(i)*2                                       8d8s14
         end do                                                          8d8s14
         if(idorel4c.ne.0)then                                          2d26s20
          nstate=1                                                      2d26s20
          isinfo(3,1)=numeminus                                         2d26s20
          isinfo(4,1)=1                                                 2d26s20
          neall=1                                                       2d26s20
         end if                                                         2d26s20
         do istat=1,nstate
          if(neall.eq.1)then
           isinfo(3,istat)=isinfo(3,istat)-nedoub
          else
          end if
         end do                                                          8d8s14
         write(6,3315)(idoub(i),i=1,nsymb)                               8d8s14
 3315    format('number of doubly occupied orbitals:',8i3)               8d8s14
         write(6,3325)(iacto(i),i=1,nsymb)                               8d8s14
 3325    format('number of active orbitals:         ',8i3)               8d8s14
         write(6,3335)(nocc2(i),i=1,nsymb)                               8d8s14
 3335    format('number of occupied orbitals:       ',8i3)               8d8s14
         nedoub=0
         do i=1,nsymb                                                   12d8s22
          nedoub=nedoub+idoub(i)*2
         end do
         casdat(1)=pthrs
         casdat(2)=econv
         casdat(3)=gconv
         casdat(4)=dfloat(maxit)
         casdat(5)=vconv                                                12d9s22
         casdat(6)=oconv                                                12d19s22
        end if                                                          11d22s19
        if(ieof.eq.0)then                                               5d19s16
         if(line(is:is).ne.'&')go to 1368                               11d21s19
         go to 254                                                      11d6s19
        else                                                            5d19s16
         sz=0d0                                                         3d23s23
         do ia=1,natom                                                  3d23s23
          sz=sz+xm(ia)                                                  3d23s23
         end do                                                         3d23s23
         if(sz.gt.1d0)then                                              3d23s23
          itmpm=ibcoff                                                    3d22s23
          ibcoff=itmpm+natom                                              3d22s23
          call enough('basisz.tmpm',bc,ibc)                               3d22s23
          jtmp=itmpm-1                                                    3d22s23
          do i=1,natom                                                   3d22s23
           bc(jtmp+i)=xm(mya(i))                                         3d22s23
          end do
          do i=1,natom                                                   3d22s23
           if(iapair(1,i).gt.0)bc(jtmp+i)=bc(jtmp+iapair(1,i))           3d22s23
          end do                                                         3d22s23
          do i=1,natom                                                   3d22s23
           xmaxs(i)=bc(jtmp+i)                                              3d22s23
          end do                                                         3d22s23
          ibcoff=itmpm                                                   3d22s23
         end if                                                         3d23s23
         return                                                         5d19s16
        end if                                                          5d19s16
c
c     system hardware parameters                                        1d31s21
c
        iddi=1                                                          9d28s22
 3253   continue                                                        1d31s21
        read(5,200,end=1008)line200                                     1d31s21
        is=1                                                            1d31s21
        write(6,*)('what we have for line200 '),line200
        call delim(line200,is,ie)                                       1d31s21
        write(6,*)('start: '),line200(is:ie)
        if(line200(is:is).eq.'&'.or.line200(is:is+1).eq.'#&')then             5d25s18
         line(1:10)=line200(1:10)                                       11d22s19
         go to 254                                                      5d25s18
        end if                                                          5d25s18
        if(line200(is:is+2).eq.'ddi')then                               9d28s22
         is=ie                                                          9d28s22
         call delim(line200,is,ie)                                      9d28s22
         if(ie.lt.is)then                                               9d28s22
          write(6,*)('there is no argument on the ddi line')            9d28s22
          write(6,*)('help!')                                           9d28s22
          irtrn=1                                                       9d28s22
          return                                                        9d28s22
         end if                                                         9d28s22
         read(line200(is:ie),*)iddi                                     9d28s22
         write(6,*)('iddi is now '),iddi                                9d28s22
         go to 3253                                                     9d28s22
        end if                                                          9d28s22
        if(line200(is:is+4).eq.'ncore')then                             1d31s21
         is=ie+1                                                        1d31s21
         call delim(line200,is,ie)                                      1d31s21
         read(line200(is:ie),*)ncore                                    1d31s21
         write(6,*)('ncore is read in to be '),ncore                    1d31s21
        end if                                                          1d31s21
        if(line200(is:is+2).eq.'max')then                               1d31s21
         is=ie+1                                                        1d31s21
         call delim(line200,is,ie)                                      1d31s21
         read(line200(is:ie),*)maxbct                                   1d31s21
         write(6,*)('max memory read in to be '),maxbct                 1d31s21
         is=ie+1                                                        1d31s21
         call delim(line200,is,ie)                                      1d31s21
         imultpr=1                                                      1d31s21
         if(ie.ge.is)then                                               1d31s21
          if(line200(is:is).eq.'k'.or.line200(is:is).eq.'K')then        1d31s21
           imultpr=1000                                                 1d31s21
          else if(line200(is:is).eq.'m'.or.line200(is:is).eq.'M')then   1d31s21
           imultpr=1 000 000                                            1d31s21
          else if(line200(is:is).eq.'g'.or.line200(is:is).eq.'G')then   1d31s21
           imultpr=1 000 000 000                                        1d31s21
          else                                                          1d31s21
           write(6,*)('unknown multiplier: '),line200(is:ie)            1d31s21
           irtrn=1                                                      1d31s21
           return                                                       1d31s21
          end if                                                        1d31s21
          write(6,*)('multiplier: '),imultpr                            1d31s21
          maxbct=maxbct*imultpr                                         1d31s21
         end if                                                         1d31s21
        end if                                                          1d31s21
        go to 3253                                                      1d31s21
c
c     mrci input
c     key words:
c     restart
c     save
c     npdiag
c     maxit
c     maxiti
c     maxds
c     maxdsi                                                            10d14s22
c     debug
c     cache
c     norb
c     2den
c     novguess
c     nowav
c     econv
c     wconv
c     toldv                                                             8d7s24
c     eiconv
c     cren
c     nointeract
c     maxopens
c     test
c     maxps
c     dyn
c     xref
c     pthrs
c     doubo
c     ref
c     0hole
c     1hole
c     2hole
c     1fill
c     2fill
c     state
c     sing
c     decomp
c     readc
c     doubx
c     itabort
c     skip4v
c     epssymci                                                          9d19s23
c     ethred
c
 2251   continue                                                        5d24s18
        include "mrci.par"                                              12d8s22
        nmn=1                                                           8d12s22
        nrxinfos0=0                                                     2d8s22
        norbci=0                                                        1d24s22
        interacts=1                                                     10d15s19
        nwavef=1                                                        3d23s21
        idomrci=1                                                       5d24s18
        iunc(1)=0                                                       3d29s19
        iunc(2)=0                                                       3d29s19
        impspin=0                                                       3d29s19
        mcache=0                                                        2d10s20
        maxopens=1000                                                   6d12s19
        nxref=0                                                         5d10s19
        novguess=0                                                      2d14s20
        nrxinfos=0                                                      4d24s21
        do i=1,8                                                        5d23s19
         irxinfo(i)=0                                                   5d23s19
        end do                                                          5d23s19
        nroot=-1                                                        4d28s21
        write(6,*)('we are going to do a mrci calculation')             5d24s18
        if(numeminus.eq.0)then                                          11d20s19
         open(unit=1,file='forbs',form='unformatted')                      11d20s19
         write(6,*)('getting input from forbs ')                           11d20s19
         call loadr(potdws,ibdat,ibstor,isstor,nsymb,istinfo,ncoreg,    1d2s20
     $        nvalg,nlzzq,iextradatd,nectradatx,0,nbasdws,element,idum, 7d26s22
     $        idum,idum,idum,1,bc,ibc,icanon)                           5d5s23
         close(unit=1)                                                     11d20s19
         first(1)=potdws                                                11d20s19
         ngrp=nsymb                                                     12d27s19
         ne=istinfo(3)                                                   1d2s20
         ismult=istinfo(2)                                               1d2s20
         isymmrci=istinfo(1)                                             1d2s20
         nroot=istinfo(4)                                               1d2s20
         lambdaci=istinfo(5)                                             1d2s20
         write(6,*)('straight copy to isymmrci '),isymmrci
         nlzzci=nlzzq                                                   1d2s20
         do ii=1,6                                                       1d2s20
          nameci(ii)=istinfo(5+ii)                                       1d2s20
         end do                                                          1d2s20
         do ii=1,5                                                      9d9s21
          iip=ii                                                        9d20s21
          if(nameci(iip).gt.0)then                                      11d19s21
           symatoms(ii:ii)=char(nameci(iip))                             9d10s21
          else                                                          11d19s21
           symatoms(ii:ii)=' '                                          11d19s21
          end if                                                        11d19s21
         end do                                                         9d9s21
         do i=1,nsymb                                                    5d24s18
          idoubo(i)=ncoreg(i)                                            1d2s20
          ne=ne+idoubo(i)*2                                             1d2s20
          irefo(i)=nvalg(i)                                              1d2s20
         end do                                                          5d24s18
        else                                                            1d2s20
         nlzzci=nlzz                                                    1d2s20
         ne=isinfo(3,1)                                                 1d2s20
         ismult=isinfo(2,1)                                             1d2s20
         isymmrci=isinfo(1,1)                                           1d2s20
         nroot=isinfo(4,1)                                              1d2s20
         lambdaci=isinfo(5,1)                                           1d2s20
         do ii=1,6                                                      1d2s20
          nameci(ii)=isinfo(5+ii,1)                                     1d2s20
         end do                                                         1d2s20
         do ii=1,5                                                      9d9s21
          iip=ii                                                        9d20s21
          if(nameci(iip).gt.0)then                                      11d19s21
           symatoms(ii:ii)=char(nameci(iip))                             9d10s21
          else                                                          11d19s21
           symatoms(ii:ii)=' '                                          11d19s21
          end if                                                        11d19s21
         end do                                                         9d9s21
         do i=1,nsymb                                                    5d24s18
          idoubo(i)=idoub(i)                                            1d2s20
          ne=ne+idoubo(i)*2                                             1d2s20
          irefo(i)=iacto(i)                                             1d2s20
         end do                                                          5d24s18
         if(nlzz.eq.2)then                                              1d8s19
          if(lambdaci.gt.0)then                                         1d8s19
           if(mod(lambdaci,2).eq.1)then                                 1d8s19
            if(isymmrci.eq.2)then                                        1d8s19
             nxref=3                                                     1d8s19
             irxinfo(2)=1                                                1d8s19
             irxinfo(3)=1                                                1d8s19
            else                                                         1d8s19
             nxref=7                                                     1d8s19
             irxinfo(6)=1                                                1d8s19
             irxinfo(7)=1                                                1d8s19
            end if                                                       1d8s19
           else
            if(isymmrci.eq.1)then                                        1d8s19
             nxref=4                                                      1d8s19
             irxinfo(1)=1                                                1d8s19
             irxinfo(4)=1                                                1d8s19
            else
             nxref=8                                                     1d8s19
             irxinfo(5)=1                                                1d8s19
             irxinfo(8)=1                                                1d8s19
            end if                                                       1d8s19
           end if                                                       1d8s19
          end if                                                        1d8s19
         else if(nlzz.eq.6)then                                         1d8s19
          if(lambdaci.eq.1)then                                         1d8s19
               if(isymmrci.eq.2)then                                    1d8s19
                nxref=8                                                 1d8s19
                irxinfo(2)=1                                            1d8s19
                irxinfo(3)=1                                            1d8s19
                irxinfo(5)=1                                            5d13s21
               else                                                     1d8s19
                nxref=7                                                 1d8s19
                irxinfo(4)=1                                            1d8s19
                irxinfo(6)=1                                            1d8s19
                irxinfo(7)=1                                            1d8s19
               end if                                                   1d8s19
          else if(lambdaci.eq.2)then                                    1d8s19
               if(isymmrci.eq.2)then                                    1d8s19
                nxref=8                                                 1d8s19
                irxinfo(2)=1                                            1d8s19
                irxinfo(3)=1                                            1d8s19
                irxinfo(5)=1                                            5d13s21
                irxinfo(8)=2                                            5d13s21
               else                                                     1d8s19
                nxref=7                                                 1d8s19
                irxinfo(1)=2                                            1d8s19
                irxinfo(4)=1                                            1d8s19
                irxinfo(6)=1                                            1d8s19
                irxinfo(7)=1                                            1d8s19
               end if                                                   1d8s19
          else if(lambdaci.eq.3)then                                    1d8s19
          else if(lambdaci.eq.4)then                                    1d8s19
          end if                                                        1d8s19
         end if                                                         1d8s19
        end if                                                          11d20s19
         ii=0                                                           5d24s18
         do isb=1,nsymb                                                 5d24s18
          do i=1,irefo(isb)                                             5d24s18
           ii=ii+1                                                      5d24s18
           if(ii.gt.ido)then                                            5d24s18
            write(6,*)('no. orbs in reference space too large ')        5d24s18
            write(6,*)('limit is ido = '),ido                           5d24s18
            stop                                                        5d24s18
           end if                                                       5d24s18
           ism(ii)=isb                                                  5d24s18
           irel(ii)=i                                                   5d24s18
          end do                                                        5d24s18
         end do                                                         1d2s20
          norb=ii
          myguesci=1                                                    4d22s21
          iunc(1)=-1                                                    5d27s19
          ethred=ethreddef                                              12d21s22
          nethred=0                                                     10s3s22
        do i=1,3                                                        5d24s18
         nhole(i)=0                                                     5d24s18
        end do                                                          5d24s18
        do i=1,4                                                        4d19s23
         nfill(i)=0                                                     5d24s18
        end do                                                          5d24s18
        nkeepl=0                                                        5d24s18
        nkeeph=0                                                        5d24s18
        iuse2den=0                                                      9d30s20
        itestmrci=0                                                     5d15s19
        ndecomp=0                                                       8d10s22
        nreadc=0                                                        8d10s22
 3251   continue                                                        5d24s18
        ieof=1                                                          5d24s18
        is=1                                                            5d24s18
        read(5,200,end=1808)line200                                     8d3s21
        leof=.false.                                                    8d3s21
        go to 1809                                                      8d3s21
 1808   continue                                                        8d3s21
        leof=.true.                                                     8d3s21
 1809   continue                                                        8d3s21
        call delim(line200,is,ie)                                          5d24s18
 6602   continue                                                        5d10s19
  200   format(a200)                                                    7d6s18
        if(line200(is:is+6).eq.'maxdiis')then                           10d16s24
         is=ie+1                                                        10d16s24
         call delim(line200,is,ie)                                      10d16s24
         if(ie.ge.is)then                                               10d16s24
          read(line200(is:ie),*)maxdiisci                               10d16s24
          write(6,*)('maxdiis read in to be '),maxdiisci                10d16s24
         else                                                           10d16s24
          write(6,*)('no argument for maxdiis!!!')                      10d16s24
          write(6,*)('help!')                                           10d16s24
          irtrn=1                                                       10d16s24
          return                                                        10d16s24
         end if                                                         10d16s24
         go to 3251                                                     10d16s24
        end if                                                          10d16s24
        if(line200(is:is+6).eq.'maxrest')then                           10d16s24
         is=ie+1                                                        10d16s24
         call delim(line200,is,ie)                                      10d16s24
         if(ie.ge.is)then                                               10d16s24
          read(line200(is:ie),*)maxrestci                               10d16s24
          write(6,*)('maxrest read in to be '),maxrestci                10d16s24
         else                                                           10d16s24
          write(6,*)('no argument for maxrest!!!')                      10d16s24
          write(6,*)('help!')                                           10d16s24
          irtrn=1                                                       10d16s24
          return                                                        10d16s24
         end if                                                         10d16s24
         go to 3251                                                     10d16s24
        end if                                                          10d16s24
        if(line200(is:is+4).eq.'toldv')then                             7d22s24
         is=ie+1                                                        9d19s23
         call delim(line200,is,ie)                                      9d19s23
         if(ie.ge.is)then                                               9d19s23
          read(line200(is:ie),*)toldv                                   7d22s24
          write(6,*)('toldv read in to be '),toldv                      7d22s24
         else                                                           7d22s24
          write(6,*)('no argument for toldv!!!')                        7d22s24
          write(6,*)('help!')                                           7d22s24
          irtrn=1                                                       7d22s24
          return                                                        7d22s24
         end if                                                         7d22s24
        end if                                                          7d22s24
        if(line200(is:is+5).eq.'epssym')then                            9d19s23
         is=ie+1                                                        9d19s23
         call delim(line200,is,ie)                                      9d19s23
         if(ie.ge.is)then                                               9d19s23
          read(line200(is:ie),*)epssymci                                9d19s23
          write(6,*)('epssymci read in to be '),epssymci                9d19s23
         else                                                           9d19s23
          write(6,*)('no argument for epssymci!!!')                     9d19s23
          write(6,*)('help!')                                           9d19s23
          irtrn=1                                                       9d19s23
          return                                                        9d19s23
         end if                                                         9d19s23
         go to 3251                                                     9d19s23
        end if                                                          9d19s23
        if(line200(is:is).eq.'&'.or.line200(is:is+1).eq.'#&'            8d3s21
     $       .or.leof)then                                              8d3s21
         write(6,*)('end of mrci input ...')
         write(6,*)('what we have for dynwci and irxinfos:'),dynwci
         if(nrxinfos.eq.0)then                                          8d3s21
          is=1                                                          8d23s21
          call delim(symatoms,is,ie)                                    8d23s21
          nrxinfos0=nrxinfos                                            3d21s22
          call parseterm(symatoms(is:ie),nroot,irxinfos,nrxinfos,       8d23s21
     $         maxxref)                                                 8d23s21
          if(nrxinfos.eq.nrxinfos0)then                                 3d21s22
           nrxinfos=nrxinfos+1                                          3d21s22
           irxinfos(1,nrxinfos)=isymmrci                                3d21s22
           irxinfos(2,nrxinfos)=nroot                                   3d21s22
           irxinfos(3,nrxinfos)=0                                       3d21s22
          end if                                                        3d21s22
         end if                                                         8d3s21
         if(nlzzci.eq.0)then                                            2d2s22
          jjjtop=2                                                      2d2s22
          write(6,*)('   #  isb nroot')                                    2d2s22
         else if(nlzzci.eq.2)then                                       2d2s22
          jjjtop=4                                                      2d2s22
          write(6,*)('   #  isb nroot nlzz lam')                           2d2s22
         else if(nlzzci.eq.6)then                                       2d2s22
          jjjtop=5                                                      2d2s22
          write(6,*)('   #  isb nroot nlzz  lam  lz')                        2d2s22
         end if                                                         2d2s22
         do jj=1,nsymb                                                  12d29s21
          irxinfo(jj)=0d0                                               12d29s21
         end do                                                         12d29s21
         ixok=ibcoff                                                    2d8s22
         ibcoff=ixok+nrxinfos                                           2d8s22
         call enough('basisz. 11',bc,ibc)
         do jj=0,nrxinfos-1                                             2d8s22
          ibc(ixok+jj)=0                                                2d8s22
         end do                                                         2d8s22
         jxok=ixok-1                                                    2d8s22
         ndup=0                                                         2d8s22
         do jj=1,nrxinfos0                                              2d8s22
          jjdup=0                                                        2d1s22
          do kk=jj+1,nrxinfos                                            2d1s22
           idel=iabs(irxinfos(1,jj)-irxinfos(1,kk))                      2d1s22
           do jjj=3,jjjtop                                               2d2s22
            idel=idel+iabs(irxinfos(jjj,jj)-irxinfos(jjj,kk))            2d2s22
           end do                                                        2d2s22
           if(idel.eq.0)jjdup=1                                          2d1s22
          end do                                                         2d1s22
          ndup=ndup+jjdup                                               2d8s22
          ibc(jxok+jj)=jjdup                                            2d9s22
         end do                                                         2d8s22
         if(ndup.gt.0)then                                              2d8s22
 2027      format(6i5)
          kk=0                                                          2d8s22
          do jj=1,nrxinfos                                              2d8s22
           if(ibc(jxok+jj).eq.0)then                                    2d8s22
            kk=kk+1                                                     2d8s22
            do jjj=1,jjjtop                                             2d8s22
             irxinfos(jjj,kk)=irxinfos(jjj,jj)                          2d8s22
            end do                                                      2d8s22
           end if                                                       2d8s22
          end do                                                        2d8s22
          nrxinfos=kk                                                   2d8s22
         end if                                                         2d8s22
         do jj=1,nrxinfos                                               8d3s21
          write(6,2027)jj,(irxinfos(jjj,jj),jjj=1,jjjtop)                  2d2s22
          irxinfo(irxinfos(1,jj))=irxinfo(irxinfos(1,jj))               12d29s21
     $        +irxinfos(2,jj)                                           12d29s21
         end do                                                         8d3s21
         if(.not.leof)then                                              8d3s21
          line(1:10)=line200(1:10)                                       11d22s19
          go to 254                                                      5d25s18
         else                                                           8d3s21
          go to 1008                                                    8d3s21
         end if                                                         8d3s21
        end if                                                          5d25s18
        if(line200(is:is+5).eq.'ethred')then                            10d5s22
         is=ie+1                                                        10d5s22
         call delim(line200,is,ie)                                      10d5s22
         if(ie.lt.is)then                                               10d5s22
          write(6,*)('sorry, I can not find an argument for ethred ')   10d5s22
          write(6,*)('help!')                                           10d5s22
          irtrn=1                                                       10d5s22
          return                                                        10d5s22
         end if                                                         10d5s22
         read(line200(is:ie),*)ethred                                   10d5s22
         write(6,*)('value of ethred read in = '),ethred                10d5s22
         nethred=1                                                      10d5s22
        end if                                                          10d5s22
        if(line200(is:is+4).eq.'itabo')then                             12d9s22
         is=ie+1                                                        12d9s22
         call delim(line200,is,ie)                                      12d9s22
         if(ie.ge.is)then                                               12d9s22
          read(line200(is:ie),*)itabort                                 12d9s22
         else                                                           12d9s22
          write(6,*)('there is no argument to the itabort card!')       12d9s22
          write(6,*)('I don''t know what to do!')                       12d9s22
          irtrn=1                                                       12d9s22
          return                                                        12d9s22
         end if                                                         12d9s22
         if(itabort.gt.0)then                                           6d26s23
          write(6,*)('MRCI calculation will be aborted at iteration '),  12d9s22
     $        itabort                                                   12d9s22
         end if                                                         6d26s23
        end if                                                          12d9s22
        if(line200(is:is+5).eq.'skip4v')then                            3d24s23
         is=ie+1                                                        3d24s23
         call delim(line200,is,ie)                                      3d24s23
         iskip4v=1                                                      3d24s23
         if(ie.ge.is)then                                               3d24s23
          if(line200(is:is+1).eq.'no')then                              3d24s23
           iskip4v=0                                                    3d24s23
          end if                                                        3d24s23
         end if                                                         3d24s23
         if(iskip4v.ne.0)write(6,*)                                     3d24s23
     $        ('We will skip 4-virtual contributions to hdd')           3d24s23
        end if                                                          3d24s23
        if(line200(is:is+6).eq.'restart')then                           8d11s22
         is=ie+1                                                        8d24s22
         call delim(line200,is,ie)                                      8d24s22
         if(ie.ge.is)then                                               8d24s22
          if(line200(is:is).eq.'y')then                                 8d24s22
           nrestart=-nrestart                                           8d24s22
          end if                                                        8d24s22
         else                                                           8d24s22
          nrestart=-nrestart                                             8d16s22
         end if                                                         8d24s22
         if(nrestart.lt.0)write(6,*)                                    8d24s22
     $       ('we will restart ci calculation from the wavef file.')    8d16s22
        end if                                                          8d11s22
        if(line200(is:is+3).eq.'save')then                              8d16s22
         is=ie+1                                                        8d16s22
         call delim(line200,is,ie)                                      8d16s22
         if(ie.lt.is)then                                               8d16s22
          nrestart=maxds                                                8d16s22
         else                                                           8d16s22
          read(line200(is:ie),*)nrestart                                8d16s22
         end if                                                         8d16s22
         write(6,*)('restart information will be saved to wavef'),      8d16s22
     $        ('every '),nrestart,(' iteration.')                       8d16s22
        end if                                                          8d16s22
        if(line200(is:is+5).eq.'npdiag')then                            8d26s22
         is=ie+1                                                        8d26s22
         call delim(line200,is,ie)                                      8d31s22
         if(ie.ge.is)then                                                8d12s22
          read(line200(is:ie),*)npdiagci                                8d31s22
          write(6,*)                                                     8d12s22
     $        ('threshold for using parallel diagonalization is now'),  8d12s22
     $        npdiagci                                                  8d26s22
         else                                                            8d12s22
          write(6,*)('there is no argument for npdiag! ')                8d12s22
          write(6,*)('please supply one!')                               8d12s22
          irtrn=1                                                        8d12s22
          return                                                         8d12s22
         end if                                                         8d26s22
        end if                                                          8d12s22
        if(line200(is:is+5).eq.'itmaxi'.or.                             6d6s23
     $       line200(is:is+5).eq.'maxiti')then                          6d6s23
         is=ie+1                                                        3d21s22
         call delim(line200,is,ie)                                      3d21s22
         if(ie.ge.is)then                                               3d21s22
          read(line200(is:ie),*)maxmiti                                 10d14s22
          write(6,*)('Maximum number of micro CI iterations is now '),  3d21s22
     $         maxmiti                                                  10d14s22
         else                                                           3d21s22
          write(6,*)('I''m sorry, but I need a value on the maxiti'),   6d6s23
     $         (' line, but I don''t see any.')                         3d21s22
          irtn=1                                                        3d21s22
          return                                                        3d21s22
         end if                                                         3d21s22
        end if                                                          3d21s22
        if(line200(is:is+4).eq.'itmax'.or.                              6d6s23
     $       line200(is:is+4).eq.'maxit')then                           6d6s23
         is=ie+1                                                        3d21s22
         call delim(line200,is,ie)                                      3d21s22
         if(ie.ge.is)then                                               3d21s22
          read(line200(is:ie),*)maxmit                                  3d21s22
          write(6,*)('Maximum number of macro CI iterations is now '),  3d21s22
     $         maxmit                                                   3d21s22
         else                                                           3d21s22
          write(6,*)('I''m sorry, but I need a value on the maxit'),    6d6s23
     $         (' line, but I don''t see any.')                         3d21s22
          irtn=1                                                        3d21s22
          return                                                        3d21s22
         end if                                                         3d21s22
        end if                                                          3d21s22
        if(line200(is:is+5).eq.'maxdsi')then                            10d14s22
         is=ie+1                                                        12d29s21
         call delim(line200,is,ie)                                      1d3s20
         if(ie.ge.is)then                                               12d29s21
          read(line200(is:ie),*)maxdsi                                   12d29s21
          write(6,*)                                                    10d14s22
     $     ('maximum no. of Davidson Vectors per internal root is now'),12d29s21
     $         maxdsi                                                   10d14s22
         end if                                                         12d29s21
        end if                                                          12d29s21
        if(line200(is:is+4).eq.'maxds')then                             12d29s21
         is=ie+1                                                        12d29s21
         call delim(line200,is,ie)                                      1d3s20
         if(ie.ge.is)then                                               12d29s21
          read(line200(is:ie),*)maxds                                   12d29s21
          write(6,*)('maximum no. of Davidson Vectors per root is now'),12d29s21
     $         maxds                                                    12d29s21
         end if                                                         12d29s21
        end if                                                          12d29s21
        if(line200(is:is+4).eq.'debug')then                             1d3s20
 6603    continue                                                        1d3s20
         is=ie+1                                                        1d3s20
         call delim(line200,is,ie)                                      1d3s20
         if(ie.ge.is)then                                               1d3s20
          write(6,*)('debug print out turned on for '),line200(is:ie)    1d3s20
          nchere=ie+1-is                                                 1d3s20
          do i=1,nprtr                                                   1d3s20
           if(line200(is:ie).eq.prtrname(i)(1:nchere))then               1d3s20
            iprtr(i)=1                                                   1d3s20
            go to 6603                                                   1d3s20
           end if                                                        1d3s20
          end do                                                         1d3s20
          write(6,*)('no matching string in prtrname found ')           1d3s20
         end if                                                         1d3s20
        end if                                                          1d3s20
        if(line200(is:is+4).eq.'cache')then                             2d10s20
         is=ie+1                                                        2d10s20
         call delim(line200,is,ie)                                      2d10s20
         if(ie.lt.is)then                                               2d10s20
          write(6,*)('need to supply an no. on cache card ')            2d10s20
          irtrn=1                                                       2d10s20
          return                                                        2d10s20
         end if                                                         2d10s20
         read(line200(is:ie),*)mcache                                   2d10s20
         write(6,*)('cache size read in is '),mcache                    2d10s20
        end if                                                          2d10s20
        if(line200(is:is+3).eq.'norb')then                              1d24s22
         is=ie+1                                                        1d24s22
         call delim(line200,is,ie)                                      1d24s22
         if(ie.lt.is)then                                               1d24s22
          write(6,*)('help! there is no input to norb!')                1d24s22
          irtrn=1                                                       1d24s22
          return                                                        1d24s22
         end if                                                         1d24s22
         if(line200(is:is+1).eq.'av')then                               1d24s22
          write(6,*)                                                    1d24s22
     $        ('ci natural orbitals will be based on average density')  1d24s22
          norbci=-1                                                     1d24s22
         else                                                           1d24s22
          read(line200(is:ie),*)norbci                                  1d24s22
          write(6,*)('ci natural orbitals will be based on root '),     1d24s22
     $         norbci                                                   1d24s22
         end if                                                         1d24s22
         go to 3251                                                     9d26s20
        end if                                                          1d24s22
        if(line200(is:is+3).eq.'2den')then                              7d13s20
         is=ie+1                                                        7d13s20
         call delim(line200,is,ie)                                      7d13s20
         if(ie.lt.is)then                                               7d13s20
          write(6,*)('no argument to 2den supplied ')                   7d13s20
          irtrn=1                                                       7d13s20
          return                                                        7d13s20
         end if                                                         7d13s20
         if(line200(is:is+1).eq.'on')then                               7d13s20
          iuse2den=1                                                    7d13s20
         else if(line200(is:is+2).eq.'tes')then                         7d13s20
          iuse2den=2                                                    7d13s20
         else if(line200(is:is+2).eq.'off')then                         7d13s20
          iuse2den=0                                                    7d13s20
         else                                                           7d13s20
          write(6,*)('unknown argument to 2den: "'),line200(is:ie),('"')7d13s20
          irtrn=1                                                       7d13s20
          return                                                        7d13s20
         end if                                                         7d13s20
         go to 3251                                                     9d26s20
        end if                                                          7d13s20
        if(line200(is:is+7).eq.'novguess')then                          2d14s20
         write(6,*)                                                     2d14s20
     $       ('we will not use old vectors as guess in iterations')     2d14s20
         novguess=1                                                     2d14s20
        end if                                                          2d14s20
        if(line200(is:is+4).eq.'nowav')then                             3d23s21
         write(6,*)('wavef file will not be created')                   3d23s21
         nwavef=0                                                       3d23s21
         nrestart=0                                                     8d16s22
         go to 3251                                                     3d23s21
        end if                                                          3d23s21
        if(line200(is:is+4).eq.'wconv')then                             6d20s24
         is=ie+1                                                        1d22s20
         call delim(line200,is,ie)                                      1d22s20
         if(ie.lt.is)then                                               1d22s20
          write(6,*)('need to supply a no. on wconv card ')             6d20s24
          irtrn=1                                                       1d22s20
          return                                                        1d22s20
         end if                                                         1d22s20
         read(line200(is:ie),*)wconvci                                  6d20s24
         write(6,*)('wavefunction convergence in mrci is '),wconvci     6d20s24
        end if                                                          1d22s20
        if(line200(is:is+4).eq.'econv')then                             1d22s20
         is=ie+1                                                        1d22s20
         call delim(line200,is,ie)                                      1d22s20
         if(ie.lt.is)then                                               1d22s20
          write(6,*)('need to supply a no. on econv card ')             1d22s20
          irtrn=1                                                       1d22s20
          return                                                        1d22s20
         end if                                                         1d22s20
         read(line200(is:ie),*)econvci                                  1d22s20
         write(6,*)('energy convergence in mrci is '),econvci           1d22s20
        end if                                                          1d22s20
        if(line200(is:is+5).eq.'eiconv')then                            2d9s21
         is=ie+1                                                        1d22s20
         call delim(line200,is,ie)                                      1d22s20
         if(ie.lt.is)then                                               1d22s20
          write(6,*)('need to supply a no. on eiconv card ')            2d9s21
          irtrn=1                                                       1d22s20
          return                                                        1d22s20
         end if                                                         1d22s20
         read(line200(is:ie),*)econvcii                                  1d22s20
         write(6,*)('energy convergence in internal ci is '),econvcii   2d9s21
        end if                                                          1d22s20
        if(line200(is:is+3).eq.'cren')then                              1d22s20
         is=ie+1                                                        1d22s20
         call delim(line200,is,ie)                                      1d22s20
         if(ie.gt.is)then                                               1d22s20
          if(line200(is:is+1).eq.'of')then                              1d22s20
           lcrenorm=0                                                   1d22s20
           write(6,*)('N-2 contractions will not be renormalized '),    1d22s20
     $          ('before orthogonalization.')                           1d22s20
          else if(line200(is:is+1).eq.'on')then                              1d22s20
           lcrenorm=1                                                   1d22s20
           write(6,*)('N-2 contractions will be renormalized '),        10d11s20
     $          ('before orthogonalization.')                           1d22s20
          end if                                                        1d22s20
         end if                                                         1d22s20
        end if                                                          1d22s20
        if(line200(is:is+9).eq.'nointeract')then                        10d15s19
         write(6,*)('we found nointeract keywork ... ')
         write(6,*)('do not make first order interacting '),            10d15s19
     $        ('space restriction')                                     10d15s19
         interacts=0                                                    10d15s19
        end if                                                          10d15s19
        if(line200(is:is+7).eq.'maxopens')then                          6d12s19
         is=ie+1                                                        6d12s19
         call delim(line200,is,ie)                                      6d12s19
         if(ie.lt.is)then                                               6d12s19
          write(6,*)('need to supply a no. on maxopens card')           6d12s19
          irtrn=1                                                       1d22s20
          return                                                        1d22s20
         end if                                                         6d12s19
         read(line200(is:ie),*)maxopens                                 6d12s19
         write(6,*)('maximum number of open shells set to '),maxopens   6d12s19
        end if                                                          6d12s19
        if(line200(is:is+3).eq.'test')then                              5d15s19
         is=ie+1                                                        5d15s19
         call delim(line200,is,ie)                                      5d15s19
         if(ie.ge.is)then                                               5d15s19
          read(line200(is:ie),*)itestmrci                               5d15s19
          write(6,*)('the argument of test is '),itestmrci              5d15s19
          if(itestmrci.eq.1)then                                        5d15s19
           write(6,*)('step 1 of testing: dump uncontracted Ham matrix')5d15s19
           pthrsci=1d10                                                 5d15s19
           maxpsci=100000                                                5d15s19
          else if(itestmrci.eq.2)then                                   5d15s19
           write(6,*)('step 2 of testing: diagonalize uncontracted Ham')5d15s19
     $          ,(' matrix and dump det information to hdata')          5d15s19
          else if(itestmrci.eq.3)then                                   5d15s19
           write(6,*)('step 3 of testing: run contracted calculation'), 5d15s19
     $          (' and dump contraction data to rdata')                 5d15s19
           write(6,*)('next run htestd')                                5d15s19
           pthrsci=1d10                                                 5d15s19
           maxpsci=100000                                                5d15s19
          else                                                          5d15s19
           write(6,*)('unknown argument to test')                       5d15s19
           irtrn=1                                                      5d15s19
           return                                                       5d15s19
          end if                                                        5d15s19
         else                                                           5d15s19
          write(6,*)('you have given the "test" keyword without an'),   5d15s19
     $         (' argument. use 1,2, or 3 please.')                     5d15s19
          irtrn=1                                                       5d15s19
          return                                                        5d15s19
         end if                                                         5d15s19
        end if                                                          5d15s19
        if(line200(is:is+4).eq.'maxps')then                             4d15s19
         is=ie+1                                                        4d15s19
         call delim(line200,is,ie)                                      4d15s19
         if(itestmrci.ne.1.and.itestmrci.ne.3)then                      9d6s19
          read(line200(is:ie),*)maxpsci                                  4d15s19
          write(6,*)('maximum size of p-space is now '),maxpsci          4d15s19
         end if                                                         5d15s19
        end if                                                          4d15s19
        if(line200(is:is+2).eq.'dyn')then                               4d28s21
         is=ie+1                                                        4d28s21
         call delim(line200,is,ie)                                      4d28s21
         read(line200(is:ie),*)dynwci(1)                                9d20s21
         write(6,*)('dynwci is read in to be '),dynwci(1)               9d20s21
         is=ie+1                                                        9d20s21
         call delim(line200,is,ie)                                      9d20s21
         if(ie.gt.is)then                                               9d20s21
          read(line200(is:ie),*)dynwci(2)                               9d20s21
          write(6,*)                                                    9d20s21
     $     ('we will be using 4/(cc+exp(2*wde)) as weighting function') 9d20s21
          write(6,*)('with absolute energy '),dynwci(2)                 9d20s21
          is=ie+1                                                       9d20s21
          call delim(line200,is,ie)                                      9d20s21
          if(ie.gt.is)then                                               9d20s21
           read(line200(is:ie),*)dynwci(3)                               9d20s21
           write(6,*)('with cc read in to be '),dynwci(3)               9d20s21
          else                                                          9d20s21
           write(6,*)('with cc default value of '),dynwci(3)            9d20s21
          end if                                                        9d20s21
         end if                                                         9d20s21
        end if                                                          4d28s21
        if(line200(is:is+3).eq.'xref')then                              5d10s19
         write(6,*)('we will have extra doubles reference fcns')        8d3s21
         write(6,*)('what we have for nrxinfos so far '),               2d8s22
     $        nrxinfos                                                  2d8s22
         is=ie+1                                                        5d10s19
         nxref=0                                                        2d11s20
c
c     make sure first is mrci function                                  4d28s21
c
         if(nroot.lt.0)then                                             4d28s21
          write(6,*)('oops!!! xref line must follow state line!')       4d28s21
          iretrn=1                                                      4d28s21
          return                                                        4d28s21
         end if                                                         4d28s21
         if(nrxinfos.eq.0)then                                          4d28s21
          write(6,*)(' stuffing in mrci function ')
          write(6,*)('what we have from symatoms: '),symatoms
          write(6,*)('parseterm the second')
          ieee=0                                                        12d6s22
          do izz=5,1,-1                                                 12d6s22
           if(symatoms(izz:izz).ne.' ')then                             12d6s22
            ieee=izz                                                    12d6s22
            go to 245                                                   12d6s22
           end if                                                       12d6s22
          end do                                                        12d6s22
  245     continue                                                      12d6s22
          write(6,*)('what we have for ieee'),ieee                      12d6s22
          if(ieee.gt.0)then                                             12d6s22
           call parseterm(symatoms(1:ieee),nroot,irxinfos,nrxinfos,     12d6s22
     $         maxxref)                                                 12d6s22
          else                                                          12d6s22
           irtrn=1                                                      12d6s22
           return                                                       12d6s22
          end if                                                        12d6s22
          nrxinfos0=nrxinfos                                            2d8s22
         end if                                                         4d28s21
 6601    continue                                                       5d10s19
          is=ie+1                                                       5d23s19
c
c     input consists of a term symbol (less spin multiplicity)          4d24s21
c     followed by an optional number of roots.                          4d24s21
c
          call delim(line200,is,ie)                                     5d10s19
          if(ie.ge.is)then                                              4d24s21
           write(6,*)('parse term for '),line200(is:ie)
           isn=ie+1                                                     4d24s21
           call delim(line200,isn,ien)                                     5d10s19
           nrooter=1                                                    4d24s21
           isnext=ie+1                                                  4d24s21
           if(ien.ge.isn)then                                           4d24s21
            do i=1,10                                                   4d24s21
             if(line200(isn:isn).eq.digit(i))then                       4d24s21
              read(line200(isn:ien),*)nrooter                           4d24s21
              isnext=ien+1                                              4d24s21
              go to 88                                                  4d24s21
             end if                                                     4d24s21
            end do                                                      4d24s21
           end if                                                       4d24s21
   88      continue                                                     4d24s21
          write(6,*)('parseterm the third')
           call parseterm(line200(is:ie),nrooter,irxinfos,nrxinfos,     4d24s21
     $          maxxref)                                                4d24s21
           ie=isnext-1                                                  4d24s21
           go to 6601                                                   4d24s21
          end if                                                        4d24s21
          write(6,*)('what we have for dynwci and irxinfos:'),dynwci
          write(6,*)(' # isb nroot nlzz lam')
          do jj=1,nrxinfos                                              4d24s21
           write(6,*)jj,(irxinfos(jjj,jj),jjj=1,4)                       4d24s21
          end do                                                        4d24s21
          go to 6602                                                    5d10s19
        end if                                                          5d10s19
        if(line200(is:is+4).eq.'pthrs')then                                5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)                                         5d24s18
         if(itestmrci.ne.1.and.itestmrci.ne.3)then                      5d15s19
          read(line200(is:ie),*)pthrsci                                     6d6s18
         end if                                                         5d15s19
        end if                                                          5d24s18
        if(line200(is:is+4).eq.'doubo')then                                 5d24s18
         write(6,*)('orbitals doubly occupied in all dets')             5d24s18
         do isb=1,nsymb                                                 10d28s20
          idoubo(isb)=0                                                 10d28s20
         end do                                                         10d28s20
         do isb=1,nsymb                                                 5d24s18
          is=ie+1                                                        5d24s18
          call delim(line200,is,ie)                                         5d24s18
          if(ie.ge.is)then                                              5d24s18
           read(line200(is:ie),*)idoubo(isb)                                      5d24s18
          else                                                          5d24s18
           go to 3251                                                   5d24s18
          end if
         end do                                                         5d24s18
        end if                                                          5d24s18
        if(line200(is:is+2).eq.'ref')then                                  5d24s18
         write(6,*)('orbitals occupied in reference CSFs')              8d16s22
         do isb=1,nsymb                                                 10d28s20
          irefo(isb)=0                                                  10d28s20
         end do                                                         10d28s20
         do isb=1,nsymb                                                 5d24s18
          is=ie+1                                                        5d24s18
          call delim(line200,is,ie)                                         5d24s18
          if(ie.ge.is)then                                              5d24s18
           read(line200(is:ie),*)irefo(isb)                                5d24s18
          else                                                          5d24s18
           go to 3252                                                   5d24s18
          end if
         end do                                                         5d24s18
 3252    continue                                                       5d24s18
         ii=0                                                           5d24s18
         write(6,*)('irel, ism: ')
         do isb=1,nsymb                                                 5d24s18
          do i=1,irefo(isb)                                             5d24s18
           ii=ii+1                                                      5d24s18
           if(ii.gt.ido)then                                            5d24s18
            write(6,*)('no. orbs in reference space too large ')        5d24s18
            write(6,*)('limit is ido = '),ido                           5d24s18
            stop                                                        5d24s18
           end if                                                       5d24s18
           ism(ii)=isb                                                  5d24s18
           irel(ii)=i                                                   5d24s18
           write(6,*)ii,i,isb
          end do                                                        5d24s18
         end do                                                         5d24s18
         norb=ii                                                        5d24s18
         if(itestmrci.eq.3)then                                         5d15s19
          open(unit=7,file='rdata')                                     5d15s19
          write(7,*)nsymb                                               8d23s19
          write(7,*)((multh(j,i),j=1,nsymb),i=1,nsymb)                  8d23s19
          write(7,*)norb                                                5d15s19
          do i=1,norb                                                   5d15s19
           write(7,*)i,irel(i),ism(i)                                   5d15s19
          end do                                                        5d15s19
         end if                                                         5d15s19
        end if                                                          5d24s18
        if(line200(is:is+4).eq.'0hole')then                                5d24s18
         write(6,*)('reading in orbs doubly occupied in all reference ')5d24s18
     $       ,('dets, but allow 1 hole in N-1 and N-2 space')           5d24s18
c
c     orbs are specified as isj where i is orbital no. and j is         5d24s18
c     symmetry block, where i points to irefo orbital, ie i=1 is first  5d24s18
c     orbital in reference space of symmetry j.                         5d24s18
c
 4251    continue                                                       5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)
         if(ie.lt.is)go to 3251                                         5d24s18
         write(6,*)('try and extra 0hole orb from '),line200(is:ie)
         if(line200(ie-1:ie-1).ne.'s')then                                 5d24s18
          write(6,*)('error in 0hole input: '),line200(is:ie)              5d24s18
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         read(line200(ie:ie),*)isb
         read(line200(is:ie-2),*)io                                        5d24s18
         if(io.gt.irefo(isb).or.isb.gt.nsymb)then                       8d27s22
          write(6,*)
     $       ('sorry, but the orbital you listed in the 0hole input "'),8d19s22
     $        line200(is:ie),('" does not exist!')                      8d19s22
          write(6,*)('double check your input ')                        8d19s22
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         8d19s22
         if(norb.eq.0)then                                              5d24s18
          write(6,*)('sorry, need to have ref card before 0hole card')  5d24s18
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         do i=1,norb                                                    5d24s18
          if(ism(i).eq.isb.and.irel(i).eq.io)then                       5d24s18
           nhole(1)=nhole(1)+1                                           5d24s18
           ihole(nhole(1),1)=i                                           5d24s18
           go to 4251                                                    5d24s18
          end if                                                        5d24s18
         end do                                                         5d24s18
        end if                                                          5d24s18
        if(line200(is:is+4).eq.'1hole')then                                5d24s18
         write(6,*)('reading in orbs having at most 1 hole in'),        5d24s18
     $       (' all reference dets')                                    5d24s18
c
c     orbs are specified as isj where i is orbital no. and j is         5d24s18
c     symmetry block, where i points to irefo orbital, ie i=1 is first  5d24s18
c     orbital in reference space of symmetry j.                         5d24s18
c
 4252    continue                                                       5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)
         if(ie.lt.is)go to 3251                                         5d24s18
         write(6,*)('try and extra 1hole orb from '),line200(is:ie)
         if(line200(ie-1:ie-1).ne.'s')then                                 5d24s18
          write(6,*)('error in 1hole input: '),line200(is:ie)              5d24s18
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         read(line200(ie:ie),*)isb
         read(line200(is:ie-2),*)io                                        5d24s18
         if(io.gt.irefo(isb).or.isb.gt.nsymb)then                       8d27s22
          write(6,*)
     $       ('sorry, but the orbital you listed in the 1hole input "'),8d19s22
     $        line200(is:ie),('" does not exist!')                      8d19s22
          write(6,*)('double check your input ')                        8d19s22
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         8d19s22
         if(norb.eq.0)then                                              5d24s18
          write(6,*)('sorry, need to have ref card before 1hole card')  5d24s18
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         do i=1,norb                                                    5d24s18
          if(ism(i).eq.isb.and.irel(i).eq.io)then                       5d24s18
           nhole(2)=nhole(2)+1                                           5d24s18
           ihole(nhole(2),2)=i                                           5d24s18
           go to 4252                                                    5d24s18
          end if                                                        5d24s18
         end do                                                         5d24s18
        end if                                                          5d24s18
        if(line200(is:is+4).eq.'2hole')then                                5d24s18
         write(6,*)('reading in orbs having at most 2 holes in'),       5d24s18
     $       (' all reference dets')                                    5d24s18
c
c     orbs are specified as isj where i is orbital no. and j is         5d24s18
c     symmetry block, where i points to irefo orbital, ie i=1 is first  5d24s18
c     orbital in reference space of symmetry j.                         5d24s18
c
 4253    continue                                                       5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)
         if(ie.lt.is)go to 3251                                         5d24s18
         write(6,*)('try and extract 2hole orb from '),line200(is:ie)
         if(line200(ie-1:ie-1).ne.'s')then                                 5d24s18
          write(6,*)('error in 2hole input: '),line200(is:ie)              5d24s18
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         read(line200(ie:ie),*)isb
         read(line200(is:ie-2),*)io                                        5d24s18
         if(io.gt.irefo(isb).or.isb.gt.nsymb)then                       8d27s22
          write(6,*)
     $       ('sorry, but the orbital you listed in the 2hole input "'),8d19s22
     $        line200(is:ie),('" does not exist!')                      8d19s22
          write(6,*)('double check your input ')                        8d19s22
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         8d19s22
         if(norb.eq.0)then                                              5d24s18
          write(6,*)('sorry, need to have ref card before 2hole card')  5d24s18
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         do i=1,norb                                                    5d24s18
          if(ism(i).eq.isb.and.irel(i).eq.io)then                       5d24s18
           nhole(3)=nhole(3)+1                                           5d24s18
           ihole(nhole(3),3)=i                                           5d24s18
           go to 4253                                                    5d24s18
          end if                                                        5d24s18
         end do                                                         5d24s18
        end if                                                          5d24s18
        if(line200(is+1:is+4).eq.'fill')then                            4d19s23
         read(line200(is:is),*)ifl                                      4d19s23
         write(6,*)('reading in orbs having at most'),ifl,              4d19s23
     $        ('e- in all reference dets')                              4d19s23
         if(ifl.lt.1.or.ifl.gt.4)then                                   4d19s23
          write(6,*)('bad number of electrons!! ')                      4d19s23
          irtrn=1                                                       4d19s23
          return                                                        4d19s23
         end if                                                         4d19s23
c
c     orbs are specified as isj where i is orbital no. and j is         5d24s18
c     symmetry block, where i points to irefo orbital, ie i=1 is first  5d24s18
c     orbital in reference space of symmetry j.                         5d24s18
c
 4254    continue                                                       5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)
         if(ie.lt.is)go to 3251                                         5d24s18
         write(6,*)('try and extract'),ifl,('fill orb from '),          4d19s23
     $        line200(is:ie)                                            4d19s23
         if(line200(is:ie).eq.'\')then                                  10d12s23
          write(6,*)('continuation encountered ...')
          read(5,200,end=1808)line200                                     8d3s21
          ie=0                                                          10d12s23
          go to 4254                                                    10d12s23
         end if                                                         10d12s23
         if(line200(ie-1:ie-1).ne.'s')then                              4d19s23
          write(6,*)('error in'),ifl,('fill input: '),line200(is:ie)    4d19s23
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         read(line200(ie:ie),*)isb
         read(line200(is:ie-2),*)io                                        5d24s18
         if(io.gt.irefo(isb).or.isb.gt.nsymb)then                       8d27s22
          write(6,*)
     $       ('sorry, but the orbital you listed in the'),ifl,          4d19s23
     $        ('fill input "'),line200(is:ie),('" does not exist!')     4d19s23
          write(6,*)('double check your input ')                        8d19s22
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         8d19s22
         if(norb.eq.0)then                                              5d24s18
          write(6,*)('sorry, need to have ref card before fill card')   4d19s23
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
         end if                                                         5d24s18
         do i=1,norb                                                    5d24s18
          if(ism(i).eq.isb.and.irel(i).eq.io)then                       5d24s18
           nfill(ifl)=nfill(ifl)+1                                      4d19s23
           ifill(nfill(ifl),ifl)=i                                      4d19s23
           go to 4254                                                    5d24s18
          end if
         end do                                                         5d24s18
         write(6,*)('bad orbital name: '),line200(is:ie)                   7d2s18
          irtrn=1                                                       8d19s22
          return                                                        8d19s22
        end if                                                          5d24s18
        if(line200(is:is+4).eq.'state')then                                5d24s18
c                                                                       5d24s18
c     ichrg: charge. If neutral, can be omitted. If cation, must be     11d18s19
c     prefaced by + sign. See cas input.                                11d18s19
c     termsym: term symbol specifying spin multiplicity and symmetry    12d24s19
c     nroot: number of roots (can be omitted if nroot=1)                12d24s19
c     dynwconv: dynamic weight to alter convergence requirements for    10d28s20
c               excited states. ignored unless nroot gt 1.              10d28s20
c                                                                       5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)                                      12d13s22
         if(ie.lt.is)then                                               12d13s22
          read(5,200,end=1808)line200                                   12d13s22
          is=1                                                          12d13s22
          call delim(line200,is,ie)                                     12d13s22
         end if                                                         12d13s22
         if(line200(is:is).eq.'+')then                                  11d18s19
          if(ie.eq.is)then                                              11d18s19
           ichrg=1                                                      11d18s19
          else                                                          11d18s19
           read(line200(is+1:ie),*)ichrg                                11d18s19
          end if                                                        11d18s19
          is=ie+1                                                       11d18s19
         else if(line200(is:is).eq.'-')then                             11d18s19
          if(ie.eq.is)then                                              11d18s19
           ichrg=-1                                                     11d18s19
          else                                                          11d18s19
           read(line200(is+1:ie),*)ichrg                                11d18s19
           ichrg=-ichrg                                                 11d18s19
          end if                                                        11d18s19
          is=ie+1                                                       11d18s19
         else                                                           11d18s19
          ichrg=0                                                       11d18s19
         end if                                                         11d18s19
         ne=numeminus-ichrg                                             11d18s19
         call delim(line200,is,ie)                                         12d24s19
         write(6,*)('getting term symbol from "'),line200(is:ie),('"')     12d24s19
            do ilook=is,ie                                              12d24s19
             itry=ichar(line200(ilook:ilook))                              12d24s19
             if(itry.lt.idiglow.or.itry.gt.idighig)then                 12d24s19
              ilm=ilook-1                                               12d24s19
              read(line200(is:ilm),*)ismult                                12d24s19
              is=ilook                                                  12d24s19
              go to 1037                                                12d24s19
             end if                                                     12d24s19
            end do                                                      12d24s19
            write(6,*)('could not find spin multiplicity in "'),        12d24s19
     $           line200(is:ie),('"')                                      12d24s19
            stop                                                        12d24s19
 1037       continue                                                    12d24s19
            is=ilook                                                    3d16s21
            call delim(line200,is,ie)                                   3d16s21
            nhere=ie+1-is                                               3d16s21
            if(nhere.gt.0)then                                          3d16s21
             symatoms='     '                                            3d16s21
             symatoms(1:nhere)=line200(is:ie)                           3d16s21
            else                                                        3d16s21
             write(6,*)('error!!!!!')
             irtn=1                                                     3d16s21
             return                                                     3d16s21
            end if                                                      3d16s21
            do i=1,6                                                    12d24s19
             nameci(i)=0                                                12d24s19
            end do                                                      12d24s19
            do i=is,ie                                                  12d24s19
             im=i+1-is                                                  12d24s19
             nameci(im)=ichar(line200(i:i))                             12d24s19
            end do                                                      12d24s19
            if(ie.eq.is.or.line200(is+1:is+1).eq.'o')then               1d2s20
             nlzzci=6                                                     12d31s19
             if(line200(is:is).eq.'S')then                              1d2s20
              lambdaci=0                                                1d2s20
              if(is.eq.ie)then                                          12d31s19
               isymmrci=1                                               1d2s20
              else                                                      12d31s19
               isymmrci=8                                               1d2s20
              end if                                                    12d31s19
             else                                                       12d31s19
              if(is.eq.ie)then                                          12d31s19
               isymmrci=6                                               1d2s20
              else                                                      12d31s19
               isymmrci=2                                               1d2s20
              end if                                                    12d31s19
              write(6,*)('setting e/o isymmrci to 2 or 6. '),
     $             ('it is: '),isymmrci
              if(line200(is:is).eq.'P')then                             1d2s20
               lambdaci=1                                               1d2s20
               if(isymmrci.eq.2)then                                    1d8s19
                isymmrci=5                                              5d14s21
                nxref=8                                                 1d8s19
                irxinfo(2)=1                                            1d8s19
                irxinfo(3)=1                                            1d8s19
                irxinfo(5)=1                                            3d31s20
               else                                                     1d8s19
                nxref=7                                                 1d8s19
                isymmrci=4                                              5d14s21
                irxinfo(4)=1                                            1d8s19
                irxinfo(6)=1                                            1d8s19
                irxinfo(7)=1                                            1d8s19
               end if                                                   1d8s19
               write(6,*)('for P, use syms '),irxinfo
               write(6,*)('reset isymmrci to '),isymmrci
              else if(line200(is:is).eq.'D')then                        1d2s20
               lambdaci=2                                               1d2s20
               if(isymmrci.eq.2)then                                    1d8s19
                isymmrci=8                                              5d14s21
                nxref=8                                                 1d8s19
                irxinfo(2)=1                                            1d8s19
                irxinfo(3)=1                                            1d8s19
                irxinfo(5)=1                                            5d13s21
                irxinfo(8)=2                                            5d13s21
               else                                                     1d8s19
                isymmrci=1                                              5d14s21
                nxref=7                                                 1d8s19
                irxinfo(1)=2                                            1d8s19
                irxinfo(4)=1                                            1d8s19
                irxinfo(6)=1                                            1d8s19
                irxinfo(7)=1                                            1d8s19
               end if                                                   1d8s19
               write(6,*)('for D, use syms '),irxinfo
               write(6,*)('reset isymmrci to '),isymmrci
              else if(line200(is:is).eq.'F')then                        1d2s20
               if(isymmrci.eq.2)then                                    5d14s21
                isymmrci=5                                              5d14s21
                nxref=8                                                 1d8s19
                irxinfo(2)=2                                            1d8s19
                irxinfo(3)=2                                            1d8s19
                irxinfo(5)=2                                            3d31s20
                irxinfo(8)=1                                            5d14s21
               else                                                     5d14s21
                isymmrci=4                                              5d14s21
                nxref=8                                                 1d8s19
                irxinfo(6)=2                                            1d8s19
                irxinfo(7)=2                                            1d8s19
                irxinfo(1)=1                                            3d31s20
                irxinfo(4)=2                                            5d14s21
               end if                                                   5d14s21
               lambdaci=3                                               1d2s20
              else if(line200(is:is).eq.'G')then                        1d2s20
               lambdaci=4                                               1d2s20
               if(isymmrci.eq.2)then                                    5d14s21
                isymmrci=8                                              5d14s21
                nxref=8                                                 1d8s19
                irxinfo(2)=2                                            1d8s19
                irxinfo(3)=2                                            1d8s19
                irxinfo(5)=2                                            3d31s20
                irxinfo(8)=3                                            5d14s21
               else                                                     5d14s21
                isymmrci=1                                              5d14s21
                nxref=8                                                 1d8s19
                irxinfo(6)=2                                            1d8s19
                irxinfo(7)=2                                            1d8s19
                irxinfo(1)=3                                            3d31s20
                irxinfo(4)=2                                            5d14s21
               end if                                                   5d14s21
              else                                                      12d31s19
               write(6,*)('unknown atomic symmetry: '),line200(is:is)   1d2s20
               irtrn=1                                                  12d31s19
               return                                                   12d31s19
              end if                                                    12d31s19
             end if                                                     12d31s19
            else if(line200(is:is+2).eq.'Sig')then                      1d2s20
             if(nsymb.eq.8.and.                                         10d3s21
     $        .not.(line200(ie:ie).eq.'g'.or.line200(ie:ie).eq.'u'))then10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line200(is:ie),('" but can not find it. Help!')       10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             iem=ie-1                                                   3d3s21
             nlzzci=2                                                   12d24s19
             lambdaci=0                                                 12d24s19
             if(line200(ie:ie).eq.'+'.or.line200(iem:iem).eq.'+')then   3d3s21
              isymmrci=1                                                12d24s19
             else                                                       12d24s19
              isymmrci=4                                                12d24s19
             end if                                                     12d24s19
             if(line200(ie:ie).eq.'u')isymmrci=isymmrci+4               3d3s21
            else if(line200(is:is+1).eq.'Pi')then                       12d24s19
             if(nsymb.eq.8.and.                                         10d3s21
     $        .not.(line200(ie:ie).eq.'g'.or.line200(ie:ie).eq.'u'))then10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line200(is:ie),('" but can not find it. Help!')       10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzzci=2                                                   12d24s19
             lambdaci=1                                                 12d24s19
             if(is+1.eq.ie)then                                         12d24s19
              isymmrci=2                                                12d24s19
              nxref=3
              irxinfo(2)=1                                              1d8s19
              irxinfo(3)=1                                              1d8s19
             else if(line200(is+2:is+2).eq.'u')then                     12d24s19
              isymmrci=2                                                12d24s19
              nxref=3                                                   1d8s19
              irxinfo(2)=1                                              1d8s19
              irxinfo(3)=1                                              1d8s19
             else                                                       12d24s19
              nxref=7                                                   1d8s19
              irxinfo(6)=1                                              1d8s19
              irxinfo(7)=1                                              1d8s19
              isymmrci=6                                                12d24s19
             end if                                                     12d24s19
            else if(line200(is:is+2).eq.'Del')then                      12d24s19
             if(nsymb.eq.8.and.                                         10d3s21
     $        .not.(line200(ie:ie).eq.'g'.or.line200(ie:ie).eq.'u'))then10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line200(is:ie),('" but can not find it. Help!')       10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzzci=2                                                     12d24s19
             lambdaci=2                                                 12d24s19
             if(is+2.eq.ie.or.(line200(ie:ie).ne.'g'.and.               3d3s21
     $            line200(ie:ie).ne.'u'))then                           3d3s21
              isymmrci=1                                                12d24s19
              nxref=4                                                   1d8s19
              irxinfo(1)=1                                              1d8s19
              irxinfo(4)=1                                              1d8s19
             else if(line200(ie:ie).eq.'g')then                         3d3s21
              nxref=4                                                   1d8s19
              isymmrci=1                                                12d24s19
              irxinfo(1)=1                                              1d8s19
              irxinfo(4)=1                                              1d8s19
             else                                                       12d24s19
              isymmrci=5                                                12d24s19
              ixref=8                                                   1d8s19
              irxinfo(5)=1                                              1d8s19
              irxinfo(8)=1                                              1d8s19
             end if                                                     12d24s19
            else if(line200(is:is+2).eq.'Phi')then                      12d24s19
             if(nsymb.eq.8.and.                                         10d3s21
     $        .not.(line200(ie:ie).eq.'g'.or.line200(ie:ie).eq.'u'))then10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line200(is:ie),('" but can not find it. Help!')       10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzzci=2                                                   12d24s19
             lambdaci=3                                                 12d24s19
             if(is+2.eq.ie)then                                         12d24s19
              isymmrci=2                                                12d24s19
              nxref=3                                                   1d8s19
              irxinfo(2)=1                                              1d8s19
              irxinfo(3)=1                                              1d8s19
             else if(line200(is+3:is+3).eq.'u')then                     12d24s19
              isymmrci=2                                                12d24s19
              nxref=3                                                   1d8s19
              irxinfo(2)=1                                              1d8s19
              irxinfo(3)=1                                              1d8s19
             else                                                       12d24s19
              isymmrci=6
              nxref=7                                                   1d8s19
              irxinfo(6)=1                                              1d8s19
              irxinfo(7)=1                                              1d8s19
             end if                                                     12d24s19
            else if(line200(is:is+2).eq.'Gam')then                      12d24s19
             if(nsymb.eq.8.and.                                         10d3s21
     $        .not.(line200(ie:ie).eq.'g'.or.line200(ie:ie).eq.'u'))then10d3s21
              write(6,*)('looking for g or u at end of "'),             10d3s21
     $            line200(is:ie),('" but can not find it. Help!')       10d3s21
              irtrn=1                                                   10d3s21
              return                                                    10d3s21
             end if                                                     10d3s21
             nlzzci=2                                                   12d24s19
             lambdaci=4                                                 12d24s19
             if(is+2.eq.ie.or.(line200(ie:ie).ne.'g'.and.               3d3s21
     $            line200(ie:ie).ne.'u'))then                           3d3s21
              isymmrci=1                                                12d24s19
              nxref=4                                                   1d8s19
              irxinfo(1)=1                                              1d8s19
              irxinfo(4)=1                                              1d8s19
             else if(line200(ie:ie).eq.'g')then                         3d3s21
              isymmrci=1                                                12d24s19
              nxref=4                                                   1d8s19
              irxinfo(1)=1                                              1d8s19
              irxinfo(4)=1                                              1d8s19
             else                                                       12d24s19
              isymmrci=5                                                12d24s19
              nxref=8                                                   1d8s19
              irxinfo(5)=1                                              1d8s19
              irxinfo(8)=1                                              1d8s19
             end if                                                     12d24s19
            else if(line200(is:is).eq.'S')then                          12d24s19
             nlzzci=0                                                   1d5s20
             isp=is+1                                                   12d24s19
             read(line200(isp:ie),*)isymmrci                               12d24s19
            else                                                        12d24s19
             nlzzci=0                                                   1d5s20
             if(ngrp.eq.1)then                                          12d24s19
              issym=1                                                   12d24s19
             else if(ngrp.eq.2)then                                     12d24s19
              if(line200(is+2:is+2).eq.''''.or.line200(is+1:is+1).eq.'"'12d24s19
     $             .or.line200(is+1:is+1).eq.'2')then                   12d24s19
               isymmrci=2                                               12d24s19
              else                                                      12d24s19
               isymmrci=1                                               12d24s19
              end if                                                    12d24s19
             else                                                       12d24s19
              write(6,*)('ngrp, nsymb '),ngrp,nsymb
              if(ngrp.lt.0.or.ngrp.gt.8)stop
              do isb=1,ngrp                                              12d24s19
               write(6,*)('compare "'),symatoms(1:3),('" with "'),
     $              stype(isb,ngrp),('"')
               if(symatoms(1:3).eq.stype(isb,ngrp))then                 3d16s21
                isymmrci=isb                                            12d24s19
                go to 2037                                              12d24s19
               end if                                                   12d24s19
              end do                                                     12d24s19
              write(6,*)('bad term symbol: '),line200(is:ie)            12d24s19
              write(6,*)('we were searching from the list ')            1d14s20
              do isb=1,ngrp                                             1d14s20
               write(6,*)stype(isb,ngrp)                                1d14s20
              end do                                                    1d14s20
              stop                                                      12d24s19
             end if                                                     12d24s19
            end if                                                      12d24s19
 2037       continue                                                    12d24s19
         is=ie+1                                                        12d24s19
         call delim(line200,is,ie)                                      12d24s19
         nroot=1                                                        4d4s22
         if(ie.ge.is)then                                               12d24s19
          ndot=0                                                        9d20s21
          do ii=is,ie                                                   9d20s21
           if(line200(ii:ii).eq.'.')ndot=ndot+1                         9d20s21
          end do                                                        9d20s21
          if(ndot.eq.0)then                                             9d20s21
           read(line200(is:ie),*)nroot                                   12d24s19
           is=ie+1                                                       10d28s20
           call delim(line200,is,ie)                                     10d28s20
          end if                                                        9d20s21
          if(ie.ge.is)then                                              10d28s20
           read(line200(is:ie),*)dynwconv(1)                               10d28s20
           write(6,*)('dynamic weight for convergence criterion is now')10d28s20
     $          ,dynwconv(1)                                               10d28s20
           is=ie+1                                                        9d20s21
           call delim(line200,is,ie)                                      9d20s21
           if(ie.gt.is)then                                               9d20s21
            read(line200(is:ie),*)dynwconv(2)                               9d20s21
            write(6,*)                                                    9d20s21
     $      ('we will be using 4/(cc+exp(2*wde)) as weighting function') 9d20s21
            write(6,*)('with absolute energy '),dynwconv(2)                 9d20s21
            is=ie+1                                                       9d20s21
            call delim(line200,is,ie)                                      9d20s21
            if(ie.gt.is)then                                               9d20s21
             read(line200(is:ie),*)dynwconv(3)                               9d20s21
             write(6,*)('with cc read in to be '),dynwconv(3)               9d20s21
            else                                                          9d20s21
             write(6,*)('with cc default value of '),dynwconv(3)            9d20s21
            end if                                                        9d20s21
           end if                                                         9d20s21
          end if                                                        10d28s20
         end if                                                         12d24s19
        end if                                                          5d24s18
        if(line200(is:is+3).eq.'sing')then                                 5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)                                         5d24s18
          if(ie.gt.is.and.line200(is:is+2).eq.'off')then                1d3s20
           iunc(1)=0                                                    1d3s20
           write(6,*)('single excitations will not be included')        1d3s20
          else                                                          1d3s20
          iunc(1)=-1                                                    5d27s19
          write(6,*)('singles will be uncontracted ')                   5d27s19
          end if                                                        1d3s20
        end if                                                          5d24s18
        if(line200(is:is+5).eq.'decomp')then                            8d10s22
         ndecomp=1                                                      8d10s22
         go to 3251                                                     8d10s22
        end if                                                          8d10s22
        if(line200(is:is+4).eq.'readc')then                             8d10s22
         nreadc=1                                                       8d10s22
         is=ie+1                                                        12d13s22
         call delim(line200,is,ie)                                         5d24s18
         ncfile=ie+1-is                                                 12d13s22
         cfile(1:ncfile)=line200(is:ie)                                 12d13s22
         go to 3251                                                     8d10s22
        end if                                                          8d10s22
        if(line200(is:is+4).eq.'doubx')then                              8d22s18
c
c     parameters for double excitations.                                8d22s18
c     ethred: 10**ethred is cutoff for eigenvalues of overlap           8d25s22
c
c                                                                       5d24s18
         is=ie+1                                                        5d24s18
         call delim(line200,is,ie)                                         5d24s18
         if(line200(is:is).eq.'u')then                                  3d29s19
          iunc(2)=-1                                                     3d29s19
          ethred=0d0                                                    1d6s20
          write(6,*)('doubles will be uncontracted')
          maxds=1                                                       8d10s22
         else if(line200(is:is+2).eq.'off')then                         1d3s20
          write(6,*)('double excitations will not be included')         1d3s20
          ethred=0d0                                                    1d3s20
          nethred=1                                                     10d6s22
         else                                                           1d3s20
          write(6,*)('doubles will be contracted ')                     1d3s20
          if(ie.ge.is)then                                              8d25s22
           read(line200(is:ie),*)ethred                                 8d25s22
           nethred=1                                                    10s3s22
           write(6,*)('ethred is now '),ethred                          8d25s22
          end if                                                        8d25s22
         end if                                                         3d29s19
        end if                                                          5d24s18
        nlzzq2(1)=nlzzci                                                12d5s22
        nlzzq2(2)=lambdaci                                              12d5s22
        go to 3251                                                      5d24s18
c
c     input for properties code.                                        11d21s19
c     just list of file name.                                           11d21s19
c     if no further input, it will be assumed that there is the single  11d21s19
c     file name wavef.                                                  11d21s19
c     sotm                                                              2d7s23
c     so4v                                                              2d7s23
c     itype
c     norb
c     genwf
c     smsz
c     geometry
c     BODC
c                                                                       11d21s19
 2252   continue                                                        11d21s19
        write(6,*)('property input')                                    11d21s19
        iprop(7)=0                                                      11d2s22
        ibodc=0                                                         9d17s24
        nmn=1                                                           8d12s22
        idoreln=0                                                       6d3s21
        nfname=0                                                          11d22s19
        nfneed=0                                                        12d5s22
        ntmso=1                                                         8d16s22
        n4vso=1                                                         2d7s23
        genwfeps=1d-10                                                  11d9s22
        smsz=1d-8                                                       11d16s23
 2225   continue
        read(5,200,end=2226)line200                                     11d22s19
        if(idoreln.ne.0)idorel=idoreln                                  6d3s21
        is=1                                                            11d22s19
        call delimb(line200,is,ie)                                      6d1s21
        if(line200(is:is).eq.'#')go to 2225                             10d26s20
        if(line200(is:is+4).eq.'debug')then                             1d3s20
 6604   continue                                                        1d3s20
         is=ie+1                                                        1d3s20
         call delim(line200,is,ie)                                      1d3s20
         if(ie.ge.is)then                                               1d3s20
          write(6,*)('debug print out turned on for '),line200(is:ie)    1d3s20
          nchere=ie+1-is                                                 1d3s20
          do i=1,nprtr                                                   1d3s20
           if(line200(is:ie).eq.prtrname(i)(1:nchere))then               1d3s20
            iprtr(i)=1                                                   1d3s20
            go to 6604                                                   1d3s20
           end if                                                        1d3s20
          end do                                                         1d3s20
          write(6,*)('no matching string in prtrname found ')           1d3s20
         end if                                                         1d3s20
         go to 2225                                                     1d3s20
        end if                                                          1d3s20
        if(line200(is:is+3).eq.'smsz')then                              11d16s23
         is=ie+1                                                        11d16s23
         call delim(line200,is,ie)                                      11d16s23
         if(ie.ge.is)then                                               11d16s23
          read(line200(is:ie),*)smsz                                    11d16s23
          if(smsz.lt.0d0)then                                           11d16s23
           smsz=1d1**smsz                                               11d16s23
          end if                                                        11d16s23
          write(6,*)('new threshold smsz is '),smsz                     11d16s23
         else                                                           11d16s23
          write(6,*)('I can not find argument to smsz keyword!!!')      11d16s23
          write(6,*)('help!')                                           11d16s23
          irtrn=1                                                       11d16s23
          return                                                        11d16s23
         end if                                                         11d16s23
         go to 2225                                                     11d16s23
        end if                                                          11d16s23
        if(line200(is:is+3).eq.'BODC')then                              9d17s24
         write(6,*)                                                     9d17s24
     $      ('we will compute the Born-Oppenheimer diagonal correction')9d17s24
         write(6,*)('using finite difference.')                         9d17s24
         ibodc=1                                                        9d17s24
        end if                                                          9d17s24
        if(line200(is:is+7).eq.'geometry')then                          9d17s24
         write(6,*)('will be using a geometry different than that')     5d15s23
         write(6,*)('used to generate the wavefunction. ')              5d15s23
         write(6,*)('did we use carts or vibrations as input?')         5d15s23
        end if                                                          5d15s23
        if(line200(is:is+4).eq.'genwf')then                             11d9s22
         is=ie+1                                                        11d9s22
         call delim(line200,is,ie)                                      11d9s22
         if(ie.ge.is)then                                               11d9s22
          read(line200(is:ie),*)genwfeps                                11d9s22
          write(6,*)('new value of genwfeps = '),genwfeps               11d9s22
         end if                                                         11d9s22
         go to 2225                                                     11d9s22
        end if                                                          11d9s22
        if(line200(is:is+3).eq.'ityp')then                              6d3s21
         is=ie+1                                                        6d3s21
         call delim(line200,is,ie)                                      6d3s21
         if(ie.lt.is)then                                               6d3s21
          write(6,*)('could not find type on ityp line')                6d3s21
          irtrn=1                                                       6d3s21
          return                                                        6d3s21
         end if                                                         6d3s21
         write(6,*)('integral type to use: '),line200(is:ie)            6d3s21
         if(line200(is:is+2).eq.'ssg')then                              6d3s21
          idoreln=-4                                                    2d18s20
         else if(line200(is:is+2).eq.'ssb')then                              6d3s21
          idoreln=-2                                                    6d3s21
         else if(line200(is:is+1).eq.'ss')then                          3d24s22
          idoreln=2                                                     3d24s22
         else if(line200(is:is+4).eq.'no so')then                       4d1s24
          write(6,*)('don''t compute spin-orbit integrtals ')           4d1s24
          idoreln=10                                                    4d1s24
         else                                                           6d3s21
          write(6,*)('unknown integral type: '),line200(is:ie)
          irtrn=1                                                       6d3s21
          return                                                        6d3s21
         end if                                                         6d3s21
         idorel=idoreln                                                 6d3s21
         go to 2225                                                     6d3s21
        end if                                                          6d3s21
        if(line200(is:is+3).eq.'norb')then                              11d2s22
         iprop(7)=1                                                     11d2s22
         write(6,*)                                                     11d2s22
     $       ('natural orbitals will be computed instead of properties')11d2s22
         write(6,*)('average over all roots present')                   11d2s22
         go to 2225                                                     11d2s22
        end if                                                          11d2s22
        if(line200(is:is+3).eq.'sotm')then                              10d6s22
         is=ie+1                                                        10d5s22
         call delim(line200,is,ie)                                      10d5s22
         if(.not.(ie.lt.is.or.line200(is:is+1).eq.'on'))then            10d6s22
          write(6,*)('spin-orbit transition moments turned off')         8d16s22
          ntmso=0                                                        8d16s22
         end if                                                         10d6s22
         go to 2225                                                     8d16s22
        end if                                                          8d16s22
        if(line200(is:is+3).eq.'so4v')then                              2d8s23
         is=ie+1                                                        10d5s22
         call delim(line200,is,ie)                                      10d5s22
         write(6,*)('is,ie '),is,ie,line200(is:is+1)
         if(.not.(ie.lt.is.or.line200(is:is+1).eq.'on'))then            2d8s23
          write(6,*)('4 virtual integrals neglected in spin-orbit'),    2d8s23
     $        (' matrix elements')                                      2d8s23
          n4vso=0                                                       2d8s23
         end if                                                         10d6s22
         go to 2225                                                     8d16s22
        end if                                                          8d16s22
        if(ie.lt.is)go to 2225                                          3d29s22
        write(6,*)('wave function file name: '),line200(is:ie)          11d22s19
        if(nfname.eq.0.and.numeminus.eq.0)then                          11d22s19
         open(unit=1,file=line200(is:ie),form='unformatted')            11d22s19
         read(1)isumlt,isymmrci,nroot,norb,nsymb                        11d25s19
         write(6,*)('other isumlt readin '),isumlt
         islow=isumlt                                                   7d26s22
         ishig=isumlt                                                   7d26s22
         iprop(1)=islow                                                 7d26s22
         iprop(2)=ishig                                                 7d26s22
         read(1)(idoubo(i),irefo(i),nbasc(i),i=1,nsymb)                 11d25s19
         read(1)(irel(i),ism(i),i=1,norb)                               11d25s19
         call loadr(potdws,ibdat,ibstor,isstor,nsymb,istinfo,ncoreg,    1d2s20
     $        nvalg,nlzzq,iextradatd,nectradatx,0,nbasdws,element,idum, 7d26s22
     $        nec,mdon,mdoop,2,bc,ibc,icanon)                           5d5s23
         nfned=1                                                        12d5s22
         if(nlzzq2(1).eq.6)then                                         12d5s22
          nfned=2*nlzzq2(2)+1                                           12d5s22
         else if(nlzzq2(1).eq.2)then                                    12d5s22
c     to be compatible with files from previous versions that did not   10d7s24
c     include lambda in nlzzq.                                          10d7s24
          nfned=2                                                       10d7s24
         end if                                                         12d5s22
         nfneed=nfneed+nfned                                            12d5s22
         iprop(3)=nec                                                   7d26s22
         iprop(4)=mdon                                                  7d26s22
         iprop(5)=mdoop                                                 7d26s22
         iprop(6)=0                                                     7d30s22
         close(unit=1)                                                     11d20s19
         first(1)=potdws                                                11d20s19
         irefdata=ibcoff                                                11d17s21
         ibcoff=irefdata+3+nsymb+3*natom                                11d17s21
         call enough('basisz. 12',bc,ibc)
         ibc(irefdata)=natom                                            11d17s21
         ibc(irefdata+1)=norb                                           11d17s21
         ibc(irefdata+2)=nsymb                                          11d17s21
         jrefd=irefdata+2
         do isb=1,nsymb
          ibc(jrefd+isb)=nbasc(isb)                                     11d17s21
         end do                                                         11d17s21
         jrefd=jrefd+nsymb                                              11d17s21
         do ia=1,natom                                                  11d17s21
          do ixyz=1,3                                                   11d17s21
           bc(jrefd+ixyz)=xcart(ixyz,ia)                                11d17s21
          end do                                                        11d17s21
          jrefd=jrefd+3                                                 11d17s21
         end do                                                         11d17s21
         if(idoreln.ne.0.and.idorel.eq.0)then                           10d20s21
          write(6,*)(' ')
          write(6,*)
     $  ('I''m sorry, you have specified a relativistic integral type,')
          write(6,*)('but the input wave function is non-relativistic.')
          write(6,*)('I can not handle that case.')                     10d20s21
          write(6,*)(' ')
          irtrn=1                                                       10d20s21
          return                                                        10d20s21
         end if                                                         10d20s21
         if(idoreln.ne.0)idorel=idoreln                                  6d3s21
         write(6,*)('back from loadr')
        else                                                            7d26s22
         open(unit=1,file=line200(is:ie),form='unformatted')            11d22s19
         read(1)isumlt,isymmrci,nroot,norb,nsymb                        11d25s19
         islow=min(islow,isumlt)                                        7d26s22
         ishig=max(ishig,isumlt)                                        7d26s22
         iprop(1)=islow                                                 7d26s22
         iprop(2)=ishig                                                 7d26s22
c     idoub
         read(1)idum                                                    7d30s22
c     irel
         read(1)idum                                                    7d30s22
c     nsymb
         nread1=8+8*8                                                   12d5s22
         nread1h=nread1/2                                               12d5s22
         nread2=6+11+1                                                  12d13s22
         nread2h=nread2/2                                               12d5s22
         ibuf2=ibcoff+nread1h                                           12d5s22
         read(1)(ibc(ibcoff+ir),ir=0,nread1h-1),dum1,dum,               12d5s22
     $        (ibc(ibuf2+ir),ir=0,nread2h-1)
         ipack8=ibc(ibuf2+nread2h-1)
         nextradata=ipack4(2)
         ipack8=ibc(ibcoff)
         write(6,*)('nsymb,idorel: '),ipack4
         ipack8=ibc(ibcoff+1)
         write(6,*)('ngauss,natom: '),ipack4
c     morb
         read(1)idum                                                    7d30s22
c     morbc
         read(1)idum                                                    7d30s22
c     morbp
         read(1)idum                                                    7d30s22
         nneed=12+2*natom*3                                                7d30s22
         ibuf=ibcoff                                                    7d30s22
         ibcoff=ibuf+nneed                                              7d30s22
         nextradata=nextradata+1                                           1d10s19
         iextradata=ibcoff                                                 1d10s19
         ibcoff=iextradata+nextradata                                      1d10s19
         call enough('basisz. 13',bc,ibc)
         jbuf=ibuf-1                                                    7d30s22
         nxtot=nneed+nextradata+nsymb
         read(1)(bc(jbuf+i),i=1,nneed),                                 12d5s22
     $        (bc(iextradata+i),i=0,nextradata-1),                      12d5s22
     $        (ibc(ibcoff+ir),ir=0,nsymb-1),nlzzq                       12d5s22
         ibcoff=iextradata                                              12d5s22
         nfned=1                                                        12d5s22
         if(nlzzq2(1).eq.6)then                                         12d5s22
          nfned=2*nlzzq2(2)+1                                           12d15s22
         else if(nlzzq2(1).eq.2)then                                    12d5s22
c     to be compatible with files from previous versions that did not   10d7s24
c     include lambda in nlzzq.                                          10d7s24
          nfned=2                                                       10d7s24
         end if                                                         12d5s22
         nfneed=nfneed+nfned                                            12d5s22
         write(6,*)('carts from this file: ')
         jbuf=ibuf+12                                                   7d30s22
         rms=0d0                                                        7d30s22
         do i=1,natom
          do ixyz=1,3
           rms=rms+(bc(jbuf)-xcart(ixyz,i))**2
           jbuf=jbuf+2
          end do
         end do
         rms=sqrt(rms/dfloat(natom*3))                                  7d30s22
         write(6,*)('rms change in geometry: '),rms                     7d30s22
         if(rms.gt.1d-10)iprop(6)=1                                     7d30s22
         close(unit=1)                                                  7d26s22
        end if                                                          11d22s19
        nfname=nfname+1                                                     11d22s19
        if(nfname.eq.1)then                                               11d22s19
         ifname=ibcoff                                                  11d22s19
         jfname=ifname                                                  11d22s19
        end if                                                          11d22s19
        nhere=ie+1-is
        ibcoff=jfname+nhere+1                                           11d22s19
        call enough('basisz. 14',bc,ibc)
        ibc(jfname)=nhere                                               11d22s19
        jfname=jfname+1                                                 11d22s19
        do i=is,ie                                                      11d22s19
         ibc(jfname)=ichar(line200(i:i))                                11d22s19
         jfname=jfname+1                                                11d22s19
        end do                                                          11d22s19
        go to 2225                                                      11d22s19
 2226   continue                                                        11d22s19
        write(6,*)('what we have for nfname: '),nfname
        if(nfname.eq.0)then                                               11d22s19
         nfname=1                                                       11d22s19
         write(6,*)('nfname is zero, numeminus = '),numeminus
         if(numeminus.eq.0)then                                         11d22s19
          open(unit=1,file='wavef',form='unformatted')                   11d22s19
          read(1)isumlt,isymmrci,nroot,norb,nsymb                        11d25s19
          write(6,*)('isumlt readin '),isumlt
          islow=isumlt                                                  7d26s22
          ishig=isumlt                                                  7d26s22
          iprop(1)=islow                                                7d26s22
          iprop(2)=ishig                                                7d26s22
          read(1)(idoubo(i),irefo(i),nbasc(i),i=1,nsymb)                 11d25s19
          read(1)(irel(i),ism(i),i=1,norb)                               11d25s19
          call loadr(potdws,ibdat,ibstor,isstor,nsymb,istinfo,ncoreg,   1d2s20
     $         nvalg,nlzzq,iextradatd,nectradatx,0,nbasdws,element,idum,7d26s22
     $         nec,mdon,mdoop,2,bc,ibc,icanon)                          5d5s23
          iprop(3)=nec                                                  7d26s22
          iprop(4)=mdon                                                 7d26s22
          iprop(5)=mdoop                                                7d26s22
          close(unit=1)                                                     11d20s19
          first(1)=potdws                                                11d20s19
         end if                                                         11d22s19
         ifname=ibcoff                                                  11d22s19
         ibcoff=ifname+6                                                11d22s19
         call enough('basisz. 15',bc,ibc)
         nfneed=1                                                       12d5s22
         if(nlzzq2(1).eq.6)then                                         12d5s22
          nfneed=2*nlzzq2(2)+1                                          12d5s22
         else if(nlzzq2(1).eq.2)then                                    12d5s22
          if(nlzzq2(2).ne.0)nfneed=2                                    12d5s22
         end if                                                         12d5s22
         ibc(ifname)=5                                                  11d22s19
         line200(1:5)='wavef'                                           11d22s19
         jfname=ifname+1                                                11d22s19
         do i=1,5                                                       11d22s19
          ibc(jfname)=ichar(line200(i:i))                               11d22s19
          jfname=jfname+1                                               11d22s19
         end do                                                         11d22s19
        end if                                                          11d22s19
         write(6,*)('return2 from basisz')
        return                                                          1d14s20
c                                                                       1d10s19
c     diabatic orbitals.                                                1d10s19
c     just input list of orbital file names, or optionally the keyword  10d8s24
c     dtest can be included to adjust the threshold for testing the     10d8s24
c     product of diagonals.                                             10d8s24
c     the first in the list is the set of orbitals to be rotated        1d10s19
c     while the remaining are reference orbitals to compare to.         1d10s19
c     force means use bcode from raw orbitals rather than ref orbitals. 7d17s23
c                                                                       1d10s19
 2253   continue                                                        1d10s19
        write(6,*)('diab input ')
        ibcodex=0                                                       2d8s20
        nfname=0                                                          11d22s19
        dtest=0.9d0                                                     10d8s24
        iforce=0                                                        7d17s23
        ifdiff=0                                                        6d24s24
 2254   continue
        read(5,200,end=2255)line200                                     11d22s19
        is=1                                                            11d22s19
        call delim(line200,is,ie)                                       11d22s19
        if(line200(is:is).eq.'#')go to 2254                             10d31s23
        if(line200(is:is+4).eq.'dtest')then                             10d8s24
         is=ie+1                                                        10d8s24
         call delim(line200,is,ie)                                      10d8s24
         if(ie.ge.is)then                                               10d8s24
          read(line200(is:ie),*)dtest                                   10d8s24
          write(6,*)('new value of dtest is '),dtest                    10d8s24
         else                                                           10d8s24
          write(6,*)('ooops there is no numerical value on the dtest'), 10d8s24
     $         ('line. help!')                                          10d8s24
          irtrn=1                                                       10d8s24
          return                                                        10d8s24
         end if                                                         10d8s24
         go to 2254                                                     10d8s24
        end if                                                          10d8s24
        if(line200(is:is+4).eq.'fdiff')then                             6d24s24
         write(6,*)('finite difference of orbitals ...')                6d24s24
         ie0=ie+1                                                       6d24s24
 2554    continue                                                       6d24s24
         is=ie+1                                                        6d24s24
         call delim(line200,is,ie)                                      6d24s24
         if(ie.ge.is)then                                               6d24s24
          write(6,*)('stepsize '),line200(is:ie)                         6d24s24
          ifdiff=ifdiff+1                                               6d24s24
          go to 2554
         end if                                                         6d24s24
         ifdstep=ibcoff                                                 6d24s24
         ibcoff=ifdstep+ifdiff                                          6d24s24
         call enough('basisz.fdstep',bc,ibc)                            6d24s24
         ifdiff=0                                                       6d24s24
         ie=ie0                                                         6d24s24
 2552    continue                                                       6d24s24
         is=ie+1                                                        6d24s24
         call delim(line200,is,ie)                                      6d24s24
         if(ie.ge.is)then                                               6d24s24
          read(line200(is:ie),*)bc(ifdstep+ifdiff)                      6d24s24
          ifdiff=ifdiff+1                                               6d24s24
          go to 2552                                                    6d24s24
         end if                                                         6d24s24
         go to 2254                                                     6d24s24
        end if                                                          6d24s24
        if(line200(is:is+3).eq.'forc')then                              7d17s23
         iforce=1                                                       7d17s23
         write(6,*)('use bcode from raw orbitals ...')                  7d17s23
         go to 2254                                                     7d17s23
        end if                                                          7d17s23
        if(line200(is:is+4).eq.'debug')then                             1d3s20
 2256   continue                                                        1d3s20
         is=ie+1                                                        1d3s20
         call delim(line200,is,ie)                                      1d3s20
         if(ie.ge.is)then                                               1d3s20
          write(6,*)('debug print out turned on for '),line200(is:ie)    1d3s20
          nchere=ie+1-is                                                 1d3s20
          do i=1,nprtr                                                   1d3s20
           if(line200(is:ie).eq.prtrname(i)(1:nchere))then               1d3s20
            iprtr(i)=1                                                   1d3s20
            go to 2256                                                  1d10s19
           end if                                                        1d3s20
          end do                                                         1d3s20
          write(6,*)('no matching string in prtrname found ')           1d3s20
         end if                                                         1d3s20
         go to 2256                                                     1d10s19
        end if                                                          1d3s20
        nfname=nfname+1                                                     11d22s19
        if(nfname.eq.1)then                                             1d10s19
         write(6,*)('raw orbital file name: '),line200(is:ie)           2d8s20
         nrawfn=ie+1-is                                                 2d9s20
         bline200(1:nrawfn)=line200(is:ie)                              2d9s20
        else                                                            1d10s19
         write(6,*)('reference orbital file name: '),line200(is:ie)     2d8s20
        end if                                                          1d10s19
        if(nfname.gt.iddx)then                                          1d10s19
         write(6,*)('no. of references exceeds dimension iddx = '),iddx 2d8s20
         stop                                                           1d10s19
        end if                                                          1d10s19
        open(unit=1,file=line200(is:ie),form='unformatted')             1d10s19
        call loadr(potdws,ibdat,ibstor,isstor,nsymb,istinfo,ncoreg,     1d10s19
     $        nvalg,nlzzq,iextradatd,nextradatx,0,nbasdws,element,idum, 7d26s22
     $       idum,idum,idum,1,bc,ibc,icanon)                            5d5s23
        call loadr2(ibcode(1,nfname),iorb(1,nfname),nsymb,nbasdws,      2d8s20
     $       iorbao(1,nfname),ibcodex,bc,ibc)                           11d9s22
        close(unit=1)                                                     11d20s19
        go to 2254                                                      1d10s19
 2255   continue                                                        1d10s19
        write(6,*)('idorel = '),idorel
        if(idorel.eq.0)then                                             2d9s20
         ncomp=1                                                        2d9s20
        else                                                            2d9s20
         ncomp=2                                                        2d9s20
        end if                                                          2d9s20
        do isb=1,nsymb                                                  2d8s20
         nbasbb=nbasb(isb)*ncomp                                        2d9s20
         if(nbasb(isb).gt.0)then                                        1d17s24
          write(6,*)('for symmetry block '),isb
          if(ifdiff.gt.0)then                                           6d24s24
           iorbxx=ibcoff                                                6d24s24
           idtmp0=iorbxx+nbasdws(isb)*nbasdws(isb)                      6d24s24
           iorbaa=idtmp0+nbasdws(isb)*nbasdws(isb)                      6d24s24
           ibcoff=iorbaa+nbasdws(isb)*nbasdws(isb)                      6d24s24
           call enough('basisz.dtmp01',bc,ibc)                          6d24s24
c
c     dv=v*dt
c
           do i=0,nbasdws(isb)-1                                        6d24s24
            do j=0,nbasdws(isb)-1                                       6d24s24
             ji=iorbxx+j+nbasdws(isb)*i                                 6d24s24
             ij=iorb(isb,1)+i+nbasdws(isb)*j                            6d24s24
             bc(ji)=bc(ij)                                              6d24s24
            end do                                                      6d24s24
           end do                                                       6d24s24
           do jfdiff=1,ifdiff                                           6d24s24
            jfdiffm=jfdiff-1                                            6d24s24
            write(6,*)('for stepsize '),bc(ifdstep+jfdiffm)             6d24s24
            if1=2+jfdiffm*2                                             6d24s24
            if2=if1+1                                                   6d24s24
            idtmp=ibcoff                                                6d24s24
            ibcoff=idtmp+nbasdws(isb)*nbasdws(isb)                      6d24s24
            ilsq=ibcoff                                                 6d24s24
            call enough('dorb.dtmp',bc,ibc)                             6d24s24
            xx=0.5d0/bc(ifdstep+jfdiffm)                                6d24s24
            do i=0,nbasdws(isb)*nbasdws(isb)-1                          6d24s24
             bc(idtmp+i)=xx*(bc(iorb(isb,if1)+i)-bc(iorb(isb,if2)+i))   6d24s24
            end do                                                      6d24s24
            call prntm2(bc(idtmp),nbasdws(isb),nbasdws(isb),            6d24s24
     $           nbasdws(isb))                                          6d24s24
            if(jfdiff.gt.1)then                                         6d24s24
             do i=0,nbasdws(isb)*nbasdws(isb)-1                         6d24s24
              xder=(4d0*bc(idtmp0+i)-bc(idtmp+i))/3d0                   6d24s24
              bc(idtmp0+i)=xder                                         6d24s24
             end do                                                     6d24s24
             write(6,*)('extrapolated ')                                6d24s24
             call prntm2(bc(idtmp0),nbasdws(isb),nbasdws(isb),          6d24s24
     $           nbasdws(isb))                                          6d24s24
             call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasdws(isb), 6d24s24
     $            1d0,bc(iorbxx),nbasdws(isb),bc(idtmp0),nbasdws(isb),  6d24s24
     $            0d0,bc(iorbaa),nbasdws(isb),'basisz.orbaa')           6d24s24
             call prntm2(bc(iorbaa),nbasdws(isb),nbasdws(isb),
     $            nbasdws(isb))
            end if                                                      6d24s24
            do i=0,nbasdws(isb)*nbasdws(isb)-1                          6d24s24
             bc(idtmp0+i)=bc(idtmp+i)                                   6d24s24
            end do                                                      6d24s24
            ibcoff=idtmp                                                6d24s24
           end do                                                       6d24s24
           ibcoff=iorbaa                                                6d24s24
          else                                                          6d24s24
           do if=2,nfname                                                 2d8s20
            do i=0,nbasdws(isb)-1                                         2d8s20
             if(ibc(ibcode(isb,1)+i).ne.ibc(ibcode(isb,if)+i))then        2d8s20
              write(6,*)('ibcode error: '),isb,if,i,                    6d24s24
     $            ibc(ibcode(isb,1)+i),ibc(ibcode(isb,if)+i)            6d24s24
              if(iforce.ne.0)then                                         7d17s23
               write(6,*)('but ignore this ')                             7d17s23
               ibc(ibcode(isb,if)+i)=ibc(ibcode(isb,1)+i)                 7d17s23
              else                                                        7d17s23
               write(6,*)('to force continue, use the force keyword')     7d27s23
               stop 'diab'                                                7d17s23
              end if                                                      7d17s23
             end if                                                       2d8s20
            end do                                                        2d8s20
           end do                                                         2d8s20
           iorbt=ibcoff                                                   2d9s20
           ibcoff=iorbt+nbasdws(isb)*nbasc(isb)                           2d9s20
           call enough('basisz. 16',bc,ibc)
           nnxh=nbasc(isb)-icanon(isb)                                    5d5s23
           do i=0,nnxh-1                                                  5d5s23
            do j=0,nbasdws(isb)-1                                         2d9s20
             ji=iorbt+j+nbasdws(isb)*i                                    2d9s20
             ij=iorb(isb,1)+i+nnxh*j                                      5d5s23
             bc(ji)=bc(ij)                                                2d9s20
            end do                                                        2d9s20
           end do                                                         2d9s20
           is=0                                                           2d8s20
 2257      continue                                                       2d8s20
           igoal=ibc(ibcode(isb,1)+is)                                    2d8s20
           do i=is,nbasdws(isb)-1                                         2d8s20
            if(ibc(ibcode(isb,1)+i).eq.igoal)ie=i                         2d8s20
           end do                                                         2d8s20
           nhere=ie+1-is                                                  2d8s20
           if(igoal.eq.1)then                                            6d20s24
            write(6,*)('no. fcns with bcode '),igoal,                    6d20s24
     $         ('(doubly occupied orbitals) is '),nhere                 6d20s24
           else if(igoal.eq.2)then                                       6d20s24
            write(6,*)('no. fcns with bcode '),igoal,                    6d20s24
     $         ('(active orbitals) is '),nhere                          6d20s24
           else if(igoal.eq.3)then                                       6d20s24
            write(6,*)('no. fcns with bcode '),igoal,                    6d20s24
     $         ('(Rydberg orbitals) is '),nhere                         6d20s24
           else
            write(6,*)('no. fcns with bcode '),igoal,                    6d20s24
     $         ('(virtual orbitals) is '),nhere                         6d20s24
           end if                                                        6d20s24
           if(nhere.gt.1)then                                             2d8s20
            sz=0d0                                                        2d8s20
            numb=0                                                        2d8s20
            iouse=ibcoff                                                  2d8s20
            ibcoff=iouse+nhere*nhere*(nfname-1)                           2d9s20
            call enough('basisz. 17',bc,ibc)
            jouse=iouse                                                   2d8s20
            do if=2,nfname                                                2d8s20
             iorb1=iorbt+is                                               2d9s20
             iorb2=iorb(isb,if)+nnxh*is                                   5d5s23
             call dgemm('n','n',nhere,nhere,nnxh,1d0,bc(iorb1),           5d5s23
     $          nbasdws(isb),bc(iorb2),nnxh,0d0,bc(jouse),nhere,        5d5s23
     d' basisz.  1')
             do i1=0,nhere-2                                              2d9s20
              do i2=i1+1,nhere-1                                          2d9s20
               i12=jouse+i1+nhere*i2                                      2d8s20
               i21=jouse+i2+nhere*i1                                      2d8s20
               sz=sz+bc(i12)**2+bc(i21)**2                                2d8s20
               numb=numb+2                                                2d8s20
              end do                                                      2d8s20
             end do                                                       2d8s20
             jouse=jouse+nhere*nhere                                      2d8s20
            end do                                                        2d8s20
            sz=sqrt(sz/dfloat(numb))                                      2d8s20
            write(6,*)('starting size: '),sz                              2d8s20
            if(igoal.eq.1.and.iforce.eq.0)go to 2258                      7d18s23
            szo=sz                                                        2d8s20
            do isweep=1,300                                               2d8s20
             do i1=0,nhere-2                                              2d9s20
              i1p=i1+is                                                   2d9s20
              do i2=i1+1,nhere-1                                          2d9s20
               i2p=i2+is                                                  2d9s20
               a=0d0                                                      2d8s20
               b=0d0                                                      2d8s20
               c=0d0                                                      2d8s20
               jouse=iouse                                                2d8s20
               do if=2,nfname                                             2d8s20
                i11=jouse+i1+nhere*i1                                     2d8s20
                i21=jouse+i2+nhere*i1                                     2d8s20
                i12=jouse+i1+nhere*i2                                     2d8s20
                i22=jouse+i2+nhere*i2                                     2d8s20
                a=a+bc(i11)*bc(i21)-bc(i12)*bc(i22)                       2d8s20
                b=b+bc(i11)**2+bc(i22)**2-bc(i12)**2-bc(i21)**2           2d8s20
                jouse=jouse+nhere*nhere                                   2d8s20
               end do                                                     2d8s20
               c=-a
               dist2=b*b-4d0*a*c
               if(dist2.lt.0d0)then
                write(6,*)('discriminant is negative!!! '),dist2
                write(6,*)isweep,i1,i2,a,b
                stop
               end if
               dist=sqrt(abs(dist2))
               q=-0.5d0*(b+sign(dist,b))
               t2=c/q
               if(a.eq.0d0)then                                           1d19s23
                t1=t2                                                     1d19s23
               else                                                       1d19s23
                t1=q/a                                                    1d19s23
               end if                                                     1d19s23
               t12=t1**2
               s12=t12/(1d0+t12)
               c1=sqrt(1d0-s12)
               s1=sign(sqrt(s12),t1)
               t22=t2**2
               s22=t22/(1d0+t22)
               c2=sqrt(1d0-s22)
               s2=sign(sqrt(s22),t2)
               phi1=0d0
               phi2=0d0
               jouse=iouse                                                2d8s20
               do if=2,nfname                                             2d8s20
                i11=jouse+i1+i1*nhere                                     2d8s20
                i21=jouse+i2+i1*nhere                                     2d8s20
                i12=jouse+i1+i2*nhere                                     2d8s20
                i22=jouse+i2+i2*nhere                                     2d8s20
                phi1=phi1+(-s1*bc(i11)+c1*bc(i21))**2                     2d8s20
     $             +(c1*bc(i12)+s1*bc(i22))**2                          2d8s20
                phi2=phi2+(-s2*bc(i11)+c2*bc(i21))**2                     2d8s20
     $             +(c2*bc(i12)+s2*bc(i22))**2                          2d8s20
                jouse=jouse+nhere*nhere                                   2d8s20
               end do                                                     2d8s20
               if(phi1.lt.phi2)then
                cu=c1
                su=s1
               else
                cu=c2
                su=s2
               end if
               jouse=iouse                                                2d9s20
               do if=2,nfname                                             2d9s20
                do n=0,nhere-1                                            2d9s20
                 i1n=jouse+i1+nhere*n                                     2d9s20
                 i2n=jouse+i2+nhere*n                                     2d9s20
                 x1=cu*bc(i1n)+su*bc(i2n)                                 2d9s20
                 y1=-su*bc(i1n)+cu*bc(i2n)
                 bc(i1n)=x1                                               2d9s20
                 bc(i2n)=y1                                               2d9s20
                end do                                                    2d9s20
                jouse=jouse+nhere*nhere                                   2d9s20
               end do                                                     2d9s20
               do n=0,nnxh-1                                              5d5s23
                in1=iorb(isb,1)+n+nnxh*i1p                                5d5s23
                in2=iorb(isb,1)+n+nnxh*i2p                                5d5s23
                x1=cu*bc(in1)+su*bc(in2)                                  2d9s20
                y1=-su*bc(in1)+cu*bc(in2)                                 2d9s20
                bc(in1)=x1                                                2d9s20
                bc(in2)=y1                                                2d9s20
               end do                                                     2d9s20
               do n=0,nbasbb-1                                            2d9s20
                in1=iorbao(isb,1)+n+nbasbb*i1p                            2d9s20
                in2=iorbao(isb,1)+n+nbasbb*i2p                            2d9s20
                x1=cu*bc(in1)+su*bc(in2)                                  2d9s20
                y1=-su*bc(in1)+cu*bc(in2)                                 2d9s20
                bc(in1)=x1                                                2d9s20
                bc(in2)=y1                                                2d9s20
               end do                                                     2d9s20
              end do                                                      2d8s20
             end do                                                       2d8s20
             sz=0d0                                                       2d9s20
             jouse=iouse                                                  2d9s20
             do if=2,nfname                                               2d9s20
              do i1=0,nhere-2                                             2d9s20
               do i2=i1+1,nhere-1                                         2d9s20
                i12=jouse+i1+nhere*i2                                     2d9s20
                i21=jouse+i2+nhere*i1                                     2d9s20
                sz=sz+bc(i12)**2+bc(i21)**2                               2d9s20
               end do                                                     2d9s20
              end do                                                      2d9s20
              jouse=jouse+nhere*nhere                                     2d9s20
             end do                                                       2d9s20
             sz=sqrt(sz/dfloat(numb))                                     2d9s20
             write(6,*)('at end of sweep we have '),sz                    2d9s20
             if(sz/szo.gt.0.9999d0)go to 2258                             2d9s20
             szo=sz                                                       2d9s20
            end do                                                        2d8s20
 2258       continue                                                      2d9s20
            det=1d0                                                       2d9s20
            do i1=0,nhere-1                                               2d9s20
             i1p=i1+is                                                    2d9s20
             i11=iouse+i1+nhere*i1                                        2d9s20
             if(bc(i11).lt.0d0)then                                       2d9s20
              do n=0,nnxh-1                                               5d5s23
               in1=iorb(isb,1)+n+nnxh*i1p                                 5d5s23
               bc(in1)=-bc(in1)                                           2d9s20
              end do                                                      2d9s20
              do n=0,nbasbb-1                                             2d9s20
               in1=iorbao(isb,1)+n+nbasbb*i1p                             2d9s20
               bc(in1)=-bc(in1)                                           2d9s20
              end do                                                      2d9s20
             end if                                                       2d9s20
             det=det*abs(bc(i11))                                         2d9s20
            end do                                                        2d9s20
            write(6,*)('product of diagonals: '),det                      2d9s20
            if(det.lt.dtest)then                                        10d8s24
             write(6,*)('yikes!')                                         2d15s22
             write(6,*)('I would say diabatization failed! ')             2d15s22
             do if=2,nfname                                                2d8s20
              iorb1=iorbt+is                                               2d9s20
              iorb2=iorb(isb,if)+nnxh*is                                  6d12s23
              call dgemm('n','n',nhere,nhere,nnxh,1d0,bc(iorb1),          6d12s23
     $          nbasdws(isb),bc(iorb2),nnxh,0d0,bc(ibcoff),nhere,       6d12s23
     d' basisz.  1.1')
              write(6,*)('this is the overlap I''m looking at: ')        3d14s24
              call prntm2(bc(ibcoff),nhere,nhere,nhere)
              rmsu=0d0                                                   3d14s24
              irmsu=0                                                    3d14s24
              do i=0,nhere*nhere-1                                       3d14s24
               if(abs(bc(ibcoff+i)).gt.1d-1)then                         3d14s24
                rmsu=rmsu+(abs(bc(ibcoff+i))-1d0)**2                     3d15s24
                irmsu=irmsu+1                                            3d14s24
               end if                                                    3d15s24
              end do                                                     3d14s24
              if(irmsu.eq.nhere)then                                     3d15s24
               rmsu=sqrt(rmsu/dfloat(irmsu))                             3d15s24
              else                                                       3d15s24
               rmsu=1d0                                                  3d15s24
              end if                                                     3d15s24
              if(rmsu.lt.1d-2)then                                       3d15s24
               write(6,*)                                                3d15s24
     $            ('this looks like a permutation of the unit matrix')  3d15s24
              write(6,*)('perhaps the orbital sets correspond to too '),3d15s24
     $             ('different geometries?')                            3d15s24
              end if                                                     3d15s24
              if(igoal.eq.1.and.iforce.eq.0)then                         3d15s24
               write(6,*)('you can also force me to unscrable this by')  3d15s24
               write(6,*)('adding the keyword "force" to the input.')    3d15s24
              end if                                                     3d15s24
             end do
             irtrn=1                                                      2d15s22
             return                                                       2d15s22
            end if                                                        2d15s22
           end if                                                         2d8s20
           is=ie+1                                                        2d8s20
           if(is.lt.nbasdws(isb))go to 2257                               2d8s20
          end if                                                        6d24s24
         end if                                                         1d17s24
        end do                                                          2d8s20
        if(ifdiff.eq.0)then                                             6d24s24
         open(unit=1,file=bline200(1:nrawfn),form='unformatted')         2d9s20
         write(6,*)('saving diabatic orbitals to file forbs')            4d15s24
         open(unit=2,file='forbs',form='unformatted')                    2d9s20
         call copyto(iorb,iorbao,1,2,bc,ibc)                             11d9s22
        end if                                                          6d24s24
      irtrn=2                                                           2d6s23
      return
   15 continue
      write(6,*)('basis file does not exist: '),                        2d18s10
     $     bpath(1:ibpath)//line(is:ie)
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
   17 continue
      write(6,*)('basis not found for atom '),atom(i)
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1001 continue
      write(6,*)('end of file 1001')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1002 continue
      write(6,*)('end of file 1002')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1003 continue
      write(6,*)('end of file 1003')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1004 continue
      write(6,*)('end of file 1004')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1005 continue
      write(6,*)('end of file 1005')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1006 continue
      write(6,*)('end of file 1006')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1007 continue
      write(6,*)('end of file 1007')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1008 continue
      write(6,*)('end of file 1008')
      return                                                            5d13s05
 1009 continue
      write(6,*)('end of file 1009')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1010 continue
      write(6,*)('end of file 1010')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1011 continue
      write(6,*)('end of file 1011')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
 1012 continue
      write(6,*)('end of file 1012')
      irtrn=1                                                           5d13s05
      return                                                            5d13s05
      end
      integer function nbr(line,is,ie,iflag)
      character*(*) line                                                1d11s23
c
c     search line, starting from is to ie,
c     for a blank delimited string.
c     if end is \, return nbr=-1, otherwise
c     reset is to start of string and
c     nbr to end of string.
c     return nbr=-2 if no more strings
c     if iflag=0, # means end of string, rest of line is a comment.                 12d2s22
c
      ienew=ie                                                          12d2s22
      if(iflag.eq.0)then                                                12d2s22
       do i=ie,is,-1                                                     12d2s22
        if(line(i:i).eq.'#')then                                         12d2s22
         ienew=i-1                                                       12d2s22
         go to 3                                                         12d2s22
        end if                                                           12d2s22
       end do                                                            12d2s22
    3  continue                                                          12d2s22
      end if                                                            12d2s22
      do i=is,ienew                                                     12d2s22
       if(line(i:i).ne.' ')then
        i0=i
        go to 1
       end if
      end do
      nbr=-2
      return
    1 continue
      do i=i0+1,ienew                                                   12d2s22
       if(line(i:i).eq.' ')then
        nbr=i-1
        go to 2
       end if
      end do
      nbr=ienew                                                         12d2s22
    2 continue
      if(line(nbr:nbr).eq.'\')then
       if(i0.eq.nbr)then
        nbr=-1
       else
        nbr=nbr-1
       end if
      end if
      is=i0
      return
      end
