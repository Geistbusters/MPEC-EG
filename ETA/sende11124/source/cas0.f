c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cas0(ih0mo,ioooo,noc4,numpro,nowpro,potdws,              8d8s14
     $                icall,iden1,jdenpt,inewl,nlzz,multh,eavg,ncomp,   8d27s19
     $                enobreit,dynw,ioconv,icasvec,lprint,eavg2,        5d7s18
     $                jmats,kmats,nvirtc,isend,irecv,ipt,ipf,idwsdebx,   12d20s19
     $     ixlzz,islz,iptr,idorb,isorb,mdonp,mdoop,ibasisp,ncsfp,         8d2s22
     $ nfcnp,nctp,ipoint2p,icsfpd,ism,irel,iorbn,morb,nbasisp,nryd,ixw1,8d2s22
     $     ixw2,iptrbitp,iorbsym,iorbsymz,nextradata,extradata,mode,     5d16s22
     $     ih0d,i4od,idvec,idenergy,tolv,ipuse,isblkder,nsblkder,       6d10s22
     $     ifirst2,ipxder,isblkderd,nsblkderd,ipxderd,dwsumi,des,bc,ibc,3d20s23
     $     cdc,nveccnt,icanon,nvect)                                    10d11s24
c
c     mode=1 diagonalize and compute densities
c     mode=2 compute derivative of vectors and densities for buildcasgrad
c     mode=3 like mode=2, except also compute densities for bodcpart
c     mode=4 compute bi-linear derivative operator 2e density
c
      implicit real*8 (a-h,o-z)
      external second
      integer*8 iarg1,iarg2,iarg3,iarg4,ma1(2,8),ma2(2,8),mb1(2,8),
     $     mb2(2,8),ma1c(2,36,2),                                       8d29s06
     $     mtmp(2,8),mct(2,8),itrial,ipack8                             8d2s22
      integer*2 ipack2(4)                                               8d2s22
      equivalence (ipack8,ipack2)                                       8d2s22
      dimension m1c(36,2),multh(8,8),mst(8,8,2),mnz(3,64,36,2),nz(36,2) 4d30s18
      parameter (irzx=32,iczx=37)
      dimension irowzx(2,irzx),icolzx(2,iczx),pszx(2),jden1i(8)         3d17s23
      data icolzx/2,108, 2,127, 1,124, 2,172, 1,169, 1,223, 2,216,      <7
     $     1,88,
     $     2,81, 1,53, 1,57, 1,188, 1,192, 1,160, 2,165, 1,145, 2,150,   <17
     $     1,92, 2,104, 1,62, 2,74, 2,16, 2,24, 1,26, 2,20, 1,1, 1,9,   <27
     $     2,11, 1,5, 1,31, 1,39, 2,41, 1,35, 1,196, 1,204, 2,206,
     $     1,200/
      data irowzx/1,130, 2,135, 1,175, 2,180, 1,212, 2,224, 1,77, 2,89, <8
     $     1,46, 1,54, 2,56, 1,50, 1,181, 1,189, 2,191, 1,185, 2,157,   <17
     $     1,154, 2,142, 1,139, 1,103, 2,96, 1,73, 2,66, 2,23, 2,27,    <26
     $     1,8, 1,12, 1,38, 1,42, 1,203, 1,207/
      data pszx/1d0,-1d0/
      data loop,loopx/0,1000000/
      integer*2 itrial2(4)                                              8d22s06
      logical bcasvec,lprint,lbail,lwrite,lprt,ldynw,ldoit              3d31s23
      character*64 ostng                                                6d30s06
      character*50 numbers                                              3d5s21
      character*2 spinm,spinroo                                         5d3s21
      equivalence (itrial,itrial2)                                      8d22s06
      dimension idata1(8),idata2(8),idatac(36),idatb1(8),idatb2(8),      8d12s14
     $     idatbc(36),igt(8),iga(8),iveco(8),igs(8),igo(8),ivecb(8),    5d7s18
     $     jmats(*),kmats(*),nvirtc(*),isend(*),irecv(*),ipt(*),ipf(*), 12d20s19
     $     ixlzz(8,*),ixlzze(8,6),ibasis(1),ncsf(1),nfcn(1),des(*),     7d28s22
     $     nct(1),icsfpd(*),ism(*),irel(*),islz(*),iorbn(*),morb(*),    12d31s19
     $     nbasisp(*),ihryd(8),nryd(8),iorbsym(*),iorbsymz(*),cdc(*),   3d20s23
     $     extradata(*),ivecfp(8),dynw(3),ivfixed(8),i4od(*),potdws(*), 6d10s22
     $     isblkder(4,*),ipxder(4,8,8,8),isblkderd(4,*),ipxderd(4,8,8,8)7d1s22
     $     ,nveccnt(*),icanon(*)                                        5d5s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      10d17s14
     $     mynnode                                                      10d17s14
      common/singcm/iuse,nff
c
c     further improvements:
c     rather than gsum to get full vec and g then generate vect and gt,
c     use gather. Also this would work for ab doubles as well.
c     then we should never need to store full anything on a single proc.
c     however, if we can hold full vectors on a single proc, Schmidt
c     orthogonalization could be done more efficiently.
c     also don't add extra vector for vectors that are converged.
c
c     cas calculation.
c     ih0mo has 1e hamiltonian ints
c     ioooo has 2e integrals
c     strategy: divide ham into p-space and q-space and explicitly
c     diagonalize in spin adapted (and perhaps Lamda^2 adpated as well)
c     p-space. q-space will be handled by Davidson diagonalization.
c     for parallelization, the trial vector c and g=hc will be
c     distributed twice: once with alpha distributed and beta not
c     and once with beta distributed and alpha not.
c     we also need a full copy of c for alpha-beta doubles              3d10s17
c     This way the only
c     communication required to build the Davidson matrix will be
c     a global sum.
c     not so. need actual full g=hc when computing updates.
c     if dynw is ne zero, use dynamic weights.                          4d6s18
c
      include "common.store"
      include "common.hf"
      include "common.cas"
      include "common.print"
      dimension noc4(8),sums(2),ih0e(8),i2e(1),m12sym(8,2,2),idva(8)    5d27s22
      dimension ihc(8),ilc(8),nherec(8),nsbeta(8),ncsym(8),ivec(8),     8d29s06
     $     ivecx(8),ig(8),ivecn(8),iden1(8,*),jdenpt(1),iooooa(idbk),   6d2s22
     $     nherect(8),ihct(8),ilct(8),ivect(8),iveca(8),ih0ed(8),       5d16s22
     $     iooood(idbk),m12symnon(8,2,2),idata1non(8),idata2non(8),     6d13s22
     $     idatb1non(8),idatb2non(8),nsbeta2(8)                         7d20s22
      dimension dbg(11),identest1(8),identest2(512)                     7d20s22
      common/timerocm/tovr,telapo(15)                                   4d26s18
      data iseed/1598526581/                                            4d20s18
      data ikall/0/
      save                                                              3d17s17
      tolbest=tolv                                                      4d27s23
      if(nsymb.eq.4)then
       igoal=1835414
      else
       igoal=2372745
      end if
      igoal=iden1(1,1)+12*3
      if(mode.eq.3)then
      end if
      if(mode.eq.4)then                                                 3d17s23
       do isb=1,nsymb                                                   3d17s23
        jden1i(isb)=iden1(isb,2)                                        3d17s23
       end do                                                           3d17s23
      end if                                                            3d17s23
      ldynw=(idynw.eq.0.and.dynw(1).ne.0d0).or.idynw.eq.1.or.           9d20s21
     $     (idynw.eq.2.and.ioconv.ne.0)                                 3d17s21
      ikall=ikall+1
      dwsumi=0d0                                                        1d17s23
      if(iprtr(8).eq.0)then                                             1d3s20
       lprt=.false.                                                     1d3s20
       idwsdeb=0                                                        5d6s21
      else                                                              1d3s20
       lprt=.true.                                                      1d3s20
       write(6,*)('entering cas0')
       write(6,*)('cas0 for kall '),ikall,nrydb,mode
       write(6,*)('idvec '),idvec
       idwsdeb=9990                                                        5d6s21
      end if                                                            1d3s20
      xnan=3.14d0
      if(nrydb.gt.0.and.mode.eq.1)then                                  5d16s22
       if(lprint)then                                                   1d18s20
        write(6,*)('making space for rydberg hamiltonians ')            1d18s20
       end if                                                           1d18s20
       do isb=1,nsymb                                                   1d18s20
        ihryd(isb)=ibcoff                                               1d18s20
        nvirtx=nbasdws(isb)-idoub(isb)-iacto(isb)                        1d18s20
        ntry=(nvirtx*(nvirtx+1))/2                                        1d18s20
        ibcoff=ibcoff+ntry                                              1d18s20
       end do                                                           1d18s20
       call enough('cas0.  1',bc,ibc)
       do i=ihryd(1),ibcoff-1                                           1d18s20
        bc(i)=0d0                                                       1d18s20
       end do                                                           1d18s20
      end if                                                            1d18s20
      bcasvec=icasvec.eq.0                                              4d18s18
      if(nusecsf.lt.0)then                                              3d17s23
       norb=0                                                           3d17s23
       do isb=1,nsymb                                                   3d17s23
        norb=norb+iacto(isb)                                            3d17s23
       end do                                                           3d17s23
       ismd=ibcoff                                                      3d17s23
       ireld=ismd+norb                                                  3d17s23
       ibcoff=ireld+norb                                                3d17s23
       call enough('cas0.ismd',bc,ibc)                                  3d17s23
       ii=0                                                             3d17s23
       do isb=1,nsymb                                                   3d17s23
        do i=1,iacto(isb)                                               3d17s23
         ibc(ismd+ii)=isb                                               3d17s23
         ibc(ireld+ii)=i                                                3d17s23
         ii=ii+1                                                        3d17s23
        end do                                                          3d17s23
       end do                                                           3d17s23
      end if                                                            3d17s23
      if(bcasvec.or.(mode.ge.2.and.mode.le.3).and.ifirst2.eq.0)then     7d7s22
       if(lprint.and.kall.eq.1)then                                     7d14s22
        write(6,*)('computing space required to store vectors ')        12d28s19
        if(nusecsf.lt.0)then                                            12d28s19
         write(6,*)('the FCI will use determinants')                    12d28s19
        else                                                            12d28s19
         if(nusecsf.ne.64)then                                          12d28s19
          write(6,*)('the CI will use CSFs with a maximum of '),         12d28s19
     $        nusecsf,(' open shells')                                  12d28s19
         else                                                           12d28s19
          write(6,*)('the FCI will use CSFs')                           12d28s19
         end if                                                         12d28s19
        end if                                                          12d28s19
       end if                                                           12d28s19
       ibctop=ibcoff                                                    4d18s18
       nvect=0                                                          4d18s18
       do i=1,nstate
        if(lprint.and.kall.eq.1)write(6,*)('for i = '),i                7d14s22
        nart=0                                                          5d19s22
        if(nusecsf.lt.0)then                                            12d28s19
         call second(time1)                                              4d26s18
         call cas1b(multh(ipuse,isinfo(1,i)),isinfo(2,i),isinfo(3,i),   6d10s22
     $        nalpha,nbeta,                                             6d10s22
     $       numpro,nowpro,ialpha,iaorb,ibeta,iborb,nconfa,              11d8s05
     $       nconfb,m1a,m2a,m1b,m2b,ma1,ma2,mb1,mb2,icall,              8d8s06
     $       m12sym,m1c,ma1c,multh,idata1,idata2,idatac,idatb1,idatb2,  8d12s14
     $       idatbc,lprint,mst,mnz,nz,0,bc,ibc)                         11d9s22
         call second(time2)                                              4d26s18
         telap=time2-time1-tovr                                          4d26s18
         telapo(6)=telapo(6)+telap                                       4d26s18
         navec=0                                                         4d18s18
         do isb=1,nsymb                                                  4d18s18
          jsb=multh(isb,multh(ipuse,isinfo(1,i)))                       6d10s22
          call ilimts(nadet(isb),1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s, 4d18s18
     $        i2e)                                                      4d18s18
          nhere=ih+1-il                                                  4d18s18
          navec=navec+nhere*nbdet(jsb)                                   4d18s18
         end do                                                          4d18s18
         navec=navec*isinfo(4,i)                                        4d20s21
        else                                                            12d28s19
         iuniq=ibc(ipoint2p+i-1)                                        8d2s22
         jct=nctp+multh(isinfo(1,i),ipuse)-1+iuniq*nsymb                8d2s22
         call ilimts(ibc(jct),isinfo(4,i),mynprocg,                     8d2s22
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          6d10s22
         navec=ih+1-il                                                  4d18s18
        end if                                                          12d28s19
        nart=nart+isinfo(4,i)                                           5d19s22
        if(lprint.and.kall.eq.-1)                                        7d14s22
     $       write(6,*)('we require '),navec,(' words ')                7d14s22
        nvect=nvect+navec                                               4d20s21
        ibcoff=ibctop                                                   4d18s18
       end do                                                           4d18s18
       if(lprint.and.kall.eq.1)write(6,*)                               7d14s22
     $      ('total no. of words to store vectors '),                   7d14s22
     $      nvect,ibcoff                                                       4d20s18
       if(mode.eq.1)then                                                5d19s22
        icasvec=ibcoff                                                   4d18s18
        ibcoff=icasvec+nvect                                             4d18s18
        call enough('cas0.  2',bc,ibc)
       else if(mode.le.3)then                                           7d7s22
        idvec=ibcoff                                                    5d19s22
        idenergy=idvec+nvect                                            5d19s22
        ibcoff=idenergy+nart                                            5d19s22
        call enough('cas0.  3',bc,ibc)
        do iz=idvec,ibcoff-1                                            5d19s22
         bc(iz)=0d0                                                     5d19s22
        end do                                                          5d19s22
        ifirst2=1                                                       5d19s22
       end if                                                           5d19s22
      end if                                                            4d18s18
      call second(timea)                                                4d23s18
      ibcoffo=ibcoff
      nroodat=0                                                         5d3s21
      mden=0                                                            3d14s23
      do i=1,nstate                                                     5d3s21
       nroodat=nroodat+isinfo(4,i)                                      5d3s21
       mden=mden+((isinfo(4,i)*(isinfo(4,i)+1))/2)                      3d14s23
      end do                                                            5d3s21
      iroodat=ibcoff                                                    5d3s21
      iroolab=iroodat+nroodat                                           5d3s21
      ibcoff=iroolab+5                                                  5d3s21
      call enough('cas0.  4',bc,ibc)
      jroolab=iroolab-6                                                 5d3s21
      spinroo='  '                                                      5d3s21
      nroolab=0                                                         5d3s21
      jroodat=iroodat                                                   5d3s21
      lwrite=.false.                                                    12d31s19
      if((icall.eq.1.or.ioconv.ne.0).and.lprint)then                    4d20s18
       lwrite=.true.                                                    12d31s19
       write(6,*)('in cas0 '),ipuse
       write(6,*)('number of symmetry/spin blocks = '),nstate
         write(6,3315)(idoub(i),i=1,nsymb)                               8d8s14
 3315    format('number of doubly occupied orbitals:',8i3)               8d8s14
         write(6,3325)(iacto(i),i=1,nsymb)                               8d8s14
 3325    format('number of active orbitals:         ',8i3)               8d8s14
      end if
      call dws_sync                                                     8d8s14
      ncdx=0                                                            11d9s06
      norb=0                                                            8d14s06
      nbodc=1                                                           6d2s22
      if(mode.eq.3)nbodc=4                                              6d2s22
      if(mode.ne.4)then                                                 7d7s22
       do isb=1,nsymb                                                    8d7s06
        nmlt=1                                                           3d15s23
        if(nbodc.ne.1)nmlt=mden                                         3d15s23
        norb=norb+iacto(isb)                                             8d14s06
        do i=0,iacto(isb)*iacto(isb)*nmlt-1                             3d23s23
         bc(iden1(isb,nbodc)+i)=0d0                                      6d15s22
        end do                                                           6d15s22
        nmlt=1                                                          3d14s23
        do im=1,nbodc-1                                                  6d15s22
         do i=1,iacto(isb)*iacto(multh(ipuse,isb))*nmlt                 3d14s23
          bc(iden1(isb,im)+i-1)=0d0                                      6d2s22
         end do                                                          6d2s22
         nmlt=mden                                                      3d14s23
        end do                                                           2d28s07
       end do                                                            8d7s06
      else                                                              7d8s22
       do isb=1,nsymb                                                   7d8s22
        norb=norb+iacto(isb)                                            7d8s22
       end do                                                           7d8s22
      end if                                                            7d7s22
      if(ipuse.eq.1.and.mode.ne.4)then                                  7d7s22
       do idws=1,nsdlk                                                      2d28s07
        n1=isblk(1,idws)                                                 8d19s14
        n2=isblk(2,idws)                                                 8d19s14
        n3=isblk(3,idws)
        n4=isblk(4,idws)                                                 8d19s14
        if(n1.eq.n2)then
         nnn=(iacto(n1)*(iacto(n1)+1))/2
        else
         nnn=iacto(n1)*iacto(n2)
        end if
        if(n3.eq.n4)then
         mmm=(iacto(n3)*(iacto(n3)+1))/2
        else
         mmm=iacto(n3)*iacto(n4)
        end if
        if(nnn*mmm.gt.0)then
         do k=0,nnn*mmm-1
          bc(jdenpt(idws)+k)=0d0
         end do
        end if
       end do                                                            2d28s07
      else                                                              6d10s22
       do idws=1,nsblkderd                                              7d1s22
        n1=isblkderd(1,idws)                                            7d1s22
        n2=isblkderd(2,idws)                                            7d1s22
        n3=isblkderd(3,idws)                                            7d1s22
        n4=isblkderd(4,idws)                                            7d1s22
        nnn=iacto(n1)*iacto(n2)                                         6d10s22
        mmm=iacto(n3)*iacto(n4)                                         6d10s22
        if(nnn*mmm.gt.0)then                                            6d10s22
         do k=0,nnn*mmm-1                                               6d10s22
          bc(jdenpt(idws)+k)=0d0                                        6d10s22
         end do                                                         6d10s22
        end if                                                          6d10s22
       end do                                                           6d10s22
      end if                                                            6d10s22
      if(mode.ge.2)then                                                 6d2s22
       if(ipuse.eq.1)then                                               6d10s22
        call core(ih0d,dshift,i4od,noc4,ih0ed,dshift1,iooood,ncomp,      5d16s22
     $     idwsdeb,.false.,bc,ibc)                                      11d14s22
        dshift=dshift+potdws(2)                                          5d16s22
       else                                                             6d10s22
        call corenota1(ipuse,ih0d,dshift,i4od,noc4,ih0ed,dshift1,       6d10s22
     $       iooood,isblkder,nsblkder,multh,idwsdeb,.false.,bc,ibc)     11d14s22
       end if                                                           6d10s22
      end if                                                            5d16s22
      call core(ih0mo,shift,ioooo,noc4,ih0e,shift1,iooooa,ncomp,idwsdeb,6d22s22
     $     .false.,bc,ibc)                                              11d14s22
      if(nlzz.ne.0)then                                                 12d22s19
       call corelz2(ixlzz,shiftlz2,noc4,ixlzze,islz,multh,nlzz,bc,ibc)  11d14s22
      end if                                                            12d22s19
      if(lprt)write(6,*)('after cores ')
      shift=shift+potdws(1)                                             5d16s22
      shift0=0d0                                                        5d17s22
      pthrs=casdat(1)
      econv=casdat(2)
      vconv=casdat(5)                                                   12d9s22
      maxit=int(casdat(4)+0.1d0)
      if((icall.eq.1.or.ioconv.ne.0).and.lprint)then                    4d20s18
       write(6,447)econv,vconv,casdat(3),maxit                          12d22s22
  447  format('energy convergence criterion: ',es9.2,/,                 12d2s22
     $      'vector convergence criterion: ',es9.2,/,                   12d9s22
     $      'gradient convergence criterion: ',es9.2,/,                 12d22s22
     $      'max. no. iters for fci diagonalization: ',i5)              12d2s22
      end if
      if(mode.eq.1)eavg=0d0                                             5d16s22
      eavg2=0d0                                                         5d1s18
       wsum=0d0
       do i=1,nstate
        if(icall.eq.1.and.lprint)then                                   5d1s18
         write(6,1132)isinfo(2,i),isinfo(1,i),pspacex(i),maxps          4d15s19
 1132   format('spin multiplicity: ',i2,', space symmetry: ',i1,        4d15s19
     $        ', p-space threshold = ',f10.2,                           4d15s19
     $        ', max no. of p-space fcns ',i8)                          4d15s19
        end if                                                          5d1s18
        do j=1,isinfo(4,i)
         wsum=wsum+eweight(j,i)                                         5d1s18
         if(icall.eq.1.and.lprint)write(6,1131)j,eweight(j,i)           5d1s18
 1131    format('root ',i2,' wgt ',es10.2)                              5d1s18
        end do
       end do
       wsumi=1d0/wsum
       do i=1,nstate                                                    5d1s18
        do j=1,isinfo(4,i)                                              5d1s18
         ewght(j,i)=eweight(j,i)*wsumi                                  5d11s21
        end do                                                          5d1s18
       end do                                                           5d1s18
       if(ldynw)then                                                    3d17s21
       wsum=0d0                                                         4d6s18
       dwsum=0d0                                                        7d20s22
      end if                                                            4d6s18
      jcasvec=icasvec                                                   4d18s18
      if(mode.ge.2)jdvec=idvec                                          6d2s22
      call second(timeb)                                                4d23s18
      telap=timeb-timea-tovr                                            4d23s18
      telapo(1)=telapo(1)+telap                                         4d23s18
      ides=1                                                            7d28s22
      mdenoff=0                                                         3d14s23
      icdc=1                                                            3d31s23
      do i=1,nstate
       call second(timea)                                               4d23s18
       ibctop=ibcoff                                                    4d18s18
       nroot=isinfo(4,i)
       if(nroot.gt.0)then                                               4d27s21
        iwgt=ibcoff
        ibcoff=iwgt+nroot
        if(mode.eq.2.or.mode.eq.3)ibcoff=ibcoff+nroot                   7d20s22
        do j=1,nroot
         bc(iwgt+j-1)=ewght(j,i)                                        5d11s21
        end do
        nlab=0                                                          12d28s19
        do j=6,11                                                       12d31s19
         if(isinfo(j,i).ne.0)nlab=j                                     12d28s19
        end do                                                          12d28s19
        if(nusecsf.ge.0)then                                            12d28s19
         jptr=iptr+(mdoo+1)*2*(isinfo(1,i)-1)                           12d28s19
         cprint=0.1d0                                                   12d28s19
         ioconvu=0                                                      12d28s19
         if((icall.eq.1.or.ioconv.ne.0).and.lprint)ioconvu=1            12d28s19
         ibcprecsf=ibcoff                                               12d28s19
         nlabu=nlab-5                                                   12d28s19
         if(mode.ge.2)then
          write(6,*)('oops, ders for csf basis not coded yet')
          call dws_synca
          call dws_finalize
          stop
         end if
         iuniq=ibc(ipoint2p+i-1)                                        8d2s22
         ipack8=ibc(ipoint2p-nstate+iuniq)                              8d2s22
         nec=ipack2(1)                                                  8d2s22
         mdon=ipack2(4)                                                 8d2s22
         mdoo=ipack2(3)                                                 8d2s22
         jbasis=ibc(ibasisp+isinfo(1,i)-1+iuniq*nsymb)                  8d2s22
         jcsf=ibc(ncsfp+iuniq)                                          8d2s22
         jfcn=nfcnp+isinfo(1,i)-1+iuniq*nsymb                           8d2s22
         jct=nctp+isinfo(1,i)-1+iuniq*nsymb                             8d2s22
         jxw1=ibc(ixw1+iuniq)                                                8d2s22
         jxw2=ibc(ixw2+iuniq)                                           8d2s22
         jptrbit=ibc(iptrbitp+iuniq)                                    8d2s22
         call intcsfcas(mdon,mdoo,ibc(jbasis),                          8d2s22
     $        ibc(jcsf),ibc(jfcn),ih0e,iooooa,shift0,                   8d2s22
     $        pspacex(i),nec,icsfpd,multh,isinfo(1,i),isinfo(4,i),      12d28s19
     $        cprint,lprint,ieigold,                                    12d28s19
     $        ivsave,maxps,ibc(jct),nvcul,nlzz,islz,ixlzze,             8d2s22
     $        ism,irel,iacto,idoub,norb,ioconvu,isinfo(6,i),nlabu,      12d28s19
     $      isinfo(2,i),isinfo(5,i),bcasvec,jcasvec,jxw1,jxw2,          8d2s22
     $        ibc(jptrbit),nroolab,spinroo,jroolab,jroodat,shift,nsymb, 8d25s22
     $        npdiag,bc,ibc)                                            11d9s22
         do i1=0,nroot-1                                                5d17s22
          bc(ieigold+i1)=bc(ieigold+i1)+shift                           5d17s22
         end do                                                         5d17s22
         if((i.eq.1.or.idynw.eq.1).and.ldynw)then                       3d17s21
          eref=bc(ieigold)                                              12d28s19
         end if                                                         12d28s19
         do i1=1,nroot                                                   5d7s18
          eavg=eavg+bc(ieigold+i1-1)*bc(iwgt+i1-1)                      12d28s19
         end do                                                          5d7s18
         if(ldynw)then                                                  3d17s21
          do i1=1,nroot                                                 12d28s19
           if(dynw(2).eq.0d0)then                                       9d20s21
            ex=exp(-dynw(1)*(bc(ieigold+i1-1)-eref))                     9d20s21
            exi=1d0/ex                                                   12d28s19
            wwww=(2d0/(ex+exi))**2                                      9d20s21
           else                                                         9d20s21
            wde=dynw(1)*(bc(ieigold+i1-1)-dynw(2))                      9d20s21
            wwww=4d0/(dynw(3)+exp(2d0*wde))                             9d20s21
           end if                                                       9d20s21
           bc(iwgt+i1-1)=eweight(i1,i)*max(1d-7,wwww)                   9d20s21
           wsum=wsum+bc(iwgt+i1-1)                                      12d28s19
           eavg2=eavg2+bc(ieigold+i1-1)*ewght(i1,i)                     5d11s21
           ewght(i1,i)=bc(iwgt+i1-1)                                    5d11s21
          end do                                                        12d28s19
         end if                                                         12d28s19
         call hccsfd(bc(ivsave),ibc(jct),ibc(jbasis),ibc(jcsf),idum,    8d2s22
     $        ibc(jfcn),                                                8d2s22
     $        iden1,jdenpt,isinfo(4,i),mdon,isorb,idorb,icsfpd,         12d28s19
     $        nec,ewght(1,i),0,ism,irel,iacto,jxw1,jxw2,ibc(jptrbit),   8d2s22
     $        mdoo,isinfo(1,i),norb,bc,ibc)                             11d10s22
         ibcoff=ibcprecsf                                               12d28s19
        else                                                            12d28s19
         do isb=1,nsymb                                                  5d11s18
          idata1(isb)=0                                                 1d18s23
          idatb1(isb)=0                                                 1d18s23
          idata2(isb)=0                                                 1d18s23
          idatb2(isb)=0                                                 1d18s23
          idatac(isb)=0                                                 1d18s23
          idatbc(isb)=0                                                 1d18s23
          do jsb=1,nsymb                                                 5d11s18
           mst(jsb,isb,1)=0                                              5d11s18
           mst(jsb,isb,2)=0                                              5d11s18
          end do                                                         5d11s18
         end do                                                          5d11s18
         do ii=1,36                                                      5d11s18
          nz(ii,1)=0                                                     5d11s18
          nz(ii,2)=0                                                     5d11s18
         end do                                                          5d11s18
         if(ifcio.ne.0)then
          call second(time1)                                              4d26s18
          call cas1b(isinfo(1,i),isinfo(2,i),isinfo(3,i),nalpha,nbeta,    11d8s05
     $       numpro,nowpro,ialpha,iaorb,ibeta,iborb,nconfa,              11d8s05
     $       nconfb,m1a,m2a,m1b,m2b,ma1,ma2,mb1,mb2,icall,              8d8s06
     $       m12sym,m1c,ma1c,multh,idata1,idata2,idatac,idatb1,idatb2,  8d12s14
     $       idatbc,lprint,mst,mnz,nz,1,bc,ibc)                         11d9s22
          if(ipuse.ne.1)then                                            6d10s22
           call cas1bnona1(ipuse,isinfo(1,i),isinfo(2,i),isinfo(3,i),   6d13s22
     $         nalpha,nbeta,numpro,nowpro,ialpha,iaorb,ibeta,iborb,     6d13s22
     $         nconfa,nconfb,m12symnon,multh,idata1non,idata2non,       6d13s22
     $         idatb1non,idatb2non,lprint,bc,ibc)                       11d9s22
          end if                                                        6d10s22
          call second(time2)                                              4d26s18
          telap=time2-time1-tovr                                          4d26s18
          telapo(6)=telapo(6)+telap                                       4d26s18
          call dws_sync                                                    12d13s06
          nconf=0                                                         8d8s06
          nconf2=0                                                      12d23s22
          ncont=0                                                         3d10s17
          do isb=1,nsymb                                                12d23s22
           nsbeta(isb)=multh(isb,isinfo(1,i))                           12d23s22
           nsbeta2(isb)=multh(ipuse,nsbeta(isb))                        12d23s22
          end do                                                        12d23s22
          do isb=1,nsymb                                                  8d8s06
           call ilimts(nadet(isb),1,numpro,nowpro,ilc(isb),ihc(isb),i1,
     $         i2,i3,i4)                                                    8d8s06
           nherec(isb)=ihc(isb)+1-ilc(isb)                                8d8s06
           call ilimts(nbdet(isb),1,numpro,nowpro,ilct(isb),ihct(isb),
     $          i1,i2,i3,i4)                                                 3d10s17
           nherect(isb)=ihct(isb)+1-ilct(isb)                             3d10s17
           ncsym(isb)=nconf                                               8d14s06
           nconf=nconf+nbdet(nsbeta(isb))*nherec(isb)                     8d8s06
           nconf2=nconf2+nbdet(nsbeta2(isb))*nherec(isb)                12d23s22
           ncont=ncont+nadet(nsbeta(isb))*nherect(isb)                    3d10s17
          end do                                                          8d8s06
          call addcomma4(nconf,numbers,0)                                 3d5s21
 2252     format(' number of dets on this processor: ',a50)               3d10s21
          xconf=dfloat(nconf)                                             10d10s14
          call dws_gsumf(xconf,1)                                         10d10s14
          nconft=nint(xconf)                                              10d10s14
          call addcomma4(nconft,numbers,0)                                3d5s21
          if(lprint.and.icall.eq.1)write(6,2251)numbers                 10d22s24
 2251     format(' total number of dets across procs: ',a50)              3d10s21
          numa=nconfa
          numb=nconfb
         else                                                             11d8s05
          write(6,*)('cas1 ')
          if(bc(132).ne.-132)then
           call dws_sync
           call dws_finalize
           stop
          end if
          call cas1(isinfo(1,i),isinfo(2,i),isinfo(3,i),nalpha,
     $      numa,iaorb,na1p,na1m,iaa1,iaa1m,na2p,na2m,iaa2,
     $      iaa2m,nbeta,numb,iborb,nb1p,nb1m,ibb1,ibb1m,
     $      nb2p,nb2m,ibb2,ibb2m,nconf,bc,ibc)                          11d9s22
          ilc=1                                                            11d8s05
          ihc=numb                                                         11d8s05
         end if                                                           11d8s05
         ihdig=ibcoff
         ibcoff=ihdig+nconf
         iihdig=ibcoff                                                    8d22s06
         ibcoff=iihdig+max(nconf2,nconf)                                12d23s22
         isdig=ibcoff                                                     3d22s17
         ibcoff=isdig+nconf                                               3d22s17
         call enough('cas0.  5',bc,ibc)
         call second(time1)                                               4d26s18
         if(lprt)write(6,*)('mode '),mode
         if(mode.ge.2)then                                              6d15s22
          idhdig=ibcoff                                                 5d16s22
          ibcoff=idhdig+nconf2                                          12d23s22
          call enough('cas0.  6',bc,ibc)
          if(ipuse.eq.1)then                                            6d15s22
           call diaghii(bc(idhdig),nconf,ibc(iaorb),nalpha,numa,         5d16s22
     $     ibc(iborb),nbeta,numb,ih0ed,iooood,0d0,ilc,ihc,1d0,          9d12s22
     $      nsymb,nsbeta,ibc(iihdig),bc,ibc)                            11d14s22
          else                                                          6d15s22
           do isb=1,nsymb                                               6d15s22
            nsbeta2(isb)=multh(ipuse,nsbeta(isb))                       6d15s22
           end do                                                       6d15s22
           call diaghii(bc(idhdig),nconf,ibc(iaorb),nalpha,numa,         5d16s22
     $     ibc(iborb),nbeta,numb,ih0e,iooooa,shift0,ilc,ihc,1d0,        6d15s22
     $      nsymb,nsbeta2,ibc(iihdig),bc,ibc)                           11d14s22
         if(lprt)then
          write(6,*)('after diag1 ')
         end if
          end if                                                        6d15s22
         end if
         call diaghii(bc(ihdig),nconf,ibc(iaorb),nalpha,numa,
     $     ibc(iborb),nbeta,numb,ih0e,iooooa,shift0,ilc,ihc,1d0,        5d17s22
     $      nsymb,nsbeta,ibc(iihdig),bc,ibc)                            11d14s22
         if(lprt)then
          write(6,*)('after diag ')
         end if
         ieconfig=ibcoff                                                8d3s22
         iconfig=ieconfig+nconft                                        8d3s22
         nclodet=iconfig+nconft                                         8d3s22
         ibcoff=nclodet+nconft                                          8d3s22
         call enough('cas0.  7',bc,ibc)
         mdoo=min(nbeta,nalpha)                                         8d3s22
         call fiddleh(bc(ihdig),nconft,ibc(iaorb),nalpha,numa,          8d3s22
     $        ibc(iborb),nbeta,numb,ilc,ihc,nsymb,nsbeta,dum,           8d3s22
     $        ibc(iconfig),bc(ieconfig),norb,mdoo,nadet,nbdet,          8d3s22
     $        ibc(nclodet),nlzz,ixlzze,islz,multh,iacto,bc,ibc)         11d14s22
         if(lprt)write(6,*)('after fiddleh ')
         ibcoff=nclodet                                                 8d3s22
         if(nlzz.ne.0)then                                                12d22s19
          ilzzdig=ibcoff                                                  12d22s19
          ibcoff=ilzzdig+nconf                                            12d22s19
          call enough('cas0.  8',bc,ibc)
          call diaglzz(bc(ilzzdig),ncoef,ibc(iaorb),nalpha,numa,          12d22s19
     $      ibc(iborb),nbeta,numb,ixlzze,nsymb,nsbeta,islz,multh,       12d22s19
     $       ilc,ihc,nlzz,bc,ibc)                                       11d14s22
         end if                                                           12d22s19
         call second(time2)                                               4d26s18
         telap=time2-time1-tovr                                           4d26s18
         telapo(7)=telapo(7)+telap                                        4d26s18
         call spindig(bc(isdig),ibc(iaorb),nalpha,ibc(iborb),nbeta,
     $        nsymb,nsbeta,nherec,ilc,bc,ibc)                           11d9s22
         if(lprt)write(6,*)('after spindig ')
         call second(time3)                                               4d26s18
         telap=time3-time2-tovr                                           4d26s18
         telapo(8)=telapo(8)+telap                                        4d26s18
         elow=bc(ieconfig)                                              8d3s22
         ehigh=elow                                                     8d3s22
         do j=1,nconft-1                                                8d3s22
          elow=min(elow,bc(ieconfig+j))                                 8d3s22
          ehigh=max(ehigh,bc(ieconfig+j))                               8d3s22
         end do                                                         8d3s22
         npsloop=0                                                        5d1s18
         ibcps=ibcoff                                                     5d1s18
 2222    continue                                                         5d1s18
         pthrs=pspacex(i)                                                 5d1s18
         npsloop=npsloop+1                                                5d1s18
         etop=elow+pthrs
         nlenth=nconf*nroot                                               4d21s06
         ncona=0                                                          3d10s17
         do isb=1,nsymb                                                   8d29s06
          ivec(isb)=ibcoff                                                8d29s06
          jsb=nsbeta(isb)                                                 3d10s17
          ivect(isb)=ivec(isb)+2*nherec(isb)*nbdet(jsb)*nroot             3d27s17
          ibcoff=ivect(isb)+2*nherect(isb)*nadet(jsb)*nroot               3d27s17
          ivecb(isb)=ibcoff                                               3d17s17
          ibcoff=ivecb(isb)+nherec(isb)*nbdet(jsb)*nroot                  3d27s17
          iveco(isb)=ibcoff                                               3d27s17
          ibcoff=iveco(isb)+nherec(isb)*nbdet(jsb)*nroot                  3d27s17
          ncona=ncona+nadet(isb)*nbdet(jsb)                               3d10s17
         end do                                                           8d29s06
         nveca=0                                                          4d18s18
         do isb=1,nsymb                                                   3d14s17
          jsb=nsbeta(isb)                                                 3d14s17
          iveca(isb)=ibcoff                                               3d14s17
          ibcoff=iveca(isb)+2*nadet(isb)*nbdet(jsb)*nroot                 3d27s17
          nveca=nveca+nadet(isb)*nbdet(jsb)*nroot                         4d18s18
         end do                                                           3d14s17
         do j=0,(nconf*2+ncont+ncona)*nroot*2-1                           3d27s17
          bc(ivec(1)+j)=0d0                                               8d29s06
         end do                                                           4d26s06
         nps=0                                                          8d3s22
         ips=ibcoff                                                     8d3s22
         do ii=0,nconft-1                                               8d3s22
          if(bc(ieconfig+ii).le.etop)then                               8d3s22
           ibc(ips+nps)=ibc(iconfig+ii)                                 8d3s22
           nps=nps+1                                                    8d3s22
          end if                                                        8d3s22
         end do                                                         8d3s22
         if(nps.gt.maxps)then                                           4d15s19
          delta=(etop-elow)/2.1d1                                        4d15s19
          deltai=1d0/delta                                               4d15s19
          igotps=ibcoff                                                  4d15s19
          ibcoff=igotps+21                                               4d15s19
          call enough('cas0.  9',bc,ibc)
          do j=0,20                                                      4d15s19
           ibc(igotps+j)=0                                               4d15s19
          end do                                                         4d15s19
          do j=0,nconft-1                                                 4d15s19
           val=deltai*(bc(ieconfig+j)-elow)                             8d3s22
           ival=val                                                      4d15s19
           ibc(igotps+ival)=ibc(igotps+ival)+1                           4d15s19
          end do                                                         4d15s19
          isum=0                                                         4d15s19
          puse=0d0                                                       4d15s19
          do j=0,20                                                      4d15s19
           nhere=ibc(igotps+j)                                          8d3s22
           isum=isum+nhere                                               4d15s19
           puse=puse+delta                                               4d15s19
           if(isum.gt.maxps)then                                         4d15s19
            if(lprint)then                                                4d15s19
             write(6,*)
     $            ('too many p-space fcns. Cut back threshold from '),  4d15s19
     $       pthrs,j                                                      4d15s19
            end if                                                         4d15s19
            if(lprint)write(6,*)('use p-space threshold '),puse-delta    4d15s19
            pspacex(i)=puse-delta                                        4d15s19
            ibcoff=ibcps                                                 10d15s20
            go to 2222                                                   4d15s19
           end if                                                        4d15s19
          end do                                                         4d15s19
         end if                                                          4d15s19
         ibcoff=ips+nps                                                 8d3s22
         if(mode.eq.4)then                                              3d17s23
          if(nstate.ge.i)then                                           3d20s23
           kcasvec=jcasvec                                              3d17s23
           lcasvec=kcasvec                                              3d17s23
           do ix=i,nstate                                               3d17s23
            kcasvec=nveccnt(ix)                                         4d18s23
            lcasvec=kcasvec                                             4d18s23
            if(isinfo(2,i).eq.isinfo(2,ix).and.                         4d18s23
     $            isinfo(3,i).eq.isinfo(3,ix))then                      3d31s23
             nxx=0                                                       3d17s23
             nyy=0
             do isb=1,nsymb                                              3d17s23
              jsb=multh(isb,isinfo(1,ix))                                3d17s23
              nxx=nxx+nadet(isb)*nbdet(jsb)                              3d17s23
              nyy=nyy+nbdet(jsb)*nherec(isb)
             end do                                                      3d17s23
             ivecxx=ibcoff
             ibcoff=ivecxx+nxx*isinfo(4,ix)
             call enough('cas0.vecxx',bc,ibc)
             do iz=ivecxx,ibcoff-1                                       3d17s23
              bc(iz)=0d0                                                 3d17s23
             end do                                                      3d17s23
             jvecxx=ivecxx                                               3d17s23
             do ir=1,isinfo(4,ix)                                       3d23s23
              jvecget=kcasvec-1                                              3d17s23
              do isb=1,nsymb                                              3d17s23
               jsb=multh(isb,isinfo(1,ix))                                3d17s23
               jad=jvecxx+nbdet(jsb)*(ilc(isb)-1)                        3d17s23
               do ia=ilc(isb),ihc(isb)                                   3d17s23
                do ibr=0,nbdet(jsb)-1                                    3d17s23
                 iqad=jvecget+ir+isinfo(4,ix)*(ia-ilc(isb)              3d23s23
     $                +nherec(isb)*ibr)                                 3d23s23
                 bc(jad+ibr)=bc(iqad)                                    3d20s23
                end do                                                   3d17s23
                jad=jad+nbdet(jsb)                                       3d17s23
               end do                                                    3d17s23
               jvecget=jvecget+isinfo(4,ix)*nherec(isb)*nbdet(jsb)      3d23s23
               jvecxx=jvecxx+nbdet(jsb)*nadet(isb)                       3d17s23
               lcasvec=lcasvec+nbdet(jsb)*(ihc(isb)+1-ilc(isb))          3d17s23
              end do                                                     3d17s23
             end do                                                      3d17s23
             call dws_gsumf(bc(ivecxx),nxx*isinfo(4,ix))                 3d24s23
             if(ix.eq.i)then
              ivecxx0=ivecxx                                             3d17s23
              nxx0=nxx                                                   3d17s23
             else                                                        3d17s23
              call tdendet(bc(ivecxx0),nxx0,isinfo(1,i),bc(ivecxx),nxx,  3d17s23
     $            isinfo(1,ix),nadet,nbdet,iacto,jden1i,ibc(iaorb),     3d17s23
     $            nalpha,ibc(iborb),nbeta,multh,nsymb,norb,ibc(ismd),   3d17s23
     $            ibc(ireld),bc,ibc,igoal,ioffdeta,ioffdetb)            3d20s23
              ibcoff=ivecxx                                              3d17s23
             end if                                                      3d17s23
            end if                                                      4d18s23
           end do                                                       3d17s23
           ibcoff=ivecxx0                                               3d17s23
          end if                                                        3d17s23
         end if                                                         3d17s23
         itmphd=ibcoff                                                  8d3s22
         ibcoff=itmphd+nps                                              8d3s22
         if(mode.ge.2)then                                              8d3s22
          itmpdhd=ibcoff                                                8d3s22
          ibcoff=itmpdhd+nps                                            8d3s22
         end if                                                         8d3s22
         if(nlzz.ne.0)then                                              8d3s22
          itmplzd=ibcoff                                                8d3s22
          ibcoff=itmplzd+nps                                            8d3s22
         end if                                                         8d3s22
         call enough('cas0. 10',bc,ibc)
         do iz=itmphd,ibcoff-1                                          8d3s22
          bc(iz)=0d0                                                    8d3s22
         end do                                                         8d3s22
         nwds=ibcoff-itmphd                                             8d3s22
         npsh=0                                                         8d3s22
         ipsh=ibcoff                                                    8d3s22
         do j=0,nps-1                                                   8d3s22
          itrial=ibc(ips+j)                                             8d3s22
          if(itrial2(1).ge.ilc(itrial2(3)).and.
     $       itrial2(1).le.ihc(itrial2(3)))then                         8d3s22
           ibc(ipsh+npsh)=ibc(ips+j)                                    8d3s22
           npsh=npsh+1                                                  8d3s22
          end if                                                        8d3s22
         end do                                                         8d3s22
         ibcoff=ipsh+npsh                                               8d3s22
         isndpsh=ibcoff                                                 8d3s22
         ibcoff=isndpsh+mynprocg                                        8d3s22
         call enough('cas0. 11',bc,ibc)
         do iz=isndpsh,ibcoff-1                                         8d3s22
          bc(iz)=0d0                                                    8d3s22
         end do                                                         8d3s22
         xnpsh=dfloat(npsh)                                             8d3s22
         bc(isndpsh+mynowprog)=xnpsh                                    8d3s22
         call dws_gsumf(bc(isndpsh),mynprocg)                           8d3s22
         npshx=0                                                        8d3s22
         ioffps=0                                                       8d3s22
         do ip=0,mynprocg-1                                             8d3s22
          npshp=nint(bc(isndpsh+ip))                                    8d3s22
          ibc(isndpsh+ip)=npshp                                         8d3s22
          npshx=max(npshx,npshp)                                        8d3s22
          if(ip.lt.mynowprog)ioffps=ioffps+npshp                        8d3s22
         end do                                                         8d3s22
         isndps=ibcoff                                                  8d3s22
         ibcoff=isndps+npshx*mynprocg                                   8d3s22
         call enough('cas0. 12',bc,ibc)
         do iz=isndps,ibcoff-1                                          8d3s22
          bc(iz)=0d0                                                    8d3s22
         end do                                                         8d3s22
         jsndps=isndps+mynowprog*npshx                                  8d3s22
         do ii=0,npsh-1                                                  8d3s22
          ibc(jsndps+ii)=ibc(ipsh+ii)                                   8d3s22
         end do                                                         8d3s22
         call dws_gbor(bc(isndps),mynprocg*npshx)                       8d3s22
         jps=ips                                                        8d3s22
         do ip=0,mynprocg-1                                             8d3s22
          jsndps=isndps+ip*npshx                                        8d3s22
          npshp=ibc(isndpsh+ip)                                         8d3s22
          do ii=0,npshp-1                                               8d3s22
           itrial=ibc(jsndps+ii)                                        8d3s22
           ibc(jps)=itrial                                              8d3s22
           jps=jps+1                                                    8d3s22
          end do                                                        8d3s22
         end do                                                         8d3s22
         ibcoff=isndpsh                                                 8d3s22
         do j=0,nps-1                                                   8d3s22
          itrial=ibc(ips+j)                                             8d3s22
          if(itrial2(1).ge.ilc(itrial2(3)).and.
     $       itrial2(1).le.ihc(itrial2(3)))then                         8d3s22
           i2=0                                                            8d22s06
           do isb=1,itrial2(3)-1                                           8d22s06
            i2=i2+nbdet(nsbeta(isb))*nherec(isb)                           8d22s06
           end do                                                          8d22s06
           i2=ihdig+i2+itrial2(2)-1+nbdet(nsbeta(itrial2(3)))              8d22s06
     $          *(itrial2(1)-ilc(itrial2(3)))                              8d22s06
           if(mode.ge.2)then                                             6d2s22
            id2=i2-ihdig+idhdig                                          5d16s22
            bc(itmpdhd+j)=bc(id2)                                       8d3s22
           end if                                                        5d16s22
           bc(itmphd+j)=bc(i2)                                           8d11s14
           if(nlzz.ne.0)then                                               12d22s19
            i2=i2+ilzzdig-ihdig                                            12d22s19
            bc(itmplzd+j)=bc(i2)                                         12d22s19
           end if                                                          12d22s19
          end if                                                        8d3s22
         end do                                                         8d3s22
         call dws_gsumf(bc(itmphd),nwds)                                8d3s22
         ipham=ibcoff
         ibcoff=ipham+nps*nps
         itmpa=ibcoff
         itmpb=itmpa+norb                                                 8d14s06
         ibcoff=itmpb+norb                                                8d14s06
         itmpc=ibcoff
         itmpd=itmpc+norb                                                 8d14s06
         itmpe=itmpd+norb                                                 8d14s06
         itmpf=itmpe+norb                                                 8d14s06
         ibcoff=itmpf+norb                                                8d14s06
         call enough('cas0. 13',bc,ibc)
         na=(iacto(1)*(iacto(1)+1))/2
         call second(time1)                                               4d26s18
         call psspin(nps,ibc(ips),bc(itmphd),bc(ipham),numb,              11d10s05
     $      ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),
     $      ibc(itmpb),ibc(itmpc),ibc(itmpd),ibc(itmpe),
     $      ibc(itmpf),bc(ihdig),nconf,norb,ncsym,nsymb,nsbeta,bc,ibc)  11d9s22
         if(lprt)write(6,*)('after psspin ')
         ieigs=ibcoff
         ivecs=ieigs+nps
         ibcoff=ivecs+nps*nps
         isymz=ibcoff                                                     8d11s14
         ibcoff=isymz+nps                                                 8d11s14
         call enough('cas0. 14',bc,ibc)
c
c     diagx tries to block diagonalize the matrix, and for
c     s**2, there are a lot of very small blocks, so this
c     is not a expensive calculation
c
         call diagx(nps,bc(ipham),bc(ieigs),bc(ivecs),ibc(isymz),bc,ibc)11d14s22
         call second(time2)                                               4d26s18
         if(lprt)write(6,*)('after diagx ')
         telap=time2-time1-tovr                                           4d26s18
         telapo(10)=telapo(10)+telap                                      4d26s18
         call dws_bcast(bc(ivecs),nps*nps)                                12d23s15
         spin=0.5d0*(dfloat(isinfo(2,i))-1d0)
         spin2=spin*(spin+1d0)                                             5d3s07
         if((icall.eq.1.or.ioconv.ne.0).and.lprint)then                 8d16s22
          write(6,446)i,spin,spin2                                      12d2s22
  446     format('for i = ',i5,', desired spin ',f5.1,                  12d2s22
     $         ', so desired S^2 eigenvalue is ',f6.2)                  12d2s22
         end if                                                         8d16s22
         ikeep=ibcoff                                                     6d30s06
         ibcoff=ikeep+nps                                                 6d30s06
         call enough('cas0. 15',bc,ibc)
         nok=0
         do i1=1,nps                                                      10d23s14
          scalc2=sqrt(1d0+4d0*bc(ieigs+i1-1))
          mtry=int(scalc2+0.001d0)
          if(abs(scalc2-dfloat(mtry)).gt.1d-10)then                       12d11s06
           write(6,*)('s**2 eigenvalue is not pure!!! '),scalc2,mtry      12d11s06
           write(6,*)('eigenvalue no. '),i1
           stop                                                           12d11s06
          end if                                                          12d11s06
          if(abs(spin2-bc(ieigs+i1-1)).lt.1d-10)then                      3d14s17
           ibc(ikeep+nok)=i1                                              6d30s06
           nok=nok+1                                                      6d30s06
          end if
          iad1=ivecs-1+nps*(i1-1)                                         6d30s06
         end do
    1    continue
         if(icall.eq.1.and.lprint)write(6,*)('no. fcns kept = '),nok      4d20s18
         itmpvec=ibcoff                                                   6d30s06
         ibcoff=itmpvec+nps*nps                                           6d30s06
         call enough('cas0. 16',bc,ibc)
         do k=1,nok                                                       6d30s06
          iad1=itmpvec-1+nps*(k-1)
          iad2=ivecs-1+nps*(ibc(ikeep+k-1)-1)                             6d30s06
          do j=1,nps                                                      6d30s06
           bc(iad1+j)=bc(iad2+j)                                          6d30s06
          end do                                                          6d30s06
         end do                                                           6d30s06
         do k=1,nok*nps                                                   6d30s06
          bc(ivecs+k-1)=bc(itmpvec+k-1)                                   6d30s06
         end do                                                           6d30s06
         if(nlzz.gt.0)then                                                6d4s07
          ill=ibcoff                                                      6d4s07
          ibcoff=ill+nps*nps                                              6d4s07
          call enough('cas0. 17',bc,ibc)
          do izzs=1,nps*nps                                               6d25s07
           bc(ill+izzs-1)=0d0                                             6d25s07
          end do                                                          6d4s07
          call psl(nps,ibc(ips),bc(itmplzd),bc(ill),numb,ibc(iaorb),     12d22s19
     $        nalpha,ibc(iborb),nbeta,ibc(itmpa),ibc(itmpb),            12d22s19
     $        ixlzze,ibc(itmpc),ibc(itmpd),ibc(itmpe),ibc(itmpf),norb,  12d22s19
     $        nsymb,nsbeta,islz,multh,bc(ivecs),nok,isinfo(5,i),nlzz,   12d31s19
     $       lwrite,bc,ibc,i)                                             11d9s22
          ibcoff=ill                                                     12d31s19
         end if                                                           6d4s07
         ibcoff=ikeep                                                     6d30s06
         iarg1=nps
         iarg2=nok
         if(nok.lt.isinfo(4,i))then                                       5d1s18
          if(lprint)                                                      5d11s18
     $      write(6,*)('there are fewer p-space configs than roots ')   5d11s18
          if(npsloop.le.15)then                                           5d1s18
           pspacex(i)=pspacex(i)+0.05d0                                   5d1s18
           ibcoff=ibcps                                                   10d15s20
           go to 2222                                                     5d1s18
          end if                                                          5d1s18
          write(6,*)('current p-space threshold = '),pspacex(i)           5d1s18
          call dws_sync                                                   5d1s18
          call dws_finalize                                               5d1s18
          stop                                                            5d1s18
         end if                                                           5d1s18
         ivecst=ibcoff
         ibcoff=ivecst+nok*nps
         call enough('cas0. 18',bc,ibc)
         do i1=1,nps
          do i2=1,nok
           ia1=ivecst+i2-1+(i1-1)*nok
           ia2=ivecs+i1-1+(i2-1)*nps
           bc(ia1)=bc(ia2)
          end do
         end do
         iphamz=ibcoff
         ibcoff=iphamz+nps*nps
         ipcopy=ibcoff
         ibcoff=ipcopy+nps*nps
         call enough('cas0. 19',bc,ibc)
         do j=1,11                                                        8d19s14
          dbg(j)=1d0                                                      8d19s14
         end do                                                           8d19s14
         call distit(nps,nhere,istart)                                    10d23s14
         iphamr=ibcoff                                                    10d23s14
         ibcoff=iphamr+nps*nhere                                          10d23s14
         call enough('cas0. 20',bc,ibc)
         do j=1,11                                                        8d19s14
          dbg(j)=1d0                                                      8d19s14
         end do                                                           8d19s14
         call second(time1)                                               4d26s18
         call pshampp(nps,ibc(ips),bc(itmphd),bc(iphamr),numb,            10d23s14
     $      ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),
     $      ibc(itmpb),ih0e,iooooa,ibc(itmpc),                          8d15s14
     $      ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $      nsymb,nsbeta,dbg,nhere,istart,.false.,bc,ibc)               11d9s22
         if(lprt)write(6,*)('after pshampp ')
         call second(time2)                                               4d26s18
         telap=time2-time1-tovr                                           4d26s18
         telapo(11)=telapo(11)+telap                                      4d26s18
         if(mode.eq.1)then                                              5d16s22
          itmp=ibcoff                                                      10d23s14
          ibcoff=itmp+nps*nps                                              3d21s17
          call enough('cas0. 21',bc,ibc)
          do j=0,nps*nps-1                                                 3d21s17
           bc(itmp+j)=0d0                                                  10d23s14
          end do                                                           10d23s14
          jtmp=itmp+istart-1                                               10d23s14
          if(nhere.gt.0.and.nps.gt.0)then                                  2d22s19
           call dgemm('n','n',nhere,nok,nps,1d0,bc(iphamr),nhere,           10d23s14
     $         bc(ivecs),nps,0d0,bc(jtmp),nps,                             10d23s14
     d' cas0.  1')
          end if                                                           2d22s19
          call dws_gsumf(bc(itmp),nps*nps)                                 3d21s17
          iarg1=nok+nroot*2                                                3d27s17
          iarg1=iarg1*iarg1+iarg1                                          10d23s14
          iarg1=iarg1/2                                                    10d23s14
          call distit8(iarg1,iarg3,iarg2)                                  10d23s14
          itmptri=ibcoff                                                   10d23s14
          ibcoff=itmptri+iarg3                                             10d23s14
          itmptric=ibcoff                                                  3d17s17
          ibcoff=itmptric+iarg3                                            3d17s17
          call enough('cas0. 22',bc,ibc)
          iarg1=0                                                          10d23s14
          iarg2=itmptri                                                    10d23s14
          do j=0,nok*nok-1                                                 3d16s17
           bc(ipham+j)=0d0                                                 3d16s17
          end do                                                           3d16s17
          do j=0,nok-1                                                     10d23s14
           do k=0,j                                                        10d23s14
            if(mod(iarg1,mynprocg).eq.mynowprog)then                       10d23s14
             dot=0d0                                                       10d23s14
             iad1=ivecs+nps*k                                              10d23s14
             iad2=itmp+nps*j                                               10d23s14
             do l=0,nps-1                                                  10d23s14
              dot=dot+bc(iad1+l)*bc(iad2+l)                                10d23s14
             end do                                                        10d23s14
             bc(iarg2)=dot                                                 10d23s14
             iad1=ipham+j+nok*k                                            3d16s17
             bc(iad1)=dot                                                  3d16s17
             iad1=ipham+k+nok*j                                            3d16s17
             bc(iad1)=dot                                                  3d16s17
             iarg2=iarg2+1                                                 10d23s14
            end if                                                         10d23s14
            iarg1=iarg1+1                                                  10d23s14
           end do                                                          10d23s14
          end do                                                           10d23s14
          call dws_gsumf(bc(ipham),nok*nok)                                3d16s17
          nroottmp=nroot                                                   10d24s14
          do k=0,iarg3-1                                                   3d17s17
           bc(itmptric+k)=bc(itmptri+k)                                     3d17s17
          end do                                                           3d17s17
          ibcsav=ibcoff                                                    3d17s17
          if(idwsdeb.ne.0)write(6,*)('phouse the first in cas0')
          if(nok.le.npdiag(1))then                                      8d12s22
           call diagdy(bc(itmptri),nok,nroottmp,ieigh,ivech,bc,ibc)     11d14s22
           isteig=1                                                     8d12s22
          else                                                          8d12s22
           call phouse(itmptri,nok,1d0,0d0,nroottmp,ieigh,ivech,0,      8d12s22
     $          isteig,bc,ibc)                                          11d9s22
          end if
        if(lprt)write(6,*)('after diagonalization ')
          if(idwsdeb.ne.0)then                                             8d27s19
           write(6,*)('eigenvalues from phouse ')
           call prntm2(bc(ieigh),1,nroottmp,1)
          end if                                                           8d27s19
          if(nroottmp.gt.0)then
           sum=0d0
           do k=0,nok-1
            sum=sum+bc(ivech+k)**2
           end do
          end if
          if(nroottmp.gt.0.and.idwsdeb.ne.0)then                           8d27s19
           write(6,*)('my vectors ')
           call prntm2(bc(ivech),nok,nroottmp,nok)
          end if                                                           8d27s19
          isteig=isteig-1                                                  10d24s14
          ieall=ibcoff                                                     10d24s14
          ivall=ieall+nroot                                                10d24s14
          ibcoff=ivall+nroot*nok                                           10d24s14
          call enough('cas0. 23',bc,ibc)
          jeall=ieall+isteig                                               10d24s14
          jvall=ivall+nok*isteig                                           10d24s14
          do j=0,nroot*(nok+1)-1                                           10d24s14
           bc(ieall+j)=0d0                                                 10d24s14
          end do                                                           10d24s14
          jvech=ivech                                                      10d24s14
          do j=0,nroottmp-1                                                10d24s14
           bc(jeall+j)=bc(ieigh+j)+shift0                               7d20s22
           do k=0,nok-1                                                    10d24s14
            bc(jvall+k)=bc(jvech+k)                                        10d24s14
           end do                                                          10d24s14
           jvall=jvall+nok                                                 10d24s14
           jvech=jvech+nok                                                 10d24s14
          end do                                                           10d24s14
          if(nok.gt.npdiag(1))then                                      8d12s22
           call dws_gsumf(bc(ieall),nroot*(nok+1))                          10d24s14
          end if                                                        8d12s22
          ieigh=ieall                                                      10d24s14
          ivech=ivall                                                      10d24s14
           ieigp=ibcsav                                                     3d17s17
          ivecp=ieigp+nroot                                                10d24s14
          ibcoff=ivecp+nroot*nok                                           10d24s14
          do j=0,nroot*(nok+1)-1                                           10d24s14
           bc(ieigp+j)=bc(ieigh+j)                                         3d17s17
          end do                                                           10d24s14
          do k=0,nroot-1
           sum=0d0
           do l=0,nok-1
            sum=sum+bc(ivecp+l+nok*k)**2
           end do
          end do
         end if                                                         5d16s22
         iarg1=1                                                          4d26s07
         iarg2=nroot                                                      4d26s07
         if(nps.eq.nconft)then                                            10d10s14
         if(lprt)write(6,*)('pspace is everything ')
          if(mode.eq.1)then                                             5d16s22
           if(ldynw.and.(i.eq.1.or.idynw.eq.1))then                        3d17s21
            eref=bc(ieigp)                                                 5d7s18
           end if                                                          5d7s18
           if(ldynw)then                                                   3d17s21
            do i1=1,nroot                                                  5d7s18
             if(dynw(2).eq.0d0)then                                      9d20s21
              ex=exp(-dynw(1)*(bc(ieigp+i1-1)-eref))                     9d20s21
              exi=1d0/ex                                                    5d7s18
              wwww=(2d0/(ex+exi))**2                                     9d20s21
             else                                                        9d20s21
              wde=dynw(1)*(bc(ieigp+i1-1)+shift-dynw(2))                5d17s22
              wwww=4d0/(dynw(3)+exp(2d0*wde))                            9d20s21
             end if                                                      9d20s21
             bc(iwgt+i1-1)=eweight(i1,i)*max(1d-7,wwww)                  9d20s21
             wsum=wsum+bc(iwgt+i1-1)                                       5d7s18
             eavg2=eavg2+(bc(ieigp+i1-1)+shift)*ewght(i1,i)             5d17s22
             ewght(i1,i)=bc(iwgt+i1-1)                                   5d11s21
            end do                                                         5d7s18
           end if                                                          5d7s18
           do i1=1,nroot                                                   5d7s18
            eavg=eavg+(bc(ieigp+i1-1)+shift)*bc(iwgt+i1-1)              7d20s22
           end do                                                          5d7s18
          else                                                          5d16s22
           if(ipuse.eq.1)then                                           6d22s22
            ivecf=ibcoff                                                 5d16s22
            idvecd=ivecf+nps*nroot                                      7d21s22
            ibcoff=idvecd+nps*nroot                                     7d21s22
            nps2=nps                                                    6d22s22
            nhere2=nhere                                                6d22s22
            istart2=istart                                              6d22s22
           else                                                         6d22s22
            nps2=0                                                      6d22s22
            do isb=1,nsymb                                              6d22s22
             jsb=nsbeta2(isb)                                           6d22s22
             nps2=nps2+nadet(isb)*nbdet(jsb)                            6d22s22
            end do                                                      6d22s22
            ivecf=ibcoff                                                 5d16s22
            idvecd=ivecf+nps*nroot                                      7d21s22
            ibcoff=idvecd+nps2*nroot                                    7d21s22
            call distit(nps2,nhere2,istart2)                            6d22s22
           end if                                                       6d22s22
           call enough('cas0. 24',bc,ibc)
           do iz=ivecf,ibcoff-1                                         5d16s22
            bc(iz)=0d0                                                  5d16s22
           end do                                                       5d16s22
           iad1=ivecf                                                      5d15s19
           iad2=jcasvec                                                    5d15s19
           idad1=idvecd                                                 5d19s22
           idad2=jdvec                                                  5d19s22
           nxn=0                                                        4d3s23
           do isb=1,nsymb                                               4d3s23
            nxn=nxn+nherec(isb)*nbdet(nsbeta(isb))                      4d3s23
           end do                                                       4d3s23
           nxn2=0                                                        4d3s23
           do isb=1,nsymb                                               4d3s23
            nxn2=nxn2+nherec(isb)*nbdet(nsbeta2(isb))                   4d13s23
           end do                                                       4d3s23
           if(ipuse.eq.1)then                                           6d22s22
            if(mode.eq.4)then                                           11d23s22
             do isb=1,nsymb                                             4d3s23
              jsb=nsbeta(isb)                                           4d3s23
              do ibr=0,nbdet(jsb)-1                                     4d3s23
               do ia=ilc(isb),ihc(isb)                                  4d3s23
                jad1=iad1+ibr+nbdet(jsb)*(ia-1)                         4d3s23
                do ir=0,nroot-1                                         4d3s23
                 bc(jad1+nps*ir)=bc(iad2+ir)                            4d9s23
                end do                                                  4d3s23
                iad2=iad2+nroot                                         4d3s23
               end do                                                   4d3s23
              end do                                                    4d3s23
              iad1=iad1+nbdet(jsb)*nadet(isb)                           4d3s23
             end do                                                     4d3s23
            else                                                        11d23s22
             do isb=1,nsymb                                             4d3s23
              jsb=nsbeta(isb)                                           4d3s23
              do ibr=0,nbdet(jsb)-1                                     4d3s23
               do ia=ilc(isb),ihc(isb)                                  4d3s23
                jad1=iad1+ibr+nbdet(jsb)*(ia-1)                         4d3s23
                do ir=0,nroot-1                                         4d3s23
                 bc(jad1+nps*ir)=bc(iad2+ir)                            4d9s23
                end do                                                  4d3s23
                iad2=iad2+nroot                                         4d3s23
               end do                                                   4d3s23
              end do                                                    4d3s23
              iad1=iad1+nbdet(jsb)*nadet(isb)                           4d3s23
              jdad1=idad1+nbdet(jsb)*(ilc(isb)-1)                       4d3s23
              do ia=ilc(isb),ihc(isb)                                        5d15s19
               do ir=0,nroot-1                                          4d3s23
                do ibr=0,nbdet(jsb)-1                                    7d20s22
                 bc(jdad1+ibr+ir*nps)=bc(idad2+ibr+ir*nxn)              4d13s23
                end do                                                  4d3s23
               end do                                                        5d15s19
               jdad1=jdad1+nbdet(jsb)                                   7d20s22
               idad2=idad2+nbdet(jsb)                                   7d20s22
              end do
              idad1=idad1+nbdet(jsb)*nadet(isb)                         7d20s22
             end do                                                     4d3s23
            end if                                                      11d23s22
            call dws_gsumf(bc(ivecf),nroot*nps*2)                       7d20s22
           else                                                         6d22s22
            do isb=1,nsymb                                              4d3s23
             jsb=nsbeta(isb)                                            4d3s23
             do ibr=0,nbdet(jsb)-1                                      4d3s23
              do ia=ilc(isb),ihc(isb)                                   4d3s23
               jad1=iad1+ibr+nbdet(jsb)*(ia-1)                          4d3s23
               do ir=0,nroot-1                                          4d3s23
                bc(jad1+nps*ir)=bc(iad2+ir)                             4d9s23
               end do                                                   4d3s23
               iad2=iad2+nroot                                          4d3s23
              end do                                                    4d3s23
             end do                                                     4d3s23
             iad1=iad1+nbdet(jsb)*nadet(isb)                            4d3s23
             jsb2=nsbeta2(isb)                                          6d22s22
             jdad1=idad1+nbdet(jsb2)*(ilc(isb)-1)                       4d3s23
             do ia=ilc(isb),ihc(isb)                                        5d15s19
              do ir=0,nroot-1                                           4d3s23
               do ibr=0,nbdet(jsb2)-1                                    7d20s22
                bc(jdad1+ibr+ir*nps2)=bc(idad2+ibr+ir*nxn2)             4d13s23
               end do                                                   4d3s23
              end do                                                        5d15s19
              jdad1=jdad1+nbdet(jsb2)                                   7d20s22
              idad2=idad2+nbdet(jsb2)                                   7d20s22
             end do
             idad1=idad1+nbdet(jsb2)*nadet(isb)                         7d20s22
            end do                                                      4d3s23
            call dws_gsumf(bc(ivecf),nroot*(nps+nps2))                  7d21s22
           end if                                                       6d22s22
           jcasvec=iad2                                                 7d20s22
           ivecfn=ibcoff                                                7d14s22
           ibcoff=ivecfn+nps*nroot                                      7d14s22
           call enough('cas0. 25',bc,ibc)
           call savetops(bc(ivecf),nps,ibc(ips),nroot,nsymb,nadet,nbdet,7d14s22
     $          nsbeta,bc(ivecfn))                                      7d14s22
           if(lprt)then                                                 7d14s22
            write(6,*)('vector in det order: ')
            call prntm2(bc(ivecf),nps,nroot,nps)                             5d16s22
            call dgemm('t','n',nroot,nroot,nps,1d0,bc(ivecf),nps,
     $           bc(ivecf),nps,0d0,bc(ibcoff),nroot,
     d'cas0.  1')
            call prntm2(bc(ibcoff),nroot,nroot,nroot)
            write(6,*)('ps order ')
            call prntm2(bc(ivecfn),nps,nroot,nps)                       7d14s22
            write(6,*)('dguess vector in standard order: ')
            call prntm2(bc(idvecd),nps2,nroot,nps2)
           end if                                                       7d14s22
           do ii=0,nps*nroot-1                                          7d14s22
            bc(ivecf+ii)=bc(ivecfn+ii)                                  7d14s22
           end do                                                       7d14s22
           ibcoff=ivecfn                                                7d14s22
           if(mode.eq.4)then                                            7d7s22
            call pshamdbl(nps,ibc(ips),bc(ivecf),                       7d7s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),jdenpt,ibc(itmpc),                                7d7s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,       6d22s22
     $       nsymb,nsbeta,nroot,dbg,nowpro,numpro,1d0,bc,ibc,mdenoff,   3d15s23
     $       iden1,igoal)                                                     3d15s23
           else                                                         7d7s22
            ieigp=ibcoff                                                7d19s22
            ibcoff=ieigp+nroot*2                                        7d20s22
            call enough('cas0. 26',bc,ibc)
            iphamdr=ibcoff                                               5d16s22
            ibcoff=iphamdr+nps*nhere2                                    6d22s22
            call enough('cas0. 27',bc,ibc)
            do iz=ieigp,ibcoff-1                                        3d31s23
             bc(iz)=0d0                                                  5d16s22
            end do                                                       5d16s22
            ivtmp=ibcoff                                                7d19s22
            ibcoff=ivtmp+nhere*nroot                                    7d19s22
            call enough('cas0. 28',bc,ibc)
            do iz=ivtmp,ibcoff-1                                        7d19s22
             bc(iz)=0d0                                                 7d19s22
            end do                                                      7d19s22
            if(nhere.gt.0)then                                          8d31s22
             call dgemm('n','n',nhere,nroot,nps,1d0,bc(iphamr),nhere,    7d19s22
     $           bc(ivecf),nps,0d0,bc(ivtmp),nhere,                     7d19s22
     d'cas0.  2')
            end if                                                      8d31s22
            do ir=0,nroot-1                                             7d19s22
             dot=0d0                                                    7d19s22
             iad1=ivtmp+nhere*ir                                        7d19s22
             iad2=ivecf+istart-1+nps*ir                                 7d19s22
             do j=0,nhere-1                                             7d19s22
              dot=dot+bc(iad1+j)*bc(iad2+j)                             7d19s22
             end do                                                     7d19s22
             bc(ieigp+ir)=dot                                           7d19s22
            end do                                                      7d19s22
            call dws_gsumf(bc(ieigp),nroot)                             7d19s22
            if(lprt)then                                                7d25s22
             do ir=1,nroot                                                7d19s22
              irm=ir-1
              write(6,*)('my energies '),ir,bc(ieigp+irm)
             end do
            end if                                                      7d25s22
            ibcoff=ivtmp
            if(ipuse.eq.1)then                                            6d22s22
             call pshampp(nps,ibc(ips),bc(itmpdhd),bc(iphamdr),numb,      5d16s22
     $      ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),
     $      ibc(itmpb),ih0ed,iooood,ibc(itmpc),                         5d16s22
     $      ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $      nsymb,nsbeta,dbg,nhere,istart,                              12d23s22
     $           .false.,bc,ibc)                                        11d9s22
             iphamr2=iphamr                                              6d23s22
            else                                                         6d22s22
             ips2=ibcoff                                                 6d23s22
             ihdigo=ips2+nps2                                            6d23s22
             ibcoff=ihdigo+nhere2                                        6d23s22
             call enough('cas0. 29',bc,ibc)
             call pshamppnona1(nps,ibc(ips),bc(iphamdr),numb,            6d22s22
     $      ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),
     $      ibc(itmpb),ih0ed,iooood,ibc(itmpc),                         5d16s22
     $      ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $      nsymb,nsbeta,nsbeta2,dbg,nhere2,istart2,ipxder,ipuse,       6d22s22
     $      multh,ibc(ips2),bc(ihdigo),bc(idhdig),bc,ibc)               11d9s22
             if(lprt)then                                               7d14s22
              call prntm2(bc(iphamdr),nhere2,nps,nhere2)                     5d16s22
             end if                                                     7d14s22
             iphamr2=ibcoff                                              6d23s22
             ibcoff=iphamr2+nhere2*nps2                                  6d23s22
             ihdigo=ibcoff                                              7d18s22
             iihdigo=ihdigo+nps2                                        7d18s22
             ibcoff=iihdigo+nps2                                        7d18s22
             call enough('cas0. 30',bc,ibc)
             call diaghii(bc(ihdigo),nps2,ibc(iaorb),nalpha,numa,       7d14s22
     $     ibc(iborb),nbeta,numb,ih0e,iooooa,shift0,ilc,ihc,1d0,        5d17s22
     $      nsymb,nsbeta2,ibc(iihdigo),bc,ibc)                          11d14s22
             do iz=iihdigo,ibcoff-1                                     7d18s22
              bc(iz)=0d0                                                7d18s22
             end do                                                     7d18s22
             jhdigo=ihdigo                                              7d18s22
             jjhdigo=iihdigo                                            7d18s22
             do isa=1,nsymb                                             7d18s22
              isb=nsbeta2(isa)                                          7d18s22
              do ia=ilc(isa),ihc(isa)                                   7d18s22
               iad1=jjhdigo+nbdet(isb)*(ia-1)                           7d18s22
               do ib=0,nbdet(isb)-1                                     7d18s22
                bc(iad1+ib)=bc(jhdigo+ib)                               7d18s22
               end do                                                   7d18s22
               jhdigo=jhdigo+nbdet(isb)                                 7d18s22
              end do                                                    7d18s22
              jjhdigo=jjhdigo+nbdet(isb)*nadet(isa)                     7d18s22
             end do                                                     7d18s22
             call dws_gsumf(bc(iihdigo),nps2)                           7d18s22
             call pshampp(nps2,ibc(ips2),bc(iihdigo),bc(iphamr2),numb,  7d18s22
     $      ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),              6d23s22
     $      ibc(itmpb),ih0e,iooooa,ibc(itmpc),                          6d23s22
     $      ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $      nsymb,nsbeta2,dbg,nhere2,istart2,ikall.eq.-25,bc,ibc)        11d28s22
             if(lprt)then                                               7d14s22
              write(6,*)('my part of ham of other symmetry ')
              call prntm2(bc(iphamr2),nhere2,nps2,nhere2)                 6d23s22
             end if                                                     7d14s22
            end if                                                       6d22s22
            if(lprt)then                                                7d14s22
             write(6,*)('my part of derivative hamiltonian: '),istart2,
     $          nhere2,nsymb,ipuse
             call prntm2(bc(iphamdr),nhere2,nps,nhere2)                     5d16s22
             if(mynprocg.gt.1)then
              iptmp=ibcoff
              iptmp2=iptmp+nhere2*nps
              ibcoff=iptmp2+nhere*nps
              call enough('cas0. 31',bc,ibc)
              do j=0,nps-1
               do k=0,nhere2-1
                kj=iphamdr+k+nhere2*j
                jk=iptmp+j+nps*k
                bc(jk)=bc(kj)
               end do
              end do
              write(6,*)('transposed ')
              call prntm2(bc(iptmp),nps,nhere2,nps)
              call savetonormal(bc(iptmp),nps,ibc(ips),nhere2,nsymb,
     $             nadet,nbdet,nsbeta,bc(iptmp2))
              write(6,*)('rows in normal order ')
              call prntm2(bc(iptmp2),nps,nhere2,nps)
              do j=0,nps-1
               do k=0,nhere2-1
                kj=iptmp+k+nhere2*j
                jk=iptmp2+j+nps*k
                bc(kj)=bc(jk)
               end do
              end do
              write(6,*)('transposed ')
              call prntm2(bc(iptmp),nhere2,nps,nhere2)
              ibcoff=iptmp
             end if
            end if                                                      7d14s22
            npsp1=nps2+1                                                 6d22s22
            irhs=ibcoff                                                  5d16s22
            ibcoff=irhs+nps2*nroot                                      7d20s22
            call enough('cas0. 32',bc,ibc)
            do iz=irhs,ibcoff-1                                          5d16s22
             bc(iz)=0d0                                                  5d16s22
            end do                                                       5d16s22
            jrhs=irhs+istart2-1                                          6d22s22
            if(nhere2.gt.0)then                                          6d22s22
             call dgemm('n','n',nhere2,nroot,nps,1d0,bc(iphamdr),nhere2,7d20s22
     $          bc(ivecf),nps,0d0,bc(jrhs),nps2,                        7d20s22
     d'cas0.  3')
            end if                                                       5d16s22
            call dws_gsumf(bc(irhs),nps2*nroot)                         7d20s22
            if(lprt)then                                                7d14s22
             write(6,*)('times vector ')
             call prntm2(bc(irhs),nps2,nroot,nps2)                      7d20s22
            end if                                                      7d14s22
            itmp=ibcoff                                                  7d18s22
            ibcoff=itmp+nps2*nroot                                       7d18s22
            call enough('cas0. 33',bc,ibc)
            if(ipuse.eq.1)then                                           7d18s22
             call savetops(bc(idvecd),nps,bc(ips),nroot,nsymb,nadet,     7d18s22
     $          nbdet,nsbeta,bc(itmp))                                  7d18s22
            else                                                         7d18s22
             call savetops(bc(idvecd),nps2,bc(ips2),nroot,nsymb,nadet,   7d18s22
     $          nbdet,nsbeta2,bc(itmp))                                 7d18s22
            end if                                                       7d18s22
            do ii=0,nps2*nroot-1                                         7d18s22
             bc(idvecd+ii)=bc(itmp+ii)                                   7d18s22
            end do                                                       7d18s22
            if(lprt)then                                                7d25s22
             write(6,*)('dguess vector in ps order: ')
             call prntm2(bc(idvecd),nps2,nroot,nps2)
            end if                                                      7d25s22
            ibcoff=itmp                                                  7d18s22
            dere=0d0                                                     5d16s22
            if(ipuse.eq.1)then                                           6d23s22
             do ir=0,nroot-1                                            7d20s22
              dere=0d0                                                  7d20s22
              jvecf=ivecf+nps*ir                                        7d20s22
              jrhs=irhs+nps*ir                                          7d20s22
              do ii=0,nps-1                                                 5d16s22
               dere=dere+bc(jvecf+ii)*bc(jrhs+ii)                          5d16s22
              end do                                                       5d16s22
              idsav=ieigp+ir+nroot                                      7d20s22
              bc(idsav)=dere                                            7d20s22
              do ii=0,nps-1                                                5d16s22
               bc(jrhs+ii)=bc(jvecf+ii)*dere-bc(jrhs+ii)                   5d16s22
              end do                                                       5d16s22
             end do                                                     7d20s22
            else                                                         6d23s22
             do ii=0,nroot*nps2-1                                       7d20s22
              bc(irhs+ii)=-bc(irhs+ii)                                   6d23s22
             end do                                                      6d23s22
            end if                                                       6d23s22
            if(lprt)write(6,*)('dere '),dere,dshift,dere+dshift         6d26s24
            if(ldynw.and.(i.eq.1.or.idynw.eq.1))then                    7d20s22
             eref=bc(ieigp)                                             7d20s22
             deref=bc(ieigp+nroot)                                      7d20s22
            end if                                                      7d20s22
            do ir=0,nroot-1                                             7d28s22
             des(ides+ir)=bc(ieigp+nroot+ir)+dshift                     12d23s22
            end do                                                      7d28s22
            ides=ides+nroot                                             7d28s22
            if(ldynw)then                                               7d20s22
             call dynwtr(nroot,bc(ieigp),bc(iwgt),eref,deref,dynw,       7d20s22
     $           eweight(1,i),wsum,dwsum,shift,dshift,1)                7d20s22
            else                                                        7d20s22
             do ir=0,nroot-1                                            7d20s22
              bc(iwgt+nroot+ir)=0d0                                     7d20s22
             end do                                                     7d20s22
            end if                                                      7d20s22
            if(lprt)then                                                7d14s22
             write(6,*)('form c dere-dh c')
             call prntm2(bc(irhs),nps2,nroot,nps2)
             if(ipuse.ne.1)then
              irhso=ibcoff
              ibcoff=irhso+nps2*nroot
              call savetonormal(bc(irhs),nps2,ibc(ips2),nroot,nsymb,
     $             nadet,nbdet,nsbeta2,bc(irhso))
              write(6,*)('in normal order ')
              call prntm2(bc(irhso),nps2,nroot,nps2)
              ibcoff=irhso
             end if
            end if                                                      7d14s22
            esub=eavg-shift
            if(lprt)write(6,*)('esub '),esub,eavg,shift                 7d14s22
            if(nps2.le.1)then                                           11d28s22
             bc(idvecd)=0d0                                             7d11s22
            else                                                        7d11s22
             call cgvec(bc(iphamr2),bc(irhs),istart2,nhere2,nps2,       7d20s22
     $            bc(ieigp),nroot,bc(ivecf),bc(idvecd),tolv,ipuse.eq.1, 7d20s22
     $            lprt,bc,ibc)                                          11d9s22
            end if                                                      7d11s22
            ivder=idvecd                                                 5d19s22
            isovec=ibcoff                                               7d14s22
            isoder=isovec+nps*nroot                                     7d14s22
            ibcoff=isoder+max(nps2,nps)*nroot                           3d16s23
            call enough('cas0. 34',bc,ibc)
            call savetonormal(bc(ivecf),nps,ibc(ips),nroot,nsymb,       7d14s22
     $           nadet,nbdet,nsbeta,bc(isovec))                         7d14s22
            if(ipuse.eq.1)then                                          7d18s22
             call savetonormal(bc(idvecd),nps,ibc(ips),nroot,nsymb,       7d14s22
     $           nadet,nbdet,nsbeta,bc(isoder))                         7d14s22
            else                                                        7d18s22
             call savetonormal(bc(idvecd),nps2,ibc(ips2),nroot,nsymb,       7d14s22
     $           nadet,nbdet,nsbeta2,bc(isoder))                         7d14s22
            end if                                                      7d18s22
            if(lprt)then                                                7d25s22
             write(6,*)('der in normal order '),loop                          7d14s22
             call prntm2(bc(isoder),nps2,nroot,nps2)                       7d14s22
             write(6,*)('nconf = '),nconf                                7d14s22
            end if                                                      7d25s22
            nshort=0                                                    3d16s23
            do isa=1,nsymb                                              3d16s23
             if(ipuse.eq.1)then                                         3d16s23
              isb=nsbeta(isa)                                           3d16s23
             else                                                       3d16s23
              isb=nsbeta2(isa)                                          3d16s23
             end if                                                     3d16s23
             nshort=nshort+nbdet(isb)*(ihc(isa)+1-ilc(isa))             3d16s23
            end do                                                      3d16s23
            ishortv=ibcoff                                              7d14s22
            ishortd=ishortv+nshort*nroot                                3d16s23
            ibcoff=ishortd+nshort*nroot                                 3d16s23
            call enough('cas0. 35',bc,ibc)
            do iz=ishortv,ibcoff-1                                      3d16s23
             bc(iz)=0d0                                                 3d16s23
            end do                                                      3d16s23
            jshortv=ishortv-1                                           7d14s22
            jshortd=ishortd-1                                           7d14s22
            jsoder=isoder-1                                             7d14s22
            jsovec=isovec-1                                             7d14s22
            do ir=1,nroot                                               7d14s22
             nshort=0                                                    7d18s22
             do isa=1,nsymb                                              7d14s22
              if(ipuse.eq.1)then                                        7d18s22
               isb=nsbeta(isa)                                          7d18s22
              else                                                      7d18s22
               isb=nsbeta2(isa)                                          7d18s22
              end if                                                    7d18s22
              jshortd0=jshortd+1
              nshort=nshort+nbdet(isb)*(ihc(isa)+1-ilc(isa))            7d18s22
              do ia=ilc(isa),ihc(isa)                                    7d14s22
               iad1=jsovec+nbdet(isb)*(ia-1)                            7d14s22
               iad2=jsoder+nbdet(isb)*(ia-1)                            7d14s22
               do ib=1,nbdet(isb)                                       7d14s22
                bc(jshortv+ib)=bc(iad1+ib)                              7d14s22
                bc(jshortd+ib)=bc(iad2+ib)                              7d14s22
               end do                                                   7d14s22
               jshortv=jshortv+nbdet(isb)                               7d14s22
               jshortd=jshortd+nbdet(isb)                               7d14s22
              end do                                                    7d14s22
              jsovec=jsovec+nadet(isa)*nbdet(isb)                       7d14s22
              jsoder=jsoder+nadet(isa)*nbdet(isb)                       7d14s22
             end do                                                     7d14s22
            end do                                                      7d14s22
            kcasvec=jcasvec                                             3d31s23
            do ix=i+1,nstate                                            3d21s23
             lcasvec=nveccnt(ix)                                        4d18s23
             ldoit=multh(ipuse,isinfo(1,i)).eq.isinfo(1,ix).and.        3d31s23
     $          isinfo(2,i).eq.isinfo(2,ix).and.                        3d31s23
     $          isinfo(3,i).eq.isinfo(3,ix)                             3d31s23
             nshortx=0                                                  3d31s23
             do isb=1,nsymb                                             3d31s23
              jsb=multh(isb,isinfo(1,ix))                               3d31s23
              nshortx=nshortx+nbdet(jsb)*nherec(isb)                    3d31s23
             end do                                                     3d31s23
             if(ldoit)then                                              3d31s23
              icdir=ibcoff                                                3d20s23
              ibcoff=icdir+isinfo(4,i)*isinfo(4,ix)                       3d20s23
              call enough('cas0.icdira',bc,ibc)                            3d20s23
              do iz=icdir,ibcoff-1                                        3d20s23
               bc(iz)=0d0                                                 3d20s23
              end do                                                      3d20s23
             end if                                                     3d21s23
             lshortd=ishortd                                            3d21s23
             do isb=1,nsymb                                               3d20s23
              jsb=multh(isb,isinfo(1,ix))                                 3d20s23
              if(jsb.eq.nsbeta2(isb).and.ldoit)then                     3d31s23
               iadc=icdir                                               3d21s23
               do irx=0,isinfo(4,ix)-1                                  3d21s23
                do ir=0,isinfo(4,i)-1                                   3d21s23
                 lshortdo=lshortd+ir*nshort                             3d24s23
                 do ia=0,nherec(isb)-1                                  3d21s23
                  do ib=0,nbdet(jsb)-1                                  3d21s23
                   iadd=lshortdo+ib+nbdet(jsb)*ia                       3d24s23
                   iadv=lcasvec+irx+isinfo(4,ix)*(ia+nherec(isb)*ib)    4d3s23
                   term=bc(iadd)*bc(iadv)                               3d21s23
                   bc(iadc+ir)=bc(iadc+ir)+bc(iadd)*bc(iadv)            3d21s23
                  end do                                                3d21s23
                 end do                                                 3d21s23
                end do                                                  3d21s23
                iadc=iadc+nroot                                         3d21s23
               end do                                                   3d21s23
              end if                                                    3d21s23
              lshortd=lshortd+nbdet(nsbeta2(isb))*nherec(isb)           3d31s23
              lcasvec=lcasvec+isinfo(4,ix)*nbdet(jsb)*nherec(isb)       4d3s23
             end do                                                     3d21s23
             kcasvec=lcasvec                                            4d3s23
             if(ldoit)then                                              3d31s23
              do ii=0,isinfo(4,i)*isinfo(4,ix)-1                        3d21s23
               cdc(icdc+ii)=bc(icdir+ii)                                3d21s23
              end do                                                    3d21s23
              ibcoff=icdir                                                3d20s23
              icdc=icdc+nroot*isinfo(4,ix)                                 3d20s23
             else                                                       3d21s23
              do ii=0,isinfo(4,i)*isinfo(4,ix)-1                        3d21s23
               cdc(icdc+ii)=0d0                                         3d21s23
              end do                                                    3d21s23
             end if                                                     3d21s23
            end do                                                      3d21s23
            do ic=0,nshort*nroot-1                                      7d18s22
             bc(jdvec+ic)=bc(ishortd+ic)                                7d18s22
            end do                                                      7d18s22
            jdvec=jdvec+nshort*nroot                                    7d20s22
            if(mode.eq.3)then                                            6d22s22
             if(ipuse.eq.1)then                                          6d23s22
              call diaghiid(bc(ishortv),bc(ishortd),nconf,ibc(iaorb),   7d14s22
     $            nalpha,                                               7d14s22
     $          numa,ibc(iborb),nbeta,numb,iden1(1,2),jdenpt,shift0,ilc,6d22s22
     $          ihc,1d0,nsymb,nsbeta,nroot,bc(iwgt),.false.,bc,ibc,     3d15s23
     $            mdenoff,igoal)                                              3d15s23
              call pshamd(nps,ibc(ips),bc(itmphd),bc(ivecf),bc(ivder),    6d22s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1(1,2),jdenpt,ibc(itmpc),                     6d22s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,       6d22s22
     $       nsymb,nsbeta,nroot,dbg,nowpro,numpro,bc(iwgt),.false.,1,   11d9s22
     $             bc,ibc,mdenoff,1)                                    3d27s23
              call diaghiid(bc(ishortd),bc(ishortv),nconf,ibc(iaorb),   7d14s22
     $             nalpha,                                              7d14s22
     $           numa,ibc(iborb),nbeta,numb,iden1(1,3),jdenpt,shift0,   6d22s22
     $           ilc,ihc,1d0,nsymb,nsbeta,nroot,bc(iwgt),.false.,bc,ibc,3d15s23
     $             mdenoff,igoal)                                             3d15s23
              call pshamd(nps,ibc(ips),bc(itmphd),bc(ivder),bc(ivecf),    6d22s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1(1,3),jdenpt,ibc(itmpc),                     6d22s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $       nsymb,nsbeta,nroot,dbg,nowpro,numpro,bc(iwgt),.false.,1,   11d9s22
     $             bc,ibc,mdenoff,1)                                    3d27s23
              call diaghiid(bc(ishortd),bc(ishortd),nconf,ibc(iaorb),   7d14s22
     $             nalpha,                                              7d14s22
     $           numa,ibc(iborb),nbeta,numb,iden1(1,4),jdenpt,shift0,   6d22s22
     $           ilc,ihc,1d0,nsymb,nsbeta,nroot,bc(iwgt),.false.,bc,ibc,3d15s23
     $             mdenoff,igoal)                                             3d15s23
              call pshamd(nps,ibc(ips),bc(itmphd),bc(ivder),bc(ivder),    6d22s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1(1,4),jdenpt,ibc(itmpc),                     6d22s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $       nsymb,nsbeta,nroot,dbg,nowpro,numpro,bc(iwgt),.false.,2,   11d9s22
     $             bc,ibc,mdenoff,1)                                    3d27s23
             else                                                        6d23s22
              call pshamdnona1(nps,ibc(ips),bc(ivecf),nsbeta,
     $            nps2,ibc(ips2),bc(ivder),nsbeta2,                     6d23s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1(1,2),jdenpt,ibc(itmpc),                     6d22s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),norb,                     6d23s22
     $       nsymb,nroot,nowpro,numpro,bc(iwgt),.false.,ipxderd,bc,ibc, 3d16s23
     $             mdenoff)                                             3d16s23
              call pshamdnona1(nps2,ibc(ips2),bc(ivder),nsbeta2,         6d23s22
     $            nps,ibc(ips),bc(ivecf),nsbeta,                        6d23s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1(1,3),jdenpt,ibc(itmpc),                     6d22s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),norb,                     6d23s22
     $       nsymb,nroot,nowpro,numpro,bc(iwgt),.false.,ipxderd,bc,ibc, 3d16s23
     $             mdenoff)                                             3d16s23
              call diaghiid(bc(ishortd),bc(ishortd),nshort,ibc(iaorb),  3d16s23
     $             nalpha,                                              7d14s22
     $           numa,ibc(iborb),nbeta,numb,iden1(1,4),jdenpt,shift0,   6d22s22
     $           ilc,ihc,1d0,nsymb,nsbeta2,nroot,bc(iwgt),.false.,bc,   11d14s22
     $             ibc,mdenoff,igoal)                                         3d15s23
              call pshamd(nps2,ibc(ips2),bc(itmphd),bc(ivder),bc(ivder),7d18s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1(1,4),jdenpt,ibc(itmpc),                     6d22s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $       nsymb,nsbeta2,nroot,dbg,nowpro,numpro,bc(iwgt),.false.,2,  11d9s22
     $             bc,ibc,mdenoff,1)                                    3d27s23
             end if                                                      6d23s22
            end if                                                       6d22s22
            if(ipuse.eq.1)then                                           6d23s22
             wuse=0d0
             call diaghiid(bc(ishortv),bc(ishortd),nconf,ibc(iaorb),    7d14s22
     $           nalpha,                                                7d14s22
     $          numa,ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc, 6d23s22
     $              1d0,nsymb,nsbeta,nroot,bc(iwgt),.true.,bc,ibc,      3d15s23
     $            mdenoff,igoal)                                              3d15s23
             call pshamd(nps,ibc(ips),bc(itmphd),bc(ivecf),bc(ivder),    6d23s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d23s22
     $     ibc(itmpb),iden1,jdenpt,ibc(itmpc),                          8d15s14
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $       nsymb,nsbeta,nroot,dbg,nowpro,numpro,bc(iwgt),.true.,1,    11d9s22
     $            bc,ibc,mdenoff,0)                                     3d27s23
             call diaghiid(bc(ishortd),bc(ishortv),nconf,ibc(iaorb),    7d14s22
     $            nalpha,                                               7d14s22
     $           numa,ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc,6d23s22
     $              1d0,nsymb,nsbeta,nroot,bc(iwgt),.true.,bc,ibc,      3d15s23
     $            mdenoff,igoal)                                              3d15s23
             call pshamd(nps,ibc(ips),bc(itmphd),bc(ivder),bc(ivecf),    6d23s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d23s22
     $     ibc(itmpb),iden1,jdenpt,ibc(itmpc),                          5d17s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,       5d17s22
     $       nsymb,nsbeta,nroot,dbg,nowpro,numpro,bc(iwgt),.true.,1,    11d9s22
     $            bc,ibc,mdenoff,0)                                     3d27s23
             if(ldynw)then                                              7d25s22
              call diaghiid(bc(ishortv),bc(ishortv),nconf,ibc(iaorb),    7d14s22
     $           nalpha,                                                7d14s22
     $         numa,ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc,  7d20s22
     $              1d0,nsymb,nsbeta,nroot,bc(iwgt+nroot),.true.,bc,ibc,3d15s23
     $             mdenoff,igoal)                                             3d15s23
              call pshamd(nps,ibc(ips),bc(itmphd),bc(ivecf),bc(ivecf),   7d20s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d23s22
     $     ibc(itmpb),iden1,jdenpt,ibc(itmpc),                          7d20s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $     nsymb,nsbeta,nroot,dbg,nowpro,numpro,bc(iwgt+nroot),.true.,1,11d9s22
     $             bc,ibc,mdenoff,0)                                    3d27s23
             end if                                                     7d25s22
            else                                                         6d23s22
             call pshamdnona1(nps,ibc(ips),bc(ivecf),nsbeta,
     $            nps2,ibc(ips2),bc(ivder),nsbeta2,                     6d23s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1,jdenpt,ibc(itmpc),                          6d23s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),norb,                     6d23s22
     $       nsymb,nroot,nowpro,numpro,bc(iwgt),.true.,ipxderd,bc,ibc)  11d9s22
             call pshamdnona1(nps2,ibc(ips2),bc(ivder),nsbeta2,          6d23s22
     $            nps,ibc(ips),bc(ivecf),nsbeta,                        6d23s22
     $           numb,ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),    6d22s22
     $     ibc(itmpb),iden1,jdenpt,ibc(itmpc),                          6d23s22
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),norb,                     6d23s22
     $       nsymb,nroot,nowpro,numpro,bc(iwgt),.true.,ipxderd,bc,ibc)  11d9s22
            end if                                                       6d23s22
           end if                                                       7d7s22
           ibcoff=ivecf                                                 5d16s22
           go to 4                                                      5d16s22
          end if                                                        5d16s22
          if((icall.eq.1..or.ioconv.ne.0).and.lprint.or.iprtr(8).ne.0)
     $         then                                                     2d13s20
           write(6,*)('p-space is everything ')
           write(6,*)('for i = '),i
           write(spinm,1066)isinfo(2,i)                                   1d5s20
 1066      format(i2)                                                     1d5s20
           write(6,*)('spin: '),isinfo(2,i),(' space symmetry: '),
     $       isinfo(1,i),('('),spinm,(char(isinfo(j,i)),j=6,nlab),(')')    12d31s19
           iroosv=0                                                     5d3s21
           do j=6,min(nroolab,nlab)                                     5d3s21
            if(isinfo(j,i).ne.ibc(jroolab+j))iroosv=1                   5d3s21
           end do                                                       5d3s21
           if(nlab.ne.nroolab.or.iroosv.ne.0.or.spinm.ne.spinroo)then   5d3s21
            spinroo=spinm                                               5d3s21
            iroosv=1                                                    5d3s21
            nroolab=nlab                                                5d3s21
            do j=6,nlab                                                 5d3s21
             ibc(jroolab+j)=isinfo(j,i)                                 5d3s21
            end do                                                      5d3s21
           end if                                                       5d3s21
           do i1=1,nroot
            write(6,1551)i1,bc(ieigp+i1-1)+shift                        5d17s22
            if(iroosv.ne.0)then                                         5d3s21
             bc(jroodat)=bc(ieigp+i1-1)+shift                           4d17s23
             jroodat=jroodat+1                                          5d3s21
            end if                                                      5d3s21
           end do
          end if
          do k=0,nroot-1
           sum=0d0
           jvecp=ivecp+nok*k
           do l=0,nok-1
            sum=sum+bc(jvecp+l)**2
           end do
          end do
          if(nrydb.gt.0)then                                            5d7s21
           call ilimts(nps,1,mynprocg,mynowprog,ilps,ihps,i1s,i1e,i2s,   5d6s21
     $         i2e)                                                     5d6s21
           nconfps=ihps+1-ilps                                           5d6s21
          else                                                          5d7s21
           nconfps=nconf                                                5d7s21
          end if                                                        5d7s21
          ivecps=ibcoff
          ibcoff=ivecps+nconfps*nroot                                   5d6s21
          ivecsc=ibcoff
          ibcoff=ivecsc+nconfps*nok                                     5d6s21
          call enough('cas0. 36',bc,ibc)
          if(nrydb.gt.0)then                                            5d7s21
           do j=0,nok-1                                                  5d6s21
            iad1=ivecs+ilps+nps*j-1                                      5d6s21
            iad2=ivecsc+nconfps*j                                        5d6s21
            do l=0,nconfps-1                                              5d6s21
             bc(iad2+l)=bc(iad1+l)                                       5d6s21
            end do                                                       5d6s21
           end do                                                        5d6s21
          else                                                          5d7s21
           irow=0                                                          10d24s14
           jvecs=ivecs                                                     12d18s15
           do isb=1,nsymb                                                  10d24s14
            do l=0,nherec(isb)-1                                           10d24s14
             lp=l+ilc(isb)-1                                               10d24s14
             do k=0,nbdet(nsbeta(isb))-1                                   10d24s14
              do j=0,nok-1                                                   10d24s14
               iad1=jvecs+k+nbdet(nsbeta(isb))*lp+nps*j                    10d24s14
               iad2=ivecsc+irow+nconf*j                                    10d24s14
               bc(iad2)=bc(iad1)                                           10d24s14
              end do                                                       10d24s14
              irow=irow+1                                                  10d24s14
             end do                                                        10d24s14
            end do                                                         10d24s14
            jvecs=jvecs+nbdet(nsbeta(isb))*nadet(isb)                      12d18s15
           end do                                                          10d24s14
          end if                                                        5d7s21
c
c     Htilde=vsT H vs
c     vhT Htilde vh=diag=vhT vsT H vs vh
c
          do k=0,nroot-1
           sumcheck=0d0
           do j=0,nok-1                                                    3d31s17
            sumcheck=sumcheck+bc(ivecp+j+nok*k)**2                         3d31s17
           end do
          end do
          if(nconfps.gt.0.and.nok.gt.0)then                             5d7s21
           call dgemm('n','n',nconfps,nroot,nok,1d0,bc(ivecsc),nconfps, 5d6s21
     $         bc(ivecp),nok,0d0,bc(ivecps),nconfps,                    5d6s21
     d' cas0.  2')
          end if                                                          2d22s19
          ibcoff=ivecsc                                                    10d24s14
          do k=0,nroot-1
           sumcheck=0d0
           do j=0,nconf-1
            sumcheck=sumcheck+bc(ivecps+j+nconf*k)**2
           end do
           call dws_gsumf(sumcheck,1)                                       12d23s15
          end do
c
c     reconstruct vectors on all procs so we can compute densities
c
          ivecf=ibcoff                                                     12d23s15
          ibcoff=ivecf+nroot*nps                                           12d23s15
          call enough('cas0. 37',bc,ibc)
          do j=0,nroot*nps-1                                               12d23s15
           bc(ivecf+j)=0d0                                                 12d23s15
          end do                                                           12d23s15
          if(nrydb.gt.0)then                                            5d7s21
           jvecf=ivecf                                                   5d6s21
           do isb=1,nsymb                                                5d6s21
            iveca(isb)=jvecf                                             5d7s21
            jvecf=jvecf+nadet(isb)*nbdet(nsbeta(isb))                    5d6s21
     $           *nroot                                                 1d31s22
           end do                                                        5d6s21
           do l=0,nconfps-1                                              5d6s21
            lp=l+ilps-1                                                  5d6s21
            itrial=ibc(ips+lp)                                           5d6s21
            j2=iveca(itrial2(3))+nroot*(itrial2(1)-1                    1d31s22
     $           +nadet(itrial2(3))*(itrial2(2)-1))                     1d31s22
            do ir=0,nroot-1                                              5d6s21
             bc(j2+ir)=bc(ivecps+l+nconfps*ir)                          1d31s22
            end do                                                       5d6s21
           end do                                                        5d6s21
          else                                                          5d7s21
           jrown=0                                                          12d23s15
           jrowo=0                                                          12d23s15
           do isb=1,nsymb                                                   12d23s15
c
c     there is a problem with how I cut down storage in ivecps to procs ...
c     although this code reassembles the full vector, the indexing of
c     ivecps does not reflect what is actually stored there.
c
            iveca(isb)=ivecf+jrown                                         1d18s20
            do l=0,nherec(isb)-1                                            12d23s15
             lp=l+ilc(isb)-1                                                12d23s15
             do k=0,nbdet(nsbeta(isb))-1                                    12d23s15
              irown=jrown+k+nbdet(nsbeta(isb))*lp                           12d23s15
              irowo=jrowo+k+nbdet(nsbeta(isb))*l                            12d23s15
              do m=0,nroot-1                                                12d23s15
               iad1=ivecps+irowo+nconf*m                                    12d23s15
               iad2=ivecf+irown+nps*m
               bc(iad2)=bc(iad1)                                            12d23s15
              end do                                                        12d23s15
             end do                                                         12d23s15
            end do                                                          12d23s15
            jrown=jrown+nbdet(nsbeta(isb))*nadet(isb)                       12d23s15
            jrowo=jrowo+nbdet(nsbeta(isb))*nherec(isb)                      12d23s15
           end do                                                           12d23s15
          end if                                                        5d7s21
          call dws_gsumf(bc(ivecf),nps*nroot)                              12d23s15
          ivecfn=ibcoff                                                 7d14s22
          ibcoff=ivecfn+nps*nroot                                       7d14s22
          call enough('cas0. 38',bc,ibc)
          call savetonormal(bc(ivecf),nps,ibc(ips),nroot,nsymb,nadet,   7d14s22
     $         nbdet,nsbeta,bc(ivecfn))                                 7d14s22
          iad1=ivecfn                                                   7d14s22
          iad2=jcasvec                                                    5d15s19
          nveccnt(i)=jcasvec                                            4d18s23
          nxn=0                                                         4d3s23
          nyy=0
          do isb=1,nsymb                                                4d3s23
           nxn=nxn+nbdet(nsbeta(isb))*nherec(isb)                       4d3s23
          end do                                                        4d3s23
          iad1=ivecfn                                                   4d3s23
          do isb=1,nsymb                                                4d3s23
           jsb=nsbeta(isb)                                              4d3s23
           do ibr=0,nbdet(jsb)-1                                        4d3s23
            do ia=ilc(isb),ihc(isb)                                      4d3s23
             jad1=iad1+ibr+nbdet(jsb)*(ia-1)                            4d3s23
             do ir=0,nroot-1                                            4d3s23
              bc(iad2+ir)=bc(jad1+nps*ir)                               4d9s23
             end do                                                     4d3s23
             iad2=iad2+nroot                                            4d3s23
            end do                                                      4d3s23
           end do                                                       4d3s23
           iad1=iad1+nbdet(jsb)*nadet(isb)                              4d3s23
          end do                                                        4d3s23
          jcasvec=iad2                                                  7d20s22
          do k=0,nroot-1                                                  5d17s18
           if((icall.eq.1.or.ioconv.ne.0).and.lprint.and.mode.eq.1)     5d16s22
     $        write(6,*)('for root no. '),k+1
           jvecf=ivecfn+nps*k                                           7d14s22
           sum=0d0                                                        5d17s18
           do isb=1,nsymb                                                  5d17s18
            jsb=nsbeta(isb)                                                5d17s18
            do ia=0,nadet(isb)-1                                          5d17s18
             iap=ia+1                                                     5d17s18
             do ib=0,nbdet(jsb)-1                                         5d17s18
              ibp=ib+1                                                    5d17s18
              if(nrydb.gt.0)then                                        1d31s22
               kvecf=iveca(isb)+k+nroot*(ia+nadet(isb)*ib)              1d31s22
              else                                                      1d31s22
               kvecf=jvecf+ib+nbdet(jsb)*ia                                5d17s18
              end if                                                    1d31s22
              sum=sum+bc(kvecf)**2
              if(abs(bc(kvecf)).gt.0.1d0.and.(icall.eq.1.or.ioconv.ne.0)
     $           .and.lprint)then
               call prtocc(bc(kvecf),isb,jsb,iap,ibp,nadet,nbdet,          5d17s18
     $           ibc(iaorb),ibc(iborb),nalpha,nbeta,norb,iacto,nsymb,   11d9s22
     $              bc,ibc)                                             11d9s22
              end if
             end do                                                       5d17s18
            end do                                                        5d17s18
            jvecf=jvecf+nadet(isb)*nbdet(jsb)                             5d17s18
           end do                                                          5d17s18
          end do                                                          5d17s18
          ibcoff=ivecfn                                                 7d14s22
          if(nrydb.gt.0.and.mode.eq.1)then                              5d16s22
           call genryd(jmats,kmats,bc(ih0mo),                           5d7s18
     $         ibc(iaorb),nalpha,ibc(iborb),nbeta,nsbeta,nconf,iveca,   1d31s22
     $         nvirtc,m12sym,idata1,idatb1,bc(ieigp),nroot,iorbn,morb,  1d30s22
     $         nbasisp,ihryd,i,nstate,lprint,nryd,ncomp,nlzz,iorbsym,   5d3s21
     $         iorbsymz,bc(iwgt),wsum,shift,bc,ibc,icanon)              5d5s23
          end if                                                        5d7s18
          nrootdbg=1
          nvecoff=ioffps                                                8d12s22
          call diaghiid(bc(ivecf+nvecoff),bc(ivecf+nvecoff),nps,        5d16s22
     $         ibc(iaorb),nalpha,numa,                                  5d16s22
     $              ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc,   8d19s14
     $              1d0,nsymb,nsbeta,nroot,ewght(1,i),.true.,bc,ibc,    3d15s23
     $         mdenoff,igoal)                                                 3d15s23
          do j=1,11
           dbg(j)=1d0
          end do
          call pshamd(nps,ibc(ips),bc(itmphd),bc(ivecf),bc(ivecf),numb, 5d16s22
     $     ibc(iaorb),nalpha,ibc(iborb),nbeta,ibc(itmpa),
     $     ibc(itmpb),iden1,jdenpt,ibc(itmpc),                          8d15s14
     $       ibc(itmpd),ibc(itmpe),ibc(itmpf),1d0,1d0,norb,ncsym,        8d15s06
     $       nsymb,nsbeta,nroot,dbg,nowpro,numpro,ewght(1,i),.true.,2,  11d9s22
     $         bc,ibc,mdenoff,0)                                        3d27s23
          go to 4
         end if
         if(lprt)write(6,*)('pspace is fractional ')
         ivecps=ibcoff
         ibcoff=ivecps+nps*nroot
         call enough('cas0. 39',bc,ibc)
         if(nps.gt.0)then                                                 2d22s19
          call dgemm('n','n',nps,nroot,nok,1d0,bc(ivecs),nps,bc(ivecp),
     $     nok,0d0,bc(ivecps),nps,
     d' cas0.  3')
         end if                                                           2d22s19
         iarg2=nps
         do isb=1,nsymb
          jsb=nsbeta(isb)
          nn=nadet(isb)*nbdet(jsb)*nroot*2-1
          do k=0,nn
           bc(iveca(isb)+k)=0d0
          end do
         end do
         if(.not.bcasvec)then                                             4d20s18
          kcasvec=jcasvec                                                 4d20s18
          do isb=1,nsymb                                                  4d20s18
           jsb=nsbeta(isb)                                                4d20s18
           iad1=kcasvec
           iad2=ivec(isb)                                                 4d20s18
           do j=0,nherec(isb)*nbdet(jsb)-1                                4d20s18
            do k=0,nroot-1                                                4d20s18
             bc(iad2+k)=bc(kcasvec+k)                                     4d20s18
            end do                                                        4d20s18
            iad2=iad2+nroot                                               4d20s18
            kcasvec=kcasvec+nroot                                         4d20s18
           end do                                                         4d20s18
           do j=0,nbdet(jsb)-1                                            4d20s18
            do k=0,nherec(isb)-1                                          4d20s18
             kp=k+ilc(isb)-1                                              4d20s18
             iad2=ivec(isb)+nroot*(k+nherec(isb)*j)                       4d20s18
             iad3=iveca(isb)+nroot*(kp+nadet(isb)*j)                      4d20s18
             do l=0,nroot-1                                               4d20s18
              bc(iad3+l)=bc(iad2+l)                                       4d20s18
             end do                                                       4d20s18
            end do                                                        4d20s18
           end do                                                         4d20s18
          end do                                                          4d20s18
         else                                                             4d20s18
          do i3=1,nroot
           do i1=1,npsh                                                    11d8s05
            itrial=ibc(ipsh+i1-1)                                          8d29s06
            isb=itrial2(3)                                                 3d14s17
            j1=i1+ioffps                                                   3d14s17
            j2=ivec(isb)+i3-1+nroot*(itrial2(1)-ilc(isb)                   3d14s17
     $        +nherec(isb)*(itrial2(2)-1))                              3d14s17
            j3=ivecps+j1-1+(i3-1)*nps                                      3d12s17
            bc(j2)=bc(j3)
            j4=iveca(isb)+i3-1+nroot*(itrial2(1)-1                         3d14s17
     $       +nadet(isb)*(itrial2(2)-1))                                3d14s17
            bc(j4)=bc(j3)                                                  3d10s17
           end do
          end do
         end if                                                           4d20s18
         call dws_gsumf(bc(iveca(1)),nroot*ncona*2)                       3d27s17
         do i3=0,nroot-1
          sum=0d0
          nsum=0
          do isb=1,nsymb
           jsb=nsbeta(isb)
           iadd=iveca(isb)+i3                                             3d14s17
           do icol=0,nbdet(jsb)-1
            do irow=0,nadet(isb)-1
             nsum=nsum+1
             iad=iveca(isb)+i3+nroot*(irow+nadet(isb)*icol)               3d14s17
             sum=sum+bc(iad)**2
            end do
           end do
          end do
         end do
         do isb=1,nsymb                                                   3d14s17
          jsb=nsbeta(isb)                                                 3d14s17
          do icol=0,nadet(jsb)-1                                          3d14s17
           do irow=0,nherect(isb)-1                                       3d14s17
            jrow=irow+ilct(isb)-1                                         3d14s17
            iad1=ivect(isb)+nroot*(irow+nherect(isb)*icol)                3d14s17
            iad2=iveca(jsb)+nroot*(icol+nadet(jsb)*jrow)                  3d14s17
            do i3=0,nroot-1                                               3d14s17
             bc(iad1+i3)=bc(iad2+i3)                                      3d14s17
            end do                                                        3d14s17
           end do                                                         3d14s17
          end do                                                          3d14s17
         end do
c
c     if iteration is multiple of kick, use random guess for next       4d9s18
c     trial vector and don't do convergence test.                       4d9s18
c                                                                       4d9s18
         kick=1000000                                                          4d9s18
         iter=0
         do isb=1,nsymb                                                   8d29s06
          ig(isb)=ibcoff                                                  8d29s06
          ibcoff=ig(isb)+2*nherec(isb)*nroot*nbdet(nsbeta(isb))           3d27s17
          igo(isb)=ibcoff                                                  8d29s06
          ibcoff=igo(isb)+2*nherec(isb)*nroot*nbdet(nsbeta(isb))           3d27s17
          igt(isb)=ibcoff                                                 3d10s17
          ibcoff=igt(isb)+2*nherect(isb)*nroot*nadet(nsbeta(isb))         3d27s17
          ivecx(isb)=ibcoff                                               8d29s06
          ibcoff=ivecx(isb)+2*nherec(isb)*nroot*nbdet(nsbeta(isb))        3d27s17
         end do                                                           8d29s06
         itmp1=ibcoff
         itmp2=itmp1+nps*nroot*2                                          5d8s20
         ibcoff=itmp2+nok
         ipxham=ibcoff
         nokr=nok+nroot*2                                                 3d27s17
         ibcoff=ipxham+nokr**2
         ieigo=ibcoff
         ibcoff=ieigo+nroot
         iconv=ibcoff
         ibcoff=iconv+nroot
         call enough('cas0. 40',bc,ibc)
         if(mode.ge.2)then                                              6d2s22
          jcasvec=kcasvec                                               7d28s22
          if(mode.eq.4)then
           call diaghiidbl(ibc(iaorb),nalpha,numa,ibc(iborb),nbeta,numb,7d21s22
     $          jdenpt,nsymb,nsbeta,ivec,ilc,ihc,bc,ibc,nroot,mdenoff,  3d15s23
     $          iden1,igoal)                                                  3d15s23
           call hc1cdbl(ivec,ivect,iveca,ilc,ihc,ilct,ihct,
     $       nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,               3d27s17
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta,m1c,idatac,idatbc,jdenpt,bc,ibc,nroot,mdenoff,      3d15s23
     $          iden1,igoal)                                                  3d15s23
           go to 4                                                      7d21s22
          end if
          ieuse=ibcoff                                                  5d27s22
          idere=ieuse+nroot                                             5d27s22
          ibcoff=idere+nroot                                            5d27s22
          call enough('cas0. 41',bc,ibc)
          do iz=ieuse,ibcoff-1                                          5d26s22
           bc(iz)=0d0                                                   5d26s22
          end do                                                        5d26s22
          do ipass=1,2                                                  6d9s22
           if(ipass.eq.1)then                                           5d27s22
            call hc1c(ivec,ivect,iveca,ig,igt,ilc,ihc,ilct,ihct,nroot,    5d26s22
     $       bc(ihdig),ih0e,nherec,nherect,ibc(iaorb),ibc(iborb),       5d26s22
     $       nalpha,nbeta,iooooa,m12sym,norb,nsymb,idata1,idata2,idatb1, 5d26s22
     $       idatb2,nsbeta,m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc)   11d10s22
            juse=ieuse                                                  5d27s22
           else                                                         5d27s22
            if(ipuse.eq.1)then                                          6d10s22
             call hc1c(ivec,ivect,iveca,ig,igt,ilc,ihc,ilct,ihct,nroot,    5d26s22
     $       bc(idhdig),ih0ed,nherec,nherect,ibc(iaorb),ibc(iborb),     5d26s22
     $       nalpha,nbeta,iooood,m12sym,norb,nsymb,idata1,idata2,idatb1,5d26s22
     $       idatb2,nsbeta,m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc)   11d10s22
            else                                                        6d10s22
             nconnona1=0                                                 6d10s22
             do isb=1,nsymb                                             6d10s22
              jsb=multh(nsbeta(isb),ipuse)                              6d10s22
              ig(isb)=ibcoff                                            6d10s22
              ibcoff=ig(isb)+2*nherec(isb)*nroot*nbdet(jsb)             6d10s22
              nconnona1=nconnona1+nadet(isb)*nbdet(jsb)                 6d10s22
              igt(isb)=ibcoff                                           6d10s22
              ibcoff=igt(isb)+2*nherect(isb)*nroot*nadet(jsb)           6d10s22
             end do                                                     6d10s22
             call enough('cas0. 42',bc,ibc)
             do iz=ig(1),ibcoff-1                                       6d13s22
              bc(iz)=0d0                                                6d13s22
             end do                                                     6d13s22
             call hc1cnona1(ivec,ivect,iveca,ig,igt,nroot,ih0ed,        6d10s22
     $            nherec,ilc,ihc,nherect,ilct,ihct,ibc(iaorb),          6d13s22
     $            ibc(iborb),nalpha,nbeta,                              6d13s22
     $            iooood,nsblkder,isblkder,norb,nsymb,ipuse,multh,      6d13s22
     $            m12symnon,idata1non,idata2non,idatb1non,idatb2non,    6d13s22
     $            ipxder,nsbeta,m12sym,idata1,idatb1,m1c,idatac,idatbc, 11d10s22
     $            bc,ibc,igoal)                                               11d10s22
            end if                                                      6d10s22
            juse=idere                                                  5d27s22
           end if                                                       5d27s22
           igfull=ibcoff                                                    3d12s17
           if(ipass.eq.1.or.ipuse.eq.1)then                             6d10s22
            ibcoff=igfull+nroot*ncona                                     3d27s17
           else                                                         6d10s22
            ibcoff=igfull+nroot*nconnona1                               6d10s22
           end if                                                       6d10s22
           call enough('cas0. 43',bc,ibc)
           nwds=ibcoff-igfull                                           6d10s22
           do iz=igfull,ibcoff-1                                        6d10s22
            bc(iz)=0d0                                                  6d10s22
           end do                                                       6d10s22
           jgfull=igfull                                                    3d12s17
           do isb=1,nsymb                                                   3d12s17
            jsb=nsbeta(isb)                                                 3d12s17
            if(ipass.eq.2)jsb=multh(jsb,ipuse)                          6d10s22
            iga(isb)=jgfull                                                 3d12s17
            do icol=0,nbdet(jsb)-1                                          3d14s17
             do irow=0,nherec(isb)-1                                        3d14s17
              iad1=ig(isb)+nroot*(irow+nherec(isb)*icol)                 3d27s17
              iad2=jgfull+nroot*(irow+ilc(isb)-1+nadet(isb)*icol)        3d27s17
              do i3=0,nroot-1                                            3d27s17
               bc(iad2+i3)=bc(iad1+i3)                                      3d14s17
              end do                                                        3d14s17
             end do                                                         3d14s17
            end do                                                          3d14s17
            jgfull=jgfull+nadet(isb)*nroot*nbdet(jsb)                    3d27s17
           end do                                                           3d12s17
           do jsb=1,nsymb                                                   3d12s17
            isb=nsbeta(jsb)                                                 3d12s17
            if(ipass.eq.2)isb=multh(isb,ipuse)                          6d10s22
            do icol=0,nadet(isb)-1                                          3d14s17
             do irow=0,nherect(jsb)-1                                       3d14s17
              iad1=igt(jsb)+nroot*(irow+nherect(jsb)*icol)               3d27s17
              iad2=iga(isb)+nroot*(icol+nadet(isb)*(irow+ilct(jsb)-1))   3d27s17
              do i3=0,nroot-1                                            3d27s17
               orig=bc(igoal)
               bc(iad2+i3)=bc(iad2+i3)+bc(iad1+i3)                          3d14s17
              end do                                                        3d14s17
             end do                                                         3d14s17
            end do                                                          3d14s17
           end do
           call dws_gsumf(bc(igfull),nwds)                                  3d12s17
           if(lprt)then                                                 7d28s22
            if(ipass.eq.1)then                                           5d27s22
             write(6,*)('what we have for gfull '),igfull
            else
             write(6,*)('what we have for dgfull '),igfull
            end if
            call prntm2(bc(igfull),nwds,1,nwds)
           end if                                                       7d28s22
           if(ipass.eq.1.or.ipuse.eq.1)then                             6d10s22
            do isb=1,nsymb                                                   3d14s17
             jsb=nsbeta(isb)                                                 3d14s17
             do icol=0,nbdet(jsb)-1                                          3d14s17
              do irow=0,nadet(isb)-1                                      5d26s22
               iad1=iga(isb)+nroot*(irow+nadet(isb)*icol)                 5d26s22
               iad2=iveca(isb)+nroot*(irow+nadet(isb)*icol)                5d26s22
               do i3=0,nroot-1                                               3d14s17
                bc(juse+i3)=bc(juse+i3)+bc(iad2+i3)*bc(iad1+i3)          5d27s22
               end do                                                        3d14s17
              end do                                                         3d14s17
             end do                                                          3d14s17
            end do
           end if                                                       6d10s22
           if(ipass.eq.1)then                                           5d27s22
            if(lprt)write(6,*)('energies ')                             7d28s22
            ibcoff=igfull                                               5d27s22
           else
            if(lprt)write(6,*)('energy derivatives '),juse                   7d28s22
           end if
           if(lprt)call prntm2(bc(juse),1,nroot,1)                      7d28s22
          end do                                                        5d27s22
          nxn=0
          if(ipuse.eq.1)then                                            6d15s22
           do isb=1,nsymb                                                5d27s22
            jsb=nsbeta(isb)                                              5d27s22
            idva(isb)=ibcoff                                             5d27s22
            ibcoff=idva(isb)+nroot*nbdet(jsb)*nadet(isb)                 5d27s22
            nxn=nxn+nbdet(jsb)*nadet(isb)
            call enough('cas0. 44',bc,ibc)
            do icol=0,nbdet(jsb)-1                                       5d27s22
             do irow=0,nadet(isb)-1                                      5d27s22
              iad1=iga(isb)+nroot*(irow+nadet(isb)*icol)                 5d27s22
              iad2=iveca(isb)+nroot*(irow+nadet(isb)*icol)               5d27s22
              iad3=idva(isb)+nroot*(irow+nadet(isb)*icol)                5d27s22
              do i3=0,nroot-1                                            5d27s22
               bc(iad1+i3)=bc(idere+i3)*bc(iad2+i3)-bc(iad1+i3)          5d27s22
               bc(iad3+i3)=0d0                                           5d27s22
              end do                                                     5d27s22
             end do                                                      5d27s22
            end do                                                       5d27s22
           end do                                                        5d27s22
          else
           do isb=1,nsymb                                                5d27s22
            jsb=multh(ipuse,nsbeta(isb))                                6d15s22
            nsbeta2(isb)=jsb                                            6d15s22
            idva(isb)=ibcoff                                             5d27s22
            ibcoff=idva(isb)+nroot*nbdet(jsb)*nadet(isb)                 5d27s22
            nxn=nxn+nbdet(jsb)*nadet(isb)
            call enough('cas0. 45',bc,ibc)
            do icol=0,nbdet(jsb)-1                                       5d27s22
             do irow=0,nadet(isb)-1                                      5d27s22
              iad1=iga(isb)+nroot*(irow+nadet(isb)*icol)                 5d27s22
              iad3=idva(isb)+nroot*(irow+nadet(isb)*icol)                5d27s22
              do i3=0,nroot-1                                            5d27s22
               bc(iad1+i3)=-bc(iad1+i3)                                 6d15s22
               bc(iad3+i3)=0d0                                           5d27s22
              end do                                                     5d27s22
             end do                                                      5d27s22
            end do                                                       5d27s22
           end do                                                        5d27s22
          end if
          kdvec=jdvec                                                   5d31s22
          do isb=1,nsymb                                                5d31s22
           jsb=multh(ipuse,nsbeta(isb))                                 6d15s22
           do icol=0,nbdet(jsb)-1                                       5d31s22
            do irow=0,nherec(isb)-1                                     5d31s22
             irowp=irow+ilc(isb)-1                                      5d31s22
             iad3=idva(isb)+nroot*(irowp+nadet(isb)*icol)               5d31s22
             iad2=kdvec+nroot*(irow+nherec(isb)*icol)                   5d31s22
             do i3=0,nroot-1                                            5d31s22
              bc(iad3+i3)=bc(iad2+i3)                                   5d31s22
             end do                                                     5d31s22
            end do                                                      5d31s22
           end do                                                       5d31s22
           kdvec=kdvec+nroot*nbdet(jsb)*nherec(isb)                     5d31s22
          end do                                                        5d31s22
          call dws_gsumf(bc(idva(1)),nxn*nroot)                         5d31s22
          if(ipuse.eq.1)then
           call cgvec2(bc(ihdig),ih0e,iooooa,iga,ieuse,iveca,idva,ig,   6d15s22
     $          igt,nroot,nadet,nbdet,nsbeta,nherec,ilc,ihc,nherect,    6d15s22
     $          ilct,ihct,ncona,iaorb,iborb,nalpha,nbeta,m12sym,norb,   6d15s22
     $          nsymb,idata1,idata2,idatb1,idatb2,m1c,idatac,idatbc,mst,6d15s22
     $          mnz,nz,tolv,1,lprt,bc,ibc,tolb)                         4d27s23
          else                                                          6d15s22
           call cgvec2(bc(idhdig),ih0e,iooooa,iga,ieuse,iveca,idva,ig,  6d15s22
     $          igt,nroot,nadet,nbdet,nsbeta2,nherec,ilc,ihc,nherect,   6d15s22
     $         ilct,ihct,nxn,iaorb,iborb,nalpha,nbeta,m12sym,norb,      4d4s23
     $          nsymb,idata1,idata2,idatb1,idatb2,m1c,idatac,idatbc,mst,6d15s22
     $          mnz,nz,tolv,0,lprt,bc,ibc,tolb)                         4d27s23
          end if                                                        6d15s22
          tolbest=max(tolbest,tolb)                                     4d27s23
          lcasvec=jcasvec                                               3d21s23
          do ix=i+1,nstate                                              3d21s23
           lcasvec=nveccnt(ix)                                          4d18s23
           ldoit=multh(ipuse,isinfo(1,i)).eq.isinfo(1,ix).and.          3d31s23
     $          isinfo(2,i).eq.isinfo(2,ix).and.                        3d31s23
     $          isinfo(3,i).eq.isinfo(3,ix)                             3d31s23
           if(ldoit)then                                                3d31s23
            icdir=ibcoff                                                3d20s23
            ibcoff=icdir+isinfo(4,i)*isinfo(4,ix)                       3d20s23
            call enough('cas0.icdir',bc,ibc)                            3d20s23
            do iz=icdir,ibcoff-1                                        3d20s23
             bc(iz)=0d0                                                 3d20s23
            end do                                                      3d20s23
           end if                                                       3d20s23
           do isb=1,nsymb                                               3d20s23
            jsb=multh(isb,isinfo(1,ix))                                 3d20s23
            if(jsb.eq.multh(ipuse,nsbeta(isb)).and.ldoit)then           3d31s23
             do icol=0,nbdet(jsb)-1                                      3d20s23
              do irow=0,nherec(isb)-1                                    3d20s23
               irowp=irow+ilc(isb)-1                                    3d20s23
               iaddc=idva(isb)+nroot*(irowp+nadet(isb)*icol)            3d20s23
               iadc=icdir                                               3d20s23
               do irx=0,isinfo(4,ix)-1                                          3d20s23
                do ir=0,isinfo(4,i)-1                                   3d20s23
                 term=bc(lcasvec+irx)*bc(iaddc+ir)
                 bc(iadc+ir)=bc(iadc+ir)+bc(lcasvec+irx)*bc(iaddc+ir)   3d20s23
                end do                                                  3d20s23
                iadc=iadc+isinfo(4,i)                                   3d20s23
               end do                                                    3d20s23
               lcasvec=lcasvec+isinfo(4,ix)                              3d20s23
              end do                                                     3d20s23
             end do                                                      3d20s23
            else                                                        3d20s23
             lcasvec=lcasvec+isinfo(4,ix)*nherec(isb)*nbdet(jsb)        3d20s23
            end if                                                      3d20s23
           end do                                                       3d20s23
           if(ldoit)then                                                3d31s23
            do ii=0,nroot*isinfo(4,ix)-1                                3d20s23
             cdc(icdc+ii)=bc(icdir+ii)                                  3d20s23
            end do                                                      3d20s23
            ibcoff=icdir                                                3d20s23
            icdc=icdc+nroot*isinfo(4,ix)                                 3d20s23
           else                                                         3d20s23
            do ii=0,nroot*isinfo(4,ix)-1                                3d20s23
             cdc(icdc+ii)=0d0                                           3d20s23
            end do                                                      3d20s23
           end if                                                       3d20s23
          end do                                                        3d20s23
          nconfder=0                                                    6d15s22
          do isb=1,nsymb                                                5d31s22
           jsb=multh(ipuse,nsbeta(isb))                                 6d15s22
           nconfder=nconfder+nbdet(jsb)*nherec(isb)                     6d15s22
           do icol=0,nbdet(jsb)-1                                       5d31s22
            do irow=0,nherec(isb)-1                                     5d31s22
             irowp=irow+ilc(isb)-1                                      5d31s22
             iad3=idva(isb)+nroot*(irowp+nadet(isb)*icol)               5d31s22
             iad2=jdvec+nroot*(irow+nherec(isb)*icol)                   5d31s22
             do i3=0,nroot-1                                            5d31s22
              bc(iad2+i3)=bc(iad3+i3)                                   5d31s22
             end do                                                     5d31s22
            end do                                                      5d31s22
           end do                                                       5d31s22
           jdvec=jdvec+nroot*nherec(isb)*nbdet(jsb)                     5d31s22
          end do                                                        5d31s22
          itmpv=ibcoff                                                  5d27s22
          itmpdv=itmpv+nconf*nroot                                      5d27s22
          ibcoff=itmpdv+nconfder*nroot                                  6d15s22
          call enough('cas0. 46',bc,ibc)
          ioff=0                                                        3d27s17
          if(ipuse.eq.1)then                                            6d15s22
           do isb=1,nsymb                                                3d27s17
            jsb=nsbeta(isb)                                              3d27s17
            nn=nherec(isb)*nbdet(jsb)                                    3d27s17
            do ia=0,nherec(isb)-1
             iap=ia+ilc(isb)-1                                           5d27s22
             do ib=0,nbdet(jsb)-1
              iad1=ivec(isb)+nroot*(ia+nherec(isb)*ib)                   3d31s17
              iad2=itmpv+ioff+ib+nbdet(jsb)*ia                            3d31s17
              iad3=idva(isb)+nroot*(iap+nadet(isb)*ib)                   5d27s22
              iad4=itmpdv+ioff+ib+nbdet(jsb)*ia                          5d27s22
              iad5=ig(isb)+nroot*(ia+nherec(isb)*ib)                     5d27s22
              do k=0,nroot-1                                             5d27s22
               bc(iad2+k*nconf)=bc(iad1+k)                               5d27s22
               bc(iad4+k*nconf)=bc(iad3+k)                               5d27s22
               bc(iad5+k)=bc(iad3+k)                                     5d27s22
              end do                                                      3d27s17
             end do                                                       3d27s17
            end do
            ioff=ioff+nn                                                 3d27s17
            do ia=0,nadet(isb)-1                                         5d27s22
             do ib=0,nherect(jsb)-1                                      5d27s22
              ibp=ib+ilct(jsb)-1                                         5d27s22
              iad1=igt(jsb)+nroot*(ib+nherect(jsb)*ia)                   5d27s22
              iad2=idva(isb)+nroot*(ia+nadet(isb)*ibp)                   5d27s22
              do k=0,nroot-1                                             5d27s22
               bc(iad1+k)=bc(iad2+k)                                     5d27s22
              end do                                                     5d27s22
             end do                                                      5d27s22
            end do                                                       5d27s22
           end do                                                        3d27s17
          else                                                          6d15s22
           ioff=0                                                       6d15s22
           ioffder=0                                                    6d15s22
           do isb=1,nsymb                                                3d27s17
            jsb=nsbeta(isb)                                              3d27s17
            nn=nherec(isb)*nbdet(jsb)                                    3d27s17
            do ia=0,nherec(isb)-1
             iap=ia+ilc(isb)-1                                           5d27s22
             do ib=0,nbdet(jsb)-1
              iad1=ivec(isb)+nroot*(ia+nherec(isb)*ib)                   3d31s17
              iad2=itmpv+ioff+ib+nbdet(jsb)*ia                            3d31s17
              do k=0,nroot-1                                             5d27s22
               bc(iad2+k*nconf)=bc(iad1+k)                               5d27s22
              end do                                                      3d27s17
             end do                                                       3d27s17
            end do
            ioff=ioff+nn                                                 3d27s17
            jsb2=nsbeta2(isb)                                           6d15s22
            nn=nherec(isb)*nbdet(jsb2)                                    3d27s17
            do ia=0,nherec(isb)-1
             iap=ia+ilc(isb)-1                                           5d27s22
             do ib=0,nbdet(jsb2)-1
              iad3=idva(isb)+nroot*(iap+nadet(isb)*ib)                   5d27s22
              iad4=itmpdv+ioffder+ib+nbdet(jsb2)*ia                          5d27s22
              iad5=ig(isb)+nroot*(ia+nherec(isb)*ib)                     5d27s22
              do k=0,nroot-1                                             5d27s22
               bc(iad4+k*nconfder)=bc(iad3+k)                               5d27s22
               bc(iad5+k)=bc(iad3+k)                                     5d27s22
              end do                                                      3d27s17
             end do                                                       3d27s17
            end do
            ioffder=ioffder+nn                                          6d15s22
            do ia=0,nadet(isb)-1                                        6d15s22
             do ib=0,nherect(jsb2)-1                                    6d15s22
              ibp=ib+ilct(jsb2)-1                                       6d15s22
              iad1=igt(jsb2)+nroot*(ib+nherect(jsb2)*ia)                6d15s22
              iad2=idva(isb)+nroot*(ia+nadet(isb)*ibp)                   5d27s22
              do k=0,nroot-1                                             5d27s22
               bc(iad1+k)=bc(iad2+k)                                     5d27s22
              end do                                                     5d27s22
             end do                                                      5d27s22
            end do                                                       5d27s22
           end do                                                       6d15s22
          end if                                                        6d15s22
          if(lprt)then                                                  7d28s22
           write(6,*)('what we have fo tmpv: ')
           call prntm2(bc(itmpv),nconf,1,nconf)
          end if                                                        7d28s22
          if(ldynw.and.(i.eq.1.or.idynw.eq.1))then                      7d27s22
           eref=bc(ieuse)                                               7d27s22
           deref=bc(ieuse+nroot)                                        7d27s22
          end if                                                        7d27s22
          do ir=0,nroot-1                                               7d28s22
           des(ides+ir)=bc(ieuse+nroot+ir)+dshift                       12d23s22
          end do                                                        7d28s22
          ides=ides+nroot                                               7d28s22
          if(ldynw)then                                                 7d27s22
           call dynwtr(nroot,bc(ieuse),bc(iwgt),eref,deref,dynw,        7d27s22
     $           eweight(1,i),wsum,dwsum,shift,dshift,1)                7d20s22
          else                                                          7d27s22
           do ir=0,nroot-1                                              7d27s22
            bc(iwgt+nroot+ir)=0d0                                       7d27s22
           end do                                                       7d27s22
          end if                                                        7d27s22
          if(mode.eq.3)then                                             6d2s22
           if(ipuse.eq.1)then                                           6d15s22
            call diaghiid(bc(itmpv),bc(itmpdv),nconf,ibc(iaorb),nalpha,   5d27s22
     $         numa,ibc(iborb),nbeta,numb,iden1(1,2),jdenpt,shift0,ilc, 6d2s22
     $         ihc,1d0,nsymb,nsbeta,nroot,1d0,.false.,bc,ibc,mdenoff,   3d16s23
     $           igoal)                                                 3d16s23
            call diaghiid(bc(itmpdv),bc(itmpv),nconf,ibc(iaorb),nalpha,   5d27s22
     $         numa,ibc(iborb),nbeta,numb,iden1(1,3),jdenpt,shift0,     6d2s22
     $          ilc,ihc,1d0,nsymb,nsbeta,nroot,1d0,.false.,bc,ibc,      3d15s23
     $           mdenoff,igoal)                                         3d16s23
            call diaghiid(bc(itmpdv),bc(itmpdv),nconf,ibc(iaorb),nalpha,   5d27s22
     $         numa,ibc(iborb),nbeta,numb,iden1(1,4),jdenpt,shift0,     6d2s22
     $          ilc,ihc,1d0,nsymb,nsbeta,nroot,1d0,.false.,bc,ibc,      3d15s23
     $           mdenoff,igoal)                                         3d16s23
            call hc1cd(ivec,ivect,iveca,ig,igt,idva,ilc,ihc,ilct,ihct,
     $         nroot,                                                   5d27s22
     $       nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,               3d27s17
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta,m1c,idatac,idatbc,iden1(1,2),jdenpt,1d0,.false.,bc, 11d10s22
     $           ibc,mdenoff,igoal)                                           3d15s23
            call hc1cd(ig,igt,idva,ivec,ivect,iveca,ilc,ihc,ilct,ihct,
     $         nroot,                                                   5d27s22
     $       nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,               3d27s17
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta,m1c,idatac,idatbc,iden1(1,3),jdenpt,1d0,.false.,bc, 11d10s22
     $           ibc,mdenoff,igoal)                                           3d15s23
            call hc1cd(ig,igt,idva,ig,igt,idva,ilc,ihc,ilct,ihct,        6d2s22
     $         nroot,                                                   5d27s22
     $       nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,               3d27s17
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta,m1c,idatac,idatbc,iden1(1,4),jdenpt,1d0,.false.,bc, 11d10s22
     $           ibc,mdenoff,igoal)                                           3d15s23
           else                                                         6d15s22
            call hc1cnona1d(ivec,ivect,iveca,ig,igt,idva,nroot,         6d15s22
     $            nherec,ilc,ihc,nherect,ilct,ihct,ibc(iaorb),          6d13s22
     $            ibc(iborb),nalpha,nbeta,                              6d13s22
     $            nsblkder,isblkder,norb,nsymb,ipuse,multh,             6d15s22
     $            m12symnon,idata1non,idata2non,idatb1non,idatb2non,    6d13s22
     $            ipxderd,nsbeta,m12sym,idata1,idatb1,iden1(1,2),jdenpt,7d1s22
     $            1d0,.false.,m1c,idatac,idatbc,bc,ibc,mdenoff)         3d15s23
             call hc1cnona1d(ig,igt,idva,ivec,ivect,iveca,nroot,        6d15s22
     $            nherec,ilc,ihc,nherect,ilct,ihct,ibc(iaorb),          6d13s22
     $            ibc(iborb),nalpha,nbeta,                              6d13s22
     $            nsblkder,isblkder,norb,nsymb,ipuse,multh,             6d15s22
     $            m12symnon,idata1non,idata2non,idatb1non,idatb2non,    6d13s22
     $           ipxderd,nsbeta2,m12sym,idata1,idatb1,iden1(1,3),jdenpt,7d1s22
     $            1d0,.false.,m1c,idatac,idatbc,bc,ibc,mdenoff)         3d15s23
            call diaghiid(bc(itmpdv),bc(itmpdv),nconfder,ibc(iaorb),    6d15s22
     $           nalpha,numa,ibc(iborb),nbeta,numb,iden1(1,4),jdenpt,   6d15s22
     $           shift0,ilc,ihc,1d0,nsymb,nsbeta2,nroot,1d0,.false.,bc, 11d14s22
     $            ibc,mdenoff,igoal)                                    3d16s23
            call hc1cd(ig,igt,idva,ig,igt,idva,ilc,ihc,ilct,ihct,        6d2s22
     $         nroot,                                                   5d27s22
     $       nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,               3d27s17
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta2,m1c,idatac,idatbc,iden1(1,4),jdenpt,1d0,.false.,   11d10s22
     $           bc,ibc,mdenoff,igoal)                                        3d15s23
           end if
          end if                                                        6d2s22
          if(ipuse.eq.1)then                                            6d15s22
           call diaghiid(bc(itmpv),bc(itmpdv),nconf,ibc(iaorb),nalpha,   5d27s22
     $         numa,ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc,  5d27s22
     $              1d0,nsymb,nsbeta,nroot,bc(iwgt),.true.,bc,ibc,      3d15s23
     $         mdenoff,igoal)                                                 3d15s23
           call diaghiid(bc(itmpdv),bc(itmpv),nconf,ibc(iaorb),nalpha,   5d27s22
     $         numa,ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc,  5d27s22
     $              1d0,nsymb,nsbeta,nroot,bc(iwgt),.true.,bc,ibc,      3d15s23
     $          mdenoff,igoal)                                                3d15s23
           call hc1cd(ivec,ivect,iveca,ig,igt,idva,ilc,ihc,ilct,ihct,
     $         nroot,                                                   5d27s22
     $       nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,               3d27s17
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta,m1c,idatac,idatbc,iden1,jdenpt,bc(iwgt),.true.,bc,  11d10s22
     $          ibc,mdenoff,igoal)                                            3d15s23
           call hc1cd(ig,igt,idva,ivec,ivect,iveca,ilc,ihc,ilct,ihct,
     $         nroot,                                                   5d27s22
     $       nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,               3d27s17
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta,m1c,idatac,idatbc,iden1,jdenpt,bc(iwgt),.true.,bc,  11d10s22
     $          ibc,mdenoff,igoal)                                            3d15s23
           if(ldynw)then                                                7d27s22
            call diaghiid(bc(itmpv),bc(itmpv),nconf,ibc(iaorb),nalpha,   7d27s22
     $         numa,ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc,  7d27s22
     $              1d0,nsymb,nsbeta,nroot,bc(iwgt+nroot),.true.,bc,ibc,3d15s23
     $          mdenoff,igoal)                                                3d15s23
            call hc1cd(ivec,ivect,iveca,ivec,ivect,iveca,ilc,ihc,ilct,  7d27s22
     $           ihct,nroot,nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,7d27s22
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,nsbeta,7d27s22
     $           m1c,idatac,idatbc,iden1,jdenpt,bc(iwgt+nroot),.true.,  11d10s22
     $           bc,ibc,mdenoff,igoal)                                        3d15s23
           end if                                                       7d27s22
          else                                                          6d16s22
            call hc1cnona1d(ivec,ivect,iveca,ig,igt,idva,nroot,         6d15s22
     $            nherec,ilc,ihc,nherect,ilct,ihct,ibc(iaorb),          6d13s22
     $            ibc(iborb),nalpha,nbeta,                              6d13s22
     $            nsblkder,isblkder,norb,nsymb,ipuse,multh,             6d15s22
     $            m12symnon,idata1non,idata2non,idatb1non,idatb2non,    6d13s22
     $            ipxderd,nsbeta,m12sym,idata1,idatb1,iden1,jdenpt,     7d1s22
     $            bc(iwgt),.true.,m1c,idatac,idatbc,bc,ibc,mdenoff)     3d15s23
             call hc1cnona1d(ig,igt,idva,ivec,ivect,iveca,nroot,        6d15s22
     $            nherec,ilc,ihc,nherect,ilct,ihct,ibc(iaorb),          6d13s22
     $            ibc(iborb),nalpha,nbeta,                              6d13s22
     $            nsblkder,isblkder,norb,nsymb,ipuse,multh,             6d15s22
     $            m12symnon,idata1non,idata2non,idatb1non,idatb2non,    6d13s22
     $            ipxderd,nsbeta2,m12sym,idata1,idatb1,iden1,jdenpt,    7d1s22
     $            bc(iwgt),.true.,m1c,idatac,idatbc,bc,ibc,mdenoff)     3d15s23
          end if
          go to 4                                                       5d27s22
         end if                                                         5d26s22
         do i1=1,nroot
          bc(ieigo+i1-1)=bc(ieigp+i1-1)
          ibc(iconv+i1-1)=0                                               4d6s18
         end do
         iarg1=numb*nroot
         iarg2=nherec(1)
         nrootu=nroot                                                     3d27s17
         timehc1c=0d0
         timemid=0d0
         timehouse=0d0
         timeend=0d0
         timeorth=0d0                                                     4d9s18
         call second(timeb)                                               4d23s18
         telap=timeb-timea-tovr                                           4d23s18
         telapo(2)=telapo(2)+telap                                        4d23s18
    2    continue
         telapo(5)=telapo(5)+1d0
         iter=iter+1
         if(iprtr(8).ne.0)write(6,*)('starting iteration no. '),iter,
     $       ibcoff,icall
         if(iter.gt.maxit)then
          write(6,*)('for i = '),i
          write(6,*)('number of iters for davidson diagonalization ')
          write(6,*)('exceeded maximum in cas0'),maxit
          write(6,*)('best results so far ')
          do i1=1,nroot
           write(6,*)bc(ieigo+i1-1),ibc(iconv+i1-1)
          end do
          write(6,*)('difx = '),difx
          write(6,*)('this is probably due to near degeneracy between'),
     $         (' the highest root computed and the lowest root not'),
     $         (' computed')                                            8d12s22
          write(6,*)('or perhaps, if this is the first iteration,')     6d7s23
          write(6,*)('your initial guess orbitals are not good enough.')6d7s23
          call dws_synca
          call dws_finalize
          stop 'cas0'                                                   8d12s22
         end if
c
c     form g=h*c
c
         call second(time1)                                              4d6s18
         call hc1c(ivec,ivect,iveca,ig,igt,ilc,ihc,ilct,ihct,nrootu,     3d27s17
     $       bc(ihdig),ih0e,nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,3d10s17
     $       nbeta,iooooa,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,3d14s17
     $       nsbeta,m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc)          11d10s22
         call second(time2)                                              4d6s18
         telap=time2-time1-tovr
         timehc1c=timehc1c+telap
         time1=time2
         nokr=nok+nrootu                                                 3d27s17
         if(iter.gt.1)then                                               4d23s18
          do i1=1,nokr*nokr
           bc(ipxham+i1-1)=0d0
          end do
c
c     get external part
c
          do isb=1,nsymb                                                   8d29s06
           do i1=0,nrootu*nherec(isb)*nbdet(nsbeta(isb))-1               3d27s17
            bc(ivecx(isb)+i1)=bc(ivec(isb)+i1)                             3d16s17
           end do                                                          8d29s06
          end do                                                           3d17s17
          igfull=ibcoff                                                    3d12s17
          ibcoff=igfull+nrootu*ncona                                     3d27s17
          call enough('cas0. 47',bc,ibc)
          do i3=0,nrootu*ncona-1                                         3d27s17
           bc(igfull+i3)=0d0                                               3d12s17
          end do                                                           3d12s17
          jgfull=igfull                                                    3d12s17
          do isb=1,nsymb                                                   3d12s17
           jsb=nsbeta(isb)                                                 3d12s17
           iga(isb)=jgfull                                                 3d12s17
           do icol=0,nbdet(jsb)-1                                          3d14s17
            do irow=0,nherec(isb)-1                                        3d14s17
             iad1=ig(isb)+nrootu*(irow+nherec(isb)*icol)                 3d27s17
             iad2=jgfull+nrootu*(irow+ilc(isb)-1+nadet(isb)*icol)        3d27s17
             do i3=0,nrootu-1                                            3d27s17
              bc(iad2+i3)=bc(iad1+i3)                                      3d14s17
             end do                                                        3d14s17
            end do                                                         3d14s17
           end do                                                          3d14s17
           jgfull=jgfull+nadet(isb)*nrootu*nbdet(jsb)                    3d27s17
          end do                                                           3d12s17
          do jsb=1,nsymb                                                   3d12s17
           isb=nsbeta(jsb)                                                 3d12s17
           do icol=0,nadet(isb)-1                                          3d14s17
            do irow=0,nherect(jsb)-1                                       3d14s17
             iad1=igt(jsb)+nrootu*(irow+nherect(jsb)*icol)               3d27s17
             iad2=iga(isb)+nrootu*(icol+nadet(isb)*(irow+ilct(jsb)-1))   3d27s17
             do i3=0,nrootu-1                                            3d27s17
              bc(iad2+i3)=bc(iad2+i3)+bc(iad1+i3)                          3d14s17
             end do                                                        3d14s17
            end do                                                         3d14s17
           end do                                                          3d14s17
          end do
          nwds=nrootu*ncona                                              3d27s17
          call dws_gsumf(bc(igfull),nwds)                                  3d12s17
          do i1=0,nps*nrootu-1                                           3d27s17
           bc(itmp1+i1)=0d0                                              3d20s17
          end do                                                          4d27s06
          do i1=0,npsh-1                                                  3d17s17
           j1=i1+ioffps                                                  3d20s17
           itrial=ibc(ipsh+i1)                                            3d17s17
           isb=itrial2(3)                                                 3d16s17
           k1=iga(isb)+nrootu*(itrial2(1)-1                              3d27s17
     $         +nadet(isb)*(itrial2(2)-1))                               3d17s17
           do k=0,nrootu-1                                               3d27s17
            bc(itmp1+j1+nps*k)=bc(k1+k)                                  3d20s17
           end do                                                         3d16s17
          end do
          iarg1=nps*nrootu                                               3d27s17
          call dws_gsumf(bc(itmp1),iarg1)                                 4d27s06
          ibcoff=igfull                                                   3d17s17
          itmp2=ibcoff
          ibcoff=itmp2+nok*nrootu                                        3d27s17
          call enough('cas0. 48',bc,ibc)
          if(nok.gt.0.and.nps.gt.0)then                                  2d22s19
           call dgemm('n','n',nok,nrootu,nps,1d0,bc(ivecst),nok,
     $         bc(itmp1),nps,0d0,bc(itmp2),nok,                               3d17s17
     d' cas0.  4')
          end if                                                         2d22s19
          jpxham=ipxham+nok*nokr                                          3d16s17
          do k=0,nrootu-1                                                3d27s17
           do l=0,nok-1                                                   3d17s17
            bc(jpxham+l+nokr*k)=bc(itmp2+l+nok*k)                         3d17s17
           end do                                                         3d17s17
          end do                                                          3d17s17
          itri=ibcoff                                                     3d16s17
          nwds=(nrootu*(nrootu+1))/2                                     3d27s17
          ibcoff=itri+nwds                                                3d16s17
          call enough('cas0. 49',bc,ibc)
          do i2=0,nwds-1                                                  3d16s17
           bc(itri+i2)=0d0                                                3d16s17
          end do                                                          3d16s17
          do isb=1,nsymb                                                  3d16s17
           nn=nbdet(nsbeta(isb))*nherec(isb)-1                            3d16s17
           do i1=0,nn                                                     3d16s17
            iadg=ig(isb)+nrootu*i1                                       3d27s17
            iadv=ivec(isb)+nrootu*i1                                     3d27s17
            jtri=itri                                                     3d16s17
            do i3=0,nrootu-1                                             3d27s17
             do i2=0,i3                                                   3d16s17
              bc(jtri)=bc(jtri)+bc(iadg+i3)*bc(iadv+i2)                   3d16s17
              jtri=jtri+1                                                 3d16s17
             end do                                                       3d16s17
            end do                                                        3d16s17
           end do                                                         3d16s17
           nn=nadet(nsbeta(isb))*nherect(isb)-1                           3d16s17
           do i1=0,nn                                                     3d16s17
            iadg=igt(isb)+nrootu*i1                                      3d27s17
            iadv=ivect(isb)+nrootu*i1                                    3d27s17
            jtri=itri                                                     3d16s17
            do i3=0,nrootu-1                                             3d27s17
             do i2=0,i3                                                   3d16s17
              bc(jtri)=bc(jtri)+bc(iadg+i3)*bc(iadv+i2)                   3d16s17
              jtri=jtri+1                                                 3d16s17
             end do                                                       3d16s17
            end do                                                        3d16s17
           end do                                                         3d16s17
          end do                                                          3d16s17
          call dws_gsumf(bc(itri),nwds)                                   3d16s17
          if(idwsdeb.ne.0)then                                           8d27s19
           write(6,*)('itri: ')
           call prntm2(bc(itri),1,nwds,1)
          end if                                                         8d27s19
          trisav=bc(itri)
          jtri=itri                                                       3d16s17
          do i3=0,nrootu-1                                               3d27s17
           do i2=0,i3                                                     3d16s17
            iad1=ipxham+nok+i2+nokr*(nok+i3)                              3d16s17
            bc(iad1)=bc(jtri)                                             3d16s17
            iad2=ipxham+nok+i3+nokr*(nok+i2)                              3d16s17
            bc(iad2)=bc(jtri)                                             3d16s17
            jtri=jtri+1                                                   3d16s17
           end do                                                         3d16s17
          end do                                                          3d16s17
          do i1=1,nok
           do i2=1,nok
            j2=ipxham+i2-1+(i1-1)*nokr
            k2=ipham+i2-1+(i1-1)*nok
            bc(j2)=bc(k2)
           end do
           do i2=1,nrootu                                                3d27s17
            j2=ipxham+i1-1+(i2+nok-1)*nokr
            k2=ipxham+i2+nok-1+(i1-1)*nokr
            bc(k2)=bc(j2)
           end do
          end do
          kl=0                                                             3d17s17
          jtmptric=itmptric                                                3d17s17
          jtmptri=itmptri                                                  3d17s17
          do k=1,nokr                                                      3d17s17
           do l=1,k                                                        3d17s17
            if(mod(kl,mynprocg).eq.mynowprog)then                          3d17s17
             if(k.gt.nok)then                                              3d17s17
              if(l.le.nok)then                                             3d17s17
               iad=itmp2+l-1+nok*(k-nok-1)                                 3d17s17
               bc(jtmptri)=bc(iad)                                          3d17s17
              else                                                         3d17s17
               kx=k-nok                                                    3d17s17
               lx=l-nok                                                    3d17s17
               iad=itri+((kx*(kx-1))/2)+lx-1                               3d17s17
               bc(jtmptri)=bc(iad)                                         3d17s17
              end if                                                       3d17s17
             else                                                          3d17s17
              bc(jtmptri)=bc(jtmptric)                                      3d17s17
             end if                                                        3d17s17
             jtmptri=jtmptri+1                                             3d17s17
             jtmptric=jtmptric+1                                           3d17s17
            end if                                                         3d17s17
            kl=kl+1                                                        3d17s17
           end do                                                          3d17s17
          end do                                                           3d17s17
          call second(time2)                                             4d6s18
          telap=time2-time1-tovr                                         4d6s18
          timemid=timemid+telap
          time1=time2
          ibcoff=itmp2                                                     3d17s17
          nroottmp=nroot
          if(idwsdeb.ne.0)write(6,*)('phouse the 2nd in cas0')
          if(nokr.le.npdiag(1))then                                      8d12s22
           call diagdy(bc(itmptri),nokr,nroottmp,ieigh,ivech,bc,ibc)    11d14s22
           isteig=1                                                     8d12s22
          else                                                          8d12s22
           call phouse(itmptri,nokr,1d0,0d0,nroottmp,ieigh,ivech,1,
     $         isteig,bc,ibc)                                           11d9s22
          end if                                                        8d12s22
          call second(time2)                                             4d6s18
          telap=time2-time1-tovr
          timehouse=timehouse+telap
          time1=time2
          if(iprtr(8).ne.0)then                                          2d13s20
           write(6,*)('eigenvalues from phouse ')
           call prntm2(bc(ieigh),1,nroottmp,1)
          end if
          ivecp=ibcoff                                                     3d17s17
          isymz=ivecp+nokr*nokr                                            3d17s17
          ibcoff=isymz+nokr                                                3d17s17
          jvecp=ivecp                                                      3d17s17
          jvech=ivech                                                      3d17s17
          do k=0,nroot-1                                                   3d17s17
           bc(ieigp+k)=bc(ieigh+k)                                         3d17s17
           do l=0,nokr-1                                                   3d17s17
            bc(jvecp+l)=bc(jvech+l)                                        3d17s17
           end do                                                          3d17s17
           jvecp=jvecp+nokr                                                3d17s17
           jvech=jvech+nokr                                                3d17s17
          end do
          itmpg=ibcoff
          ibcoff=itmpg+nokr*nroot
          call enough('cas0. 50',bc,ibc)
          if(nokr.gt.0)then                                              2d22s19
           call dgemm('n','n',nokr,nroot,nokr,1d0,bc(ipxham),nokr,
     $        bc(ivecp),nokr,0d0,bc(itmpg),nokr,
     d' cas0.  5')
          end if                                                         2d22s19
          difx=0d0
          nconv=0
          do i1=1,nroot
           diff=bc(ieigp+i1-1)+shift0-bc(ieigo+i1-1)                     5d17s22
           difx=max(difx,abs(diff))                                      4d30s18
           if((abs(diff).lt.econv.and.ibc(iconv+i1-1).eq.1).or.          4d9s18
     $         bc(iwgt+i1-1).lt.1d-8.or.nokr.eq.nconft)nconv=nconv+1    5d8s18
           vx=0d0
           do i2=nok+1,nokr
            vx=vx+bc(ivecp+i2-1+(i1-1)*nokr)**2
           end do
           call dws_bcast(vx,1)                                          5d1s18
 1551      format(i3,f19.12,1p2e10.2)
           if(vx.gt.0.3d0)then                                           10d16s20
            if(lprint)then                                               5d11s18
             write(6,*)('warning: q-space part of vector for root '),i1
             write(6,*)('is '),vx
             write(6,*)('increasing p-space threshold should improve'),
     $          (' convergence')
             write(6,*)('current p-space threshold: '),pspacex(i)         5d1s18
            end if                                                       5d11s18
            if(npsloop.le.1)then                                         2d9s20
             pspacex(i)=pspacex(i)+0.05d0                                5d1s18
             if(lprint)then                                              5d11s18
              write(6,*)('increase p-space threshold to '),pspacex(i)     5d1s18
              write(6,*)('and try again ')                                5d1s18
             end if                                                      5d11s18
             ibcoff=ibcps                                                 10d15s20
             go to 2222                                                  5d1s18
            end if                                                       5d1s18
            if(vx.gt.0.9d0)then                                            5d30s06
             if(lprint)then                                              5d11s18
              write(6,*)(' weight in external space is too large ')         5d30s06
              write(6,*)('to be meaningful')                                5d30s06
             end if                                                      5d11s18
             call dws_synca
             call dws_finalize                                             5d30s06
             stop                                                          5d30s06
            end if                                                         5d30s06
           end if
          end do
          if(mod(iter,kick).eq.1)nconv=0                                 4d9s18
          if(nconv.eq.nroot)then                                         4d9s18
           if(ldynw.and.(i.eq.1.or.idynw.eq.1))then                       3d17s21
            eref=bc(ieigp)+shift                                        5d17s22
           end if                                                        4d6s18
           if(ldynw)then                                                 3d17s21
            do i1=1,nroot                                                4d6s18
             if(dynw(2).eq.0)then                                       9d20s21
              ex=exp(-dynw(1)*(bc(ieigp+i1-1)+shift-eref))              5d17s22
              exi=1d0/ex                                                  4d6s18
              wwww=(2d0/(ex+exi))**2                                    9d20s21
             else                                                       9d20s21
              wde=dynw(1)*(bc(ieigp+i1-1)+shift-dynw(2))                5d17s22
              wwww=4/(dynw(3)+exp(2d0*wde))                             9d20s21
             end if                                                     9d20s21
             bc(iwgt+i1-1)=eweight(i1,i)*max(1d-7,wwww)                 9d20s21
             wsum=wsum+bc(iwgt+i1-1)                                     4d6s18
             eavg2=eavg2+(bc(ieigp+i1-1)+shift)*ewght(i1,i)             5d17s22
c     ? 4d30s21
            ewght(i1,i)=bc(iwgt+i1-1)                                   5d11s21
            end do                                                       4d6s18
           end if                                                        4d6s18
           do i1=0,nroot-1                                              1d18s20
            eavg=eavg+(bc(ieigp+i1)+shift)*bc(iwgt+i1)                  5d17s22
           end do                                                       1d18s20
           if((icall.eq.1.or.ioconv.ne.0).and.lprint)then                         5d2s18
            write(spinm,1066)isinfo(2,i)                                 1d5s20
            write(6,*)('spin: '),isinfo(2,i),(' space symmetry: '),
     $        isinfo(1,i),('('),spinm,(char(isinfo(j,i)),j=6,nlab),(')')1d5s20
            write(6,*)('diagonalization converged at iteration '),iter
            iroosv=0                                                     5d3s21
            do j=6,min(nroolab,nlab)                                     5d3s21
             if(isinfo(j,i).ne.ibc(jroolab+j))iroosv=1                   5d3s21
            end do                                                       5d3s21
            if(nlab.ne.nroolab.or.iroosv.ne.0.or.spinm.ne.spinroo)then   5d3s21
             spinroo=spinm                                               5d3s21
             iroosv=1                                                    5d3s21
             nroolab=nlab                                                5d3s21
             do j=6,nlab                                                 5d3s21
              ibc(jroolab+j)=isinfo(j,i)                                 5d3s21
             end do                                                      5d3s21
            end if                                                       5d3s21
            ivecpsx=ibcoff                                               4d23s18
            npsx=nps+nrootu                                              4d23s18
            ibcoff=ivecpsx+npsx*nroot                                    4d23s18
            call enough('cas0. 51',bc,ibc)
            if(nps.gt.0.and.nokr.gt.0.and.npsx.gt.0)then                 2d22s19
             call dgemm('n','n',nps,nroot,nok,1d0,bc(ivecs),nps,          4d23s18
     $        bc(ivecp),nokr,0d0,bc(ivecpsx),npsx,                      4d23s18
     d' cas0.  6')
            end if                                                       2d22s19
            do i1=1,nroot
             write(ostng(1:22),1551)i1,bc(ieigp+i1-1)+shift             5d17s22
             if(iroosv.ne.0)then                                        5d3s21
              bc(jroodat)=bc(ieigp+i1-1)+shift                          5d17s22
              jroodat=jroodat+1                                          5d3s21
             end if                                                      5d3s21
             diff=bc(ieigp+i1-1)+shift0-bc(ieigo+i1-1)                  7d20s22
             if(abs(diff).gt.1d-13.and.nokr.lt.nconft)then               5d8s18
              ndig=nint(-log10(abs(diff)))+1                                 3d17s17
              ist=10+ndig
              do k=ist,22
               ostng(k:k)=' '
              end do
             end if
             write(6,1552)ostng(1:22)
 1552        format(a22)
             sum=0d0
             isrt=ibcoff                                                 4d23s18
             jsrt=isrt+nps                                               4d23s18
             ibcoff=jsrt+nps                                             4d23s18
             call enough('cas0. 52',bc,ibc)
             do j=0,nps-1                                                4d23s18
              iad=ivecpsx+j+npsx*(i1-1)                                  4d23s18
              bc(isrt+j)=-abs(bc(iad))                                   4d23s18
             end do                                                      4d23s18
             call dsortdws(bc(isrt),ibc(jsrt),nps)                      1d18s23
             do j=0,nps-1                                                4d23s18
              ju=ibc(jsrt+j)-1                                           4d23s18
              iad=ivecpsx+ju+npsx*(i1-1)                                 4d23s18
              if(abs(bc(iad)).gt.0.1d0)then                              4d23s18
               sum=sum+bc(iad)**2                                        4d23s18
               itrial=ibc(ips+ju)                                        5d7s18
               isb=itrial2(3)                                            4d23s18
               jsb=nsbeta(isb)                                           4d23s18
               kp=itrial2(1)                                             4d23s18
               jp=itrial2(2)                                             4d23s18
               call prtocc(bc(iad),isb,jsb,kp,jp,nadet,nbdet,            4d23s18
     $              ibc(iaorb),ibc(iborb),nalpha,nbeta,norb,iacto,nsymb,11d9s22
     $              bc,ibc)                                             11d9s22
              end if                                                     4d23s18
             end do                                                      4d23s18
             ibcoff=isrt                                                 4d23s18
            end do
           end if
          end if
          do k=0,nroot-1                                                  3d17s17
           bc(ieigo+k)=bc(ieigp+k)
          end do                                                          3d17s17
          ivecpsx=ibcoff                                                   3d17s17
          npsx=nps+nrootu                                                3d27s17
          ibcoff=ivecpsx+npsx*nroot                                        3d17s17
          call enough('cas0. 53',bc,ibc)
          if(nps.gt.0.and.nokr.gt.0.and.npsx.gt.0)then                   2d22s19
           call dgemm('n','n',nps,nroot,nok,1d0,bc(ivecs),nps,              3d17s17
     $        bc(ivecp),nokr,0d0,bc(ivecpsx),npsx,                        3d17s17
     d' cas0.  7')
          end if                                                         2d22s19
          do k=0,nroot-1                                                   3d17s17
           iad1=ivecpsx+nps+npsx*k                                         3d17s17
           iad2=ivecp+nok+nokr*k                                           3d17s17
           do l=0,nroot-1                                                  3d17s17
            bc(iad1+l)=bc(iad2+l)                                          3d17s17
           end do                                                          3d17s17
          end do                                                           3d17s17
          if(nroot.ne.nrootu)then                                        3d27s17
           do k=0,nroot-1                                                   3d17s17
            iad1=ivecpsx+nps+nroot+npsx*k                                3d27s17
            iad2=ivecp+nok+nroot+nokr*k                                  3d27s17
            do l=0,nroot-1                                                  3d17s17
             bc(iad1+l)=bc(iad2+l)                                          3d17s17
            end do                                                          3d17s17
           end do                                                           3d17s17
          end if                                                         3d27s17
          do isb=1,nsymb                                                   8d29s06
           nn=nroot*nbdet(nsbeta(isb))*nherec(isb)                       3d27s17
           do i1=1,nn                                                      8d29s06
            bc(ivec(isb)+i1-1)=0d0                                         8d29s06
           end do                                                          8d29s06
          end do                                                           8d29s06
          do i1=0,npsh-1                                                  3d17s17
           itrial=ibc(ipsh+i1)                                            3d17s17
           j2=ivec(itrial2(3))+nroot*(itrial2(1)-ilc(itrial2(3))+        3d27s17
     $         nherec(itrial2(3))*(itrial2(2)-1))                        3d16s17
           do i3=0,nroot-1                                               3d27s17
            j3=ivecpsx+i1+ioffps+i3*npsx                                  3d17s17
            bc(j2+i3)=bc(j3)                                              3d16s17
           end do
          end do
          ibcoff=idot
          ivecpxt=ibcoff                                                   3d17s17
          ibcoff=ivecpxt+nrootu*nroot                                    3d27s17
          call enough('cas0. 54',bc,ibc)
          do k=0,nroot-1                                                  3d17s17
           do l=0,nroot-1                                                3d27s17
            iad1=ivecpsx+nps+l+npsx*k                                      3d17s17
            iad2=ivecpxt+k+nroot*l                                        3d17s17
            bc(iad2)=bc(iad1)                                             3d17s17
           end do                                                         3d17s17
          end do                                                          3d17s17
          if(nrootu.ne.nroot)then                                        3d27s17
           do k=0,nroot-1                                                  3d17s17
            do l=0,nroot-1                                                3d27s17
             iad1=ivecpsx+nps+l+nroot+npsx*k                                      3d17s17
             iad2=ivecpxt+k+nroot*(l+nroot)                              3d27s17
             bc(iad2)=bc(iad1)                                             3d17s17
            end do                                                         3d17s17
           end do                                                          3d17s17
          end if                                                         3d27s17
          do isb=1,nsymb                                                   8d29s06
           nn=nbdet(nsbeta(isb))*nherec(isb)                               8d29s06
           if(nn.gt.0.and.nroot.gt.0.and.nrootu.gt.0)then                2d22s19
            call dgemm('n','n',nroot,nn,nrootu,1d0,bc(ivecpxt),nroot,    3d27s17
     $         bc(ivecx(isb)),nrootu,1d0,bc(ivec(isb)),nroot,           3d27s17
     d' cas0.  8')
           end if                                                          1d31s07
          end do                                                           8d29s06
          ibcoff=ivecp                                                     3d17s17
          do i3=0,ncona*nroot*2-1                                        3d27s17
           bc(iveca(1)+i3)=0d0                                            3d16s17
          end do                                                          3d16s17
          do isb=1,nsymb                                                  3d16s17
           jsb=nsbeta(isb)                                                3d16s17
           do ib=0,nbdet(jsb)-1                                           3d16s17
            do ia=0,nherec(isb)-1                                         3d16s17
             iap=ia+ilc(isb)-1                                            3d16s17
             iad1=ivec(isb)+nroot*(ia+nherec(isb)*ib)                    3d27s17
             iad2=iveca(isb)+nroot*(iap+nadet(isb)*ib)                   3d27s17
             do k=0,nroot-1                                               3d16s17
              bc(iad2+k)=bc(iad1+k)                                       3d16s17
             end do                                                       3d16s17
            end do                                                        3d16s17
           end do                                                         3d16s17
          end do                                                          3d16s17
          call dws_gsumf(bc(iveca(1)),ncona*nroot*2)                     3d27s17
          do isb=1,nsymb                                                  3d16s17
           jsb=nsbeta(isb)                                                3d16s17
           do ia=0,nadet(jsb)-1                                           3d16s17
            do ib=0,nherect(isb)-1                                        3d16s17
             ibp=ib+ilct(isb)-1                                           3d16s17
             iad1=iveca(jsb)+nroot*(ia+nadet(jsb)*ibp)                   3d27s17
             iad2=ivect(isb)+nroot*(ib+nherect(isb)*ia)                  3d27s17
             do k=0,nroot-1                                               3d16s17
              bc(iad2+k)=bc(iad1+k)                                       3d16s17
             end do                                                       3d16s17
            end do                                                        3d16s17
           end do                                                         3d16s17
          end do                                                          3d16s17
          call second(time2)
          telap=time2-time1-tovr
          timeend=timeend+telap
          if(nconv.eq.nroot)then                                          3d20s17
           call second(timec)                                             4d23s18
           telap=timec-timeb-tovr                                         4d23s18
           telapo(3)=telapo(3)+telap                                      4d23s18
c
c     diaghiid expects vector to be dimensioned nconf,nroot with
c     beta inner loop and alpha outer loop.                             3d31s17
c
           itmp=ibcoff                                                   3d27s17
           ibcoff=itmp+nconf*nroot                                       3d31s17
           call enough('cas0. 55',bc,ibc)
           ioff=0                                                        3d27s17
           do isb=1,nsymb                                                3d27s17
            jsb=nsbeta(isb)                                              3d27s17
            nn=nherec(isb)*nbdet(jsb)                                    3d27s17
            do ia=0,nherec(isb)-1
             do ib=0,nbdet(jsb)-1
              iad1=ivec(isb)+nroot*(ia+nherec(isb)*ib)                   3d31s17
              iad2=itmp+ioff+ib+nbdet(jsb)*ia                            3d31s17
              do k=0,nroot-1                                              3d27s17
               bc(iad2+k*nconf)=bc(iad1+k)                                3d27s17
              end do                                                      3d27s17
             end do                                                       3d27s17
            end do
            ioff=ioff+nn                                                 3d27s17
           end do                                                        3d27s17
           jsav=jcasvec
           nsav=0
           do isb=1,nsymb                                                   4d18s18
            jsb=nsbeta(isb)                                                 4d18s18
            iad1=jcasvec                                                    4d18s18
            iad2=ivec(isb)                                                  4d18s18
            nsav=nsav+nherec(isb)*nbdet(jsb)
            do j=0,nherec(isb)*nbdet(jsb)-1                                 4d18s18
             do k=0,nroot-1                                                 4d18s18
              bc(iad1+k)=bc(iad2+k)                                         4d18s18
             end do                                                         4d18s18
             iad1=iad1+nroot                                                4d18s18
             iad2=iad2+nroot
            end do                                                          4d18s18
            jcasvec=iad1                                                    4d18s18
           end do                                                           4d18s18
           nveccnt(i)=jsav                                              4d18s23
           if(nrydb.gt.0)then                                            5d7s18
            call genryd(jmats,kmats,bc(ih0mo),                           5d7s18
     $         ibc(iaorb),nalpha,ibc(iborb),nbeta,nsbeta,nconf,iveca,   5d7s18
     $         nvirtc,m12sym,idata1,idatb1,bc(ieigp),nroot,iorbn,morb,  12d31s19
     $         nbasisp,ihryd,i,nstate,lprint,nryd,ncomp,nlzz,iorbsym,   5d3s21
     $          iorbsymz,bc(iwgt),wsum,shift,bc,ibc,icanon)             5d5s23
           end if                                                        5d7s18
           call second(time1)                                            4d26s18
           telap=time1-timec-tovr
           telapo(15)=telapo(15)+telap
           call diaghiid(bc(itmp),bc(itmp),nconf,ibc(iaorb),nalpha,numa,5d16s22
     $              ibc(iborb),nbeta,numb,iden1,jdenpt,shift0,ilc,ihc,  5d17s22
     $              1d0,nsymb,nsbeta,nroot,ewght(1,i),.true.,bc,ibc,    3d15s23
     $          mdenoff,igoal)                                                3d15s23
           call second(time2)                                            4d26s18
           telap=time2-time1-tovr                                        4d26s18
           telapo(12)=telapo(12)+telap                                   4d26s18
           ibcoff=itmp                                                   3d27s17
           call second(time1)                                            4d26s18
           call hc1cd(ivec,ivect,iveca,ivec,ivect,iveca,ilc,ihc,ilct,   5d27s22
     $          ihct,nroot,nherec,nherect,ibc(iaorb),ibc(iborb),nalpha, 5d27s22
     $       nbeta,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,       3d27s17
     $       nsbeta,m1c,idatac,idatbc,iden1,jdenpt,ewght(1,i),.true.,   11d10s22
     $          bc,ibc,mdenoff,igoal)                                         3d15s23
           call second(timed)                                            4d23s18
           telap=timed-time1-tovr                                        4d26s18
           telapo(13)=telapo(13)+telap                                   4d26s18
           telap=timed-timec-tovr                                        4d23s18
           telapo(4)=telapo(4)+telap                                     4d23s18
           go to 4                                                        3d20s17
          end if                                                          3d20s17
          call second(time1)                                             4d9s18
          call hc1c(ivec,ivect,iveca,ig,igt,ilc,ihc,ilct,ihct,nroot,      3d10s17
     $       bc(ihdig),ih0e,nherec,nherect,ibc(iaorb),ibc(iborb),nalpha,3d10s17
     $       nbeta,iooooa,m12sym,norb,nsymb,idata1,idata2,idatb1,idatb2,3d14s17
     $       nsbeta,m1c,idatac,idatbc,mst,mnz,nz,lbail,bc,ibc)          1d6s23
          call second(time2)                                             4d9s18
          telap=time2-time1-tovr                                         4d9s18
          timehc1c=timehc1c+telap
          time1=time2                                                    4d9s18
          if(mod(iter,kick).ne.1)then                                    4d9s18
c
c     if degenerate roots, vectors will not neccessiarly be unique      10d28s20
c                                                                       10d28s20
           if(nroot.gt.0)then                                            10d28s20
            is=0                                                         10d28s20
  311       continue                                                     10d28s20
            e000=bc(ieigp+is)                                           10d28s20
            do itry=is,nroot-1                                          10d28s20
             if(abs(e000-bc(ieigp+itry)).gt.1d-10)then                  12d7s20
              ndggg=itry-is                                             10d28s20
              go to 312                                                 10d28s20
             end if                                                     10d28s20
            end do                                                      10d28s20
            ndggg=nroot-is                                              10d28s20
  312       continue                                                    10d28s20
            if(ndggg.gt.1)then                                          10d28s20
             ndggg2=ndggg**2                                            10d28s20
             iovlpr=ibcoff                                              10d28s20
             ibcoff=iovlpr+ndggg2                                       10d28s20
             call enough('cas0. 56',bc,ibc)
             itest=0                                                    10d28s20
  315        continue                                                   10d28s20
             itest=itest+1                                              10d28s20
             do ii=iovlpr,ibcoff-1                                       10d28s20
              bc(ii)=0d0                                                 10d28s20
             end do                                                     10d28s20
             do isb=1,nsymb                                             10d28s20
              nn=nherec(isb)*nbdet(nsbeta(isb))-1                       10d28s20
              do iab=0,nn                                               10d28s20
               iadk=ivec(isb)+is+nroot*iab                              10d28s20
               iadb=ivecb(isb)+is+nroot*iab                             10d28s20
               do ik=0,ndggg-1                                           10d28s20
                do ib=0,ndggg-1                                         10d28s20
                 ibb=iovlpr+ib+ndggg*ik                                 10d28s20
                 bc(ibb)=bc(ibb)+bc(iadb+ib)*bc(iadk+ik)                10d28s20
                end do                                                  10d28s20
               end do                                                   10d28s20
              end do                                                    10d28s20
             end do                                                     10d28s20
             call dws_gsumf(bc(iovlpr),ndggg2)                          10d28s20
             nsweep=0                                                   10d28s20
             isweep=0                                                   10d28s20
  313        continue                                                   10d28s20
             isweep=isweep+1                                            10d28s20
             if(isweep.gt.10)then
              write(6,*)('tooooooo many sweeps!!! '),offxsz             12d7s20
              go to 314                                                 12d7s20
             end if
             offxsz=0d0                                                 10d28s20
             do ik=0,ndggg-1                                             10d28s20
              do ib=0,ik-1                                              10d28s20
               ibk=iovlpr+ib+ndggg*ik                                   10d28s20
               ikb=iovlpr+ik+ndggg*ib                                   10d28s20
               offxsz=max(offxsz,abs(bc(ibk)),abs(bc(ikb)))              10d28s20
              end do                                                    10d28s20
             end do                                                     10d28s20
             call dws_bcast(offxsz,1)                                   10d28s20
             if(offxsz.lt.1d-7)go to 314                                10d31s20
             do ik=0,ndggg-1                                            10d28s20
              do ib=0,ik-1                                              10d28s20
               ibk=iovlpr+ib+ndggg*ik                                   10d28s20
               ikb=iovlpr+ik+ndggg*ib                                    10d28s20
               if(max(abs(bc(ikb)),abs(bc(ibk))).gt.thrsox)then         10d28s20
                ikk=iovlpr+ik+ndggg*ik                                  10d28s20
                ibb=iovlpr+ib+ndggg*ib                                  10d28s20
                t=(bc(ibb)*bc(ibk)-bc(ikk)*bc(ikb))/                     10d28s20
     $                 (bc(ibb)**2+bc(ikk)**2)                          10d28s20
                nsweep=nsweep+1                                         10d28s20
                ccx=1d0/sqrt(1d0+t*t)                                   10d28s20
                ssx=ccx*t                                               10d28s20
                crot11=ccx                                              10d28s20
                crot12=ssx                                              10d28s20
                crot21=-ssx                                             10d28s20
                crot22=ccx                                              10d28s20
                trybb=crot11*bc(ibb)+crot12*bc(ibk)                     10d28s20
                if(trybb.lt.0d0)then                                    10d28s20
                 crot11=-crot11                                         10d28s20
                 crot12=-crot12                                         10d28s20
                end if                                                  10d28s20
                trykk=crot21*bc(ikb)+crot22*bc(ikk)                     10d28s20
                if(trykk.lt.0d0)then                                    10d28s20
                 crot21=-crot21                                         10d28s20
                 crot22=-crot22                                         10d28s20
                end if                                                  10d28s20
                isendto=ibcoff                                          10d28s20
                ibcoff=isendto+3                                        10d28s20
                call enough('cas0. 57',bc,ibc)
                bc(isendto)=crot11                                      10d28s20
                bc(isendto+1)=crot21                                    10d28s20
                bc(isendto+2)=crot12                                    10d28s20
                bc(isendto+3)=crot22                                    10d28s20
                call dws_bcast(bc(isendto),4)                           10d28s20
                crot11=bc(isendto)                                      10d28s20
                crot21=bc(isendto+1)                                    10d28s20
                crot12=bc(isendto+2)                                    10d28s20
                crot22=bc(isendto+3)                                    10d28s20
                ibcoff=isendto                                          10d28s20
                do io=0,ndggg-1                                         10d28s20
                 iadob=iovlpr+io+ndggg*ib                               10d28s20
                 iadok=iovlpr+io+ndggg*ik                               10d28s20
                 xnewb=crot11*bc(iadob)+crot12*bc(iadok)                10d28s20
                 xnewk=crot21*bc(iadob)+crot22*bc(iadok)                10d28s20
                 bc(iadob)=xnewb                                        10d28s20
                 bc(iadok)=xnewk                                        10d28s20
                end do                                                  10d28s20
                call dws_bcast(bc(iovlpr),ndggg2)                       10d28s20
                do isb=1,nsymb                                          10d28s20
                 nn=nherec(isb)*nbdet(nsbeta(isb))-1                       10d28s20
                 do iab=0,nn                                               10d28s20
                  iadk=ivec(isb)+is+ik+nroot*iab                        10d28s20
                  iadb=ivec(isb)+is+ib+nroot*iab                        10d28s20
                  xnewb=crot11*bc(iadb)+crot12*bc(iadk)                 10d28s20
                  xnewk=crot21*bc(iadb)+crot22*bc(iadk)                 10d28s20
                  bc(iadb)=xnewb                                        10d28s20
                  bc(iadk)=xnewk                                        10d28s20
                 end do                                                 10d28s20
                end do                                                  10d28s20
               end if                                                   10d28s20
              end do                                                    10d28s20
             end do                                                     10d28s20
             go to 313                                                  10d28s20
  314        continue                                                   10d28s20
             ibcoff=iovlpr
            end if                                                      10d28s20
            is=is+ndggg                                                 10d28s20
            if(is.lt.nroot)go to 311                                     10d28s20
           end if                                                        10d28s20
           idot=ibcoff                                                     3d17s17
           ibcoff=idot+nroot*2                                           5d11s18
           call enough('cas0. 58',bc,ibc)
           do k=0,2*nroot-1                                              5d11s18
            bc(idot+k)=0d0                                                 3d17s17
           end do                                                          3d17s17
           jdot=idot+nroot                                               5d11s18
           do isb=1,nsymb                                                  3d17s17
            nn=nherec(isb)*nbdet(nsbeta(isb))-1                            3d17s17
            nsumt=nsumt+nn+1
            do iab=0,nn                                                    3d17s17
             iad1=ivec(isb)+nroot*iab                                      3d17s17
             iad2=ivecb(isb)+nroot*iab                                     3d17s17
             do k=0,nroot-1                                                3d17s17
              bc(idot+k)=bc(idot+k)+bc(iad1+k)*bc(iad2+k)                  3d17s17
             end do                                                        3d17s17
             do k=0,nroot-2                                              5d11s18
              bc(jdot+k)=bc(jdot+k)+bc(iad1+k)*bc(iad2+k+1)              5d11s18
             end do                                                      5d11s18
            end do                                                         3d17s17
           end do                                                          3d17s17
           call dws_gsumf(bc(idot),nroot*2)                              5d11s18
           do k=0,nroot-1                                                  3d17s17
            bc(idot+k)=bc(idot+k)-1d0                                      3d17s17
            if(bc(idot+k).gt.1d-7)then                                   3d3s21
             write(6,*)('bad dot!!! '),bc(idot+k)+1d0
             call dws_synca                                              3d3s21
             call dws_finalize                                           3d3s21
             stop                                                        3d3s21
            end if                                                       3d3s21
            if(abs(bc(idot+k)).lt.vconv)then                               3d17s17
             ibc(iconv+k)=1                                                3d17s17
            else                                                           3d17s17
             ibc(iconv+k)=0                                                 3d17s17
            end if                                                          3d17s17
           end do                                                          3d17s17
           ibcoff=idot                                                     3d17s17
          end if                                                         4d9s18
c
c     this is best vec so far
c
          do isb=1,nsymb                                                   3d17s17
           nn=nherec(isb)*nbdet(nsbeta(isb))*nroot-1                       3d17s17
           do k=0,nn                                                       3d17s17
            bc(ivecb(isb)+k)=bc(ivec(isb)+k)                               3d17s17
           end do                                                          3d17s17
          end do                                                           3d17s17
         end if                                                          3d20s17
c
c     update
c
         igfull=ibcoff                                                    3d12s17
         ibcoff=igfull+nroot*ncona                                        3d12s17
         call enough('cas0. 59',bc,ibc)
         do i3=0,nroot*ncona-1                                            3d12s17
          bc(igfull+i3)=0d0                                               3d12s17
         end do                                                           3d12s17
         jgfull=igfull                                                    3d12s17
         do isb=1,nsymb                                                   3d12s17
          jsb=nsbeta(isb)                                                 3d12s17
          iga(isb)=jgfull                                                 3d12s17
          do icol=0,nbdet(jsb)-1                                          3d14s17
           do irow=0,nherec(isb)-1                                        3d14s17
            iad1=ig(isb)+nroot*(irow+nherec(isb)*icol)                    3d14s17
            iad2=jgfull+nroot*(irow+ilc(isb)-1+nadet(isb)*icol)           3d14s17
            do i3=0,nroot-1                                               3d14s17
             bc(iad2+i3)=bc(iad1+i3)                                      3d14s17
            end do                                                        3d14s17
           end do                                                         3d14s17
          end do                                                          3d14s17
          jgfull=jgfull+nadet(isb)*nroot*nbdet(jsb)                       3d12s17
         end do                                                           3d12s17
         do jsb=1,nsymb                                                   3d12s17
          isb=nsbeta(jsb)                                                 3d12s17
          do icol=0,nadet(isb)-1                                          3d14s17
           do irow=0,nherect(jsb)-1                                       3d14s17
            iad1=igt(jsb)+nroot*(irow+nherect(jsb)*icol)                  3d14s17
            iad2=iga(isb)+nroot*(icol+nadet(isb)*(irow+ilct(jsb)-1))      3d14s17
            do i3=0,nroot-1                                               3d14s17
             bc(iad2+i3)=bc(iad2+i3)+bc(iad1+i3)                          3d14s17
            end do                                                        3d14s17
           end do                                                         3d14s17
          end do                                                          3d14s17
         end do
         nwds=nroot*ncona                                                 3d12s17
         call dws_gsumf(bc(igfull),nwds)                                  3d12s17
         if(iter.eq.1)then
         isumq=ibcoff
         ibcoff=isumq+nrootu
         call enough('cas0. 60',bc,ibc)
         do iz=isumq,ibcoff-1
          bc(iz)=0d0
         end do
         do isb=1,nsymb
          jsb=nsbeta(isb)
          iad1=iga(isb)
          iad2=iveca(isb)
          do iiii=0,nadet(isb)*nbdet(jsb)-1
           do ir=0,nrootu-1
            bc(isumq+ir)=bc(isumq+ir)+bc(iad1+ir)*bc(iad2+ir)
           end do
           iad1=iad1+nrootu
           iad2=iad2+nrootu
          end do
         end do
        end if
         idot=ibcoff
         ibcoff=idot+nroot
         call enough('cas0. 61',bc,ibc)
         do k=0,nroot-1
          bc(idot+k)=0d0
         end do
         do isb=1,nsymb                                                   3d17s17
          jsb=nsbeta(isb)                                                 3d17s17
          if(nadet(isb)*nbdet(jsb).gt.0.and.idwsdeb.ne.0)then
           write(6,*)('g & v for symmetry block '),isb
           do k=0,nroot-1
            write(6,*)('for root no. '),k
            call prntm2(bc(iga(isb)+k),1,nadet(isb)*nbdet(jsb),
     $         nroot)
            call prntm2(bc(iveca(isb)+k),1,nadet(isb)*nbdet(jsb),
     $         nroot)
           end do
          end if
          do ib=0,nbdet(jsb)-1                                            3d17s17
           do ia=0,nadet(isb)-1                                          3d20s17
            iad1=iveca(isb)+nroot*(ia+nadet(isb)*ib)
            iad2=iga(isb)+nroot*(ia+nadet(isb)*ib)
            do k=0,nroot-1                                               3d20s17
             bc(idot+k)=bc(idot+k)+bc(iad1+k)*bc(iad2+k)
            end do
           end do
           do ia=0,nherec(isb)-1                                          3d17s17
            iap=ia+ilc(isb)-1                                             3d17s17
            iad1=iga(isb)+nroot*(iap+nadet(isb)*ib)                       3d17s17
            iad2=ig(isb)+nroot*(ia+nherec(isb)*ib)                        3d17s17
            do k=0,nroot-1                                                3d17s17
             bc(iad2+k)=bc(iad1+k)                                        3d17s17
            end do                                                        3d17s17
           end do                                                         3d17s17
          end do                                                          3d17s17
         end do                                                           3d17s17
         if(iter.eq.1)then
          if(.not.bcasvec)then                                           4d23s18
           do i1=0,nroot-1                                               4d23s18
            bc(ieigp+i1)=bc(idot+i1)                                     4d23s18
           end do                                                        4d23s18
          end if                                                         4d23s18
         end if
         if(iter.gt.1)then                                               3d27s17
          nrootu=nroot*2                                                 3d27s17
          do isb=1,nsymb                                                 3d27s17
           jsb=nsbeta(isb)                                               3d27s17
           nn=nroot*nherec(isb)*nbdet(jsb)-1                             3d27s17
           do k=0,nn                                                     3d27s17
            bc(iveco(isb)+k)=bc(ivec(isb)+k)                             3d27s17
           end do                                                        3d27s17
           nn=nherec(isb)*nbdet(jsb)-1                                   3d27s17
           do iab=0,nn                                                   3d27s17
            iad1=iveco(isb)+nroot*iab                                    3d27s17
            iad2=ivec(isb)+nrootu*iab                                    3d27s17
            do k=0,nroot-1                                               3d27s17
             bc(iad2+k)=bc(iad1+k)                                       3d27s17
            end do                                                       3d27s17
           end do                                                        3d27s17
          end do
         end if
         ibcoff=igfull                                                    3d17s17
         do i3=0,nroot-1                                                 3d12s17
          energy=bc(ieigp+i3)                                            3d12s17
          if(iprtr(8).ne.0)write(6,*)('energy = '),energy,bc(idot+i3)
          dot=0d0                                                        4d21s06
          dot2=0d0
          ioff=0                                                         8d29s06
          do isb=1,nsymb                                                 8d29s06
           if(iprtr(8).ne.0)write(6,*)('isb = '),isb
           jsb=nsbeta(isb)                                               3d12s17
           nn=nbdet(jsb)*nherec(isb)                                     3d16s17
           iad1=ig(isb)+i3                                               3d17s17
           iad2=ivec(isb)+i3                                             3d17s17
           do ib=0,nbdet(jsb)-1                                          3d17s17
            do ia=0,nherec(isb)-1                                        3d17s17
             iab=ia+nherec(isb)*ib                                       3d17s17
             iba=ib+nbdet(jsb)*ia                                        3d17s17
             j2=iad1+iab*nroot                                             3d16s17
             k2=iad2+iab*nrootu                                          3d27s17
             iad3=ips2+ioff+iba                                          3d17s17
             sav=bc(k2)
             sav2=bc(j2)
             l2=ihdig+ioff+iba                                           3d17s17
             resid=bc(j2)-energy*bc(k2)                                   4d26s06
             bot=bc(l2)-energy                                            4d26s06
             if(resid.eq.0d0)then                                       5d2s18
              update=0d0                                                5d2s18
             else                                                       5d2s18
              update=-resid/bot                                            4d26s06
             end if                                                     5d2s18
             if(idwsdeb.ne.0.or.iprtr(8).ne.0)
     $         write(6,3352)ia,ib,bc(j2),energy*bc(k2),resid,bc(l2),bot,
     $             update,bc(k2),bc(k2)+update
 3352        format(2i5,8es22.14)
             if(iter.gt.1)then                                          3d27s17
              bc(k2+nroot)=bc(k2)+update                                3d27s17
              dot=dot+bc(k2+nroot)**2                                             4d26s06
             else                                                       3d27s17
              bc(k2)=bc(k2)+update                                         4d21s06
              dot=dot+bc(k2)**2                                             4d26s06
             end if                                                     3d27s17
32312        format(i5,i2,1p7e15.7,i8,5i5)
            end do                                                       3d17s17
           end do                                                         4d21s06
           ioff=ioff+nn                                                  8d29s06
          end do                                                         8d29s06
          call dws_gsumf(dot,1)                                          3d17s17
          if(iprtr(8).ne.0)write(6,*)('norm of update for root '),i3,
     $        ('is '),
     $        dot
          if(dot.eq.0d0.or.dot.ne.dot)dot=1d0                            5d2s18
           doti=1d0/sqrt(dot)
           do isb=1,nsymb                                                 8d29s06
            nn=nbdet(nsbeta(isb))*nherec(isb)                             8d29s06
            iad1=ivec(isb)+i3-nrootu                                      3d27s17
            if(iter.gt.1)iad1=iad1+nroot                                  3d27s17
            do i2=1,nn                                                    8d29s06
             k2=iad1+i2*nrootu                                            3d27s17
             bc(k2)=bc(k2)*doti                                            4d26s06
            end do
           end do                                                         8d29s06
          end do
          do i1=0,npsh-1
           itrial=ibc(ipsh+i1)
           j2=ivec(itrial2(3))+nrootu*(itrial2(1)-ilc(itrial2(3))+        3d27s17
     $        nherec(itrial2(3))*(itrial2(2)-1))
           do k=0,nrootu-1                                                3d27s17
            bc(j2+k)=0d0
           end do
          end do
c
c     orthogonalize
c
c     p-space configs knocked up by ips2 test above
          do isb=1,nsymb                                                  8d29s06
           numbx=nbdet(nsbeta(isb))                                       8d29s06
           do l=1,numbx                                                   8d29s06
            do j=1,nherec(isb)                                             8d29s06
             do k=1,nrootu                                                   4d26s06
              iad1=ivec(isb)+k-1+nrootu*(j-1+nherec(isb)*(l-1))           3d27s17
              iad2=ig(isb)+l-1+numbx*(j-1+nherec(isb)*(k-1))              3d16s17
              bc(iad2)=bc(iad1)                                            4d26s06
             end do                                                        4d26s06
            end do                                                         4d26s06
           end do                                                          4d26s06
          end do                                                          3d16s17
          do i3=1,nrootu                                                  3d27s17
           i3m=i3-1                                                       3d16s17
           if(iprtr(8).ne.0)write(6,*)('orthogonalize root '),i3,nconf    5d12s20
           iloopit=0                                                      4d9s18
 8892      continue                                                       2d17s07
           iloopit=iloopit+1                                              4d9s18
           if(iloopit.gt.10)then                                          4d9s18
            write(6,*)('iloopit is too large in cas0 ')                   4d9s18
            write(6,*)('for root no. '),i3                                4d9s18
            call dws_synca                                                 4d9s18
            call dws_finalize                                             4d9s18
            stop                                                          4d9s18
           end if                                                         4d9s18
           do i2=1,i3m                                                    3d16s17
            i2m=i2-1                                                      3d16s17
            dot=0d0                                                       4d21s06
            do isb=1,nsymb                                                8d29s06
             numbx=nbdet(nsbeta(isb))                                     3d17s17
             nn=numbx*nherec(isb)                                         3d17s17
             nnm=nn-1                                                     3d16s17
             iad3=ig(isb)+nn*i3m                                          3d16s17
             iad2=ig(isb)+nn*i2m                                          3d16s17
             do i4=0,nnm                                                  3d16s17
              dot=dot+bc(iad3+i4)*bc(iad2+i4)                             3d16s17
             end do                                                        4d21s06
            end do                                                        8d29s06
            iarg1=1                                                       4d26s06
            call dws_gsumf(dot,iarg1)                                     3d16s17
            if(iprtr(8).ne.0)write(6,*)('dot after gsumf '),i3,i2,dot     5d12s20
            do isb=1,nsymb                                                8d29s06
             numbx=nbdet(nsbeta(isb))                                     8d29s06
             nn=numbx*nherec(isb)                                         3d16s17
             nnm=nn-1                                                     3d16s17
             iad3=ig(isb)+nn*i3m                                          3d16s17
             iad2=ig(isb)+nn*i2m                                          3d16s17
             do i4=0,nnm                                                  3d16s17
              bc(iad3+i4)=bc(iad3+i4)-dot*bc(iad2+i4)                     3d16s17
             end do                                                       3d16s17
            end do                                                        8d29s06
           end do                                                         4d21s06
           sum=0d0                                                        4d21s06
           do isb=1,nsymb                                                 8d29s06
            numbx=nbdet(nsbeta(isb))                                      8d29s06
            nn=numbx*nherec(isb)                                          3d16s17
            nnm=nn-1                                                      3d16s17
            iad3=ig(isb)+nn*i3m                                           3d16s17
            do i4=0,nnm                                                   3d16s17
             sum=sum+bc(iad3+i4)**2                                       3d16s17
            end do                                                        3d16s17
           end do                                                         4d21s06
           iarg1=1                                                        5d30s06
           call dws_gsumf(sum,iarg1)                                      3d16s17
           if(iprtr(8).ne.0)                                              5d12s20
     $        write(6,*)('sum after projecting out best so far '),sum
           if(sum.lt.1d-20)then                                           4d6s18
            if(iprtr(8).ne.0)                                             5d12s20
     $        write(6,*)('are we already an eigenvector? '),sum         2d14s20
            sum=0d0                                                       2d14s20
            ibc(iconv+i3-1)=1                                             4d6s18
           end if                                                         4d6s18
           if(ibc(iconv+i3-1).eq.0.and.(sum.lt.1d-10.or.                  5d8s18
     $        (i3.gt.nroot.and.mod(iter,kick).eq.0.and.                 5d8s18
     $        iloopit.eq.1)))then                                       5d8s18
            if(iprtr(8).ne.0)
     $        write(6,*)('in random block ')                            5d12s20
            if(iseed.gt.0)then                                            4d20s18
             if(lprint)then                                               4d20s18
              write(6,*)('oops: for root '),i3                          4d27s23
              write(6,*)('norm of orthogonalized guess = '),sum,        4d27s23
     $             ('thus use random number generator instead')         4d27s23
             end if                                                       4d20s18
             call randset(iseed)                                          3d17s17
             iseed=0                                                      3d17s17
            end if                                                        3d17s17
            do isb=1,nsymb                                                2d17s07
             numbx=nbdet(nsbeta(isb))                                     2d17s07
             do ia=0,nherec(isb)-1                                        3d27s17
              do ib=0,numbx-1                                             3d27s17
               iad=ig(isb)+ib+numbx*(ia+nherec(isb)*i3m)                 4d9s18
               bc(iad)=randdws()-0.5d0                                   3d27s17
              end do                                                      2d17s07
             end do                                                       2d17s07
            end do                                                        2d17s07
            do i1=0,npsh-1                                                  4d13s20
             itrial=ibc(ipsh+i1)                                            4d13s20
             numbx=nbdet(nsbeta(itrial2(3)))                                10d16s20
             j2=ig(itrial2(3))+itrial2(2)-1+numbx*(                         10d16s20
     $        itrial2(1)-ilc(itrial2(3))+                               10d16s20
     $        nherec(itrial2(3))*i3m)                                   10d16s20
             bc(j2)=0d0                                                     10d16s20
            end do                                                          4d13s20
            go to 8892                                                    2d17s07
           end if
           if(iprtr(8).ne.0)
     $        write(6,*)('normalizing with sum = '),sum                 5d12s20
           if(sum.eq.0d0)then                                             9d16s19
            sumi=0d0                                                      9d16s19
           else                                                           9d16s19
            sumi=1d0/sqrt(sum)
           end if                                                         9d16s19
           do isb=1,nsymb                                                 8d29s06
            numbx=nbdet(nsbeta(isb))                                      8d29s06
            nn=numbx*nherec(isb)                                          3d16s17
            nnm=nn-1                                                      3d16s17
            iad3=ig(isb)+nn*i3m                                           3d16s17
            do i4=0,nnm                                                   3d16s17
             bc(iad3+i4)=bc(iad3+i4)*sumi                                 3d16s17
            end do                                                        3d16s17
           end do                                                         8d29s06
          end do                                                          4d21s06
          idot=ibcoff
          ibcoff=idot+nrootu
          call enough('cas0. 62',bc,ibc)
          do ii=idot,ibcoff-1
           bc(ii)=0d0
          end do
          jdot=idot-1                                                     10d16s20
          do isb=1,nsymb                                                  8d29s06
           numbx=nbdet(nsbeta(isb))                                       8d29s06
           do l=1,numbx                                                   8d29s06
            do j=1,nherec(isb)                                             8d29s06
             do k=1,nrootu                                                   4d26s06
              iad1=ivec(isb)+k-1+nrootu*(j-1+nherec(isb)*(l-1))           3d27s17
              iad2=ig(isb)+l-1+numbx*(j-1+nherec(isb)*(k-1))              3d16s17
              bc(jdot+k)=bc(jdot+k)+bc(iad2)**2
              bc(iad1)=bc(iad2)                                            4d26s06
             end do                                                        4d26s06
            end do                                                         4d26s06
           end do                                                          4d26s06
          end do                                                          3d16s17
          call dws_gsumf(bc(idot),nrootu)
          do i3=0,ncona*nroot*2-1                                         3d27s17
           bc(iveca(1)+i3)=0d0                                            3d16s17
          end do                                                          3d16s17
          do isb=1,nsymb                                                  3d16s17
           jsb=nsbeta(isb)                                                3d16s17
           do ib=0,nbdet(jsb)-1                                           3d16s17
            do ia=0,nherec(isb)-1                                         3d16s17
             iap=ia+ilc(isb)-1                                            3d16s17
             iad1=ivec(isb)+nrootu*(ia+nherec(isb)*ib)                    3d27s17
             iad2=iveca(isb)+nrootu*(iap+nadet(isb)*ib)                   3d27s17
             do k=0,nrootu-1                                               3d16s17
              bc(iad2+k)=bc(iad1+k)                                       3d16s17
             end do                                                       3d16s17
            end do                                                        3d16s17
           end do                                                         3d16s17
          end do                                                          3d16s17
          call dws_gsumf(bc(iveca(1)),ncona*nroot*2)                      3d27s17
          do isb=1,nsymb                                                  3d16s17
           jsb=nsbeta(isb)                                                3d16s17
           do ia=0,nadet(jsb)-1                                           3d16s17
            do ib=0,nherect(isb)-1                                        3d16s17
             ibp=ib+ilct(isb)-1                                           3d16s17
             iad1=iveca(jsb)+nrootu*(ia+nadet(jsb)*ibp)                   3d27s17
             iad2=ivect(isb)+nrootu*(ib+nherect(isb)*ia)                  3d27s17
             do k=0,nrootu-1                                              3d27s17
              bc(iad2+k)=bc(iad1+k)                                       3d16s17
             end do                                                       3d16s17
            end do                                                        3d16s17
           end do                                                         3d16s17
          end do                                                          3d16s17
          call second(time2)                                              4d9s18
          telap=time2-time1-tovr                                          4d9s18
          timeorth=timeorth+telap                                         4d9s18
          go to 2
    4     continue
          do j=1,nroot                                                     3d16s21
           ewght(j,i)=bc(iwgt+j-1)                                      5d11s21
          end do                                                           3d16s21
          ibcoff=ibctop                                                    4d18s18
         end if                                                           12d28s19
       end if                                                           4d27s21
       mdenoff=mdenoff+((nroot*(nroot+1))/2)                            3d14s23
      end do
      sum=0d0
      nbodc=1                                                           6d2s22
      if(mode.eq.3)nbodc=4                                              6d2s22
      if(mode.eq.4)nbodc=0                                              7d7s22
      do im=1,nbodc                                                     7d18s22
       do isb=1,nsymb                                                    12d28s19
        if(im.ne.4)then                                                 7d18s22
         jsb=multh(isb,ipuse)                                             7d18s22
        else                                                            7d18s22
         jsb=isb                                                        7d18s22
        end if                                                          7d18s22
        if(min(iacto(isb),iacto(jsb)).gt.0)then                          7d18s22
         nsum=iacto(isb)*iacto(jsb)                                     3d16s23
         if(im.ne.1)then                                                3d16s23
          nsumu=nsum*mden                                               3d16s23
         else                                                           3d16s23
          nsumu=nsum                                                    3d16s23
         end if                                                         3d16s23
         call dws_gsumf(bc(iden1(isb,im)),nsumu)                        3d16s23
        end if                                                           12d28s19
       end do                                                           6d2s22
      end do                                                            12d28s19
      do isb=1,nsymb
       if(iacto(isb).gt.0.and.mode.ne.4)then                            7d8s22
        if(iacto(isb).gt.1.and.ipuse.eq.1)then                           7d8s22
         do j=0,iacto(isb)-1                                             8d19s14
          do k=0,j-1                                                     8d19s14
           iad1=iden1(isb,1)+k+iacto(isb)*j                              6d2s22
           iad2=iden1(isb,1)+j+iacto(isb)*k                              6d2s22
           avg=0.5d0*(bc(iad1)+bc(iad2))                                 8d19s14
           bc(iad1)=avg                                                  8d19s14
           bc(iad2)=avg                                                  8d19s14
          end do                                                         8d19s14
         end do                                                          8d19s14
  310    format(10es26.16)
        end if
        if(lprt)then                                                    1d3s20
         write(6,*)('one particle density for symmetry block '),isb     4d20s18
         call prntm2(bc(iden1(isb,1)),iacto(isb),iacto(isb),iacto(isb)) 6d2s22
        end if                                                          1d3s20
        sz=0d0                                                           5d20s22
        do j=0,iacto(isb)*iacto(isb)-1                                   8d18s14
         sum=sum+bc(iden1(isb,1)+j)*bc(ih0e(isb)+j)                      6d2s22
         sz=sz+bc(iden1(isb,1)+j)**2                                     6d2s22
        end do                                                           8d18s14
        sz=sqrt(sz/dfloat(iacto(isb)*iacto(isb)))                        5d20s22
        if(sz.lt.1d-14.and.mode.eq.1)then                                5d20s22
         write(6,*)('one particle density for symmetry '),isb            5d20s22
         write(6,*)('is the null matrix')                                5d20s22
         write(6,*)('I can not continue in this situation.')             5d20s22
         write(6,*)
     $        ('I suggest you eliminate orbitals of this symmetry ')
     $       ,('from you active space and try again.')                  5d20s22
         call dws_synca                                                  5d20s22
         call dws_finalize                                               5d20s22
         stop 'cas0'                                                     5d20s22
        end if                                                           5d20s22
       end if                                                            8d18s14
      end do
      if(ipuse.eq.1)then                                                6d16s22
       do ib=1,nsdlk                                                     12d28s19
        n1=isblk(1,ib)                                                   12d28s19
        n2=isblk(2,ib)                                                   12d28s19
        n3=isblk(3,ib)                                                   12d28s19
        n4=isblk(4,ib)                                                   12d28s19
        if(n1.eq.n2.and.mode.ne.4)then                                  7d14s22
         nnn=(iacto(n1)*(iacto(n1)+1))/2                                 12d28s19
        else                                                             12d28s19
         nnn=iacto(n1)*iacto(n2)                                         12d28s19
        end if                                                           12d28s19
        if(n3.eq.n4.and.mode.ne.4)then                                  7d14s22
         mmm=(iacto(n3)*(iacto(n3)+1))/2                                 12d28s19
        else                                                             12d28s19
         mmm=iacto(n3)*iacto(n4)                                         12d28s19
        end if                                                           12d28s19
        if(nnn*mmm.gt.0)then                                             12d28s19
         if(mode.eq.4)then                                              3d16s23
          call dws_gsumf(bc(jdenpt(ib)),nnn*mmm*mden)                   3d16s23
         else                                                           3d16s23
          call dws_gsumf(bc(jdenpt(ib)),nnn*mmm)                          12d28s19
         end if                                                         3d16s23
         if(lprt)then                                                    1d3s20
          write(6,*)('starting value of '),(isblk(j,ib),j=1,4),ib,
     $       jdenpt(ib)
          call prntm2(bc(jdenpt(ib)),nnn,mmm,nnn)
         end if                                                          1d3s20
        end if                                                           12d28s19
       end do                                                            12d28s19
       if(mode.ne.4)then                                                7d8s22
        sum2=0d0
        do ib=1,nsdlk                                                     8d19s14
         n1=isblk(1,ib)                                                   8d19s14
         n2=isblk(2,ib)                                                   8d19s14
         n3=isblk(3,ib)
         n4=isblk(4,ib)                                                   8d19s14
         if(n1.eq.n2)then
          nnn=(iacto(n1)*(iacto(n1)+1))/2
         else
          nnn=iacto(n1)*iacto(n2)
         end if
         if(n3.eq.n4)then
          mmm=(iacto(n3)*(iacto(n3)+1))/2
         else
          mmm=iacto(n3)*iacto(n4)
         end if
         if(nnn*mmm.gt.0)then
          if(min(n1,n2).eq.min(n3,n4).and.max(n1,n2).eq.max(n3,n4))then   8d20s14
           do j=0,nnn-1                                                    8d20s14
            do k=0,j-1                                                     8d20s14
             iad1=jdenpt(ib)+k+nnn*j                                       8d20s14
             iad2=jdenpt(ib)+j+nnn*k                                       8d20s14
             avg=0.5d0*(bc(iad1)+bc(iad2))                                 8d20s14
             bc(iad1)=avg                                                  8d20s14
             bc(iad2)=avg                                                  8d20s14
            end do                                                         8d20s14
           end do                                                          8d20s14
          end if                                                           8d20s14
          if(.not.(n1.eq.n2.and.n3.eq.n4.and.n1.eq.n3))then                9d18s14
           do io=ib+1,nsdlk                                                 9d18s14
            if(isblk(1,io).eq.n3.and.isblk(2,io).eq.n4.and.                 9d18s14
     $       isblk(3,io).eq.n1.and.isblk(4,io).eq.n2)then               9d18s14
             do j=0,nnn-1                                                    9d18s14
              do k=0,mmm-1                                                   9d18s14
               iad1=jdenpt(ib)+j+nnn*k                                       9d18s14
               iad2=jdenpt(io)+k+mmm*j                                       9d18s14
               avg=0.5d0*(bc(iad1)+bc(iad2))                                 9d18s14
               bc(iad1)=avg                                                  9d18s14
               bc(iad2)=avg                                                  9d18s14
              end do                                                         9d18s14
             end do                                                          9d18s14
            end if                                                           9d18s14
           end do                                                            9d18s14
          end if                                                            9d18s14
         end if                                                            9d18s14
        end do                                                            9d18s14
        do ib=1,nsdlk                                                     8d19s14
         n1=isblk(1,ib)                                                   8d19s14
         n2=isblk(2,ib)                                                   8d19s14
         n3=isblk(3,ib)
         n4=isblk(4,ib)                                                   8d19s14
         if(n1.eq.n2)then
          nnn=(iacto(n1)*(iacto(n1)+1))/2
         else
          nnn=iacto(n1)*iacto(n2)
         end if
         if(n3.eq.n4)then
          mmm=(iacto(n3)*(iacto(n3)+1))/2
         else
          mmm=iacto(n3)*iacto(n4)
         end if
         if(nnn*mmm.gt.0)then
          orig=sum2
          do j=0,nnn*mmm-1                                                8d19s14
           sum2=sum2+bc(iooooa(ib)+j)*bc(jdenpt(ib)+j)                    8d19s14
          end do                                                          8d19s14
          if(n1.eq.n2)then                                                6d16s22
           do i2=0,iacto(n3)-1                                            4d10s17
            do i1=0,i2                                                    4d10s17
             fact=0.5d0                                                   4d10s17
             icol=jdenpt(ib)+nnn*(((i2*(i2+1))/2)+i1)                     4d10s17
             if(i1.eq.i2)fact=1d0                                         4d10s17
             facth=fact*0.5d0                                             4d10s17
             do i4=0,iacto(n1)-1                                          4d10s17
              do i3=0,i4-1                                                4d10s17
               bc(icol)=bc(icol)*facth                                    4d10s17
               icol=icol+1                                                4d10s17
              end do                                                      4d10s17
              bc(icol)=bc(icol)*fact                                      4d10s17
              icol=icol+1                                                 4d10s17
             end do                                                       4d10s17
            end do                                                        4d10s17
           end do                                                         4d10s17
          else
           itmp=ibcoff
           ibcoff=itmp+nnn*mmm
           call enough('cas0. 63',bc,ibc)
           do i4=0,iacto(n4)-1
            do i3=0,iacto(n3)-1
             do i2=0,iacto(n2)-1
              do i1=0,iacto(n1)-1
               iad1=jdenpt(ib)+i1+iacto(n1)*(i2+iacto(n2)
     $            *(i3+iacto(n3)*i4))
               iad2=itmp+i2+iacto(n2)*(i1+iacto(n1)*(i4+iacto(n4)*i3))
               bc(iad2)=bc(iad1)*0.25d0
              end do
             end do
            end do
           end do
           ibcoff=itmp
          end if                                                          4d10s17
          if(lprt)then
           write(6,*)('2 part density for symmetry '),
     $          (isblk(j,ib),j=1,4),jdenpt(ib)
           call prntm2(bc(jdenpt(ib)),nnn,mmm,nnn)
          end if                                                          1d3s20
         end if                                                           8d19s14
        end do
       end if                                                           7d8s22
      else                                                              7d8s22
       do is=1,nsblkderd                                                7d1s22
        ncol=iacto(isblkderd(3,is))*iacto(isblkderd(4,is))              7d1s22
        nrow=iacto(isblkderd(1,is))*iacto(isblkderd(2,is))              7d1s22
        call dws_gsumf(bc(jdenpt(is)),ncol*nrow)                        7d18s22
       end do                                                           7d1s22
       do is=1,nsblkderd                                                7d1s22
        if(isblkderd(1,is).eq.isblkderd(2,is).and.                      7d1s22
     $       iacto(isblkderd(1,is)).gt.0)then                            7d1s22
         ncol=iacto(isblkderd(3,is))*iacto(isblkderd(4,is))             7d1s22
         nrow=iacto(isblkderd(1,is))**2                                 7d1s22
         do icol=0,ncol-1                                               6d17s22
          j2=jdenpt(is)+nrow*icol                                       6d17s22
          do i2=0,iacto(isblkderd(1,is))-1                              7d1s22
           do i1=0,i2-1                                                 6d17s22
            i12=i1+iacto(isblkderd(1,is))*i2                            7d1s22
            i21=i2+iacto(isblkderd(1,is))*i1                            7d1s22
            o1=bc(j2+i12)
            o2=bc(j2+i21)
            avg=0.5d0*(bc(j2+i12)+bc(j2+i21))                           6d17s22
            bc(j2+i12)=avg                                              6d17s22
            bc(j2+i21)=avg                                              6d17s22
           end do                                                       6d17s22
          end do                                                        6d17s22
         end do                                                         6d17s22
        end if                                                          6d17s22
        if(isblkderd(3,is).eq.isblkderd(4,is).and.                      7d1s22
     $       iacto(isblkderd(1,is)).gt.0)then                           7d1s22
         nrow=iacto(isblkderd(1,is))*iacto(isblkderd(2,is))             7d1s22
         do i4=0,iacto(isblkderd(4,is))-1                               7d1s22
          do i3=0,i4-1                                                  6d17s22
           i34=jdenpt(is)+nrow*(i3+iacto(isblkderd(3,is))*i4)           7d1s22
           i43=jdenpt(is)+nrow*(i4+iacto(isblkderd(3,is))*i3)           7d1s22
           do irow=0,nrow-1                                             6d17s22
            o1=bc(i34+irow)
            o2=bc(i43+irow)
            avg=0.5d0*(bc(i34+irow)+bc(i43+irow))                       6d17s22
            bc(i34+irow)=avg                                            6d17s22
            bc(i43+irow)=avg                                            6d17s22
           end do                                                       6d17s22
          end do                                                        6d17s22
         end do                                                         6d17s22
        end if                                                          6d17s22
       end do                                                           6d17s22
      end if                                                            6d16s22
      ibcoff=ibcoffo
      icall=icall+1                                                     7d14s06
c
c     if dynamic weights, we now need to normalize ...                  4d6s18
c
      if(ldynw)then                                                     3d17s21
       if(wsum.ne.0d0)then                                              4d4s23
        wsumi=1d0/wsum                                                   1d18s20
       else                                                             4d4s23
        wsumi=0d0                                                       4d4s23
       end if                                                           4d4s23
       dwsumi=-wsumi*dwsum                                              7d20s22
       eavg=eavg*wsumi                                                  1d18s20
       if((icall.eq.1.or.ioconv.ne.0).and.lprint)                       5d1s18
     $     write(6,*)('new dynamic weights ')
       do i=1,nstate                                                    5d1s18
        do j=1,isinfo(4,i)                                              5d1s18
         ewght(j,i)=ewght(j,i)*wsumi                                    5d11s21
        end do                                                          5d1s18
        if((icall.eq.1.or.ioconv.ne.0).and.lprint)                      5d1s18
     $       write(6,1351)(ewght(j,i),j=1,isinfo(4,i))                  5d11s21
 1351   format(20es14.5)                                                5d1s18
       end do                                                           5d1s18
       if(idynw.ne.2)then                                               4d27s21
        nbodc=1                                                         6d2s22
        if(mode.eq.3)nbodc=4                                            6d2s22
        if(mode.eq.4)nbodc=0                                            7d7s22
        do im=1,nbodc                                                   6d2s22
         do isb=1,nsymb                                                   4d6s18
          if(im.eq.4)then                                               7d20s22
           jsb=isb                                                      7d20s22
          else                                                          7d20s22
           jsb=multh(isb,ipuse)                                         7d20s22
          end if                                                        7d20s22
          do i=0,iacto(isb)*iacto(jsb)-1                                6d16s22
           bc(iden1(isb,im)+i)=bc(iden1(isb,im)+i)*wsumi                6d2s22
          end do                                                        6d2s22
         end do                                                          4d6s18
        end do                                                           4d6s18
        if(ipuse.eq.1.and.mode.ne.4)then                                7d7s22
         do idws=1,nsdlk                                                  4d6s18
          n1=isblk(1,idws)                                                4d6s18
          n2=isblk(2,idws)                                                4d6s18
          n3=isblk(3,idws)                                                4d6s18
          n4=isblk(4,idws)                                                4d6s18
          if(n1.eq.n2)then                                                4d6s18
           nnn=(iacto(n1)*(iacto(n1)+1))/2                                4d6s18
          else                                                            4d6s18
           nnn=iacto(n1)*iacto(n2)                                        4d6s18
          end if                                                          4d6s18
          if(n3.eq.n4)then                                                4d6s18
           mmm=(iacto(n3)*(iacto(n3)+1))/2                                4d6s18
          else                                                            4d6s18
           mmm=iacto(n3)*iacto(n4)                                        4d6s18
          end if                                                          4d6s18
          if(nnn*mmm.gt.0)then                                            4d6s18
           do k=0,nnn*mmm-1                                               4d6s18
            bc(jdenpt(idws)+k)=bc(jdenpt(idws)+k)*wsumi                   4d6s18
           end do                                                         4d6s18
          end if                                                          4d6s18
         end do                                                           4d6s18
        else                                                            6d16s22
         do idws=1,nsblkderd                                            4d4s23
          n1=isblkderd(1,idws)                                           6d16s22
          n2=isblkderd(2,idws)                                           6d16s22
          n3=isblkderd(3,idws)                                           6d16s22
          n4=isblkderd(4,idws)                                           6d16s22
          nnn=iacto(n1)*iacto(n2)                                       6d16s22
          mmm=iacto(n3)*iacto(n4)                                       6d16s22
          if(min(nnn,mmm).gt.0)then                                     6d16s22
           do k=0,nnn*mmm-1                                             6d16s22
            bc(jdenpt(idws)+k)=bc(jdenpt(idws)+k)*wsumi                   4d6s18
           end do                                                         4d6s18
          end if                                                        6d16s22
         end do                                                         6d16s22
        end if                                                          6d16s22
       end if                                                           4d27s21
      end if                                                            4d6s18
      if(nsumcas.ne.0.and.ioconv.ne.0.and.mynowprog.eq.0)then           5d3s21
       write(42,42)(extradata(i),i=1,nextradata),                       5d4s21
     $     (bc(i),i=iroodat,jroodat-1)                                  4d17s23
   42  format(100f16.9)                                                 5d4s21
      end if                                                            5d3s21
      if(mode.eq.2.or.mode.eq.3)tolv=tolbest                            4d27s23
      return
      end
