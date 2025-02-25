c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine intcsfcas(mdon,mdoo,ibasis,ncsf,nfcn,                  8d2s22
     $     ih0a,i2e,shift,pthresin,nec,icsfpd,multhd,mysym,             12d28s19
     $     nrootw,cprint,lprint,ieigold,ivsave,maxpsci,ncsft,           12d28s19
     $     nvcul,nlzzu,islz,ixlzz,ism,irel,irefo,idoubo,norb,ioconv,    12d28s19
     $     namesym,nlab,ismult,lambdacas,bcasvec,icasvec,ixw1,ixw2,     4d12s21
     $     iptrbit,nroolab,spinroo,jroolab,jroodat,shiftr,nsymb,npdiag, 11d9s22
     $     bc,ibc)                                                      11d9s22
      implicit real*8 (a-h,o-z)
      external second                                                   8d2s19
      logical lprint,lwrite,bcasvec                                     1d17s20
      integer*8 ipack(2)                                                7d5s19
      integer*4 ipack4(4)                                               7d5s19
      character*50 numbers                                              4d13s21
      equivalence (ipack,ipack4)                                        7d5s19
      character*2 spinm,spinroo                                         5d3s21
      dimension ncsf(*),ih0a(*),i2e(*),icsfpd(*),ibasis(*),             8d2s22
     $     multhd(8,8),timers(5),nvcul(*),ixlzz(8,*),islz(*),           12d31s19
     $     ixlzzu(8,4),ism(*),irel(*),irefo(*),idoubo(*),namesym(*)     12d28s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      data icall/0/                                                     12d28s19
      save icall                                                        12d28s19
      include "common.basis"                                               8d2s19
      include "common.input"                                            8d2s19
      include "common.store"                                            6d10s19
      include "common.print"                                            1d23s19
      icall=icall+1                                                     12d28s19
      pthresp=pthresin                                                  7d15s19
      ibcoffo=ibcoff                                                    7d12s19
      lwrite=lprint.and.cprint.lt.0.99d0.and.ioconv.eq.1                1d23s19
     $     .or.iprtr(12).ne.0                                           1d23s19
c
c     do ci using csfs for active space ...
c
      isbasis=ibcoff                                                    7d5s19
      call r2sb(ibasis,ibc(isbasis),nfcn,myfcn)                         1d25s21
      need=3*myfcn                                                      7d5s19
      need8=need/2                                                      7d5s19
      if(need8*2.lt.need)need8=need8+1                                  7d5s19
      ibcoff=isbasis+need8                                              7d5s19
      ihdig=ibcoff                                                      6d10s19
      ibcoff=ihdig+myfcn                                                7d5s19
      call enough('intcsfcas.  1',bc,ibc)
      call hiicsfcas(bc(ihdig),myfcn,mdon,mdoo,ibc(isbasis),            8d2s22
     $     ih0a,i2e,shift,nec,ism,irel,irefo,iptrbit,mysym,norb,bc,ibc) 11d10s22
      ihdigps=ibcoff                                                    8d24s22
      ibcoff=ihdigps+myfcn                                              8d24s22
      call enough('intcsfcas.  2',bc,ibc)
      ilzzdig=ibcoff                                                    12d24s19
      if(nlzzu.ne.0)then                                                12d24s19
       ibcoff=ilzzdig+myfcn                                             12d24s19
       if(nlzzu.eq.6)ibcoff=ibcoff+myfcn                                5d14s21
       call enough('intcsfcas.  3',bc,ibc)
       npass=nlzzu/2                                                    12d31s19
       call lzziicsfcas(bc(ilzzdig),myfcn,mdon,mdoo,                    8d2s22
     $     ibc(isbasis),ixlzz,shift,nec,islz,multhd,ism,irel,irefo,     8d2s22
     $      idoubo,npass,iptrbit,mysym,norb,bc,ibc)                     11d10s22
      end if                                                            12d24s19
      hlow=bc(ihdig)                                                    6d11s19
      bc(ihdigps)=hlow                                                  8d24s22
      do i=1,myfcn-1                                                    7d5s19
       hlow=min(hlow,bc(ihdig+i))
       bc(ihdigps+i)=bc(ihdig+i)                                        8d24s22
      end do
      if(lwrite)write(6,*)('hlow = '),hlow,pthresp
      hlow0=hlow                                                        7d15s19
      loopit=0                                                          7d15s19
 1676 continue                                                          7d15s19
      hlow=hlow0+pthresp                                                7d15s19
      loopit=loopit+1                                                   7d15s19
      jbasis=0
      nps=0
      npsf=0
      infcn=0                                                           7d1s19
      ipack(1)=ibc(isbasis)                                             7d5s19
      jsbasis=isbasis+1                                                 7d5s19
      ipack(2)=ibc(jsbasis)                                             7d5s19
      j0=1                                                              7d5s19
      ncsft=0                                                           7d11s19
      do i=1,myfcn                                                      7d5s19
       im=i-1
       edel=bc(ihdigps+im)-hlow                                         8d24s22
       nclo=ipack4(j0)                                                  7d5s19
       j1=j0+1                                                          7d5s19
       j2=j1+1                                                          7d5s19
       nclop=nclo+1
       nopen=nec-2*nclo
       iarg=nclop-mdon
       ncsft=ncsft+ncsf(iarg)                                           7d11s19
       if(edel.le.0d0)then
        npsf=npsf+1
        nps=nps+ncsf(iarg)
       end if
       if(edel.le.0d0)then                                              6d17s19
       infcn=infcn+1                                                    7d1s19
       end if                                                           6d17s19
       if(j0.eq.1)then                                                  7d5s19
        ipack(1)=ibc(jsbasis)                                           7d5s19
        jsbasis=jsbasis+1                                               7d5s19
        ipack(2)=ibc(jsbasis)                                           7d5s19
        jsbasis=jsbasis+1                                               7d5s19
        j0=2                                                            7d5s19
       else                                                             7d5s19
        ipack(1)=ibc(jsbasis)                                           7d5s19
        jsbasis=jsbasis+1                                               7d5s19
        ipack(2)=ibc(jsbasis)                                           7d5s19
        j0=1                                                            7d5s19
       end if                                                           7d5s19
      end do
      if(lwrite)write(6,*)('no. of ps fcns '),npsf,('csfs'),nps
      if(nps.gt.maxpsci)then                                            7d15s19
       if(lwrite)write(6,*)('nps exceeds maxpsci ... '),maxpsci         8d4s22
       ihdcopy=ibcoff                                                   2d9s20
       isort=ihdcopy+myfcn                                              2d9s20
       ibcoff=isort+myfcn                                               2d9s20
       call enough('intcsfcas.  4',bc,ibc)
       do i=0,myfcn-1                                                   2d9s20
        bc(ihdcopy+i)=bc(ihdig+i)                                       2d9s20
       end do                                                           2d9s20
       call dsortdws(bc(ihdcopy),ibc(isort),myfcn)                      1d18s23
       nnew=0                                                           2d9s20
       do i=0,myfcn-1                                                   2d9s20
        ifcn=ibc(isort+i)                                               2d9s20
        ioff=3*(ifcn-1)                                                 2d9s20
        if(mod(ioff,2).eq.0)then                                        2d9s20
         ioff=ioff/2                                                    2d9s20
         j0=1                                                           2d9s20
        else                                                            2d9s20
         ioff=(ioff-1)/2                                                2d9s20
         j0=2                                                           2d9s20
        end if                                                          2d9s20
        ipack(1)=ibc(isbasis+ioff)                                      2d9s20
        ipack(2)=ibc(isbasis+ioff+1)                                    2d9s20
        nclo=ipack4(j0)                                                 2d9s20
        nclop=nclo+1
        iarg=nclop-mdon
        nnew=nnew+ncsf(iarg)                                            2d9s20
        if(nnew.gt.maxpsci)then                                         2d9s20
         pthresp=bc(ihdcopy+i-2)-hlow0                                  2d9s20
         ijump=i
         go to 1166                                                     2d9s20
        end if                                                          2d9s20
       end do                                                           2d9s20
 1166  continue                                                         2d9s20
       ibcoff=ihdcopy                                                   2d9s20
       if(loopit.gt.4)then                                              7d15s19
        write(6,*)('but this would have been our 5th try ...')          7d15s19
        call dws_sync                                                   7d15s19
        call dws_finalize                                               7d15s19
        stop 'intcsf'                                                   7d15s19
       end if                                                           7d15s19
       if(lwrite)write(6,*)('reduce p space threshold to  '),pthresp              7d15s19
       go to 1676                                                       7d15s19
      end if                                                            7d15s19
      if(lwrite)then                                                    2d28s21
       call addcomma4(ncsft,numbers,0)                                  2d28s21
       write(6,*)('total number of csfs: '),numbers                     2d28s21
      end if                                                            2d28s21
      igcode=0                                                          7d19s22
 1414 continue                                                          7d19s22
      need=3*npsf*4                                                     6d12s19
      need8=need/8                                                      6d11s19
      if(8*need8.lt.need)need8=need8+1                                  6d11s19
      ipsbase=ibcoff
      ihdps=ipsbase+need8                                               6d11s19
      ibcoff=ihdps+npsf                                                 6d11s19
      call enough('intcsfcas.  5',bc,ibc)
      ipointp=ibcoff                                                    4d12s21
      ibcoff=ipointp+npsf                                               4d12s21
      ipointf=ibcoff                                                    7d11s19
      ibcoff=ipointf+nps                                                7d11s19
      ihivs=1                                                           12d28s19
      ilzzpsdig=ibcoff                                                  12d24s19
      ndigs=0                                                           5d14s21
      if(nlzzu.ne.0)then                                                12d24s19
       ibcoff=ilzzpsdig+npsf                                            12d24s19
       ndigs=1                                                          5d14s21
       call enough('intcsfcas.  6',bc,ibc)
       call gcsfpslz(myfcn,bc(ihdig),bc(ihdps),hlow,ibc(isbasis),          7d5s19
     $     ibc(ipsbase),ncsf,ibc(ipointf),mdon,0,bc(ihivs),hiv,ncsft,   8d2s22
     $     nlzzu,bc(ilzzdig),bc(ilzzpsdig),ibc(ipointp),ndigs,nps,npsf,     7d19s22
     $      iptrbit,mysym,nec,mdoo,nsymb,ixlzz,islz,multh,igcode,lprint,8d9s22
     $      norb,irefo,ism,irel,bc(ihdigps),bc,ibc)                     11d14s22
       if(igcode.eq.0)then                                              7d19s22
        igcode=1                                                        7d19s22
        ibcoff=ipsbase                                                  7d19s22
        go to 1414                                                      7d19s22
       end if                                                           7d19s22
      else                                                              8d2s22
       call gcsfps(myfcn,bc(ihdig),bc(ihdps),hlow,ibc(isbasis),          7d5s19
     $     ibc(ipsbase),ncsf,ibc(ipointf),mdon,0,bc(ihivs),hiv,ncsft,   12d28s19
     $     nlzzu,bc(ilzzdig),bc(ilzzpsdig),ibc(ipointp),ndigs,npsf)     5d14s21
      end if                                                            8d2s22
      npsx=nps                                                          12d28s19
      ntri=(npsx*(npsx+1))/2                                            12d28s19
      nhere=ntri/mynprocg                                               7d11s19
      nleft=ntri-nhere*mynprocg                                         7d12s19
      if(mynowprog.lt.nleft)nhere=nhere+1                               7d12s19
      ilzps=ibcoff                                                      12d24s19
      npsk=nps                                                          12d24s19
      if(nlzzu.ne.0)then                                                12d24s19
       ibcoff=ilzps+nps*nps                                             12d24s19
       ntri=nps*nps                                                     12d24s19
       call enough('intcsfcas.  7',bc,ibc)
       npass=nlzzu/2                                                    12d31s19
       call pslzzcsfcas(npsf,bc(ilzzpsdig),ibc(ipsbase),nec,bc(ilzps),
     $      nps,ncsf,mdon,ixlzz,multh,islz,npsk,                        8d1s22
     $      ism,irel,irefo,idoubo,lambdacas,npass,lwrite,mdoo,iptrbit,  4d12s21
     $      ixw1,ixw2,mysym,norb,bc,ibc)                                11d10s22
       ibcoff=ilzps+nps*npsk                                            12d24s19
      end if                                                            12d24s19
      ihps=ibcoff                                                       7d11s19
      ibcoff=ihps+ntri                                                  7d11s19
      call pshamcsfcas(bc,ibc,npsf,bc(ihdps),ibc(ipsbase),nec,bc(ihps),
     $     nps,ncsf,mdon,i2e,ih0a,multh,                                    8d1s22
     $     0,bc(ihivs),hvv,nlzzu,npsk,bc(ilzps),ism,irel,irefo,mdoo,    4d12s21
     $     mysym,iptrbit,ixw1,ixw2,norb,bc,ibc)                         11d10s22
      if(nlzzu.ne.0)then                                                12d24s19
       npsxk=npsk                                                       12d28s19
       ntri=(npsxk*(npsxk+1))/2                                         12d24s19
       nhere=ntri/mynprocg                                              12d24s19
       nleft=ntri-nhere*mynprocg                                        12d24s19
       if(mynowprog.lt.nleft)nhere=nhere+1                              12d24s19
      end if                                                            12d24s19
      ibcoff=ihps+nhere                                                 7d11s19
      ibctop=ibcoff                                                     6d25s18
      iter=0                                                            7d11s19
      ntail=0                                                           6d27s18
      ihtail=1
      nrootz=nrootw                                                     7d13s18
      isto=0                                                            6d22s18
      nokv=nps                                                          12d28s19
      nokkv=npsk                                                        12d28s19
      ibasisv=1                                                         7d11s19
    1 continue                                                          7d11s19
       iter=iter+1                                                      7d11s19
       if(iprtr(12).ne.0)write(6,*)('starting iteration no. '),iter,
     $      icall
       if(iter.gt.100)then                                              7d11s19
        write(6,*)('exceeded 100 iterations in intcsf')
        call dws_sync
        call dws_finalize
        stop
       end if
       ihcpy=ibcoff                                                     6d27s18
       ibcoff=ihcpy+nhere+ntail                                         6d27s18
       call enough('intcsfcas.  8',bc,ibc)
       do i=0,nhere-1                                                   6d22s18
        bc(ihcpy+i)=bc(ihps+i)                                           6d22s18
       end do                                                            6d22s18
       if(iprtr(12).ne.0)then                                           8d2s22
        write(6,*)('hcpy ')
        call prntm2(bc(ihcpy),1,nhere,1)
        if(ntail.gt.0)then
         write(6,*)('htail ')
         call prntm2(bc(ihtail),1,ntail,1)
        end if                                                          8d2s22
       end if                                                           8d2s22
       if(ntail.gt.0)then                                               6d27s18
        jhcpy=ihcpy+nhere                                               6d27s18
        do i=0,ntail-1                                                  6d27s18
         bc(jhcpy+i)=bc(ihtail+i)                                       6d27s18
        end do                                                          6d27s18
       end if                                                           6d27s18
       if(nlzzu.eq.0)then                                               12d27s19
        nuse=nokv+isto                                                   11d13s18
       else                                                             12d27s19
        nuse=nokkv+isto                                                 12d27s19
       end if                                                           12d27s19
       ntri=(nuse*(nuse+1))/2                                           6d27s18
       nhx=ntri/mynprocg                                                6d27s18
       nleft=ntri-nhx*mynprocg                                          6d27s18
       if(mynowprog.lt.nleft)nhx=nhx+1                                  6d27s18
        if(nuse.le.npdiag)then                                          8d25s22
         call diagdy(bc(ihcpy),nuse,nrootz,ieig,ivec,bc,ibc)            11d14s22
         idum=1                                                         8d25s22
        else                                                            8d25s22
         call phouse(ihcpy,nuse,1d0,0d0,nrootz,ieig,ivec,1,idum,bc,ibc) 11d9s22
        end if                                                          8d25s22
       if(iprtr(12).ne.0)then                                           2d14s20
        write(6,*)('eigenvalues from phouse '),ivec,ieig,ibcoff
        call prntm2(bc(ieig),1,nrootz,1)
        call prntm2(bc(ivec),nuse,nrootz,nuse)
        if(bc(ieig).ne.bc(ieig))then
         call dws_synca                                                 8d2s22
         call dws_finalize                                              8d2s22
         stop                                                           8d2s22
        end if
       end if                                                           2d14s20
       ivecx=ibcoff                                                     7d11s19
       igx=ivecx+ncsft*nrootz                                           7d11s19
       ibcoff=igx+ncsft*nrootz                                          7d11s19
       call enough('intcsfcas.  9',bc,ibc)
       nx=isto                                                          7d11s19
       call epdvec(nrootz,nuse,bc(ivec),bc(ivecx),ncsft,ibc(ipointf),   7d11s19
     $      bc(ibasisv),nx,nps,0,nlzzu,npsk,bc(ilzps),bc,ibc)           11d14s22
       if(iprtr(12).ne.0)write(6,*)('after epdvec, eig: '),
     $      (bc(ieig+ix),ix=0,nrootz-1),ieig
  310  continue                                                         1d17s20
       ivecxsave=0                                                      4d21s21
       if(iter.eq.1.and..not.bcasvec.and.nuse.ne.ncsft)then             2d14s20
c                                                                       4d21s21
c     save copy of vectors just in case isto comes out 0 from update.   4d21s21
c     then p-space is actually everything and what is currently under   4d21s21
c     ivecx is actually the converged vector.                           4d21s21
c                                                                       4d21s21
        ivecxsave=ibcoff                                                4d21s21
        ibcoff=ivecxsave+ncsft*nrootz                                   4d21s21
        call enough('intcsfcas. 10',bc,ibc)
        do i=0,ncsft*nrootz-1                                           1d17s20
         bc(ivecxsave+i)=bc(ivecx+i)                                    4d21s21
         bc(ivecx+i)=0d0                                                1d17s20
        end do                                                          1d17s20
        ncsfr=ncsft*nrootw                                               1d17s20
        call ilimts(ncsfr,1,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s,m2e)    1d17s20
        jcasvec=icasvec                                                   1d17s20
        jvsave=ivecx-1                                                   1d17s20
        nherex=mh+1-ml
        do i=ml,mh                                                        1d17s20
         bc(jvsave+i)=bc(jcasvec)                                        1d17s20
         jcasvec=jcasvec+1                                                1d17s20
        end do                                                            1d17s20
        call dws_gsumf(bc(ivecx),ncsfr)                                  1d17s20
       end if                                                           1d17s20
       if(iprtr(12).ne.0)write(6,*)('after bcasvec, eig: '),
     $      (bc(ieig+ix),ix=0,nrootz-1),ieig
       ieigp=ibcoff                                                     6d22s18
       ibcoff=ieigp+nrootz                                              7d13s18
       call enough('intcsfcas. 11',bc,ibc)
       call hccsfcas(bc(ivecx),bc(igx),ncsft,ibc(isbasis),ncsf,         8d1s22
     $      myfcn,ih0a,i2e,bc(ihdig),nrootz,mdon,nec,                   8d1s22
     $      bc(ieigp),ism,irel,irefo,mysym,iptrbit,ixw1,ixw2,mdoo,norb, 11d10s22
     $      bc,ibc)                                                     11d10s22
       if(iprtr(12).ne.0)write(6,*)('after hccsfcas, eig: '),
     $      (bc(ieig+ix),ix=0,nrootz-1),ieig
       if(iter.eq.1.and..not.bcasvec.and.nuse.ne.ncsft)then             2d14s20
        do i=0,nrootz-1                                                 1d17s20
         dot=0d0                                                        1d17s20
         sz=0d0
         iadv=ivecx+ncsft*i                                             1d17s20
         iadg=igx+ncsft*i                                               1d17s20
         do j=0,ncsft-1                                                 1d17s20
          dot=dot+bc(iadv+j)*bc(iadg+j)                                 1d17s20
          sz=sz+bc(iadv+j)**2
         end do                                                         1d17s20
         bc(ieig+i)=dot                                                 1d17s20
        end do                                                          1d17s20
       end if                                                           1d17s20
       telap=time2-time1-tovr                                           8d2s19
       timers(3)=timers(3)+telap                                        8d2s19
       if(iter.eq.1)then                                                6d27s18
        ieigold=ibcoff                                                  6d27s18
        ibcoff=ieigold+nrootz                                           7d13s18
        nrooto=nrootz                                                   7d13s18
        call enough('intcsfcas. 12',bc,ibc)
       else                                                             6d27s18
        if(nrootz.gt.nrooto)then                                        7d13s18
         iegn=ibcoff                                                    6d27s18
         ibcoff=iegn+nrootz                                             7d13s18
         call enough('intcsfcas. 13',bc,ibc)
         do i=0,nrooto-1                                                6d27s18
          bc(iegn+i)=bc(ieigold+i)                                      6d27s18
         end do                                                         6d27s18
         do i=nrooto,nrootz-1                                           7d13s18
          bc(iegn+i)=0d0                                                6d27s18
         end do                                                         6d27s18
         nrooto=nrootz                                                  7d13s18
         ieigold=iegn                                                   6d27s18
        end if                                                          6d27s18
       end if                                                           6d27s18
       diffx=0d0                                                        6d27s18
       do i=0,nrootz-1                                                  7d13s18
        bc(ieigp+i)=bc(ieig+i)                                          6d26s18
        if(iprtr(12).ne.0)write(6,*)('setting ieigp to '),bc(ieig+i),
     $       ieig+i
        if(iter.eq.1)then                                               6d25s18
         diffx=1d0                                                      6d27s18
        else                                                            6d25s18
         diffx=max(diffx,abs(bc(ieigp+i)-bc(ieigold+i)))                6d27s18
        end if                                                          6d25s18
        bc(ieigold+i)=bc(ieigp+i)                                       6d25s18
        if(iprtr(12).ne.0)write(6,*)('setting ieigold to '),bc(ieigp+i),
     $       ieigp+i,ieigold+i
       end do                                                           6d25s18
c
c     at this point, under ivecx are the best vectors so far.
c
       ivsave=ibcoff                                                    6d27s18
       ibcoff=ivsave+nrootz*ncsft                                       12d28s19
       call enough('intcsfcas. 14',bc,ibc)
       do i=0,nrootz*ncsft-1                                            7d12s19
        bc(ivsave+i)=bc(ivecx+i)                                        7d12s19
       end do                                                           6d27s18
       call updateiicsf(bc(ivecx),bc(ieigp),nrootz,ncsft,bc(ihdig),     7d12s19
     $      myfcn,ibc(isbasis),ncsf,mdon,nps,ibc(ipointf),isto,         7d12s19
     $      bc(ivecx),bc,ibc)                                           11d10s22
       if(ivecxsave.ne.0.and.isto.eq.0)then                             4d21s21
c
c     oops: p-space is everything but we overwrote vectors ...          4d21s21
c     so reload correct vectors.
c
        do i=0,ncsft*nrootz-1                                           4d21s21
         bc(ivecx+i)=bc(ivecxsave+i)                                    4d21s21
         bc(ivsave+i)=bc(ivecx+i)                                       4d21s21
        end do                                                          4d21s21
        call hccsfcas(bc(ivecx),bc(igx),ncsft,ibc(isbasis),ncsf,        8d1s22
     $      myfcn,ih0a,i2e,bc(ihdig),nrootz,mdon,nec,                   8d1s22
     $      bc(ieigp),ism,irel,irefo,mysym,iptrbit,ixw1,ixw2,mdoo,norb, 11d10s22
     $       bc,ibc)                                                    11d10s22
        do i=0,nrootz-1                                                 1d17s20
         dot=0d0                                                        1d17s20
         sz=0d0
         iadv=ivecx+ncsft*i                                             1d17s20
         iadg=igx+ncsft*i                                               1d17s20
         do j=0,ncsft-1                                                 1d17s20
          dot=dot+bc(iadv+j)*bc(iadg+j)                                 1d17s20
          sz=sz+bc(iadv+j)**2
         end do                                                         1d17s20
         bc(ieig+i)=dot                                                 1d17s20
         bc(ieigp+i)=dot                                                4d21s21
         bc(ieigold+i)=dot                                              4d21s21
        end do                                                          1d17s20
       end if                                                           4d21s21
       tcov=smallest*1d4                                                8d2s19
       if(isto.eq.0.or.((diffx.lt.tcov.or.isto.eq.nrootz)               8d2s19
     $      .and.iter.gt.1))then                                        6d27s18
        if(nrootz.eq.0)then                                             7d13s18
         if(lprint)                                                     2d19s19
     $       write(6,*)('no roots found satisfying energy criterion')   2d19s19
        else                                                            6d26s18
         if(lwrite)then                                                 9d5s19
          write(6,*)('calculations converged after iteration '),iter
          write(spinm,1066)ismult                                       1d5s20
 1066     format(i2)                                                    1d5s20
          write(6,*)('spin multiplicity: '),ismult,                     12d28s19
     $ ('symmetry '),mysym,('('),spinm,(char(namesym(j)),j=1,nlab),(')')1d5s20
          iroosv=0                                                      5d3s21
          do j=6,min(nroolab,nlab)                                      5d3s21
           if(namesym(j).ne.ibc(jroolab+j))iroosv=1                     5d3s21
          end do                                                        5d3s21
          if(nlab.ne.nroolab.or.iroosv.ne.0.or.spinm.ne.spinroo)then    5d3s21
           spinroo=spinm                                                5d3s21
           iroosv=1                                                     5d3s21
           nroolab=nlab                                                 5d3s21
           do j=6,nlab                                                  5d3s21
            ibc(jroolab+j)=namesym(j)                                   5d3s21
           end do                                                       5d3s21
           do i=0,nrootz-1
            bc(jroodat)=bc(ieigold+i)                                   5d3s21
            jroodat=jroodat+1                                           5d3s21
           end do                                                       5d3s21
          end if                                                        5d3s21
          do i=0,nrootz-1                                                7d13s18
           ip=i+1                                                       8d10s22
           write(6,*)ip,bc(ieigold+i)+shiftr                            8d10s22
          end do
         end if                                                         2d19s19
        end if                                                          6d26s18
        go to 2
       end if
       ibasisv=ivecx                                                    7d12s19
       ieigp=ibcoff                                                     6d22s18
       igx=ieigp+isto                                                   7d12s19
       ibcoff=igx+ncsft*isto                                            7d12s19
       call enough('intcsfcas. 15',bc,ibc)
       call hccsfcas(bc(ivecx),bc(igx),ncsft,ibc(isbasis),ncsf,         8d1s22
     $      myfcn,ih0a,i2e,bc(ihdig),isto,mdon,nec,                     8d1s22
     $      bc(ieigp),ism,irel,irefo,mysym,iptrbit,ixw1,ixw2,mdoo,norb, 11d10s22
     $      bc,ibc)                                                     11d10s22
       call stuffintohiicsf(bc(igx),isto,ncsft,ihcpy,bc(ivecx),nps,          6d27s18
     $      ibc(ipointf),icall,iter,ntail,0,hiv,nlzzu,npsk,bc(ilzps),   4d27s23
     $      bc,ibc)                                                     4d27s23
c
c     close up memory                                                   6d27s18
c
       jsto=ibctop                                                      6d27s18
       do i=0,nrootz-1                                                  7d13s18
        bc(jsto+i)=bc(ieigold+i)                                        6d27s18
       end do                                                           6d27s18
       ieigold=jsto                                                     6d27s18
       jsto=jsto+nrootz                                                 7d13s18
       do i=0,isto*ncsft                                                6d27s18
        bc(jsto+i)=bc(ibasisv+i)                                        6d27s18
       end do                                                           6d27s18
       ibasisv=jsto                                                     6d27s18
       jsto=jsto+isto*ncsft                                             6d27s18
       do i=0,ntail-1                                                   6d27s18
        bc(jsto+i)=bc(ihcpy+i)                                          6d27s18
       end do                                                           6d27s18
       ihtail=jsto                                                      6d27s18
       ibcoff=jsto+ntail                                                6d27s18
      go to 1                                                           6d22s18
    2  continue                                                         7d12s19
c
c     close up memory                                                   6d27s18
c
      if(lwrite)then                                                    12d28s19
       call prtcsfveccas(bc(ivsave),ncsft,nrootw,myfcn,ibc(isbasis),    8d2s22
     $     ncsf,mdon,cprint,xnorm,nec,6,norb,irefo,iptrbit,mysym,mdoo,  11d10s22
     $     bc,ibc)                                                      11d10s22
      end if                                                            7d12s19
      nrootw=nrootz                                                     7d13s18
      jsto=ibcoffo                                                      6d27s18
      do i=0,nrootz-1                                                   7d13s18
       bc(jsto+i)=bc(ieigold+i)                                         6d27s18
      end do                                                            6d27s18
      ieigold=jsto                                                      6d27s18
      jsto=jsto+nrootz                                                  7d13s18
      sz=0d0                                                            1d17s20
      do i=0,nrootz*ncsft-1                                             12d28s19
       bc(jsto+i)=bc(ivsave+i)                                          6d27s18
       sz=sz+bc(jsto+i)**2
      end do                                                            6d27s18
      ivsave=jsto                                                       6d27s18
      ibcoff=ivsave+ncsft*nrootw                                        12d28s19
      ii=ivsave
      ncsfr=ncsft*nrootw                                                1d17s20
      call ilimts(ncsfr,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,j2e)     1d24s20
      jcasvec=icasvec                                                   1d17s20
      jvsave=ivsave-1                                                   1d17s20
      do i=il,ih                                                        1d17s20
       bc(jcasvec)=bc(jvsave+i)                                         1d17s20
       jcasvec=jcasvec+1                                                1d17s20
      end do                                                            1d17s20
      icasvec=jcasvec                                                   4d20s21
      if(nrootw.gt.0)then                                               2d16s19
      xnorm=0d0
      end if                                                            2d16s19
      return
      end
