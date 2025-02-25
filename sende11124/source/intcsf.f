c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine intcsf(idorb,isorb,mdon,mdoo,ibasis,iptr,ncsf,nfcn,    6d10s19
     $     ih0a,i2e,shift,pthresin,nec,icsfpd,multhd,mysym,nvv,hiv,hvv,  7d24s19
     $     nrootw,ethresx,cprint,lprint,ieigold,ivsave,maxpsciq,ncsft,  8d2s19
     $     nvcul,nlzzu,islz,ixlzz,igoon,bintvec,intvec,nfall,ipall,     4d2s20
     $     iball,ih0ae,iptrbit,ixw1,ixw2,nball,ifcnpx,nsymb,iptribit,   10d15s20
     $     idoi,isoi,nfcnc,ibasisc,iptrcb,vidotvi,vihvi,nameciu,lz2,bc, 11d9s22
     $     ibc,isaved,ostng,idroot)                                     7d13s23
      implicit real*8 (a-h,o-z)
      external second                                                   8d2s19
      logical lprint,lwrite,bintvec,ltest,ltest2,lconv                  6d20s24
      character*2 spinm                                                 1d9s18
      character*50 numbers                                              2d28s21
      character*(*) ostng(idroot)                                       7d13s23
      integer*1 idorb(*),isorb(*)
      integer*8 ipack(2)                                                7d5s19
      integer*4 ipack4(4)                                               7d5s19
      equivalence (ipack,ipack4)                                        7d5s19
      dimension ncsf(*),iptr(4,*),ih0a(*),i2e(*),icsfpd(*),ibasis(*),   6d18s19
     $     multhd(8,8),hvv(*),hiv(*),timers(5),nvcul(*),ixlzz(8,5),     5d14s21
     $     ixlzzu(8,5),islz(*),iball(*),ipall(*),ih0ae(*),              5d14s21
     $     iptrbit(2,mdoo+1,*),vidotvi(*),vihvi(*),                     3d31s21
     $     nball(*),ifcnpx(*),iptribit(*),nameciu(*)                    4d25s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      common/hcdscm/timds(12)
      include "common.basis"                                               8d2s19
      include "common.input"                                            8d2s19
      include "common.store"                                            6d10s19
      include "common.mrci"                                             7d17s19
      include "common.print"                                            1d13s20
      data icall/0/                                                     4d8s20
      save icall                                                        4d8s20
      icall=icall+1                                                     4d8s20
      pthresp=pthresin                                                  7d15s19
      ibcoffo=ibcoff                                                    7d12s19
      diffwf=0d0                                                        7d15s24
      shift0=0d0                                                        2d9s21
      isaved=0                                                          2d3s23
      lwrite=(lprint.and.cprint.lt.0.99d0).or.iprtr(12).ne.0            1d13s20
c
c     do ci using csfs for active space ...
c
      if(lwrite)write(6,*)('in intcsf '),nfcn,mysym,nec,nrootw,         2d3s23
     $     bintvec                                                      2d3s23
      if(iprtr(12).ne.0.and.nvv.gt.0)then                               1d13s20
       write(6,*)('input hiv: ')                                        1d13s20
       call prntm2(hiv,nrootw,nvv,nrootw)
       write(6,*)('input hvv: ')
       call prntm2(hvv,nrootw,nrootw,nrootw)
      end if                                                            1d13s20
      isbasis=ibcoff                                                    7d5s19
      call r2sb(ibasis,ibc(isbasis),nfcn,myfcn)                         1d25s21
      need=3*myfcn                                                      7d5s19
      need8=need/2                                                      7d5s19
      if(need8*2.lt.need)need8=need8+1                                  7d5s19
      ibcoff=isbasis+need8                                              7d5s19
      nhist=10                                                          4d22s22
      ihist=ibcoff                                                      4d22s22
      ibcoff=ihist+nhist                                                4d22s22
      ihdig=ibcoff                                                      6d10s19
      ibcoff=ihdig+myfcn                                                7d5s19
      call enough('intcsf.  1',bc,ibc)
      call hiicsf(bc(ihdig),myfcn,idorb,isorb,mdon,mdoo,ibc(isbasis),   7d5s19
     $     iptr,ih0a,i2e,shift0,nec,bc,ibc)                             11d10s22
      ihdigps=ibcoff                                                    8d24s22
      ibcoff=ihdigps+myfcn                                              8d24s22
      call enough('intcsf.  2',bc,ibc)
      ilzzdig=ibcoff                                                    12d24s19
      if(nlzzu.ne.0)then                                                12d24s19
       ibcoff=ilzzdig+myfcn                                             12d24s19
       if(nlzzu.eq.6)ibcoff=ibcoff+myfcn                                5d14s21
       call enough('intcsf.  3',bc,ibc)
       npass=nlzzu/2                                                    1d2s20
       call lzziicsf(bc(ilzzdig),myfcn,idorb,isorb,mdon,mdoo,           12d24s19
     $      ibc(isbasis),iptr,ixlzz,shift0,nec,islz,multhd,npass,bc,ibc)11d10s22
      end if                                                            12d24s19
      hlow=bc(ihdig)                                                    6d11s19
      bc(ihdigps)=hlow                                                  8d24s22
      do i=1,myfcn-1                                                    7d5s19
       hlow=min(hlow,bc(ihdig+i))
       bc(ihdigps+i)=bc(ihdig+i)                                        8d24s22
      end do
      hlow0=hlow                                                        7d15s19
      loopit=0                                                          7d15s19
 1676 continue                                                          7d15s19
      do i=0,nhist-1                                                    4d22s22
       ibc(ihist+i)=0                                                   4d22s22
      end do                                                            4d22s22
      frach=dfloat(nhist)/pthresp                                       4d22s22
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
       edel=bc(ihdigps+im)-hlow                                         8d25s22
       nclo=ipack4(j0)                                                  7d5s19
       j1=j0+1                                                          7d5s19
       j2=j1+1                                                          7d5s19
       nclop=nclo+1
       nopen=nec-2*nclo
       iarg=nclop-mdon
       ncsft=ncsft+ncsf(iarg)                                           7d11s19
       if(edel.le.0d0)then
        ddd=(bc(ihdigps+im)-hlow0)*frach                                8d25s22
        iddd=ddd                                                        4d22s22
        ibc(ihist+iddd)=ibc(ihist+iddd)+ncsf(iarg)                      8d4s22
        npsf=npsf+1
        nps=nps+ncsf(iarg)
       end if
       jc=iptr(2,nclop)+nclo*(ipack4(j1)-1)-1                           7d5s19
       jo=iptr(4,nclop)+nopen*(ipack4(j2)-1)-1                          7d5s19
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
       if(lwrite)write(6,*)('nps exceeds maxpsci ... '),maxpsci                            7d15s19
       ntoth=0                                                           4d22s22
       do i=0,nhist-1                                                   4d22s22
        ntoth=ntoth+ibc(ihist+i)                                          4d22s22
        if(ntoth.gt.maxpsci)then                                        4d22s22
         facth=dfloat(max(i,1))/dfloat(nhist)                           4d22s22
         pthresp=pthresp*facth                                          4d22s22
         go to 12321                                                    4d22s22
        end if                                                          4d22s22
       end do                                                           4d22s22
       pthresp=pthresp*dfloat(nhist-1)/dfloat(nhist)                    8d4s22
12321  continue                                                         4d22s22
       if(loopit.gt.4)then                                              7d15s19
        write(6,*)('but this would have been our 5th try ...')          7d15s19
        call dws_sync                                                   7d15s19
        call dws_finalize                                               7d15s19
        stop 'intcsf'                                                   7d15s19
       end if                                                           7d15s19
       if(lwrite)write(6,*)('reduce p space threshold to  '),pthresp              7d15s19
       go to 1676                                                       7d15s19
      end if                                                            7d15s19
      ratio=dfloat(nps)/dfloat(ncsft)                                   5d25s21
      if(ratio.gt.0.9d0.and.ncsft.lt.maxpsci.and.nps.ne.ncsft)then      5d25s21
       if(lwrite)
     $ write(6,*)('p-space would essentially include everything so'),
     $     (' force it to be everything.')                              5d25s21
       pthresp=1d10                                                     5d25s21
       if(loopit.gt.4)then                                              5d25s21
        call dws_synca                                                  5d25s21
        call dws_finalize                                               5d25s21
        stop 'loopit'                                                            5d25s21
       end if                                                           5d25s21
       go to 1676                                                       5d25s21
      end if                                                            5d25s21
      if(lwrite)write(6,*)('total number of csfs: '),ncsft
      igcode=0                                                          7d19s22
 1414 continue                                                          7d19s22
      need=3*npsf*4                                                     6d12s19
      need8=need/8                                                      6d11s19
      if(8*need8.lt.need)need8=need8+1                                  6d11s19
      ipsbase=ibcoff
      ihdps=ipsbase+need8                                               6d11s19
      ibcoff=ihdps+npsf                                                 6d11s19
      call enough('intcsf.  4',bc,ibc)
      ipointp=ibcoff                                                    3d24s21
      ibcoff=ipointp+npsf                                               3d25s21
      ipointf=ibcoff                                                    7d11s19
      ibcoff=ipointf+nps                                                7d11s19
      if(nvv.gt.0)then                                                  9d4s19
       ihivs=ibcoff                                                     9d4s19
       ibcoff=ihivs+ncsft*nvv                                           9d4s19
       call enough('intcsf.  5',bc,ibc)
      else                                                              9d4s19
       ihivs=1                                                          9d4s19
      end if                                                            9d4s19
      ilzzpsdig=ibcoff                                                  12d24s19
      ndigs=0                                                           5d14s21
      if(nlzzu.ne.0)then                                                12d24s19
       ibcoff=ilzzpsdig+npsf                                            12d24s19
       ndigs=1                                                          5d14s21
       if(nlzzu.eq.6)then                                               5d14s21
        ibcoff=ibcoff+npsf                                              5d14s21
        ndigs=2                                                         5d14s21
       end if                                                           5d14s21
       call enough('intcsf.  6',bc,ibc)
       call gcsfpslz(myfcn,bc(ihdig),bc(ihdps),hlow,ibc(isbasis),          7d5s19
     $     ibc(ipsbase),ncsf,ibc(ipointf),mdon,nvv,bc(ihivs),hiv,ncsft, 12d24s19
     $     nlzzu,bc(ilzzdig),bc(ilzzpsdig),ibc(ipointp),ndigs,nps,npsf,     7d19s22
     $      iptrbit,mysym,nec,mdoo,nsymb,ixlzz,islz,multh,igcode,lprint,8d9s22
     $      norb,irefo,ism,irel,bc(ihdigps),bc,ibc)                     11d14s22
       if(igcode.eq.0)then                                              7d19s22
        igcode=1                                                        7d19s22
        ibcoff=ipsbase                                                  7d19s22
        go to 1414                                                      7d19s22
       end if                                                           7d19s22
      else                                                              7d19s22
       call gcsfps(myfcn,bc(ihdig),bc(ihdps),hlow,ibc(isbasis),          7d5s19
     $     ibc(ipsbase),ncsf,ibc(ipointf),mdon,nvv,bc(ihivs),hiv,ncsft, 12d24s19
     $     nlzzu,bc(ilzzdig),bc(ilzzpsdig),ibc(ipointp),ndigs,npsf)     5d14s21
      end if                                                            7d19s22
      mdoop=mdoo+1                                                      3d23s21
      if(iuse2den.ne.2)then                                             9d5s24
       ipsbasen=ibcoff                                                   3d24s21
       iptrps=ipsbasen+myfcn*2                                           3d24s21
       ibcoff=iptrps+mdoop                                               3d24s21
       call enough('intcsf.  7',bc,ibc)
       call int1tobit2(ibc(ipsbase),ibc(ipsbasen),iptr,ibc(iptrps),      3d24s21
     $     idorb,isorb,npsf,nec,mdoop,bc,ibc)                           11d10s22
       ipointq=ibcoff                                                   3d24s21
       ipointq2=ipointq+myfcn                                           3d25s21
       ipsbaseq=ipointq2+myfcn                                          3d25s21
       iptrqs=ipsbaseq+myfcn*2                                           3d24s21
       ibcoff=iptrqs+mdoop                                              3d24s21
       call enough('intcsf.  8',bc,ibc)
       call int1tobit3(ibc(isbasis),ibc(ipsbaseq),iptr,ibc(ipointp),    3d24s21
     $      npsf,ibc(ipointq),ibc(iptrqs),idorb,isorb,myfcn,nec,mdoop,  3d24s21
     $      nfcnq,ibc(ipointq2),bc,ibc)                                 11d10s22
       ikeep=ibcoff                                                     3d24s21
       ikeepa=ikeep+nfcnc                                               3d24s21
       ibcoff=ikeepa+nfcnc                                              3d25s21
       call enough('intcsf.  9',bc,ibc)
       call strip2ps(ibc(iptrps),ibc(ipsbasen),npsf,                    3d24s21
     $     nfcnc,ibasisc,iptrcb,ibc(ikeep),nkeep,norb,                  3d24s21
     $      ibc(iptrqs),ibc(ipsbaseq),nfcnq,ibc(ikeepa),bc,ibc)         11d10s22
      end if                                                            3d24s21
      npsx=nps+nvv                                                      9d4s19
      if(itestmrci.eq.1)then                                            8d1s19
       ntri=0                                                           8d1s19
      else                                                              8d1s19
       ntri=(npsx*(npsx+1))/2                                           9d4s19
       nhere=ntri/mynprocg                                               7d11s19
       nleft=ntri-nhere*mynprocg                                         7d12s19
       if(mynowprog.lt.nleft)nhere=nhere+1                               7d12s19
      end if                                                            8d1s19
      ilzps=ibcoff                                                      12d24s19
      npsk=nps                                                          12d24s19
      if(nlzzu.ne.0)then                                                12d24s19
       ibcoff=ilzps+nps*nps                                             12d24s19
       ntri=max(ntri,nps*nps)                                           10d14s20
       call enough('intcsf. 10',bc,ibc)
       call pslzzcsf(npsf,bc(ilzzpsdig),ibc(ipsbase),nec,bc(ilzps),nps, 12d24s19
     $      ncsf,mdon,icsfpd,iptr,ixlzz,idorb,isorb,multh,islz,npsk,    1d2s20
     $      npass,lwrite,mdoo,iptrbit,ixw1,ixw2,mysym,lz2,bc,ibc)       11d10s22
       ibcoff=ilzps+nps*npsk                                            12d24s19
      end if                                                            12d24s19
      ihps=ibcoff                                                       7d11s19
      ibcoff=ihps+ntri                                                  7d11s19
      call pshamcsf(npsf,bc(ihdps),ibc(ipsbase),nec,bc(ihps),nps,       7d11s19
     $     ncsf,mdon,icsfpd,iptr,i2e,idorb,isorb,ih0a,multh,            9d4s19
     $     nvv,bc(ihivs),hvv,nlzzu,npsk,bc(ilzps),nfall,ipall,          9d20s20
     $     ibc(iball(1)),                                               9d20s20
     $     mdoo,mysym,iptrbit,ih0ae,ixw1,ixw2,nball,ifcnpx,nsymb,shift0,11d10s22
     $     bc,ibc)                                                      11d10s22
      if(nlzzu.ne.0)then                                                12d24s19
       npsxk=npsk+nvv                                                   12d24s19
       ntri=(npsxk*(npsxk+1))/2                                         12d24s19
       nhere=ntri/mynprocg                                              12d24s19
       nleft=ntri-nhere*mynprocg                                        12d24s19
       if(mynowprog.lt.nleft)nhere=nhere+1                              12d24s19
      end if                                                            12d24s19
      ibcoff=ihps+nhere                                                 7d11s19
      iter=0                                                            7d11s19
      ntail=0                                                           6d27s18
      ihtail=1
      nrootz=nrootw                                                     7d13s18
      if(wconvci.ne.0)then                                              6d20s24
       iveclast=ibcoff                                                  6d20s24
       ibcoff=iveclast+npsf*nrootw                                      6d20s24
       call enough('intcsf.veclast',bc,ibc)                             6d20s24
      end if                                                            6d20s24
      isto=0                                                            6d22s18
      nokv=nps+nvv                                                      7d11s19
      nokkv=npsk+nvv                                                    12d27s19
      ibasisv=1                                                         7d11s19
      ivecps=ibcoff                                                     3d25s21
      ivecqs=ivecps+nps*nrootw                                          3d25s21
      ncsfq=ncsft-nps                                                   3d25s21
      ibcoff=ivecqs+ncsfq*nrootw                                        3d25s21
      mxxds=maxdsi                                                      10d14s22
      isavqv=ibcoff                                                     3d25s21
      isavqg=isavqv+ncsfq*nrootw*2*mxxds                                      3d25s21
      isavpg=isavqg+ncsfq*nrootw*2*mxxds                                      3d25s21
      ibcoff=isavpg+nps*nrootw*2*mxxds                                        3d25s21
      ihivq=ibcoff                                                      3d25s21
      ibcoff=ihivq+nvv*ncsfq                                            3d25s21
      ngot=0                                                            3d25s21
      call enough('intcsf. 11',bc,ibc)
      do i=ivecqs,ibcoff-1                                              3d25s21
       bc(i)=0d0                                                        3d25s21
      end do                                                            3d25s21
      if(nvv.gt.0)then                                                  3d25s21
       call gotoqs(hiv,bc(ihivq),ncsft,ncsfq,nvv,ibc(ipointq),myfcn,    3d25s21
     $     ibc(isbasis),ncsf,mdon)                                      3d25s21
      end if                                                            3d25s21
      ibctop=ibcoff                                                     6d25s18
      ltest=.false.                                                     3d22s22
      ltest2=.true.                                                     4d6s22
      iterd=0                                                           7d18s22
    1 continue                                                          7d11s19
       iter=iter+1                                                      7d11s19
       iterd=iterd+1                                                    7d18s22
       if(iprtr(12).ne.0)write(6,*)('starting iteration no. '),iter
       if(iter.gt.maxmiti)then                                              7d11s19
        write(6,*)('exceeded '),maxmiti,(' iterations in intcsf')       10d14s22
        call dws_sync
        call dws_finalize
        stop 'intcsf'                                                   10d14s22
       end if
       ihcpy=ibcoff                                                     6d27s18
       ibcoff=ihcpy+nhere+ntail                                         6d27s18
       call enough('intcsf. 12',bc,ibc)
       do i=0,nhere-1                                                   6d22s18
        bc(ihcpy+i)=bc(ihps+i)                                           6d22s18
       end do                                                            6d22s18
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
       if(ethresx.eq.0d0)then                                           6d26s18
        if(nuse.lt.npdiagci)then                                        8d26s22
         call diagdy(bc(ihcpy),nuse,nrootz,ieig,ivec,bc,ibc)            11d14s22
         idum=1                                                         8d26s22
        else                                                            8d26s22
         call phouse(ihcpy,nuse,1d0,0d0,nrootz,ieig,ivec,1,idum,bc,ibc) 11d9s22
        end if                                                          8d26s22
       else                                                             6d26s18
        call phouse(ihcpy,nuse,elow,ethresx,nrootz,ieig,ivec,1,idum,bc, 11d9s22
     $       ibc)                                                       11d9s22
       end if                                                           6d26s18
       if(igoon.eq.0.and.mynowprog.eq.0.and.iter.gt.1.and.              7d13s23
     $      (mod(iter-2,mxxds).ne.0.or.iter.eq.2))then                  7d14s23
        if(nrootz.le.idroot)then                                        7d13s23
         do ir=1,nrootz                                                 7d13s23
          irm=ir-1                                                      7d13s23
          write(ostng(ir),1551)bc(ieig+irm)+shift                       7d13s23
 1551     format(f19.12)                                                7d13s23
          change=bc(ieig+irm)-bc(ieigold+irm)                            7d13s23
          if(change.ne.0d0)then                                            2d26s19
           ndig=nint(-log10(abs(change)))+1                              2d26s19
           ist=min(19,7+ndig)                                            1d12s23
           if(ostng(ir)(ist:ist).ne.' ')then                              1d6s20
            if(ichar(ostng(ir)(ist:ist)).lt.ichar('0').or.                1d21s20
     $        ichar(ostng(ir)(ist:ist)).gt.ichar('9'))then               1d21s20
             if(ist.lt.20)then                                           10d20s20
              write(6,*)('bad last digit: '),ostng(ir)(ist:ist)
              write(6,*)('ndig,ist '),ndig,ist
             end if                                                      10d20s20
             left=0                                                      1d21s20
            else                                                         1d21s20
             read(ostng(ir)(ist:ist),*)left                                 11d17s19
            end if                                                       1d21s20
            if(left.gt.5)then                                             11d17s19
             iadd=1                                                       11d17s19
            else if(left.lt.5)then                                        11d17s19
             iadd=0                                                       11d17s19
            else                                                          11d17s19
             if(ist.lt.19)then                                            11d17s19
              istp=ist+1                                                  11d17s19
              read(ostng(ir)(istp:istp),*)left2                            11d17s19
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
             write(ostng(ir),1551)bc(ieig+irm)+shift-fact                7d13s23
            end if                                                        11d17s19
           end if                                                        1d6s20
           do k=ist,19                                                   2d26s19
            ostng(ir)(k:k)=' '                                            2d26s19
           end do                                                        2d26s19
          end if                                                        7d13s23
         end do                                                          7d13s23
         write(6,1216)iter,(ostng(ir),ir=1,nrootz)                      7d13s23
 1216    format(i5,(a19))                                               7d13s23
        else                                                            7d13s23
         write(6,1215)iter,(bc(ieig+ir)+shift,ir=0,nrootz-1)             6d23s23
 1215    format(i5,20f18.8)                                              6d23s23
        end if                                                          7d13s23
       else if(igoon.eq.0.and.mynowprog.eq.0.and.iter.gt.1)then         7d14s23
        write(6,*)('contract down to best Davidson vector')             7d14s23
       end if                                                           6d23s23
       if(iprtr(12).ne.0)then
        write(6,*)('eigenvalues from phouse '),ivec,loc(bc(ivec))
        call prntm2(bc(ieig),1,nrootz,1)
        call prntm2(bc(ivec),nuse,nrootz,nuse)
       end if
       if(wconvci.ne.0d0)then                                           6d20s24
        if(iter.gt.1)then                                               6d20s24
         diffwf=0d0                                                     6d20s24
         do i=0,nrootz-1                                                6d20s24
          jvec=ivec+nuse*i                                              6d20s24
          jveclast=iveclast+npsf*i                                      6d20s24
          do j=0,npsf-1                                                  6d20s24
           diffwf=diffwf+(bc(jvec+j)-bc(jveclast+j))**2                 6d20s24
          end do                                                        6d20s24
         end do                                                         6d20s24
         diffwf=sqrt(diffwf/dfloat(npsf*nrootz))                        6d20s24
        end if                                                          6d20s24
        do i=0,nrootz-1                                                 6d20s24
         jvec=ivec+nuse*i                                               6d20s24
         jveclast=iveclast+npsf*i                                       6d20s24
         do j=0,npsf-1                                                  6d20s24
          bc(jveclast+j)=bc(jvec+j)                                     6d20s24
         end do                                                         6d20s24
        end do                                                          6d20s24
       end if                                                           6d20s24
c
c     for Davidson correction, we need Vi*Vi and Vi*Hii*Vi.             3d30s21
c     the structure of ivec is npsk p-space contributions,
c     nvv external contributions, and nroot or 2*nroot q-space
c     contributions. We only want p and q-space contributions           6d20s24
c
       if(nvv.gt.0)then                                                 3d31s21
        do ir=1,nrootw                                                  3d31s21
         vidotvi(ir)=0d0                                                3d31s21
         lvec=ivec+nuse*(ir-1)                                          3d31s21
         do i=0,npsk-1                                                    3d31s21
          vidotvi(ir)=vidotvi(ir)+bc(lvec+i)**2                                   3d31s21
         end do                                                           3d31s21
         jvec=lvec+npsk+nvv                                               3d31s21
         if(iter.eq.2)then                                                3d31s21
          do i=0,isto-1                                                 4d30s21
           vidotvi(ir)=vidotvi(ir)+bc(jvec+i)**2                                  3d31s21
          end do                                                          3d31s21
         else if(iter.gt.2)then                                         3d31s21
          do i=0,isto-1                                                 4d14s21
           vidotvi(ir)=vidotvi(ir)+bc(jvec+i)**2                        3d31s21
          end do                                                          3d31s21
         end if                                                           3d31s21
        end do                                                          3d31s21
        if(iprtr(12).ne.0)then
         write(6,*)('vidotvi ')
         call prntm2(vidotvi,nrootw,1,nrootw)
        end if
        imxt=ibcoff                                                     3d31s21
        ibcoff=imxt+nuse*nrootw                                         3d31s21
        call enough('intcsf. 13',bc,ibc)
        do i=imxt,ibcoff-1                                              3d31s21
         bc(i)=0d0                                                      3d31s21
        end do                                                          3d31s21
        jhcpy=ihps                                                      9d8s22
        jhcpytop=ihps+nhere                                             5d2s23
        ihigh=npsk-1                                                     3d31s21
        ilow=npsk+nvv                                                   3d31s21
        do i=0,nuse-1                                                   3d31s21
         do j=0,i                                                       3d31s21
          ii=((i*(i+1))/2)+j                                            3d31s21
          if(mod(ii,mynprocg).eq.mynowprog)then                         3d31s21
           if((i.le.ihigh.or.i.ge.ilow).and.(j.le.ihigh.or.j.ge.ilow))  3d31s21
     $         then                                                     3d31s21
            do ir=0,nrootw-1                                             3d31s21
             iadx=imxt+i+nuse*ir                                         3d31s21
             iadv=ivec+j+nuse*ir                                         3d31s21
             bc(iadx)=bc(iadx)+bc(iadv)*bc(jhcpy)                       3d31s21
            end do                                                       3d31s21
            if(i.ne.j)then                                              3d31s21
             do ir=0,nrootw-1                                             3d31s21
              iadx=imxt+j+nuse*ir                                         3d31s21
              iadv=ivec+i+nuse*ir                                         3d31s21
              bc(iadx)=bc(iadx)+bc(iadv)*bc(jhcpy)                       3d31s21
             end do                                                       3d31s21
            end if                                                      3d31s21
           end if                                                       3d31s21
           jhcpy=jhcpy+1                                                3d31s21
           if(jhcpy.eq.jhcpytop)then                                    5d2s23
            jhcpytop=0                                                  5d2s23
            jhcpy=ihtail                                                5d2s23
           end if                                                       5d2s23
          end if                                                        3d31s21
         end do                                                         3d31s21
        end do                                                          3d31s21
        if(iprtr(12).ne.0)then
         write(6,*)('my part of imxt')
         call prntm2(bc(imxt),nuse,nrootw,nuse)
         do i=0,nuse*nrootw-1
          bc(ibcoff+i)=bc(imxt+i)
         end do
         call dws_gsumf(bc(ibcoff),nuse*nrootw)
         write(6,*)('global summed')
         call prntm2(bc(ibcoff),nuse,nrootw,nuse)
        end if
        jvec=ivec                                                       3d31s21
        jmxt=imxt                                                       3d31s21
        do ir=1,nrootw                                                  3d31s21
         vihvi(ir)=0d0                                                  3d31s21
         do i=0,nuse-1                                                  3d31s21
          vihvi(ir)=vihvi(ir)+bc(jvec+i)*bc(jmxt+i)                     3d31s21
         end do                                                         3d31s21
         jvec=jvec+nuse                                                 3d31s21
         jmxt=jmxt+nuse                                                 3d31s21
        end do                                                          3d31s21
        if(iprtr(12).ne.0)then                                          3d31s21
         write(6,*)('my part of vihvi ')
         call prntm2(vihvi,nrootw,1,nrootw)
         do i=1,nrootw
          bc(ibcoff+i)=vihvi(i)
         end do
         call dws_gsumf(bc(ibcoff+1),nrootw)
         write(6,*)('global summed ')
         call prntm2(bc(ibcoff+1),nrootw,1,nrootw)
        end if                                                          3d31s21
        ibcoff=imxt                                                     3d31s21
       end if                                                           3d31s21
       ivecx=ibcoff                                                     7d11s19
       igx=ivecx+ncsft*nrootz                                           7d11s19
       ibcoff=igx+ncsft*nrootz                                          7d11s19
       igq=ibcoff                                                       3d25s21
       ibcoff=igq+ncsfq*nrootz                                          3d25s21
       call enough('intcsf. 14',bc,ibc)
       nx=isto                                                          7d11s19
       if(iuse2den.ne.2)then                                            9d5s24
        call epdvec(nrootz,nuse,bc(ivec),bc(ivecx),ncsft,ibc(ipointf),   7d11s19
     $      bc(ibasisv),0,nps,nvv,nlzzu,npsk,bc(ilzps),bc,ibc)          11d14s22
       else                                                             3d25s21
        call epdvec(nrootz,nuse,bc(ivec),bc(ivecx),ncsft,ibc(ipointf),   7d11s19
     $      bc(ibasisv),nx,nps,nvv,nlzzu,npsk,bc(ilzps),bc,ibc)         11d14s22
       end if                                                           3d25s21
       if(iter.eq.1.and..not.bintvec.and.ncsft.gt.nps.and.              2d14s20
     $      novguess.eq.0.and.iuse2den.eq.2)then                        9d5s24
        write(6,*)('loading vectors '),intvec                                  1d17s20
        write(6,*)('what I have for external part of vectors: ')
        call prntm2(bc(intvec),nrootz*2,nrootz,nrootz*2)                1d17s20
        do i=0,ncsft*nrootz-1                                           1d17s20
         bc(ivecx+i)=0d0                                                1d17s20
        end do                                                          1d17s20
        call ilimts(nrootz,ncsft,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s,   1d17s20
     $      m2e)                                                        1d17s20
        i10=m1s                                                          1d17s20
        i1n=nrootz                                                       1d17s20
        jntvec=intvec+nrootz*2*nrootz                                   1d17s20
        do i2=m2s,m2e                                                    1d17s20
         if(i2.eq.m2e)i1n=m1e                                            1d17s20
         do i1=i10,i1n                                                   1d17s20
          iad=ivecx+i2-1+ncsft*(i1-1)                                   1d17s20
          bc(iad)=bc(jntvec)                                            1d17s20
          jntvec=jntvec+1                                                1d17s20
         end do                                                          1d17s20
         i10=1                                                           1d17s20
        end do                                                           1d17s20
        call dws_gsumf(bc(ivecx),ncsft*nrootz)                          1d17s20
        call enough('intcsf. 15',bc,ibc)
        call checkps(bc(ivecx),ncsft,nrootz,nps,ibc(ipointf),iok,bc,ibc)11d14s22
        if(iok.ne.nrootz)then                                           2d25s21
         novguess=1                                                     2d25s21
         go to 1491                                                     2d25s21
        end if                                                          2d25s21
        if(iprtr(12).ne.0)then                                          10d2s20
         write(6,*)('reconstituted vecx: ')                              1d17s20
         call prntm2(bc(ivecx),ncsft,nrootz,ncsft)                       1d17s20
        end if                                                          10d2s20
        if(nvv.gt.0)then                                                1d17s20
         if(nlzzu.ne.0)then                                             1d17s20
          ivecp=ivec+npsk                                               1d17s20
         else                                                           1d17s20
          ivecp=ivec+nps                                                1d17s20
         end if                                                         1d17s20
         do i=0,nrootz-1                                                1d17s20
          do j=0,nvv-1                                                  1d17s20
           iadv=ivecp+j+nuse*i                                          1d17s20
           iado=intvec+j+nrootw*2*i                                     1d17s20
           bc(iadv)=bc(iado)                                            1d17s20
          end do                                                        1d17s20
         end do                                                         1d17s20
        end if                                                          1d17s20
       end if                                                           1d17s20
 1491  continue                                                         2d25s21
       ieigp=ibcoff                                                     6d22s18
       ibcoff=ieigp+nrootz                                              7d13s18
       call enough('intcsf. 16',bc,ibc)
c
c     here we need q space part of h*vfull.
c     gq=hqp*Vp+hqq*Vq. The second part we make from what we have.
c
       if(iuse2den.ne.2)then                                            10d8s20
        call mov2ps(bc(ivecx),bc(ivecps),ncsft,ibc(ipointf),nps,nrootw) 3d25s21
        if(iuse2den.eq.0)then                                           9d5s24
         call hccsfrilr(bc(igq),ncsfq,bc(ivecps),nps,ibc(ipsbaseq),      3d24s21
     $       ibc(ipsbasen),ncsf,nfcnq,npsf,ih0ae,i2e,nrootw,mdon,nec,   3d25s21
     $       multh,nkeep,ibc(ikeep),ibasisc,mdoo,ibc(iptrqs),           3d24s21
     $       ibc(iptrps),ixw1,ixw2,iptrcb,bc,ibc)                       11d10s22
        else                                                            9d5s24
         call hccsflr(bc(igq),ncsfq,bc(ivecps),nps,ibc(ipsbaseq),        9d4s24
     $       ibc(ipsbasen),ncsf,nec,nfcnq,npsf,ih0a,i2e,nrootw,         9d4s24
     $       multh,mdon,mdoo,ibc(iptrqs),ibc(iptrps),ixw1,ixw2,bc,ibc)  9d4s24
        end if                                                          9d5s24
         if(iter.gt.1)then                                              3d25s21
          jvec=ivec+nuse-ngot                                           3d22s22
          if(iter.eq.2)then                                             3d25s21
           ngotu=0                                                      3d25s21
          else                                                          3d25s21
           ngotu=ngot                                                   3d25s21
          end if                                                        3d25s21
          call dgemm('n','n',ncsfq,nrootw,ngot,1d0,bc(isavqg),ncsfq,    3d22s22
     $         bc(jvec),nuse,1d0,bc(igq),ncsfq,                         3d23s22
     d' intcsf.  1')
          if(iter.gt.2)ngot=ngotu                                       5d5s21
          ltest=.true.                                                  3d22s22
          call dgemm('n','n',ncsfq,nrootw,ngot,1d0,bc(isavqv),ncsfq,    3d23s22
     $         bc(jvec),nuse,0d0,bc(ivecqs),ncsfq,                      3d23s22
     d' intcsf.  2')
          if(mod(iter-1,mxxds).eq.0)then                                3d22s22
           iterd=0                                                      7d18s22
           ivtmp=ibcoff                                                 4d8s22
           ibcoff=ivtmp+nrootw*ncsfq                                    4d8s22
           call enough('intcsf. 17',bc,ibc)
           isto=0                                                       4d8s22
           do i=0,nrootw-1                                              3d22s22
            iadi=ivecqs+ncsfq*i                                         3d22s22
            do j=0,isto-1                                               4d8s22
             iadj=ivtmp+ncsfq*j                                         4d8s22
             dot=0d0                                                    3d22s22
             do k=0,ncsfq-1                                             3d22s22
              dot=dot+bc(iadi+k)*bc(iadj+k)                             3d22s22
             end do                                                     3d22s22
             do k=0,ncsfq-1                                             3d22s22
              bc(iadi+k)=bc(iadi+k)-dot*bc(iadj+k)                      3d22s22
             end do                                                     3d22s22
            end do                                                      3d22s22
            dot=0d0                                                     3d22s22
            do k=0,ncsfq-1                                              3d22s22
             dot=dot+bc(iadi+k)**2                                      3d22s22
            end do                                                      3d22s22
            if(dot.gt.1d-12)then                                        4d8s22
             dot=1d0/sqrt(dot)                                           3d22s22
             jvtmp=ivtmp+ncsfq*isto                                     4d8s22
             do k=0,ncsfq-1                                              3d22s22
              bc(jvtmp+k)=bc(iadi+k)*dot                                4d8s22
             end do                                                      3d22s22
             isto=isto+1                                                4d8s22
            end if                                                      4d8s22
           end do                                                       3d22s22
           ngot=0                                                       3d22s22
           do i=0,isto*ncsfq-1                                          4d8s22
            bc(ivecqs+i)=bc(ivtmp+i)                                    4d8s22
           end do                                                       4d8s22
           ibcoff=ivtmp                                                 4d8s22
           ltest=.false.                                                3d22s22
          end if                                                        3d22s22
          ltest2=mod(iter-2,mxxds).eq.0                                      3d22s22
          if(.not.ltest)ltest2=.false.                                  3d23s22
         end if                                                         3d25s21
        if(iprtr(12).ne.0)then                                          3d25s21
         write(6,*)('gq after hccsfri ')                                3d25s21
         call prntm2(bc(igq),ncsfq,nrootz,ncsfq)                        3d25s21
        end if                                                          3d25s21
       else                                                             10d8s20
        call hccsf(bc(ivecx),bc(igx),ncsft,ibc(isbasis),ncsf,iptr,myfcn, 10d2s20
     $      ih0a,i2e,bc(ihdig),nrootz,mdon,isorb,idorb,icsfpd,nec,      7d12s19
     $       bc(ieigp),mysym,iptrbit,ixw1,ixw2,mdoo,bc,ibc)             11d10s22
        if(iprtr(12).ne.0)then
         write(6,*)('gx after hccsfri ')
         call prntm2(bc(igx),ncsft,nrootz,ncsft)
        end if                                                          3d25s21
       end if                                                           10d8s20
       if(iprtr(12).ne.0)then
        write(6,*)('vectors after hccsf')
        call prntm2(bc(ivecx),ncsft,nrootz,ncsft)
       end if
       call addqtop(bc(ivecx),ncsft,nrootz,bc(ivecqs),ncsfq,            3d23s22
     $      ibc(ipointq),myfcn,ibc(isbasis),mdon,ncsf,bc,ibc)           11d14s22
       if(iprtr(12).ne.0)then
        write(6,*)('vectors after addqtop')
        call prntm2(bc(ivecx),ncsft,nrootz,ncsft)
       end if
       if(iter.eq.1.and..not.bintvec.and.ncsft.gt.nps.and.novguess.eq.0.
     $      and.iuse2den.eq.2)then                                      9d5s24
        do i=0,nrootz-1                                                 1d17s20
         dot=0d0                                                        1d17s20
         sz=0d0                                                         1d17s20
         do j=0,ncsft-1                                                 1d17s20
          iadv=ivecx+j+ncsft*i                                          1d17s20
          iadg=igx+j+ncsft*i                                            1d17s20
          dot=dot+bc(iadv)*bc(iadg)                                     1d17s20
          sz=sz+bc(iadv)**2                                             1d17s20
         end do                                                         1d17s20
         bc(ieig+i)=dot                                                 1d17s20
        end do                                                          1d17s20
       end if                                                           1d17s20
       if(iter.gt.1.and.intvec.gt.0)then                                3d25s21
        if(iprtr(12).ne.0)then                                          3d25s21
         write(6,*)('saving isavqv')                                    3d25s21
         call prntm2(bc(isavqv),ncsfq,nrootw,ncsfq)                     3d25s21
        end if                                                          3d25s21
        call ilimts(nrootw,ncsfq,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s,  3d25s21
     $      m2e)                                                        1d17s20
        mhere=mh+1-ml                                                   1d17s20
        jntvec=intvec                                                   3d25s21
        isaved=1                                                        2d3s23
        i10=m1s                                                         1d17s20
        i1n=nrootw                                                      1d17s20
        do i2=m2s,m2e                                                   1d17s20
         if(i2.eq.m2e)i1n=m1e                                           1d17s20
         do i1=i10,i1n                                                  1d17s20
          iadg=isavqv+i2-1+ncsfq*(i1-1)                                 3d25s21
          bc(jntvec)=bc(iadg)                                           1d17s20
          jntvec=jntvec+1                                               1d17s20
         end do                                                         1d17s20
         i10=1                                                          1d17s20
        end do                                                          1d17s20
       end if                                                           3d25s21
       if(.not.bintvec)then                                             1d17s20
        if(iuse2den.ne.2)then                                           9d5s24
        else                                                            3d25s21
         call ilimts(nrootw,ncsft,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s,  1d17s20
     $      m2e)                                                        1d17s20
         mhere=mh+1-ml                                                   1d17s20
         jntvec=intvec+mhere+nrootw*nrootw*2                             1d17s20
         write(6,*)('stuff igx into intvec ')
         call prntm2(bc(igx),ncsft,nrootw,ncsft)
         isaved=1                                                       2d3s23
         i10=m1s                                                         1d17s20
         i1n=nrootw                                                      1d17s20
         do i2=m2s,m2e                                                   1d17s20
          if(i2.eq.m2e)i1n=m1e                                           1d17s20
          do i1=i10,i1n                                                  1d17s20
           iadg=igx+i2-1+ncsft*(i1-1)                                    1d17s20
           bc(jntvec)=bc(iadg)                                           1d17s20
           jntvec=jntvec+1                                               1d17s20
          end do                                                         1d17s20
          i10=1                                                          1d17s20
         end do                                                          1d17s20
        end if                                                          3d25s21
       end if                                                           1d17s20
       if(nvv.gt.0)then                                                 9d4s19
        if(nlzzu.ne.0)then                                              12d27s19
         ivecp=ivec+npsk                                                12d27s19
        else                                                            12d27s19
         ivecp=ivec+nps                                                  9d4s19
        end if                                                          12d27s19
        if(ncsfq.gt.0)then                                              5d7s21
         call dgemm('n','n',ncsfq,nrootw,nvv,1d0,bc(ihivq),ncsfq,        3d25s21
     $       bc(ivecp),nuse,1d0,bc(igq),ncsfq,                          3d25s21
     d' intcsf.  3')
        end if                                                          5d7s21
       end if                                                           9d4s19
       if(iter.eq.1)then                                                6d27s18
        ieigold=ibcoff                                                  6d27s18
        ibcoff=ieigold+nrootz                                           7d13s18
        nrooto=nrootz                                                   7d13s18
        call enough('intcsf. 18',bc,ibc)
       else                                                             6d27s18
        if(nrootz.gt.nrooto)then                                        7d13s18
         iegn=ibcoff                                                    6d27s18
         ibcoff=iegn+nrootz                                             7d13s18
         call enough('intcsf. 19',bc,ibc)
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
        if(iter.eq.1)then                                               6d25s18
         diffx=1d0                                                      6d27s18
        else                                                            6d25s18
         if(iprtr(12).ne.0)                                             3d22s22
     $        write(6,*)i,bc(ieigp+i),bc(ieigp+i)-bc(ieigold+i)         3d22s22
         diffx=max(diffx,abs(bc(ieigp+i)-bc(ieigold+i)))                6d27s18
        end if                                                          6d25s18
        bc(ieigold+i)=bc(ieigp+i)                                       6d25s18
       end do                                                           6d25s18
       if(wconvci.ne.0d0)then                                           6d20s24
        bc(ibcoff)=diffx                                                6d20s24
        bc(ibcoff+1)=diffwf                                             6d20s24
        call dws_bcast(bc(ibcoff),2)                                    6d20s24
        diffx=bc(ibcoff)                                                6d20s24
        diffwf=bc(ibcoff+1)                                             6d20s24
       else                                                             6d20s24
        call dws_bcast(diffx,1)                                          3d11s21
       end if                                                           6d20s24
c
c     at this point, under ivecx are the best vectors so far.
c
       ivsave=ibcoff                                                    6d27s18
       ibcoff=ivsave+nrootz*(ncsft+nvv)                                 9d4s19
       call enough('intcsf. 20',bc,ibc)
       do i=0,nrootz*ncsft-1                                            7d12s19
        bc(ivsave+i)=bc(ivecx+i)                                        7d12s19
       end do                                                           6d27s18
       if(nvv.gt.0)then                                                 9d4s19
        if(nlzzu.ne.0)then                                              12d27s19
         ivecp=ivec+npsk                                                12d27s19
        else                                                            12d27s19
         ivecp=ivec+nps                                                  9d4s19
        end if                                                          12d27s19
        ivsave2=ivsave+nrootz*ncsft                                     9d4s19
        do i=0,nrootz-1                                                 9d4s19
         do j=0,nvv-1                                                   9d4s19
          iadv=ivecp+j+nuse*i                                           9d4s19
          bc(ivsave2)=bc(iadv)                                          9d4s19
          ivsave2=ivsave2+1                                             9d4s19
         end do                                                         9d4s19
        end do                                                          9d4s19
       end if                                                           9d4s19
       if(iuse2den.eq.2)then                                            9d5s24
        call updateiicsf(bc(ivecx),bc(ieigp),nrootz,ncsft,bc(ihdig),     7d12s19
     $      myfcn,ibc(isbasis),ncsf,mdon,nps,ibc(ipointf),isto,         7d12s19
     $      bc(ivecx),bc,ibc)                                           11d10s22
       else                                                             3d25s21
        if(ncsfq.gt.0.and.(iter.eq.1.or.ltest))then                     3d22s22
        call updateq(bc(ivecqs),bc(igq),bc(ieigp),nrootz,ncsfq,          3d25s21
     $      bc(ihdig),ibc(ipointq2),nfcnq,ibc(ipsbaseq),ncsf,mdon,isto, 3d25s21
     $      bc(isavqv),ngot,iprtr(12))                                  3d25s21
        end if                                                          5d7s21
        if(iter.eq.1.and..not.bintvec.and.ncsfq.gt.0.and.novguess.eq.0) 3d25s21
     $       then                                                       3d25s21
         if(iprtr(12).ne.0)                                             3d25s21
     $       write(6,*)('overwrite qs guess vectors with saved ones ')  3d25s21
         do i=0,nrootz*ncsfq-1                                          3d25s21
          bc(ivecqs+i)=0d0                                              3d25s21
         end do                                                         3d25s21
         call ilimts(nrootw,ncsfq,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s, 3d25s21
     $      m2e)                                                        1d17s20
         mhere=mh+1-ml                                                   1d17s20
         jntvec=intvec                                                  3d25s21
         i10=m1s                                                         1d17s20
         i1n=nrootw                                                      1d17s20
         do i2=m2s,m2e                                                   1d17s20
          if(i2.eq.m2e)i1n=m1e                                           1d17s20
          do i1=i10,i1n                                                  1d17s20
           iadg=ivecqs+i2-1+ncsfq*(i1-1)                                3d25s21
           bc(iadg)=bc(jntvec)                                          3d25s21
           jntvec=jntvec+1                                               1d17s20
          end do                                                         1d17s20
          i10=1                                                          1d17s20
         end do                                                          1d17s20
         call dws_gsumf(bc(ivecqs),ncsfq*nrootw)                        3d25s21
         if(iprtr(12).ne.0)then                                         3d25s21
          call prntm2(bc(ivecqs),ncsfq,nrootw,ncsfq)                     3d25s21
          call dgemm('t','n',nrootw,nrootw,ncsfq,1d0,bc(ivecqs),
     $         ncsfq,bc(ivecqs),ncsfq,0d0,bc(ibcoff),nrootw)
          call prntm2(bc(ibcoff),nrootw,nrootw,nrootw)
         end if                                                         3d25s21
        end if                                                          3d25s21
       end if                                                           3d25s21
       tcov=econvcii                                                    2d9s21
       if(iprtr(12).ne.0)then                                           7d19s22
        write(6,*)('ltest2 '),ltest2,isto,diffx,ngot,nrootz,iter
        write(6,*)isto,diffx.lt.tcov,isto+ngot.eq.nrootz
        write(6,*)('iterd '),iterd,tcov
       end if                                                           7d19s22
       lconv=(diffx.lt.tcov.or.isto+ngot.eq.nrootz).and.iter.gt.1       6d20s24
       if(wconvci.ne.0d0)then                                           6d20s24
        if(diffwf.gt.wconvci)lconv=.false.                              6d20s24
       end if                                                           6d20s24
       if((isto.eq.0.or.lconv).and.(iterd.gt.1.or.iterd.eq.iter))then   6d20s24
        if(nrootz.eq.0)then                                             7d13s18
         if(lprint)                                                     2d19s19
     $       write(6,*)('no roots found satisfying energy criterion')   2d19s19
        else                                                            6d26s18
         if(lwrite)then                                                 9d5s19
          ntop=0                                                        1d9s18
          do i=1,6                                                      1d9s18
           if(nameciu(i).ne.0)ntop=i                                    4d25s21
          end do                                                        1d9s18
          write(spinm,37)ismult                                         1d9s18
   37     format(i2)                                                    1d9s18
          write(6,*)('calculations converged after iteration '),iter
          write(6,*)('state '),spinm,(char(nameciu(i)),i=1,ntop)        4d25s21
          if(igoon.eq.0)then                                            1d10s19
           do i=0,nrootz-1                                              1d10s19
            ip=i+1                                                      8d10s22
            write(6,*)spinm,(char(nameciu(j)),j=1,ntop),('>'),ip,       8d10s22
     $           bc(ieigold+i)+shift                                    2d9s21
           end do                                                       1d10s19
          else                                                          1d10s19
           do i=0,nrootz-1                                              1d10s19
            ip=i+1                                                      8d10s22
            write(6,*)ip,bc(ieigold+i)+shift                            8d10s22
           end do                                                       1d10s19
          end if                                                        1d10s19
         end if                                                         2d19s19
        end if                                                          6d26s18
        go to 2
       end if
       ibasisv=ivecx                                                    7d12s19
       ieigp=ibcoff                                                     6d22s18
       igx=ieigp+isto                                                   7d12s19
       igx2=igx+ncsft*isto                                              3d25s21
       ibcoff=igx2+ncsft*isto                                           3d25s21
       call enough('intcsf. 21',bc,ibc)
c
c     here vecx only has q-space non-zero elements, but we need
c     full space h*vecx
c
       if(iuse2den.ne.2)then                                            9d5s24
        if(iuse2den.ne.0)then                                           9d5s24
         call hccsflr(bc(igx2),ncsft,bc(ivecqs),ncsfq,ibc(isbasis),        9d4s24
     $       ibc(ipsbaseq),ncsf,nec,myfcn,nfcnq,ih0a,i2e,nrootw,         9d4s24
     $       multh,mdon,mdoo,iptrbit(1,1,mysym),ibc(iptrqs),ixw1,ixw2,  9d4s24
     $       bc,ibc)                                                    9d4s24
        else                                                            9d5s24
         call hccsfrilr(bc(igx2),ncsft,bc(ivecqs),ncsfq,ibc(isbasis),    3d25s21
     $       ibc(ipsbaseq),ncsf,myfcn,nfcnq,ih0ae,i2e,isto,mdon,nec,    3d25s21
     $       multh,nfcnc,ibc(ikeepa),ibasisc,mdoo,iptrbit(1,1,mysym),   3d31s21
     $       ibc(iptrqs),ixw1,ixw2,iptrcb,bc,ibc)                       11d10s22
        end if                                                          9d5s24
        igx3=ibcoff                                                     3d25s21
        ibcoff=igx3+ncsfq*isto                                          3d25s21
        call enough('intcsf. 22',bc,ibc)
        call gotoqs(bc(igx2),bc(igx3),ncsft,ncsfq,isto,ibc(ipointq),    3d25s21
     $       myfcn,ibc(isbasis),ncsf,mdon)                              3d25s21
        call stuffintohii2(bc(igx2),ncsft,bc(ivecqs),ihcpy,ncsfq,isto,  3d25s21
     $      bc(isavqv),bc(isavqg),bc(isavpg),ibc(ipointf),nps,          3d25s21
     $      ibc(ipointq),ntail,nvv,bc(ihivq),nlzzu,npsk,bc(ilzps),ngot, 3d25s21
     $      bc(igx3),nigot,iprtr(12),bc,ibc)                            11d10s22
        isto=ngot                                                       3d22s22
        isto0=isto                                                      3d25s21
       else                                                             10d8s20
        call hccsf
     $      (bc(ivecx),bc(igx),ncsft,ibc(isbasis),ncsf,iptr,myfcn,      5d6s20
     $      ih0a,i2e,bc(ihdig),isto,mdon,isorb,idorb,icsfpd,nec,        7d12s19
     $      bc(ieigp),mysym,iptrbit,ixw1,ixw2,mdoo,bc,ibc)              11d10s22
        call stuffintohiicsf(bc(igx),isto,ncsft,ihcpy,bc(ivecx),nps,          6d27s18
     $      ibc(ipointf),icall,iter,ntail,nvv,hiv,nlzzu,npsk,bc(ilzps), 11d10s22
     $       bc,ibc)                                                    11d10s22
       end if                                                           10d8s20
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
      if(lprint)then                                                    7d12s19
       call prtcsfvec(bc(ivsave),ncsft,nrootw,myfcn,ibc(isbasis),ncsf,   7d12s19
     $     mdon,cprint,xnorm,iptr,idorb,isorb,nec,6)                    7d17s19
       if(itestmrci.eq.2.or.itestmrci.eq.3)then                         7d17s19
      write(7,*)ncsft,nrootw
       call prtcsfvec(bc(ivsave),ncsft,nrootw,myfcn,ibc(isbasis),ncsf,   7d12s19
     $     mdon,0d0,xnorm,iptr,idorb,isorb,nec,7)                       7d17s19
       end if                                                           7d17s19
      end if                                                            7d12s19
      nrootw=nrootz                                                     7d13s18
      jsto=ibcoffo                                                      6d27s18
      do i=0,nrootz-1                                                   7d13s18
       bc(jsto+i)=bc(ieigold+i)                                         6d27s18
      end do                                                            6d27s18
      ieigold=jsto                                                      6d27s18
      jsto=jsto+nrootz                                                  7d13s18
      do i=0,nrootz*(ncsft+nvv)-1                                       9d4s19
       bc(jsto+i)=bc(ivsave+i)                                          6d27s18
      end do                                                            6d27s18
      ivsave=jsto                                                       6d27s18
      ibcoff=ivsave+(ncsft+nvv)*nrootw                                  11d13s18
      jvsave=ivsave
      if(intvec.ne.0.and.iuse2den.eq.2)then                             9d5s24
       call ilimts(nrootw,ncsft,mynprocg,mynowprog,ml,mh,m1s,m1e,m2s,   1d17s20
     $      m2e)                                                        1d17s20
       do i=0,nrootw*nrootw*2-1                                         1d17s20
        bc(intvec+i)=0d0                                                1d17s20
       end do                                                           1d17s20
       write(6,*)('saving intvec ')
       isaved=1                                                         2d3s23
       write(6,*)('ivsave: ')
       call prntm2(bc(ivsave),ncsft+nvv,nrootw,ncsft+nvv)
       call dgemm('t','n',nrootw,nrootw,ncsft+nvv,1d0,
     $      bc(ivsave),ncsft+nvv,bc(ivsave),ncsft+nvv,0d0,
     $      bc(ibcoff),nrootw)
       call prntm2(bc(ibcoff),nrootw,nrootw,nrootw)
       do i=0,nrootw-1                                                  1d17s20
        do j=0,nvv-1                                                    1d17s20
         iad1=intvec+j+nrootw*2*i                                       1d17s20
         iad2=ivsave+j+ncsft+(ncsft+nvv)*i                              1d17s20
         bc(iad1)=bc(iad2)                                              1d17s20
        end do                                                          1d17s20
       end do                                                           1d17s20
       jntvec=intvec+nrootw*nrootw*2                                    1d17s20
       i10=m1s                                                          1d17s20
       i1n=nrootw                                                       1d17s20
       do i2=m2s,m2e                                                    1d17s20
        if(i2.eq.m2e)i1n=m1e                                            1d17s20
        do i1=i10,i1n                                                   1d17s20
         iad=ivsave+i2-1+(ncsft+nvv)*(i1-1)                             1d17s20
         bc(jntvec)=bc(iad)                                             1d17s20
         jntvec=jntvec+1                                                1d17s20
        end do                                                          1d17s20
        i10=1                                                           1d17s20
       end do                                                           1d17s20
       jntvec=intvec+nrootw*nrootw*2                                    1d17s20
       mhere=mh+1-ml
      end if                                                            1d17s20
      ii=ivsave
      if(nrootw.gt.0)then                                               2d16s19
      xnorm=0d0
      if(lprint.and.nvv.ne.0.and.cprint.lt.0.99d0)                      2d26s19
     $     write(6,*)('internal part of first root: '),xnorm            2d19s19
      end if                                                            2d16s19
      if(itestmrci.eq.2)then                                            5d15s19
       if(lprint)write(6,*)('we have completed test 2 process')         5d15s19
       call dws_sync                                                    5d15s19
       call dws_finalize                                                5d15s19
       stop                                                             5d15s19
      end if                                                            5d15s19
      return
      end
