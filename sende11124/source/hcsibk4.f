c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcsibk4(nhand,nff1,iff1,i2sb,i2smb,mdoobp,ncsfs,       10d20s21
     $     nff0,iff0,vec,ncsfv,i2sk,i2smk,mdookp,nec,                   10d20s21
     $     mdon,nsymb,multh,irw1,irw2,nvirt,isymc,irorip,isopt,         10d20s21
     $     maxbx,ism,irel,                                              10d20s21
     $    irefo,norb,nroot,ih0n,nh0,ionexr,iifmx,ntype,isymmrci,l2e,    1d26s23
     $     nwiacc,bc,ibc)                                               1d26s23
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2),icode1(64),imap1(64),ipackc1(8)         10d20s21
      integer*8 nhand(mdoobp,nsymb),i1,in,i1c,i1o,i0c,i0o,itesta,        8d21s21
     $     itestb,ipackc,ipack,gandcc,gandco,gandcb                     2d7s23
      integer*2 ipack2(4)                                               10d22s21
      logical ldebug                                                    1d25s21
      equivalence (ipackc,ipackc1),(ipack,ipack2)                       10d22s21
      dimension nff1(mdoobp,nsymb,2),iff1(*),nff0(mdookp,3),            10d22s21
     $     iff0(*),vec(ncsfv,*),multh(8,8),nab4(2,3),idorb1(32),        10d20s21
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),ih0n(*),          10d20s21
     $     nvirt(*),iden1x(8,8),nden1x(8,8),nh0(*),ncsfs(*),            10d20s21
     $     ism(*),irel(*),irefo(*),ionexr(8,8,*),iifmx(*),              10d22s21
     $     mcsf(2),igya(4),imy(4),isopt(*),isy(4),iwpb(4),iwpk(4),      10d22s21
     $     iwpb1(4),iwpk1(4),iwpb2(4),iwpk2(4),ioxx(2)                  2d7s23
      include "common.store"                                            11d25s20
      data loopx/1000/
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/paddcm/npadddi                                             6d7s22
      loop=0
      irori=irorip-1                                                    10d20s21
      if(isopt(3).ne.0)then
       ieoro=1-irori                                                      9d6s21
      else                                                              10d14s21
       ieoro=irori                                                      9d6s21
      end if                                                            9d6s21
      ldebug=.false.
      if(ldebug)write(6,*)('hi, I am the new and improved hcsibk4!')       1d25s21
      mdoob=mdoobp-1                                                    10d20s21
      mdook=mdookp-1                                                    10d20s21
      ibcoffo=ibcoff                                                    12d8s20
c
c     maxbx is for input wavefcns ...
c     but here bra has ket roots, so maxbx may not be large enough.     6d7s22
c
      maxbxb=0                                                          6d7s22
      do isb=1,nsymb                                                    6d7s22
       isbv=multh(isb,isymmrci)                                         6d7s22
       do ii=mdon+1,mdoobp                                              6d30s23
        iarg=ii-mdon                                                    6d7s22
        maxbxb=max(maxbxb,ncsfs(iarg)*nff1(ii,isb,1)*nvirt(isbv)*nroot) 6d8s22
       end do                                                           6d7s22
      end do                                                            6d7s22
      maxbxb=maxbxb+npadddi                                             6d7s22
      ighere=ibcoff                                                     3d22s21
      iacc=ighere+maxbxb                                                6d7s22
      ibcoff=iacc+mynprocg                                              3d22s21
      call enough('hcsibk4.  1',bc,ibc)
      nfirst=0                                                          7d14s21
      nlast=0                                                           7d14s21
      nacc=0                                                            3d22s21
       do isb=1,nsymb                                                    11d26s20
        do ii=mdon+1,mdoobp                                              7d9s21
         if(nhand(ii,isb).ge.0)then                                     8d21s21
          call ddi_zero(bc,ibc,nhand(ii,isb))                           11d15s22
         end if                                                          11d26s20
        end do                                                           11d26s20
       end do                                                            11d26s20
       call dws_synca                                                   12s15s21
      norbx=norb+1                                                      11d25s20
      i1=1                                                              11d25s20
      ifirsttime=0                                                      11d26s20
      do isb=1,nsymb                                                    11d25s20
       isbv=multh(isb,isymmrci)
       do nclo1p=mdon+1,mdoobp                                          10d20s21
        if(min(nvirt(isbv),nff1(nclo1p,isb,1)).gt.0)then                  11d25s20
         nclo1=nclo1p-1                                                  11d25s20
         nopen1=nec-nclo1*2                                             11d25s20
         iarg=nclo1p-mdon                                                11d25s20
         ncolt=nff1(nclo1p,isb,1)*ncsfs(iarg)                           10d20s21
         call ilimts(1,ncolt,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  12d7s20
         ncolh=ih+1-il                                                  12d7s20
         igflag=0                                                       3d22s21
         nwds=nvirt(isbv)*ncolh*nroot                                   3d22s21
         nn=ncsfs(iarg)*nroot                                           10d20s21
         nnm=nn-1                                                       11d26s20
         isbo=multh(isbv,isymc)                                         10d20s21
         idenh=ibcoff                                                   11d25s20
         ndenh=idenh+nn*irefo(isbo)                                     10d20s21
         ibcoff=ndenh+irefo(isbo)                                       10d20s21
         do isa=1,nsymb                                                 11d26s20
          isav=multh(isa,isbo)                                          7d12s21
          do isc=1,nsymb                                                7d14s21
           nnv=irefo(isa)*irefo(isc)                                    7d12s21
           iscav=multh(isc,isav)                                        11d26s20
           iden1x(isc,isa)=ibcoff                                       10d20s21
           nden1x(isc,isa)=iden1x(isc,isa)+nn*                          10d20s21
     $           nnv*irefo(iscav)*ntype                                 10d20s21
           ibcoff=nden1x(isc,isa)+nnv*irefo(iscav)*ntype                10d20s21
          end do                                                        11d26s20
         end do                                                         11d26s20
         call enough('hcsibk4.  2',bc,ibc)
         ibctop=ibcoff-1                                                11d26s20
         jff1=nff1(nclo1p,isb,2)                                        11d25s20
         ifcol=1                                                        12d7s20
         do if1=1,nff1(nclo1p,isb,1)                                    11d26s20
          itcol=ifcol+ncsfs(iarg)-1                                     10d20s21
          if(ifcol.le.ih.and.itcol.ge.il)then                           12d7s20
           nnlow=max(0,il-ifcol)                                        12d7s20
           nnhgh=min(ncsfs(iarg)-1,ih-ifcol)                            10d20s21
           nnuse=nnhgh+1-nnlow                                          12d7s20
           nnuser=nnuse*nroot                                           11d10s20
           do i=idenh,ibctop                                             11d26s20
            bc(i)=0d0                                                    11d26s20
           end do                                                        11d26s20
           i1c=iff1(jff1)                                                11d25s20
           ntest=popcnt(i1c)
           do i=1,norb                                                   11d25s20
            itest(i,1)=0                                                 11d25s20
           end do                                                        11d25s20
           ii=1                                                          11d25s20
           do i=1,norb                                                   11d25s20
            if(btest(i1c,i))then                                         11d25s20
             idorb1(ii)=i                                                11d25s20
             itest(i,1)=2                                                11d25s20
             ii=ii+1                                                     11d25s20
            end if                                                       11d25s20
           end do                                                        11d25s20
           jff1=jff1+1                                                   11d25s20
           i1o=iff1(jff1)                                                11d25s20
           i1o=ibset(i1o,norbx)                                          11d25s20
           ii=1                                                          11d25s20
           do i=1,norb                                                   11d25s20
            if(btest(i1o,i))then                                         11d25s20
             isorb1(ii)=i                                                11d25s20
             itest(i,1)=1                                                11d25s20
             ii=ii+1                                                     11d25s20
            end if                                                       11d25s20
           end do                                                        11d25s20
           jff1=jff1+1                                                   11d25s20
           ngot=0
           do nclop=max(mdon+1,nclo1p-2),min(mdookp,nclo1p+2)           10d20s21
            nclo=nclop-1                                                  11d25s20
            nopen=nec-nclo*2                                             11d25s20
            jarg=nclop-mdon                                               11d25s20
            jff0=nff0(nclop,2)                                           11d25s20
            jvs=nff0(nclop,3)                                            11d26s20
            call wfetch(nopen,mdon,idum,i2sk,i2smk,ndet,ncsfk,iad1,iad2,10d22s21
     $           idum,bc,ibc)                                           11d9s22
            do jf=1,nff0(nclop,1)                                        11d26s20
             i0c=iff0(jff0)                                              11d25s20
             jff0=jff0+1                                                 11d25s20
             i0o=iff0(jff0)                                              11d25s20
             jff0=jff0+1                                                 11d25s20
             gandcc=ieor(i1c,i0c)                                        2d6s23
             gandco=ieor(i1o,i0o)                                        2d6s23
             gandcb=ior(gandcc,gandco)                                     2d6s23
             ndifb=popcnt(gandcb)                                          2d6s23
             if(ndifb.le.4)then                                         2d7s23
              ndifs=popcnt(gandco)                                        2d6s23
              ndifd=popcnt(gandcc)                                        2d6s23
              if(ndifs.eq.2.and.ndifb.eq.2)then                         2d6s23
               do i=1,norbx                                             2d6s23
                if(btest(gandco,i))then                                   2d6s23
                 if((btest(i1o,i).and..not.btest(i0c,i)).or.            2d6s23
     $             (btest(i1c,i).and.btest(i0o,i)))then                 2d6s23
                  nab4(1,1)=i                                             2d6s23
                 else                                                     2d6s23
                  nab4(2,1)=i                                             2d6s23
                 end if                                                   2d6s23
                end if                                                    2d6s23
               end do                                                     2d6s23
               call gandcr(i1c,i1o,i0c,i0o,nopen1,nopen,                 10d20s21
     $            norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,iwpb,iwpk,11d14s22
     $            bc,ibc)                                               11d14s22
               call gencup(i2sb,i2smb,i2sk,i2smk,nopen1,nopen,nab1,iwpb, 9d17s21
     $            iwpk,iouta,imata,ntypea,mcsf,vec(jvs,1),ncsfv,        10d20s21
     $            nroot,bc,ibc)                                         11d14s22
               ii=1                                                          11d25s20
               do i=1,norb                                                   11d25s20
                itest(i,2)=0                                              11d25s20
                if(btest(i0c,i))then                                         11d25s20
                 idorb(ii)=i                                                11d25s20
                 ii=ii+1                                                     11d25s20
                 itest(i,2)=2                                             11d25s20
                end if                                                       11d25s20
               end do                                                     11d25s20
               ii=1                                                          11d25s20
               do i=1,norb                                                   11d25s20
                if(btest(i0o,i))then                                         11d25s20
                 isorb(ii)=i                                                11d25s20
                 ii=ii+1                                                     11d25s20
                 itest(i,2)=1                                             11d25s20
                end if                                                       11d25s20
               end do                                                        11d25s20
               nok=0                                                      11d25s20
               do i=1,norb                                                11d25s20
                ixn=min(itest(i,1),itest(i,2))                            11d25s20
                if(ixn.gt.0)then                                          11d25s20
                 nok=nok+1                                                11d25s20
                 itest(nok,3)=ixn                                         11d25s20
                 itest(nok,4)=i                                           11d25s20
                end if                                                    11d25s20
               end do                                                     11d25s20
               do ii=0,ntypea-1                                           10d20s21
                ipackc=ibc(iouta+ii)                                      9d10s21
                imat=imata+mcsf(1)*nroot*ii                               10d20s21
                if(ipackc1(1).gt.0)then                                   9d10s21
                 imy(3)=0                                                 10d12s21
                else
                 imy(3)=1                                                 10d12s21
                end if
                if(ipackc1(2).gt.0)then                                   9d10s21
                 lsk=ism(ipackc1(2))                                      10d12s21
                 lgk=irel(ipackc1(2))-1                                   10d12s21
                 igya(4)=lgk                                              10d12s21
                 imy(4)=0                                                 10d12s21
                else
                 lsk=ism(-ipackc1(2))                                     10d12s21
                 igya(4)=irel(-ipackc1(2))-1                              10d12s21
                 imy(4)=1                                                 10d12s21
                end if
                if(lsk.ne.isbo)then                                       7d12s21
                 write(6,*)('lsk vs. isbo: '),lsk,isbo
                 stop 'hcsibk4'
                end if
                jdenh=idenh+nn*igya(4)                                    10d20s21
                phs=1d0                                                   10d20s21
                if(imy(3).eq.1.and.isopt(4).ne.0)phs=-1d0                 10d20s21
                ibc(ndenh+igya(4))=1                                         7d19s21
                do i=0,nnm                                                 11d26s20
                 bc(jdenh+i)=bc(jdenh+i)+phs*bc(imat+i)                   10d20s21
                end do                                                     11d26s20
                do i=1,norb                                               9d10s21
                 if(btest(i1c,i).and.btest(i0c,i))then                    10d20s21
                  js=ism(i)                                               9d10s21
                  jg=irel(i)-1                                            9d10s21
                  do jpass=0,1                                            9d10s21
c
c     integral is (jg,jg|v imy(4))
c     onex: b,a,d,c: here imy(3) is c and imy(4) is d, and b=a=jg.
c
c     in non-rel code, den1(c,a) went with integrals c,a,savc,v
c     with savc=isav*c=a*v*c
c     den1(c,a) has addressing
c
                   ltest=jpass+2*(jpass+2*(imy(3)+2*imy(4)))              10d20s21
                   itestp=ltest+1                                         10d12s21
                   jtest=1-jpass+2*(1-jpass+2*(1-imy(3)+2*(1-imy(4))))    10d20s21
                   jtestp=jtest+1                                         10d12s21
                   iuse=-1                                                10d22s21
                   phs=1d0                                               10d20s21
                   if(iifmx(itestp).ge.0)then                             10d20s21
                    iuse=iifmx(itestp)                                    10d20s21
                   else if(iifmx(jtestp).ge.0)then                        10d20s21
                    iuse=iifmx(jtestp)                                    10d20s21
                    if(isopt(4).ne.0)phs=-phs                             10d20s21
                   end if                                                 10d20s21
                   if(iuse.ge.0)then                                      10d20s21
                    iicolj=jg+irefo(js)*(jg+irefo(js)*(igya(4)            10d22s21
     $                 +irefo(isbo)*iuse))                               10d22s21
                    jden=iden1x(js,js)+nn*iicolj                          10d22s21
                    ibc(nden1x(js,js)+iicolj)=1                           10d22s21
                    do iii=0,nnm                                          10d22s21
                     bc(jden+iii)=bc(jden+iii)+bc(imat+iii)*phs           10d20s21
                    end do                                                10d22s21
                   end if                                                 10d20s21
c
c     integral is (jg,imy(4)|v jg)
c     onex: b,a,d,c: here imy(3) is c and imy(4) is b, and a=d=jg.
c
                   ltest=jpass+2*(imy(4)+2*(imy(3)+2*jpass))              10d22s21
                   itestp=ltest+1                                         10d12s21
                   jtest=1-jpass+2*(1-imy(4)+2*(1-imy(3)+2*(1-jpass)))    10d22s21
                   jtestp=jtest+1                                         10d12s21
                   iuse=-1                                               11d12s21
                   phs=1d0                                               10d20s21
                   if(iifmx(itestp).ge.0)then                             10d20s21
                    iuse=iifmx(itestp)                                    10d20s21
                    kchoice=iuse
                   else if(iifmx(jtestp).ge.0)then                        10d20s21
                    iuse=iifmx(jtestp)                                    10d20s21
                    if(isopt(4).ne.0)phs=-phs                             10d20s21
                    kchoice=100+iuse
                   end if                                                 10d20s21
                   if(iuse.ge.0)then                                      10d20s21
                    iicolj=igya(4)+irefo(isbo)*(jg+irefo(js)*(jg          10d20s21
     $                 +irefo(js)*iuse))                                 10d20s21
                    jden=iden1x(isbo,js)+nn*iicolj                        10d20s21
                    ibc(nden1x(isbo,js)+iicolj)=1                         10d20s21
                    do iii=0,nnm                                               11d26s20
                     bc(jden+iii)=bc(jden+iii)-bc(imat+iii)*phs           10d22s21
                    end do                                                   11d26s20
                   end if                                                 10d20s21
                  end do                                                  10d22s21
                 end if                                                   10d22s21
                end do                                                    10d22s21
                igya(2)=igya(4)                                           10d22s21
                if(ipackc1(2).gt.0)then                                   10d12s21
                 im=ipackc1(2)                                            10d12s21
                 imy(2)=1                                                 10d13s21
                else                                                      10d12s21
                 im=-ipackc1(2)                                           10d12s21
                 imy(2)=0                                                 10d13s21
                end if                                                    10d12s21
c
c     isbo*isbo*isbo*isbv=isbo*isbv=isbv*isymc*isbv=isymc,
c     so symmetry is always correct.
c
                if(btest(i0c,im))then                                    10d24s21
c
c     "y" term.
c     integral is (imy(2)imy(2)|v imy(4))
c     onex: b,a,d,c: here imy(3) is c and imy(4) is d, and b=a=imy(2)
c
                 ltest=imy(2)+2*(imy(2)+2*(imy(3)+2*imy(4)))              10d12s21
                 itestp=ltest+1                                           10d12s21
                 jtest=1-imy(2)+2*(1-imy(2)+2*(1-imy(3)+2*(1-imy(4))))    10d12s21
                 jtestp=jtest+1                                           10d12s21
                 iuse=-1                                                  10d22s21
                 phs=1d0                                                 10d22s21
                 if(iifmx(itestp).ge.0)then                               10d12s21
                  iuse=iifmx(itestp)                                      10d22s21
                 else if(iifmx(jtestp).ge.0)then                          10d12s21
                  iuse=iifmx(jtestp)                                      10d22s21
                  if(isopt(4).ne.0)phs=-phs                               10d22s21
                 end if                                                   10d12s21
                 if(iuse.ge.0)then                                        10d22s21
                  iicolj=igya(2)+irefo(isbo)*(igya(2)+irefo(isbo)*       10d22s21
     $               (igya(4)+irefo(isbo)*iuse))                        10d22s21
                  jden=iden1x(isbo,isbo)+nn*iicolj                        10d22s21
                  ibc(nden1x(isbo,isbo)+iicolj)=1                         10d22s21
                  do iii=0,nnm                                            10d22s21
                   bc(jden+iii)=bc(jden+iii)+bc(imat+iii)*phs             10d22s21
                  end do                                                  10d22s21
                 end if                                                   10d22s21
c     integral is (imy(2)imy(4)|v imy(2))
c     onex: b,a,d,c: here imy(3) is c and imy(2) is d=a, and b=imy(4)
c
                 ltest=imy(2)+2*(imy(4)+2*(imy(3)+2*imy(2)))              10d12s21
                 itestp=ltest+1                                           10d12s21
                 jtest=1-imy(2)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(2))))    10d12s21
                 jtestp=jtest+1                                           10d12s21
                 iuse=-1                                                  10d22s21
                 phs=1d0                                                 10d22s21
                 if(iifmx(itestp).ge.0)then                               10d12s21
                  iuse=iifmx(itestp)                                      10d22s21
                 else if(iifmx(jtestp).ge.0)then                          10d12s21
                  iuse=iifmx(jtestp)                                      10d22s21
                  if(isopt(4).ne.0)phs=-phs                               10d22s21
                 end if                                                   10d12s21
                 if(iuse.ge.0)then                                        10d22s21
                  iicolj=igya(2)+irefo(isbo)*(igya(4)+irefo(isbo)*(      10d22s21
     $               igya(2)+irefo(isbo)*iuse))                         10d22s21
                  jden=iden1x(isbo,isbo)+nn*iicolj                        10d22s21
                  ibc(nden1x(isbo,isbo)+iicolj)=1                         10d22s21
                  do iii=0,nnm                                            10d22s21
                   bc(jden+iii)=bc(jden+iii)-bc(imat+iii)*phs             10d22s21
                  end do                                                  10d22s21
                 end if                                                   10d22s21
                end if                                                    9d13s21
               end do                                                     10d22s21
               ibcoff=iouta                                               10d22s21
               do i=1,nok                                                    11d13s20
                if(itest(i,3).eq.1)then                                      11d13s20
                 itesta=i1c                                               11d26s20
                 itestb=i1o                                               11d26s20
                 nopenk=nopen1                                                11d13s20
c
c     anihilate common
c
                 if(btest(itesta,itest(i,4)))then                             11d13s20
                  itesta=ibclr(itesta,itest(i,4))                             11d13s20
                  itestb=ibset(itestb,itest(i,4))                             11d13s20
                  karg=iarg-1                                                11d13s20
                  nopenk=nopenk+1                                             11d13s20
                 else                                                         11d13s20
                  itestb=ibclr(itestb,itest(i,4))                             11d13s20
                  karg=iarg                                                  11d13s20
                  nopenk=nopenk-1                                             11d13s20
                 end if                                                       11d13s20
c
c     create ket
c
                 if(btest(itestb,nab4(2,1)))then                          11d27s20
                  itesta=ibset(itesta,nab4(2,1))                          11d27s20
                  itestb=ibclr(itestb,nab4(2,1))                          11d27s20
                  karg=karg+1                                                11d13s20
                  nopenk=nopenk-1                                             11d13s20
                 else                                                         11d13s20
                  itestb=ibset(itestb,nab4(2,1))                          11d27s20
                  nopenk=nopenk+1                                             11d13s20
                 end if                                                       11d13s20
                 call gandcr(i1c,i1o,itesta,itestb,nopen1,nopenk,        10d22s21
     $          norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,iwpb1,iwpk1,11d14s22
     $               bc,ibc)                                            11d14s22
c     nab2(1) is virt
                 call gandcr(itesta,itestb,i0c,i0o,nopenk,nopen,         10d22s21
     $          norbx,nnot2,nab2,icode1,imap1,nx1,irw1,irw2,iwpb2,iwpk2,11d14s22
     $               bc,ibc)                                            11d14s22
                 if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                  call spinloop(i2sb,i2smb,i2sk,i2smk,nopen1,nopen,      10d22s21
     $               nopenk,ncsfs(iarg),nroot,itype,imatx,ntypex,nab1,  10d24s21
     $               iwpb1,iwpk1,nab2,iwpb2,iwpk2,vec(jvs,1),ncsfv,     10d22s21
     $               ieoro,bc,ibc)                                      11d14s22
                  do ini=0,ntypex-1                                      10d24s21
                   ipack=ibc(itype+ini)
                   imat=imatx+nn*ini                                       10d22s21
                   if(ipack2(3).eq.norbx.or.-ipack2(3).eq.norbx)then     10d22s21
                   else
                    write(6,*)
     $                   ('oh no, virt is not were I expected it!!')
                    write(6,*)('ini = '),ini,itype+ini,ipack
                    write(6,*)ipack2
                    stop 'hcsibk4'
                   end if
                   isy(3)=isbv                                           10d22s21
                   if(ipack2(3).gt.0)then                                10d22s21
                    imy(3)=0                                             10d22s21
                   else                                                  10d22s21
                    imy(3)=1                                             10d22s21
                   end if                                                10d22s21
                   do j=1,4                                               9d9s21
                    if(j.ne.3)then                                       10d22s21
                     if(ipack2(j).gt.0)then                                9d9s21
                      isy(j)=ism(ipack2(j))                                9d9s21
                      igya(j)=irel(ipack2(j))-1                             9d9s21
                      imy(j)=0                                             10d13s21
                     else                                                  9d9s21
                      isy(j)=ism(-ipack2(j))                                9d9s21
                      igya(j)=irel(-ipack2(j))-1                           10d13s21
                      imy(j)=1                                             10d13s21
                     end if                                              10d22s21
                    end if                                                9d9s21
                   end do                                                 9d9s21
c
c                  ab cd
c     integral is (12|v4)
c     onex: b,a,d,c: here v(=3) is c, b=2, a=1, and d=4.
c
                   ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))            10d22s21
                   itestp=ltest+1                                         10d13s21
                   jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*(1-imy(4))))  10d13s21
                   jtestp=jtest+1                                         10d13s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolj=igya(2)+irefo(isy(2))*(igya(1)               2d7s23
     $                  +irefo(isy(1))*(                                2d7s23
     $                  igya(4)+irefo(isy(4))*iuse))                       10d22s21
                    jden=iden1x(isy(2),isy(1))+nn*iicolj                  10d22s21
                    ibc(nden1x(isy(2),isy(1))+iicolj)=1                   10d22s21
                    do iii=0,nnm                                            10d22s21
                     bc(jden+iii)=bc(jden+iii)+bc(imat+iii)*phs             10d22s21
                    end do                                                  10d22s21
                   end if                                                   10d22s21
c
c                  ab cd
c     integral is (14|v2)
c     onex: b,a,d,c: here v(=3) is c, b=1, a=4, and d=2.
c
                   ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            10d22s21
                   itestp=ltest+1                                         10d13s21
                   jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(2))))  10d13s21
                   jtestp=jtest+1                                         10d13s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolj=igya(4)+irefo(isy(4))*(igya(1)               2d7s23
     $                  +irefo(isy(1))*(                                2d7s23
     $                  igya(2)+irefo(isy(2))*iuse))                       10d22s21
                    jden=iden1x(isy(4),isy(1))+nn*iicolj                  10d22s21
                    ibc(nden1x(isy(4),isy(1))+iicolj)=1                   10d22s21
                    do iii=0,nnm                                            10d22s21
                     bc(jden+iii)=bc(jden+iii)-bc(imat+iii)*phs             10d22s21
                    end do                                                  10d22s21
                   end if                                                   10d22s21
                  end do                                                 10d22s21
                  ibcoff=itype                                           10d22s21
                 end if                                                    11d13s20
                end if                                                   10d22s21
               end do                                                        11d13s20
              else if(l2e.eq.0)then                                     2d7s23
               nnot=0
               if(ndifs.eq.4.and.ndifb.eq.4)then                         2d6s23
                nnot=4                                                   2d6s23
                ioxx(1)=1                                                2d6s23
                ioxx(2)=1                                                2d6s23
                do i=1,norbx                                             2d6s23
                 if(btest(gandcb,i))then                                 2d6s23
                  if((btest(i0c,i).and.btest(i1o,i)).or.                 2d6s23
     $               (btest(i0o,i).and..not.btest(i1c,i)))then          2d6s23
                   nab4(2,ioxx(2))=i                                     2d6s23
                   ioxx(2)=ioxx(2)+1                                     2d6s23
                  else                                                   2d6s23
                   nab4(1,ioxx(1))=i                                     2d6s23
                   ioxx(1)=ioxx(1)+1                                     2d6s23
                  end if                                                 2d6s23
                 end if                                                  2d6s23
                end do                                                   2d6s23
               else if(ndifb.eq.3)then                                   2d6s23
                nnot=3                                                   2d6s23
                ioxx(1)=1                                                2d6s23
                ioxx(2)=1                                                2d6s23
                iswap=0                                                  2d6s23
                do i=1,norbx                                             2d6s23
                 if(btest(gandcb,i))then                                 2d6s23
                  if(btest(gandcc,i).and.                                2d6s23
     $               ((btest(i1c,i).and..not.btest(i0o,i)).or.          2d6s23
     $               (btest(i0c,i).and..not.btest(i1o,i))))then         2d6s23
                   if(btest(i0c,i))iswap=1                               2d6s23
                   nab4(1,1)=i                                           2d6s23
                   nab4(1,2)=i                                           2d6s23
                  else                                                   2d6s23
                   nab4(2,ioxx(2))=i                                     2d6s23
                   ioxx(2)=ioxx(2)+1                                     2d6s23
                  end if                                                 2d6s23
                 end if                                                  2d6s23
                end do                                                   2d6s23
                if(iswap.ne.0)then                                       2d6s23
                 icpy=nab4(1,1)                                          2d6s23
                 nab4(1,1)=nab4(2,1)                                     2d6s23
                 nab4(2,1)=icpy                                          2d6s23
                 icpy=nab4(1,2)                                          2d6s23
                 nab4(1,2)=nab4(2,2)                                     2d6s23
                 nab4(2,2)=icpy                                          2d6s23
                 nbt=0                                                   2d6s23
                 if(btest(i0c,nab4(2,2)).and.                            2d6s23
     $                   .not.btest(i0c,nab4(2,1)))nbt=1                2d6s23
                else                                                     2d6s23
                 nbt=0                                                   2d6s23
                 if(btest(i1c,nab4(1,2)).and.                            2d6s23
     $                .not.btest(i1c,nab4(1,1)))nbt=1                    2d6s23
                end if                                                   2d6s23
                if(nbt.ne.0)then                                         2d6s23
                 nab4(1,1)=nab4(1,2)                                     2d6s23
                 nab4(2,1)=nab4(2,2)                                     2d6s23
                end if                                                   2d6s23
               else if(ndifs.eq.0.and.ndifd.eq.2)then                    2d6s23
                nnot=3                                                   2d6s23
                do i=1,norbx                                             2d6s23
                 if(btest(gandcb,i))then                                 2d6s23
                  if(btest(i1c,i))then                                   2d6s23
                   nab4(1,1)=i                                           2d6s23
                   nab4(1,2)=i                                           2d6s23
                  else                                                   2d6s23
                   nab4(2,1)=i                                           2d6s23
                   nab4(2,2)=i                                           2d6s23
                  end if                                                 2d6s23
                 end if                                                  2d6s23
                end do                                                   2d6s23
               end if                                                    2d6s23
               if(nnot.ne.0)then                                        2d7s23
                iu1=1                                                     1d25s21
                iu2=1                                                     1d25s21
c     sing is bra
                itesta=i1c                                                 11d26s20
                itestb=i1o                                                 11d26s20
                if(btest(itesta,nab4(1,iu1)))then                         1d25s21
                 itesta=ibclr(itesta,nab4(1,iu1))                         1d25s21
                 itestb=ibset(itestb,nab4(1,iu1))                         1d25s21
                 nopenk=nopen1+1                                              11d13s20
                 karg=iarg-1                                                  11d13s20
                else if(btest(itestb,nab4(1,iu1)))then                    1d25s21
                 itestb=ibclr(itestb,nab4(1,iu1))                         1d25s21
                 nopenk=nopen1-1                                              11d13s20
                 karg=iarg                                                  11d13s20
                end if                                                        11d13s20
                if(btest(itestb,nab4(2,iu2)))then                         1d25s21
                 itesta=ibset(itesta,nab4(2,iu2))                         1d25s21
                 itestb=ibclr(itestb,nab4(2,iu2))                         1d25s21
                 nopenk=nopenk-1                                              11d13s20
                 karg=karg+1                                                  11d13s20
                else                                                          11d13s20
                 itestb=ibset(itestb,nab4(2,iu2))                         1d25s21
                 nopenk=nopenk+1                                              11d13s20
                end if                                                        11d13s20
                call gandcr(i1c,i1o,itesta,itestb,nopen1,nopenk,            11d26s20
     $         norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,iwpb1,iwpk1, 11d14s22
     $             bc,ibc)                                              11d14s22
                call gandcr(itesta,itestb,i0c,i0o,nopenk,nopen,             11d26s20
     $         norbx,nnot2,nab2,icode1,imap1,nx1,irw1,irw2,iwpb2,iwpk2, 11d14s22
     $             bc,ibc)                                              11d14s22
                if(nnot1.eq.2.and.nnot2.eq.2)then                         1d25s21
                 call spinloop(i2sb,i2smb,i2sk,i2smk,nopen1,nopen,        10d22s21
     $               nopenk,ncsfs(iarg),nroot,itype,imatx,ntypex,nab1,  10d24s21
     $             iwpb1,iwpk1,nab2,iwpb2,iwpk2,vec(jvs,1),ncsfv,ieoro, 11d14s22
     $             bc,ibc)                                              11d14s22
                 do i=0,ntypex-1                                          10d24s21
                  ipack=ibc(itype+i)
                  imat=imatx+nn*i                                         10d22s21
                  do j=1,4                                                9d9s21
                   ltest=ipack2(j)                                        10d22s21
                   if(iabs(ltest).eq.norbx)then                           10d22s21
                    isy(j)=isbv                                           10d22s21
                    if(ltest.gt.0)then                                    10d22s21
                     imy(j)=0                                             10d22s21
                    else                                                  10d22s21
                     imy(j)=1                                             10d22s21
                    end if                                                10d22s21
                   else                                                   10d22s21
                    if(ipack2(j).gt.0)then                                9d9s21
                     isy(j)=ism(ipack2(j))                                9d9s21
                     igya(j)=irel(ipack2(j))-1                             9d9s21
                     imy(j)=0                                             10d12s21
                    else                                                  9d9s21
                     isy(j)=ism(-ipack2(j))                                9d9s21
                     igya(j)=irel(-ipack2(j))-1                           10d12s21
                     imy(j)=1                                             10d12s21
                    end if                                                9d9s21
                   end if                                                 10d22s21
                  end do                                                  10d22s21
c
c     v is 1 or 3.
c
                  if(ipack2(1).eq.norbx.or.-ipack2(1).eq.norbx)then       10d22s21
c
c     if v is 1, integral is (v2|34)=(34|v2)
c     onex: b,a,d,c: here v(=1) is c, b=4, a=3, and d=2.
c
                   ltest=imy(3)+2*(imy(4)+2*(imy(1)+2*imy(2)))            10d22s21
                   itestp=ltest+1                                         10d22s21
                   jtest=1-imy(3)+2*((1-imy(4))+2*((1-imy(1))             10d22s21
     $                 +2*(1-imy(2))))                                  10d22s21
                   jtestp=jtest+1                                         10d22s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolj=igya(4)+irefo(isy(4))*(igya(3)               2d7s23
     $                  +irefo(isy(3))*(                                2d7s23
     $               igya(2)+irefo(isy(2))*iuse))                       10d22s21
                    jden=iden1x(isy(4),isy(3))+nn*iicolj                  10d22s21
                    ibc(nden1x(isy(4),isy(3))+iicolj)=1                   10d22s21
                    do iii=0,nnm                                            10d22s21
                     bc(jden+iii)=bc(jden+iii)+bc(imat+iii)*phs             10d22s21
                    end do                                                  10d22s21
                   end if                                                   10d22s21
                   if(nnot.eq.4)then                                      10d22s21
c
c     integral is (v4|32)=(32|v4)
c     onex: b,a,d,c: here v(=1) is c, b=2, a=3, and d=4.
c
                    ltest=imy(3)+2*(imy(2)+2*(imy(1)+2*imy(4)))           10d22s21
                    itestp=ltest+1                                        10d22s21
                    jtest=1-imy(3)+2*((1-imy(2))+2*((1-imy(1))            10d22s21
     $                  +2*(1-imy(4))))                                  10d22s21
                    jtestp=jtest+1                                        10d22s21
                    iuse=-1                                                  10d22s21
                    phs=1d0                                                 10d22s21
                    if(iifmx(itestp).ge.0)then                               10d12s21
                     iuse=iifmx(itestp)                                      10d22s21
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iuse=iifmx(jtestp)                                      10d22s21
                     if(isopt(4).ne.0)phs=-phs                               10d22s21
                    end if                                                   10d12s21
                    if(iuse.ge.0)then                                        10d22s21
                     iicolj=igya(2)+irefo(isy(2))*(igya(3)               2d7s23
     $                  +irefo(isy(3))*(                                2d7s23
     $                   igya(4)+irefo(isy(4))*iuse))                       10d22s21
                     jden=iden1x(isy(2),isy(3))+nn*iicolj                  10d22s21
                     ibc(nden1x(isy(2),isy(3))+iicolj)=1                   10d22s21
                     do iii=0,nnm                                            10d22s21
                      bc(jden+iii)=bc(jden+iii)-bc(imat+iii)*phs             10d22s21
                     end do                                                  10d22s21
                    end if                                                   10d22s21
                   end if                                                 10d22s21
                  else if(ipack2(3).eq.norbx.or.-ipack2(3).eq.norbx)then  10d22s21
c
c     if v is 3, integral is (12|v4)
c     onex: b,a,d,c: here v(=3) is c, b=2, a=1, and d=4.
c
                   ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))            10d22s21
                   itestp=ltest+1                                         10d22s21
                   jtest=1-imy(1)+2*((1-imy(2))+2*((1-imy(3))             10d22s21
     $                  +2*(1-imy(4))))                                  10d22s21
                   jtestp=jtest+1                                         10d22s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolj=igya(2)+irefo(isy(2))*(igya(1)               2d7s23
     $                  +irefo(isy(1))*(                                2d7s23
     $               igya(4)+irefo(isy(4))*iuse))                       10d22s21
                    jden=iden1x(isy(2),isy(1))+nn*iicolj                  10d22s21
                    ibc(nden1x(isy(2),isy(1))+iicolj)=1                   10d22s21
                    do iii=0,nnm                                            10d22s21
                     bc(jden+iii)=bc(jden+iii)+bc(imat+iii)*phs             10d22s21
                    end do                                                  10d22s21
                   end if                                                   10d22s21
                   if(nnot.eq.4)then                                      10d22s21
c
c     integral is (14|v2)
c     onex: b,a,d,c: here v(=3) is c, b=4, a=1, and d=2.
c
                    ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))           10d22s21
                    itestp=ltest+1                                        10d22s21
                    jtest=1-imy(1)+2*((1-imy(4))+2*((1-imy(3))            10d22s21
     $                  +2*(1-imy(2))))                                  10d22s21
                    jtestp=jtest+1                                        10d22s21
                    iuse=-1                                                  10d22s21
                    phs=1d0                                                 10d22s21
                    if(iifmx(itestp).ge.0)then                               10d12s21
                     iuse=iifmx(itestp)                                      10d22s21
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iuse=iifmx(jtestp)                                      10d22s21
                     if(isopt(4).ne.0)phs=-phs                               10d22s21
                    end if                                                   10d12s21
                    if(iuse.ge.0)then                                        10d22s21
                     iicolj=igya(4)+irefo(isy(4))*(igya(1)               2d7s23
     $                  +irefo(isy(1))*(                                2d7s23
     $               igya(2)+irefo(isy(2))*iuse))                       10d22s21
                     jden=iden1x(isy(4),isy(1))+nn*iicolj                  10d22s21
                     ibc(nden1x(isy(4),isy(1))+iicolj)=1                   10d22s21
                     do iii=0,nnm                                            10d22s21
                      bc(jden+iii)=bc(jden+iii)-bc(imat+iii)*phs             10d22s21
                     end do                                                  10d22s21
                    end if                                                   10d22s21
                   end if                                                 10d22s21
                  else                                                    10d22s21
                   write(6,*)('oh no!! v is in wrong place!!')            10d22s21
                   stop 'hcsibk4'
                  end if                                                  10d22s21
                 end do                                                   10d22s21
                 ibcoff=itype                                             10d22s21
                end if                                                    10d22s21
               end if                                                   2d7s23
              end if                                                    2d7s23
             end if                                                      11d25s20
             jvs=jvs+ncsfk                                              10d22s21
            end do                                                       11d25s20
           end do                                                        11d25s20
           itmpg=ibcoff                                                 12d7s20
           ibcoff=itmpg+nnuser*nvirt(isbv)                              11d10s20
           call enough('hcsibk4.  3',bc,ibc)
           do iz=itmpg,ibcoff-1                                         10d22s21
            bc(iz)=0d0                                                  10d22s21
           end do                                                       10d22s21
           factg=0d0                                                    12d7s20
           if(irefo(isbo).gt.0)then                                     7d12s21
            nok=0                                                        12d2s20
            jdenh=idenh                                                  12d2s20
            ikeep=ibcoff                                                 12d2s20
            ibcoff=ikeep+irefo(isbo)                                    7d12s21
            call enough('hcsibk4.  4',bc,ibc)
            do i=0,irefo(isbo)-1                                        7d12s21
             if(ibc(ndenh+i).gt.0)then                                  3d19s21
              iad=idenh+nn*i
              ibc(ikeep+nok)=i                                          10d22s21
              do ir=0,nroot-1                                           12d12s20
               do j=nnlow,nnhgh                                          12d12s20
                bc(jdenh)=bc(iad+j)                                      12d12s20
                jdenh=jdenh+1                                            12d12s20
               end do                                                   12d12s20
               iad=iad+ncsfs(iarg)                                       12d12s20
              end do                                                     12d2s20
              nok=nok+1                                                  12d2s20
             end if                                                      12d2s20
            end do
            if(nok.gt.0)then                                             12d2s20
             itmp=ibcoff                                                 12d2s20
             ibcoff=itmp+nvirt(isbv)*nok                                 12d2s20
             call enough('hcsibk4.  5',bc,ibc)
             if(ih0n(isbv).gt.0)then                                    10d22s21
              do iv=0,nvirt(isbv)-1                                       12d2s20
               iad1=itmp+nok*iv                                           12d2s20
               iad2=ih0n(isbv)+irefo(isbv)+iv                            10d22s21
               do i=0,nok-1                                               12d2s20
                ip=ibc(ikeep+i)                                          7d19s21
                bc(iad1+i)=bc(iad2+ip*nh0(isbv))                         10d22s21
               end do                                                     12d2s20
              end do                                                      12d2s20
             else if(ih0n(isbo).gt.0)then                               10d22s21
              phs=1d0                                                   10d22s21
              if(isopt(2).ne.0)phs=-phs                                 10d25s21
              if(isopt(4).ne.0.and.isopt(3).ne.0)phs=-phs               10d25s21
              do iv=0,nvirt(isbv)-1                                       12d2s20
               iad1=itmp+nok*iv                                           12d2s20
               do i=0,nok-1                                               12d2s20
                ip=ibc(ikeep+i)                                          7d19s21
                iad2=ih0n(isbo)+ip+nh0(isbo)*(irefo(isbv)+iv)           10d22s21
                bc(iad1+i)=bc(iad2)*phs                                 10d24s21
               end do                                                     12d2s20
              end do                                                      12d2s20
             else                                                       10d22s21
              write(6,*)('ih0n(isbv='),isbv,(') = '),ih0n(isbv)
              write(6,*)('and ih0n(isbo='),isbo,(') = '),ih0n(isbo)
              stop 'hcsibk4'
             end if                                                     10d22s21
             jdenh=idenh+nnlow                                          12d12s20
             call dgemm('n','n',nnuser,nvirt(isbv),nok,1d0,bc(idenh),   10d22s21
     $            nnuser,bc(itmp),nok,factg,bc(itmpg),nnuser,           12d12s20
     d' hcsibk4.  1')
             factg=1d0                                                  12d7s20
            end if                                                       12d2s20
            ibcoff=ikeep                                                 12d2s20
           end if
           do isa=1,nsymb                                                11d26s20
            isav=multh(isa,isbo)                                        7d12s21
            do isc=1,nsymb                                              7d14s21
             isavc=multh(isav,isc)                                       11d26s20
             nnv=irefo(isc)*irefo(isa)
             nnn=nnv*irefo(isavc)                                       10d22s21
             nnvn=nnv*ntype                                             10d22s21
             ncol=nnvn*irefo(isavc)                                     10d22s21
             if(ncol.gt.0.and.l2e.eq.0)then                             2d17s22
              nok=0                                                     12d7s20
              do i=0,ncol-1
               if(ibc(nden1x(isc,isa)+i).ne.0)then                      10d22s21
                ibc(nden1x(isc,isa)+nok)=i                               10d22s21
                nok=nok+1                                               3d22s21
               end if
              end do
              if(nok.gt.0)then                                          10d22s21
               itmpi=ibcoff                                             12d7s20
               ibcoff=itmpi+nok*nvirt(isbv)                             12d7s20
               call enough('hcsibk4.  6',bc,ibc)
               jden=iden1x(isc,isa)                                     10d22s21
               jtmpi=itmpi                                              12d7s20
               do iz=itmpi,ibcoff-1                                     10d22s21
                bc(iz)=0d0                                              10d22s21
               end do                                                   10d22s21
               do i=0,nok-1                                             3d22s21
                icol=ibc(nden1x(isc,isa)+i)                             10d22s21
                ifrm=iden1x(isc,isa)+nn*icol                            10d22s21
                do ir=0,nroot-1                                         12d12s20
                 do j=nnlow,nnhgh                                       12d12s20
                  bc(jden)=bc(ifrm+j)                                   12d12s20
                  jden=jden+1                                           12d12s20
                 end do                                                 12d12s20
                 ifrm=ifrm+ncsfs(iarg)                                   12d12s20
                end do                                                  12d12s20
                iq=icol/nnn                                             10d24s21
                ileft=icol-iq*nnn                                       10d24s21
                id=ileft/nnv                                            10d24s21
                ileft=ileft-id*nnv                                      10d24s21
                ia=ileft/irefo(isc)                                     1d25s21
                ic=ileft-ia*irefo(isc)                                  1d25s21
                iad=ionexr(isc,isa,isavc)+ic+irefo(isc)*(ia+irefo(isa)  10d22s21
     $               *(id+irefo(isavc)*nvirt(isbv)*iq))                 10d22s21
                ktmpi=jtmpi                                             7d9s21
                do iv=0,nvirt(isbv)-1                                   7d9s21
                 bc(ktmpi)=bc(ktmpi)+bc(iad+iv*nnn)                     10d22s21
                 ktmpi=ktmpi+nok                                        7d9s21
                end do                                                  7d9s21
                jtmpi=jtmpi+1
               end do                                                   12d7s20
               call dgemm('n','n',nnuser,nvirt(isbv),nok,1d0,           10d22s21
     $              bc(iden1x(isc,isa)),nnuser,bc(itmpi),nok,factg,     10d22s21
     $              bc(itmpg),nnuser,                                   12d12s20
     d' hcsibk4.  2')
               factg=1d0                                                12d7s20
              end if                                                    12d7s20
             end if
            end do
           end do                                                        11d26s20
           if(abs(factg-1d0).lt.1d-3)then                               12d7s20
            if(igflag.eq.0)then                                         3d22s21
             call ddi_done(ibc(iacc),nacc)                                     3d22s21
             do i=0,nwds-1                                                  3d22s21
              bc(ighere+i)=0d0                                          3d22s21
             end do                                                         3d22s21
             igflag=1                                                   3d22s21
            end if                                                      3d22s21
            do iv=0,nvirt(isbv)-1                                       12d7s20
             do i=nnlow,nnhgh                                           12d7s20
              ip=ifcol+i-il                                             12d7s20
              im=i-nnlow                                                12d7s20
              do ir=0,nroot-1
               iad1=itmpg+im+nnuse*(ir+nroot*iv)                        11d10s20
               iad2=ighere+ir+nroot*(iv+nvirt(isbv)*ip)                 12d14s20
               bc(iad2)=bc(iad2)+bc(iad1)                                12d7s20
              end do                                                    11d10s20
             end do                                                     12d7s20
            end do                                                      12d7s20
           end if                                                       12d7s20
          else                                                          12d8s20
           jff1=jff1+2                                                  12d8s20
          end if                                                        12d7s20
          ifcol=itcol+1                                                 12d7s20
         end do                                                         11d25s20
         if(ifirsttime.eq.0)call dws_synca                              11d26s20
         ifirsttime=1                                                   11d26s20
         i1c=nvirt(isbv)*nroot                                          11d10s20
         i1o=il                                                         12d7s20
         in=ih                                                          12d7s20
         if(igflag.ne.0)then                                            3d22s21
          call ddi_iacc(bc,ibc,nhand(nclo1p,isb),i1,i1c,i1o,in,         11d15s22
     $        bc(ighere),ibc(iacc),nacc)                                11d15s22
          nwiacc=nwiacc+i1c*(ih+1-il)                                   1d26s23
         end if                                                         3d22s21
         ibcoff=idenh                                                   3d22s21
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      call ddi_done(ibc(iacc),nacc)                                     3d22s21
      ibcoff=ibcoffo                                                    3d22s21
      if(ldebug)then                                                    1d25s21
       write(6,*)('all done in hcsibk')
       call dws_synca                                                    1d24s21
       write(6,*)('what we have for hcsibk:')
       do isb=1,nsymb                                                    12d7s20
        write(6,*)('for isb = '),isb
        isbv=multh(isb,isymmrci)                                         12d7s20
        i1c=nvirt(isbv)*nroot                                            11d10s20
        nhere=0                                                          12d7s20
        do ii=1,mdoobp                                                   7d9s21
         iarg=ii-mdon                                                    12d7s20
         ndblock=ncsfs(iarg)*nff1(ii,isb,1)                               11d10s20
         nhere=nhere+ndblock
        end do
        mdblock=nhere                                                    11d10s20
        nsing=mdblock*nvirt(isbv)                                        12d7s20
        igs=ibcoff
        ibcoff=igs+nhere*nvirt(isbv)*nroot                               11d10s20
        call enough('hcsibk4.  7',bc,ibc)
        jgs=igs
        do ii=1,mdoobp                                                   7d9s21
         iarg=ii-mdon                                                    12d7s20
         ndblock=ncsfs(iarg)*nff1(ii,isb,1)                               11d10s20
         if(ndblock.gt.0)then
          itmp=ibcoff                                                    12d7s20
          ibcoff=itmp+ndblock*nvirt(isbv)*nroot                          11d10s20
          call enough('hcsibk4.  8',bc,ibc)
          in=ndblock                                                     12d7s20
          write(6,*)('getting '),ii,isb,nhand(ii,isb),i1,i1c,i1,in
          call ddi_get(bc,ibc,nhand(ii,isb),i1,i1c,i1,in,bc(itmp))      11d15s22
          jtmp=itmp                                                      12d7s20
          do iff=1,nff1(ii,isb,1)                                        12d7s20
           do m=0,ncsfs(iarg)-1
            do iv=0,nvirt(isbv)-1
             do ir=0,nroot-1
              kgs=jgs+nsing*ir+iv                                        12d14s20
              bc(kgs)=bc(jtmp+ir)                                        12d14s20
             end do
             jtmp=jtmp+nroot                                             12d14s20
            end do
            jgs=jgs+nvirt(isbv)                                          11d10s20
           end do
          end do                                                         12d7s20
          ibcoff=itmp                                                    12d7s20
         end if
        end do                                                           12d7s20
        write(6,*)('all together now')
        call prntm2(bc(igs),nsing,nroot,nsing)                           11d10s20
        ibcoff=igs
       end do                                                            12d7s20
      end if                                                            1d25s21
      return
      end
