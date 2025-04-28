c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcisbk4(ihsdiag,nff1,iff1,i2sk,i2smk,mdookp,ncsfs,
     $     nff0,iff0,gi,ncsfv,i2sb,i2smb,mdoobp,nec,                    10d25s21
     $     mdon,nsymb,multh,irw1,irw2,nvirt,nrootu,                     10d25s21
     $     lprt,isymc,irorip,isopt,ism,irel,irefo,norb,                       10d25s21
     $     ih0n,nh0,ionexr,iifmx,ntype,isymmrci,l2e,bc,ibc)             11d9s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2),imap1(64),icode1(64),ipackc1(8)         10d26s21
      integer*2 ipack2(4)                                               10d26s21
      logical lprt                                                      1d25s21
      integer*8 i1,in,i1c,i1o,i0c,i0o,itesta,itestb,i28,i38,i48,        7d9s21
     $     ihsdiag(mdookp,nsymb,2),ipackc,ipack,gandcc,gandco,gandcb    2d7s23
      equivalence (ipackc1,ipackc),(ipack2,ipack)                       10d26s21
      dimension nff1(mdookp,nsymb,2),iff1(*),nff0(mdoobp,3),gi(ncsfv,*),  7d9s21
     $     iff0(*),ncsfs(*),multh(8,8),nab4(2,3),idorb1(32),ioxx(2),    2d7s23
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),iwpb(4),iwpk(4),  10d26s21
     $     nvirt(*),iden1x(8,8),nden1x(8,8),ism(*),irel(*),irefo(*),    7d9s21
     $     ih0n(*),nh0(*),ionexr(8,8,8),iifmx(*),iwpb1(4),iwpk1(4),     10d26s21
     $     iwpb2(4),iwpk2(4),mcsf(2),isy(4),imy(4),igya(4),isopt(*)     10d26s21
      include "common.store"                                            11d25s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      irori=irorip-1                                                    10d25s21
      if(lprt)write(6,*)('hi, I am the new and improved hcisbk4!')
      if(isopt(3).ne.0)then                                             10d28s21
       ieoro=1-irori                                                    10d28s21
      else                                                              10d28s21
       ieoro=irori                                                      10d28s21
      end if                                                            10d28s21
      nrootm=nrootu-1                                                   7d12s21
      i1=1                                                              11d10s20
      ig=ibcoff                                                         12d12s20
      ibcoff=ig+ncsfv*nrootu                                            12d12s20
      call enough('hcisbk4.  1',bc,ibc)
      do i=ig,ibcoff-1                                                  12d12s20
       bc(i)=0d0                                                        12d12s20
      end do                                                            12d12s20
      ibcoffo=ibcoff                                                    12d8s20
      norbx=norb+1                                                      11d25s20
      i1=1                                                              11d25s20
      ivs=1
      do isb=1,nsymb                                                    11d25s20
       isbv=multh(isb,isymmrci)
       isbo=multh(isbv,isymc)                                           10d26s21
       do nclo1p=mdon+1,mdookp                                          10d26s21
        if(min(nvirt(isbv),nff1(nclo1p,isb,1)).gt.0)then                  11d25s20
         nclo1=nclo1p-1                                                  11d25s20
         nopen1=nec-nclo1*2                                             11d25s20
         iarg=nclo1p-mdon                                                11d25s20
         ncolt=nff1(nclo1p,isb,1)*ncsfs(iarg)                           10d26s21
         call ilimts(1,ncolt,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  12d7s20
         ncolh=ih+1-il                                                  12d7s20
         ivst=ibcoff                                                    12d14s20
         ibcoff=ivst+nvirt(isbv)*nrootu*ncolh                           12d14s20
         ivtmp=ibcoff                                                   7d9s21
         ibcoff=ivtmp+nvirt(isbv)*nrootu*ncolh                          7d9s21
         call enough('hcisbk4.  2',bc,ibc)
         i28=nvirt(isbv)*nrootu                                         7d9s21
         i38=il                                                         7d9s21
         i48=ih                                                         7d9s21
         call ddi_get(bc,ibc,ihsdiag(nclo1p,isb,1),i1,i28,i38,i48,      11d15s22
     $        bc(ivtmp))                                                11d15s22
         do i=0,ncolh-1                                                 12d14s20
          do iv=0,nvirt(isbv)-1                                         12d14s20
           do ir=0,nrootu-1                                             12d14s20
            iad1=ivtmp+ir+nrootu*(iv+nvirt(isbv)*i)                     7d9s21
            iad2=ivst+iv+nvirt(isbv)*(ir+nrootu*i)                      12d14s20
            bc(iad2)=bc(iad1)                                           7d9s21
           end do                                                       12d14s20
          end do                                                        12d14s20
         end do                                                         12d14s20
         ibcoff=ivtmp                                                   7d9s21
         nn=ncsfs(iarg)*nrootu                                          10d26s21
         nnm=nn-1                                                       11d26s20
         jff1=nff1(nclo1p,isb,2)                                        11d25s20
         ifcol=1                                                        12d7s20
         jvst=ivst                                                      7d9s21
         do if1=1,nff1(nclo1p,isb,1)                                    11d26s20
          itcol=ifcol+ncsfs(iarg)-1                                     10d26s21
          if(ifcol.le.ih.and.itcol.ge.il)then                           12d7s20
           nnlow=max(0,il-ifcol)                                        12d7s20
           nnhgh=min(ncsfs(iarg)-1,ih-ifcol)                            10d26s21
           nnuse=nnhgh+1-nnlow                                          12d7s20
           nnuser=nnuse*nrootu                                          12d11s20
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
           do nclop=max(mdon+1,nclo1p-2),min(mdoobp,nclo1p+2)           10d26s21
            nclo=nclop-1                                                  11d25s20
            nopen=nec-nclo*2                                             11d25s20
            jarg=nclop-mdon                                               11d25s20
            call wfetch(nopen,mdon,idum,i2sb,i2smb,ndetb,ncsfb,ivecb,   10d26s21
     $           iaorbb,idum,bc,ibc)                                    11d9s22
            jff0=nff0(nclop,2)                                           11d25s20
            jvs=nff0(nclop,3)                                            11d26s20
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
                 if((btest(i0o,i).and..not.btest(i1c,i)).or.            2d6s23
     $             (btest(i0c,i).and.btest(i1o,i)))then                 2d6s23
                  nab4(1,1)=i                                             2d6s23
                 else                                                     2d6s23
                  nab4(2,1)=i                                             2d6s23
                 end if                                                   2d6s23
                end if                                                    2d6s23
               end do                                                     2d6s23
               call gandcr(i0c,i0o,i1c,i1o,nopen,nopen1,norbx,nnot1,    2d7s23
     $             nab1,icode1,imap1,nx1,irw1,irw2,iwpb,iwpk,bc,ibc)    2d7s23
               itmpp=ibcoff                                              10d26s21
               ibcoff=itmpp+nrootu*ncsfb                                 10d26s21
               call enough('hcisbk4.  3',bc,ibc)
               do i=itmpp,ibcoff-1                                       3d22s21
                bc(i)=0d0                                                3d22s21
               end do                                                    3d22s21
               iint=ibcoff                                               3d22s21
               ibcoff=iint+nvirt(isbv)                                   3d22s21
               call enough('hcisbk4.  4',bc,ibc)
               call gencup(i2sb,i2smb,i2sk,i2smk,nopen,nopen1,nab1,iwpb, 10d28s21
     $            iwpk,ioutg,imatg,ntypeg,mcsf,dum,1,0,bc,ibc)          11d14s22
               iprod=ibcoff                                              10d28s21
               ibcoff=iprod+ncsfb*ncsfs(iarg)                            10d28s21
               call enough('hcisbk4.  5',bc,ibc)
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
               do ii=0,ntypeg-1                                          10d26s21
                ipackc=ibc(ioutg+ii)                                      9d10s21
                if(ipackc1(2).gt.0)then                                   9d10s21
                 imy(4)=0                                                 10d12s21
                else
                 imy(4)=1                                                 10d12s21
                end if
                if(ipackc1(1).gt.0)then                                   9d10s21
                 lsk=ism(ipackc1(1))                                      10d12s21
                 igya(3)=irel(ipackc1(1))-1                                   10d12s21
                 imy(3)=0                                                 10d12s21
                else
                 lsk=ism(-ipackc1(1))                                     10d12s21
                 igya(3)=irel(-ipackc1(1))-1                              10d12s21
                 imy(3)=1                                                 10d12s21
                end if
                if(lsk.ne.isbo)then                                       7d12s21
                 write(6,*)('lsk vs. isbo: '),lsk,isbo
                 stop 'hcisbk4'
                end if                                                    7d12s21
                phs=1d0                                                  10d26s21
                if(ih0n(isbo).gt.0)then                                  2d24s22
                 if(imy(3).ne.0.and.isopt(4).ne.0)phs=-phs               10d26s21
                 ihq=ih0n(isbo)+igya(3)+nh0(isbo)*irefo(isbv)            10d26s21
                 do iv=0,nvirt(isbv)-1
                  bc(iint+iv)=bc(ihq+iv*nh0(isbo))*phs                   10d26s21
                 end do                                                  10d26s21
                else if(ih0n(isbv).gt.0)then                             2d24s22
                 if(imy(4).ne.0.and.isopt(4).ne.0)phs=-1d0               10d28s21
                 if(isopt(2).ne.0)phs=-phs                               10d26s21
                 ihq=ih0n(isbv)+irefo(isbv)+nh0(isbv)*igya(3)            10d28s21
                 do iv=0,nvirt(isbv)-1
                  bc(iint+iv)=bc(ihq+iv)*phs                             10d26s21
                 end do                                                  10d26s21
                else                                                     10d26s21
                 do iv=0,nvirt(isbv)-1                                   10d26s21
                  bc(iint+iv)=0d0                                        10d26s21
                 end do                                                  10d26s21
                end if                                                   10d26s21
                do i=1,norb                                               9d10s21
                 if(btest(i1c,i).and.btest(i0c,i).and.l2e.eq.0)then      2d17s22
                  js=ism(i)                                               9d10s21
                  jg=irel(i)-1                                            9d10s21
                  do jpass=0,1                                            9d10s21
c
c     integral is (jg,jg|imy(3) v)=+/- (jg,jg|v imy(3)
c     onex: b,a,d,c: here imy(3) is c and imy(4) is d, and b=a=jg.
c
                   ltest=jpass+2*(jpass+2*(imy(4)+2*imy(3)))              10d20s21
                   itestp=ltest+1                                         10d12s21
                   jtest=1-jpass+2*(1-jpass+2*(1-imy(4)+2*(1-imy(3))))    10d20s21
                   jtestp=jtest+1                                         10d12s21
                   iuse=-1                                                10d22s21
                   phs=1d0                                               10d20s21
                   if(isopt(2).ne.0)phs=-phs                             10d28s21
                   if(iifmx(itestp).ge.0)then                             10d20s21
                    iuse=iifmx(itestp)                                    10d20s21
                   else if(iifmx(jtestp).ge.0)then                        10d20s21
                    iuse=iifmx(jtestp)                                    10d20s21
                    if(isopt(4).ne.0)phs=-phs                             10d20s21
                   end if                                                 10d20s21
                   if(iuse.ge.0)then                                      10d20s21
                    iicolj=ionexr(js,js,isbo)+jg+irefo(js)*(jg           10d26s21
     $                +irefo(js)*(igya(3)+irefo(isbo)*nvirt(isbv)*iuse))10d26s21
                    nbad=irefo(js)*irefo(js)*irefo(isbo)                 10d26s21
                    do iv=0,nvirt(isbv)-1                                10d26s21
                     bc(iint+iv)=bc(iint+iv)+phs*bc(iicolj+iv*nbad)      10d26s21
                    end do                                               10d26s21
                   end if                                                 10d20s21
c
c     integral is (jg,v|imy(3)jg)=+/-(jg,imy(3)|v,jg)
c
                   ltest=jpass+2*(imy(3)+2*(imy(4)+2*jpass))              10d22s21
                   itestp=ltest+1                                         10d12s21
                   jtest=1-jpass+2*(1-imy(3)+2*(1-imy(4)+2*(1-jpass)))    10d22s21
                   jtestp=jtest+1                                         10d12s21
                   iuse=-1                                               11d12s21
                   phs=1d0                                               10d20s21
                   if(isopt(2).ne.0)phs=-phs                             10d28s21
                   if(iifmx(itestp).ge.0)then                             10d20s21
                    iuse=iifmx(itestp)                                    10d20s21
                   else if(iifmx(jtestp).ge.0)then                        10d20s21
                    iuse=iifmx(jtestp)                                    10d20s21
                    if(isopt(4).ne.0)phs=-phs                             10d20s21
                   end if                                                 10d20s21
                   if(iuse.ge.0)then                                      10d20s21
                    iicolk=ionexr(isbo,js,js)+igya(3)+irefo(isbo)*(jg    10d28s21
     $                 +irefo(js)*(jg+irefo(js)*nvirt(isbv)*iuse))      10d26s21
                    nbad=irefo(isbo)*irefo(js)*irefo(js)                 10d26s21
                    do iv=0,nvirt(isbv)-1                                10d26s21
                     bc(iint+iv)=bc(iint+iv)-phs*bc(iicolk+iv*nbad)      10d26s21
                    end do                                               10d26s21
                   end if                                                 10d20s21
                  end do                                                  10d22s21
                 end if                                                   10d22s21
                end do                                                    10d22s21
                igya(2)=igya(3)                                          10d28s21
                if(ipackc1(1).gt.0)then                                  10d28s21
                 im=ipackc1(1)                                           10d28s21
                 imy(2)=1                                                 10d13s21
                else                                                      10d12s21
                 im=-ipackc1(1)                                          10d28s21
                 imy(2)=0                                                 10d13s21
                end if                                                    10d12s21
c
c     isbo*isbo*isbo*isbv=isbo*isbv=isbv*isymc*isbv=isymc,
c     so symmetry is always correct.
c
                if(btest(i0c,im).and.l2e.eq.0)then                       2d17s22
c
c     "x" term.                                                         10d28s21
c                                            a     b    c    d          10d28s21
c     integral is (v imy(3)|imy(2)imy(2))=(imy(2)imy(2)|v imy(3)        10d28s21
c     onex: b,a,d,c: here imy(3) is c and imy(4) is d, and b=a=imy(2)
c
                 ltest=imy(2)+2*(imy(2)+2*(imy(4)+2*imy(3)))             10d28s21
                 itestp=ltest+1                                           10d12s21
                 jtest=1-imy(2)+2*(1-imy(2)+2*(1-imy(4)+2*(1-imy(3))))   10d28s21
                 jtestp=jtest+1                                           10d12s21
                 iuse=-1                                                  10d22s21
                 phs=1d0                                                 10d22s21
                 if(isopt(2).ne.0)phs=-phs                               10d28s21
                 if(iifmx(itestp).ge.0)then                               10d12s21
                  iuse=iifmx(itestp)                                      10d22s21
                 else if(iifmx(jtestp).ge.0)then                          10d12s21
                  iuse=iifmx(jtestp)                                      10d22s21
                  if(isopt(4).ne.0)phs=-phs                               10d22s21
                 end if                                                   10d12s21
                 if(iuse.ge.0)then                                        10d22s21
                  iicolj=ionexr(isbo,isbo,isbo)+igya(2)+irefo(isbo)*(
     $       igya(2)+irefo(isbo)*(igya(3)+irefo(isbo)*nvirt(isbv)*iuse))10d28s21
                  nbad=irefo(isbo)*irefo(isbo)*irefo(isbo)               10d26s21
                  do iv=0,nvirt(isbv)-1                                  10d26s21
                   bc(iint+iv)=bc(iint+iv)+phs*bc(iicolj+iv*nbad)        10d26s21
                  end do                                                 10d26s21
                 end if                                                   10d22s21
c                                             a      b    c   d         10d28s21
c     integral is (v,imy(2)|imy(2) imy(3))=(imy(2),imy(3)|v imy(2))     10d28s21
c     onex: b,a,d,c: here imy(3) is c and imy(2) is d=a, and b=imy(4)
c
                 ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(2)))             10d28s21
                 itestp=ltest+1                                           10d12s21
                 jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)+2*(1-imy(2))))   10d28s21
                 jtestp=jtest+1                                           10d12s21
                 iuse=-1                                                  10d22s21
                 phs=1d0                                                 10d22s21
                 if(isopt(2).ne.0)phs=-phs                               10d28s21
                 if(iifmx(itestp).ge.0)then                               10d12s21
                  iuse=iifmx(itestp)                                      10d22s21
                 else if(iifmx(jtestp).ge.0)then                          10d12s21
                  iuse=iifmx(jtestp)                                      10d22s21
                  if(isopt(4).ne.0)phs=-phs                               10d22s21
                 end if                                                   10d12s21
                 if(iuse.ge.0)then                                        10d22s21
                  iicolk=ionexr(isbo,isbo,isbo)+igya(2)+irefo(isbo)*(    10d26s21
     $       igya(3)+irefo(isbo)*(igya(2)+irefo(isbo)*nvirt(isbv)*iuse))10d28s21
                  nbad=irefo(isbo)*irefo(isbo)*irefo(isbo)               10d26s21
                  do iv=0,nvirt(isbv)-1                                  10d26s21
                   bc(iint+iv)=bc(iint+iv)-phs*bc(iicolk+iv*nbad)        10d26s21
                  end do                                                 10d26s21
                 end if                                                   10d22s21
                end if                                                    9d13s21
                imat=imatg+ncsfs(iarg)*ncsfb*ii                          10d28s21
                do ib=0,ncsfb-1                                          10d28s21
                 do ik=0,ncsfs(iarg)-1                                   10d28s21
                  ibk=imat+ib+ncsfb*ik                                   10d28s21
                  ikb=iprod+ik+ncsfs(iarg)*ib                            10d28s21
                  bc(ikb)=bc(ibk)                                        10d28s21
                 end do                                                  10d28s21
                end do                                                   10d28s21
                jprod=iprod+nnlow                                        10d28s21
                if(nvirt(isbv).gt.nnuse)then                              3d22s21
                 nrow=nvirt(isbv)*nrootu                                  3d22s21
                 ivdprod=ibcoff                                           3d22s21
                 ibcoff=ivdprod+nrow*ncsfb                                10d26s21
                 call enough('hcisbk4.  6',bc,ibc)
                 call dgemm('n','n',nrow,ncsfb,nnuse,1d0,                 10d26s21
     $              bc(jvst),nrow,bc(jprod),ncsfs(iarg),0d0,            10d26s21
     $              bc(ivdprod),nrow,                                   3d22s21
     d' hcisbk4.  1')
                 jvdprod=ivdprod                                          3d22s21
                 do i=0,nrootu*ncsfb-1                                    10d26s21
                  do iv=0,nvirt(isbv)-1                                   3d22s21
                   bc(itmpp+i)=bc(itmpp+i)+bc(jvdprod+iv)*bc(iint+iv)     3d22s21
                  end do                                                  3d22s21
                  jvdprod=jvdprod+nvirt(isbv)                             3d22s21
                 end do                                                   3d22s21
                 ibcoff=ivdprod                                           3d22s21
                else                                                      3d22s21
                 iivprod=ibcoff                                           3d22s21
                 ibcoff=iivprod+nnuser                                    3d22s21
                 call enough('hcisbk4.  7',bc,ibc)
                 do i=iivprod,ibcoff-1                                    3d22s21
                  bc(i)=0d0                                               3d22s21
                 end do                                                   3d22s21
                 iix=jvst                                                 7d9s21
                 do i=0,nnuser-1                                          3d22s21
                  do iv=0,nvirt(isbv)-1                                   3d22s21
                   bc(iivprod+i)=bc(iivprod+i)+bc(iint+iv)*bc(iix+iv)     7d9s21
                  end do                                                  3d22s21
                  iix=iix+nvirt(isbv)                                     3d22s21
                 end do                                                   3d22s21
                 call dgemm('n','n',nrootu,ncsfb,nnuse,1d0,               10d26s21
     $               bc(iivprod),nrootu,bc(jprod),ncsfs(iarg),1d0,       10d26s21
     $               bc(itmpp),nrootu,                                   3d22s21
     d' hcisbk4.  2')
                 ibcoff=iivprod
                end if                                                    3d22s21
               end do                                                     10d22s21
               ibcoff=ioutg                                               10d22s21
               do i=1,nok                                                    11d13s20
                if(itest(i,3).eq.1.and.l2e.eq.0)then                     2d17s22
                 itesta=i0c                                              10d28s21
                 itestb=i0o                                              10d28s21
                 nopenk=nopen                                            10d28s21
c
c     anihilate common
c
                 if(btest(itesta,itest(i,4)))then                             11d13s20
                  itesta=ibclr(itesta,itest(i,4))                             11d13s20
                  itestb=ibset(itestb,itest(i,4))                             11d13s20
                  nopenk=nopenk+1                                             11d13s20
                 else                                                         11d13s20
                  itestb=ibclr(itestb,itest(i,4))                             11d13s20
                  nopenk=nopenk-1                                             11d13s20
                 end if                                                       11d13s20
c
c     create ket
c
                 if(btest(itestb,nab4(2,1)))then                          11d27s20
                  itesta=ibset(itesta,nab4(2,1))                          11d27s20
                  itestb=ibclr(itestb,nab4(2,1))                          11d27s20
                  nopenk=nopenk-1                                             11d13s20
                 else                                                         11d13s20
                  itestb=ibset(itestb,nab4(2,1))                          11d27s20
                  nopenk=nopenk+1                                             11d13s20
                 end if                                                       11d13s20
                 call gandcr(i0c,i0o,itesta,itestb,nopen,nopenk,         10d28s21
     $          norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,iwpb1,iwpk1,11d14s22
     $               bc,ibc)                                            11d14s22
c     nab2(1) is virt
                 call gandcr(itesta,itestb,i1c,i1o,nopenk,nopen1,        10d28s21
     $          norbx,nnot2,nab2,icode1,imap1,nx1,irw1,irw2,iwpb2,iwpk2,11d14s22
     $               bc,ibc)                                            11d14s22
                 if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                  iunit=ibcoff
                  ibcoff=iunit+ncsfs(iarg)*ncsfs(iarg)                   10d28s21
                  call enough('hcisbk4.  8',bc,ibc)
                  do iz=iunit,ibcoff-1
                   bc(iz)=0d0
                  end do
                  do i1i=0,ncsfs(iarg)-1                                 10d28s21
                   iad=iunit+i1i*(ncsfs(iarg)+1)                         10d28s21
                   bc(iad)=1d0                                           10d26s21
                  end do                                                 10d26s21
                  call spinloop(i2sb,i2smb,i2sk,i2smk,nopen,nopen1,      10d28s21
     $               nopenk,ncsfb,ncsfs(iarg),itype,imatx,ntypex,nab1,  10d24s21
     $               iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),ncsfs(iarg),10d28s21
     $               ieoro,bc,ibc)                                      11d14s22
                  iprod=ibcoff                                           10d28s21
                  ibcoff=iprod+ncsfs(iarg)*ncsfb                         10d28s21
                  call enough('hcisbk4.  9',bc,ibc)
                  do ini=0,ntypex-1                                      10d24s21
                   do iv=0,nvirt(isbv)-1                                 10d26s21
                    bc(iint+iv)=0d0                                      10d26s21
                   end do                                                10d26s21
                   ipack=ibc(itype+ini)
                   imat=imatx+ncsfs(iarg)*ncsfb*ini                      10d26s21
                   do ib=0,ncsfb-1                                       10d28s21
                    do ik=0,ncsfs(iarg)-1                                10d28s21
                     ikb=iprod+ik+ncsfs(iarg)*ib                         10d28s21
                     ibk=imat+ib+ncsfb*ik                                10d28s21
                     bc(ikb)=bc(ibk)                                     10d28s21
                    end do                                               10d28s21
                   end do                                                10d28s21
                   if(ipack2(2).eq.norbx.or.-ipack2(2).eq.norbx)then     10d28s21
                   else
                    write(6,*)
     $                   ('oh no, virt is not were I expected it!!')
                    write(6,*)ipack2
                    stop 'hcsibk4'
                   end if
                   isy(2)=isbv                                           10d28s21
                   if(ipack2(2).gt.0)then                                10d28s21
                    imy(2)=0                                             10d28s21
                   else                                                  10d28s21
                    imy(2)=1                                             10d28s21
                   end if                                                10d28s21
                   do j=1,4                                               9d9s21
                    if(j.ne.2)then                                       10d28s21
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
c                                     ab cd
c     integral is (1v|34)=(34|1v)=+/-(43|v1)
c     onex: b,a,d,c: here v(=3) is c, b=2, a=1, and d=4.
c
                   ltest=imy(4)+2*(imy(3)+2*(imy(2)+2*imy(1)))           10d28s21
                   itestp=ltest+1                                         10d13s21
                   jtest=1-imy(4)+2*(1-imy(3)+2*(1-imy(2)+2*(1-imy(1)))) 10d28s21
                   jtestp=jtest+1                                         10d13s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(isopt(2).ne.0)phs=-phs                             10d28s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolj=ionexr(isy(3),isy(4),isy(1))+igya(3)          10d28s21
     $                 +irefo(isy(3))*(igya(4)+irefo(isy(4))*(          10d28s21
     $               igya(1)+irefo(isy(1))*nvirt(isbv)*iuse))           10d28s21
                    nbad=irefo(isy(3))*irefo(isy(4))*irefo(isy(1))       10d28s21
                    do iv=0,nvirt(isbv)-1                                  10d26s21
                     bc(iint+iv)=bc(iint+iv)+phs*bc(iicolj+iv*nbad)        10d26s21
                    end do                                                 10d26s21
                   end if                                                   10d22s21
c
c                             ab cd
c     integral is (14|3v)=+/-(41|v3)
c     onex: b,a,d,c: here v(=3) is c, b=1, a=4, and d=2.
c
                   ltest=imy(4)+2*(imy(1)+2*(imy(2)+2*imy(3)))           10d28s21
                   itestp=ltest+1                                         10d13s21
                   jtest=1-imy(4)+2*(1-imy(1)+2*(1-imy(2)+2*(1-imy(3)))) 10d28s21
                   jtestp=jtest+1                                         10d13s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(isopt(2).ne.0)phs=-phs                             10d28s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolk=ionexr(isy(1),isy(4),isy(3))+igya(1)          10d28s21
     $                 +irefo(isy(1))*(igya(4)+irefo(isy(4))*(          10d28s21
     $               igya(3)+irefo(isy(3))*nvirt(isbv)*iuse))           10d28s21
                    nbad=irefo(isy(1))*irefo(isy(4))*irefo(isy(3))       10d28s21
                    do iv=0,nvirt(isbv)-1                                  10d26s21
                     bc(iint+iv)=bc(iint+iv)-phs*bc(iicolk+iv*nbad)        10d26s21
                    end do                                                 10d26s21
                   end if                                                   10d22s21
                   jprod=iprod+nnlow                                     10d28s21
                   if(nvirt(isbv).gt.nnuse)then                              3d22s21
                    nrow=nvirt(isbv)*nrootu                                  3d22s21
                    ivdprod=ibcoff                                           3d22s21
                    ibcoff=ivdprod+nrow*ncsfb                                10d26s21
                    call enough('hcisbk4. 10',bc,ibc)
                    call dgemm('n','n',nrow,ncsfb,nnuse,1d0,                 10d26s21
     $                  bc(jvst),nrow,bc(jprod),ncsfs(iarg),0d0,            10d26s21
     $                  bc(ivdprod),nrow,                                   3d22s21
     d' hcisbk4.  3')
                    jvdprod=ivdprod                                          3d22s21
                    do ii=0,nrootu*ncsfb-1                                    10d26s21
                     do iv=0,nvirt(isbv)-1                                   3d22s21
                      bc(itmpp+ii)=bc(itmpp+ii)                          10d26s21
     $                     +bc(jvdprod+iv)*bc(iint+iv)                   10d26s21
                     end do                                                  3d22s21
                     jvdprod=jvdprod+nvirt(isbv)                             3d22s21
                    end do                                                   3d22s21
                    ibcoff=ivdprod                                           3d22s21
                   else                                                      3d22s21
                    iivprod=ibcoff                                           3d22s21
                    ibcoff=iivprod+nnuser                                    3d22s21
                    call enough('hcisbk4. 11',bc,ibc)
                    do iz=iivprod,ibcoff-1                                    3d22s21
                     bc(iz)=0d0                                               3d22s21
                    end do                                                   3d22s21
                    iix=jvst                                                 7d9s21
                    do ii=0,nnuser-1                                          3d22s21
                     do iv=0,nvirt(isbv)-1                                   3d22s21
                      bc(iivprod+ii)=bc(iivprod+ii)                      10d26s21
     $                     +bc(iint+iv)*bc(iix+iv)                       10d26s21
                     end do                                                  3d22s21
                     iix=iix+nvirt(isbv)                                     3d22s21
                    end do                                                   3d22s21
                    call dgemm('n','n',nrootu,ncsfb,nnuse,1d0,               10d26s21
     $               bc(iivprod),nrootu,bc(jprod),ncsfs(iarg),1d0,       10d26s21
     $               bc(itmpp),nrootu,                                   3d22s21
     d' hcisbk4.  4')
                    ibcoff=iivprod
                   end if                                                    3d22s21
                  end do                                                 10d22s21
                  ibcoff=itype                                           10d22s21
                 end if                                                    11d13s20
                end if                                                   10d22s21
               end do                                                        11d13s20
               do ir=0,nrootm                                            3d22s21
                iadg=ig+jvs-1+ncsfv*ir                                   3d22s21
                do j=0,ncsfb-1                                           10d26s21
                 jad=itmpp+ir+nrootu*j                                   3d22s21
                 bc(iadg+j)=bc(iadg+j)+bc(jad)                           3d22s21
                end do                                                   3d22s21
               end do                                                    3d22s21
               ibcoff=itmpp                                              10d26s21
              else if(l2e.eq.0)then                                     2d7s23
               nnot=0                                                   2d6s23
               if(ndifs.eq.4.and.ndifb.eq.4)then                        2d6s23
                nnot=4                                                  2d6s23
                ioxx(1)=1                                               2d6s23
                ioxx(2)=1                                               2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if((btest(i1c,i).and.btest(i0o,i)).or.                2d6s23
     $                   (btest(i1o,i).and..not.btest(i0c,i)))then      2d6s23
                   nab4(2,ioxx(2))=i                                    2d6s23
                   ioxx(2)=ioxx(2)+1                                    2d6s23
                  else                                                  2d6s23
                   nab4(1,ioxx(1))=i                                    2d6s23
                   ioxx(1)=ioxx(1)+1                                    2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
               else if(ndifb.eq.3)then                                  2d6s23
                nnot=3                                                  2d6s23
                ioxx(1)=1                                               2d6s23
                ioxx(2)=1                                               2d6s23
                iswap=0                                                 2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if(btest(gandcc,i).and.                               2d6s23
     $                   ((btest(i0c,i).and..not.btest(i1o,i)).or.      2d6s23
     $                   (btest(i1c,i).and..not.btest(i0o,i))))then     2d6s23
                   if(btest(i1c,i))iswap=1                              2d6s23
                   nab4(1,1)=i                                          2d6s23
                   nab4(1,2)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,ioxx(2))=i                                    2d6s23
                   ioxx(2)=ioxx(2)+1                                    2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
                if(iswap.ne.0)then                                      2d6s23
                 icpy=nab4(1,1)                                         2d6s23
                 nab4(1,1)=nab4(2,1)                                    2d6s23
                 nab4(2,1)=icpy                                         2d6s23
                 icpy=nab4(1,2)                                         2d6s23
                 nab4(1,2)=nab4(2,2)                                    2d6s23
                 nab4(2,2)=icpy                                         2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(i0c,nab4(1,2)).and.                           2d6s23
     $                   .not.btest(i0c,nab4(1,1)))nbt=1                2d6s23
                else                                                    2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(i1c,nab4(2,2)).and.                           2d6s23
     $                   .not.btest(i1c,nab4(2,1)))nbt=1                2d6s23
                end if                                                  2d6s23
                if(nbt.ne.0)then                                        2d6s23
                 nab4(1,1)=nab4(1,2)                                    2d6s23
                 nab4(2,1)=nab4(2,2)                                    2d6s23
                end if                                                  2d6s23
               else if(ndifs.eq.0.and.ndifd.eq.2)then                   2d6s23
                nnot=3                                                  2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if(btest(i0c,i))then                                  2d6s23
                   nab4(1,1)=i                                          2d6s23
                   nab4(1,2)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,1)=i                                          2d6s23
                   nab4(2,2)=i                                          2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
               end if                                                   2d6s23
               if(nnot.ne.0)then                                        2d7s23
                iint=ibcoff                                               10d26s21
                itmpp=iint+nvirt(isbv)                                    10d26s21
                ibcoff=itmpp+nrootu*ncsfb                                 10d26s21
                call enough('hcisbk4. 12',bc,ibc)
                do iz=iint,ibcoff-1                                       10d26s21
                 bc(iz)=0d0                                               10d26s21
                end do                                                    10d26s21
                iu1=1                                                     1d25s21
                iu2=1                                                     1d25s21
c     sing is bra
                itesta=i0c                                                10d28s21
                itestb=i0o                                                10d28s21
                if(btest(itesta,nab4(1,iu1)))then                         1d25s21
                 itesta=ibclr(itesta,nab4(1,iu1))                         1d25s21
                 itestb=ibset(itestb,nab4(1,iu1))                         1d25s21
                 nopenk=nopen+1                                           10d28s21
                else if(btest(itestb,nab4(1,iu1)))then                    1d25s21
                 itestb=ibclr(itestb,nab4(1,iu1))                         1d25s21
                 nopenk=nopen-1                                           10d28s21
                end if                                                        11d13s20
                if(btest(itestb,nab4(2,iu2)))then                         1d25s21
                 itesta=ibset(itesta,nab4(2,iu2))                         1d25s21
                 itestb=ibclr(itestb,nab4(2,iu2))                         1d25s21
                 nopenk=nopenk-1                                              11d13s20
                else                                                          11d13s20
                 itestb=ibset(itestb,nab4(2,iu2))                         1d25s21
                 nopenk=nopenk+1                                              11d13s20
                end if                                                        11d13s20
                call gandcr(i0c,i0o,itesta,itestb,nopen,nopenk,           10d28s21
     $         norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,iwpb1,iwpk1, 11d14s22
     $             bc,ibc)                                              11d14s22
                call gandcr(itesta,itestb,i1c,i1o,nopenk,nopen1,          10d28s21
     $         norbx,nnot2,nab2,icode1,imap1,nx1,irw1,irw2,iwpb2,iwpk2, 11d14s22
     $             bc,ibc)                                              11d14s22
                if(nnot1.eq.2.and.nnot2.eq.2)then                         1d25s21
                 iunit=ibcoff
                 ibcoff=iunit+ncsfs(iarg)*ncsfs(iarg)                   10d28s21
                 call enough('hcisbk4. 13',bc,ibc)
                 do iz=iunit,ibcoff-1
                  bc(iz)=0d0
                 end do
                 do i1i=0,ncsfs(iarg)-1                                 10d28s21
                  iad=iunit+i1i*(ncsfs(iarg)+1)                         10d28s21
                  bc(iad)=1d0                                           10d26s21
                 end do                                                 10d26s21
                 call spinloop(i2sb,i2smb,i2sk,i2smk,nopen,nopen1,        10d28s21
     $               nopenk,ncsfb,ncsfs(iarg),itype,imatx,ntypex,nab1,  10d28s21
     $             iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),ncsfs(iarg),  10d28s21
     $                ieoro,bc,ibc)                                     11d14s22
                 iprod=ibcoff                                             10d28s21
                 ibcoff=iprod+ncsfb*ncsfs(iarg)                           10d28s21
                 call enough('hcisbk4. 14',bc,ibc)
                 do i=0,ntypex-1                                          10d24s21
                  ipack=ibc(itype+i)
                  do iz=0,nvirt(isbv)-1                                   10d26s21
                   bc(iint+iz)=0d0                                        10d26s21
                  end do                                                  10d26s21
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
                  if(ipack2(2).eq.norbx.or.-ipack2(2).eq.norbx)then       10d28s21
c
c                                                ab cd
c     if v is 2, integral is (1v|34)=(34|1v)=+/-(43|v1)
c     onex: b,a,d,c: here v(=1) is c, b=4, a=3, and d=2.
c
                   ltest=imy(4)+2*(imy(3)+2*(imy(2)+2*imy(1)))            10d28s21
                   itestp=ltest+1                                         10d22s21
                   jtest=1-imy(4)+2*((1-imy(3))+2*((1-imy(2))             10d28s21
     $                  +2*(1-imy(1))))                                  10d28s21
                   jtestp=jtest+1                                         10d22s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(isopt(2).ne.0)phs=-phs                              10d28s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolj=ionexr(isy(3),isy(4),isy(1))+igya(3)           10d28s21
     $                  +irefo(isy(3))*(igya(4)+irefo(isy(4))*(           10d28s21
     $               igya(1)+irefo(isy(1))*nvirt(isbv)*iuse))           10d28s21
                    nbad=irefo(isy(4))*irefo(isy(3))*irefo(isy(1))        10d28s21
                    do iv=0,nvirt(isbv)-1                                  10d26s21
                     bc(iint+iv)=bc(iint+iv)+phs*bc(iicolj+iv*nbad)        10d26s21
                    end do                                                 10d26s21
                   end if                                                   10d22s21
                   if(nnot.eq.4)then                                      10d22s21
c
c                             ab cd
c     integral is (14|3v)=+/-(41|v3)
c     onex: b,a,d,c: here v(=1) is c, b=2, a=3, and d=4.
c
                    ltest=imy(4)+2*(imy(1)+2*(imy(2)+2*imy(3)))           10d28s21
                    itestp=ltest+1                                        10d22s21
                    jtest=1-imy(4)+2*((1-imy(1))+2*((1-imy(2))            10d28s21
     $                   +2*(1-imy(3))))                                 10d28s21
                    jtestp=jtest+1                                        10d22s21
                    iuse=-1                                                  10d22s21
                    phs=1d0                                                 10d22s21
                    if(isopt(2).ne.0)phs=-phs                             10d28s21
                    if(iifmx(itestp).ge.0)then                               10d12s21
                     iuse=iifmx(itestp)                                      10d22s21
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iuse=iifmx(jtestp)                                      10d22s21
                     if(isopt(4).ne.0)phs=-phs                               10d22s21
                    end if                                                   10d12s21
                    if(iuse.ge.0)then                                        10d22s21
                     iicolj=ionexr(isy(1),isy(4),isy(3))+igya(1)+         10d28s21
     $                 irefo(isy(1))*(igya(4)+irefo(isy(4))*(           10d28s21
     $               igya(3)+irefo(isy(3))*nvirt(isbv)*iuse))           10d28s21
                     nbad=irefo(isy(1))*irefo(isy(3))*irefo(isy(4))       10d28s21
                     do iv=0,nvirt(isbv)-1                                  10d26s21
                      bc(iint+iv)=bc(iint+iv)-phs*bc(iicolj+iv*nbad)        10d26s21
                     end do                                                 10d26s21
                    end if                                                   10d22s21
                   end if                                                 10d22s21
                  else if(ipack2(4).eq.norbx.or.-ipack2(4).eq.norbx)then  10d28s21
c
c                             ab cd
c     integral is (12|3v)=+/-(21|v3)
c     onex: b,a,d,c: here v(=3) is c, b=2, a=1, and d=4.
c
                   ltest=imy(2)+2*(imy(1)+2*(imy(4)+2*imy(3)))            10d28s21
                   itestp=ltest+1                                         10d22s21
                   jtest=1-imy(2)+2*((1-imy(1))+2*((1-imy(4))             10d28s21
     $                  +2*(1-imy(3))))                                  10d28s21
                   jtestp=jtest+1                                         10d22s21
                   iuse=-1                                                  10d22s21
                   phs=1d0                                                 10d22s21
                   if(isopt(2).ne.0)phs=-phs                              10d28s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iuse=iifmx(itestp)                                      10d22s21
                    jchoice=iuse
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iuse=iifmx(jtestp)                                      10d22s21
                    if(isopt(4).ne.0)phs=-phs                               10d22s21
                    jchoice=iuse+100
                   end if                                                   10d12s21
                   if(iuse.ge.0)then                                        10d22s21
                    iicolj=ionexr(isy(1),isy(2),isy(3))+igya(1)           10d28s21
     $                +irefo(isy(1))*(igya(2)+irefo(isy(2))*(           10d28s21
     $               igya(3)+irefo(isy(3))*nvirt(isbv)*iuse))           10d28s21
                    nbad=irefo(isy(2))*irefo(isy(1))*irefo(isy(3))        10d26s21
                    do iv=0,nvirt(isbv)-1                                  10d26s21
                     bc(iint+iv)=bc(iint+iv)+phs*bc(iicolj+iv*nbad)        10d26s21
                    end do                                                 10d26s21
                   end if                                                   10d22s21
                   if(nnot.eq.4)then                                      10d22s21
c
c                                     ab cd
c     integral is (1v|32)=(32|1v)=+/-(23|v1)
c     onex: b,a,d,c:
c
                    ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))           10d28s21
                    itestp=ltest+1                                        10d22s21
                    jtest=1-imy(2)+2*((1-imy(3))+2*((1-imy(4))            10d28s21
     $                  +2*(1-imy(1))))                                 10d28s21
                    jtestp=jtest+1                                        10d22s21
                    iuse=-1                                                  10d22s21
                    phs=1d0                                                 10d22s21
                    if(isopt(2).ne.0)phs=-phs                             10d28s21
                    if(iifmx(itestp).ge.0)then                               10d12s21
                     iuse=iifmx(itestp)                                      10d22s21
                     kchoice=iuse
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iuse=iifmx(jtestp)                                      10d22s21
                     if(isopt(4).ne.0)phs=-phs                               10d22s21
                     kchoice=iuse+100
                    end if                                                   10d12s21
                    if(iuse.ge.0)then                                        10d22s21
                     iicolj=ionexr(isy(3),isy(2),isy(1))+igya(3)          10d28s21
     $                 +irefo(isy(3))*(igya(2)+irefo(isy(2))*(          10d28s21
     $               igya(1)+irefo(isy(1))*nvirt(isbv)*iuse))           10d28s21
                     nbad=irefo(isy(3))*irefo(isy(1))*irefo(isy(2))       10d28s21
                     do iv=0,nvirt(isbv)-1                                  10d26s21
                      bc(iint+iv)=bc(iint+iv)-phs*bc(iicolj+iv*nbad)        10d26s21
                     end do                                                 10d26s21
                    end if                                                   10d22s21
                   end if                                                 10d22s21
                  else                                                    10d22s21
                   write(6,*)('oh no!! v is in wrong place!!')            10d22s21
                   stop 'hcisbk4'
                  end if                                                  10d22s21
                  imat=imatx+ncsfs(iarg)*ncsfb*i                          10d28s21
                  do ib=0,ncsfb-1                                         10d28s21
                   do ik=0,ncsfs(iarg)-1                                  10d28s21
                    ibk=imat+ib+ncsfb*ik                                  10d28s21
                    ikb=iprod+ik+ncsfs(iarg)*ib                           10d28s21
                    bc(ikb)=bc(ibk)                                       10d28s21
                   end do                                                 10d28s21
                  end do                                                  10d28s21
                  jprod=iprod+nnlow                                       10d28s21
                  if(nvirt(isbv).gt.nnuse)then                              3d22s21
                   nrow=nvirt(isbv)*nrootu                                  3d22s21
                   ivdprod=ibcoff                                           3d22s21
                   ibcoff=ivdprod+nrow*ncsfb                                10d26s21
                   call enough('hcisbk4. 15',bc,ibc)
                   call dgemm('n','n',nrow,ncsfb,nnuse,1d0,                 10d26s21
     $               bc(jvst),nrow,bc(jprod),ncsfs(iarg),0d0,            10d26s21
     $               bc(ivdprod),nrow,                                   3d22s21
     d' hcisbk4.  5')
                   jvdprod=ivdprod                                          3d22s21
                   do ii=0,nrootu*ncsfb-1                                    10d26s21
                    do iv=0,nvirt(isbv)-1                                   3d22s21
                     bc(itmpp+ii)=bc(itmpp+ii)                          2d7s23
     $                    +bc(jvdprod+iv)*bc(iint+iv)                   2d7s23
                    end do                                                  3d22s21
                    jvdprod=jvdprod+nvirt(isbv)                             3d22s21
                   end do                                                   3d22s21
                   ibcoff=ivdprod                                           3d22s21
                  else                                                      3d22s21
                   iivprod=ibcoff                                           3d22s21
                   ibcoff=iivprod+nnuser                                    3d22s21
                   call enough('hcisbk4. 16',bc,ibc)
                   do ii=iivprod,ibcoff-1                                    3d22s21
                    bc(ii)=0d0                                               3d22s21
                   end do                                                   3d22s21
                   iix=jvst                                                 7d9s21
                   do ii=0,nnuser-1                                          3d22s21
                    do iv=0,nvirt(isbv)-1                                   3d22s21
                     bc(iivprod+ii)=bc(iivprod+ii)                      2d7s23
     $                    +bc(iint+iv)*bc(iix+iv)                       2d7s23
                    end do                                                  3d22s21
                    iix=iix+nvirt(isbv)                                     3d22s21
                   end do                                                   3d22s21
                   call dgemm('n','n',nrootu,ncsfb,nnuse,1d0,               10d26s21
     $               bc(iivprod),nrootu,bc(jprod),ncsfs(iarg),1d0,       10d26s21
     $               bc(itmpp),nrootu,                                   3d22s21
     d' hcisbk4.  6')
                   ibcoff=iivprod
                  end if                                                    3d22s21
                 end do                                                   10d22s21
                 ibcoff=itype                                             10d22s21
                 do ir=0,nrootm                                            3d22s21
                  iadg=ig+jvs-1+ncsfv*ir                                   3d22s21
                  do j=0,ncsfb-1                                           10d26s21
                   jad=itmpp+ir+nrootu*j                                   3d22s21
                   bc(iadg+j)=bc(iadg+j)+bc(jad)                           3d22s21
                  end do                                                   3d22s21
                 end do                                                    3d22s21
                 ibcoff=iint                                              10d26s21
                end if                                                    10d22s21
               end if                                                   2d7s23
              end if                                                      11d25s20
             end if                                                     2d7s23
             jvs=jvs+ncsfb                                              10d26s21
            end do                                                       11d25s20
           end do                                                        11d25s20
           jvst=jvst+nvirt(isbv)*nnuser                                 7d9s21
          else                                                          12d8s20
           jff1=jff1+2                                                  12d8s20
          end if                                                        12d7s20
          ifcol=itcol+1                                                 12d7s20
         end do                                                         11d25s20
         ibcoff=ivst                                                    3d22s21
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      if(lprt)then
       write(6,*)('all done in hcisbk')
       write(6,*)('in hcisbk, my part of gi ')
       call prntm2(bc(ig),ncsfv,nrootu,ncsfv)
      end if
      jg=ig-1                                                           7d9s21
      do ir=1,nrootu                                                    7d9s21
       do i=1,ncsfv                                                     7d9s21
        gi(i,ir)=gi(i,ir)+bc(jg+i)                                      7d9s21
       end do                                                           7d9s21
       jg=jg+ncsfv                                                      7d9s21
      end do                                                            7d9s21
      ibcoff=ig                                                         12d12s20
      return
      end
