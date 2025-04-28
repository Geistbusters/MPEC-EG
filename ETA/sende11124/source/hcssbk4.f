c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcssbk4(ihsdiagb,nff1b,iff1b,i2sb,i2smb,mdoobp,ncsfb,  11d5s21
     $     isymbra,ihsdiagk,nff1k,iff1k,i2sk,i2smk,mdookp,ncsfk,isymket,11d5s21
     $     nec,mdon,nsymb,multh,irw0,irw1,irw2,nvirt,nrootu,lprt,isymc, 11d5s21
     $     irorip,isopt,ism,irel,irefo,norb,ih0n,nh0,i4or,jmatsr,kmatsr,11d5s21
     $     iifmx,ntype,maxbx,l2e,nwiacc,bc,ibc)                         1d26s23
      implicit real*8 (a-h,o-z)
      external second                                                   2d18s21
      integer*8 ihsdiagb(mdoobp,nsymb),i18,i28,i38,i48,i1c,i1o,j1c,j1o,  8d21s21
     $     itestc,itesto,last8(2),ihsdiagk(mdookp,nsymb,2),ipackc,ipack,2d7s23
     $     gandcc,gandco,gandcb                                         2d7s23
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2),ipackc1(8),icode1(64),11d5s21
     $     imap1(64),icode2(64),imap2(64)                               1d18s23
      integer*2 ipack2(4)                                               11d5s21
      logical lkeep,lprt,lchoice,lkey                                        3d17s21
      equivalence (ipackc1,ipackc),(ipack2,ipack)                       11d5s21
      dimension nff1b(mdoobp,nsymb,2),ncsfb(*),ncsfk(*),multh(8,8),
     $     iff1k(*),nff1k(mdookp,nsymb,2),nvirt(*),ism(*),irel(*),      11d5s21
     $     itest(32,4),isorb1(32),idorb1(32),jsorb1(32),jdorb1(32),     11d5s21
     $     nab4(2,3),idenj(8),idenk(8),irefo(*),iff1b(*),ndenj(8),      11d5s21
     $     ndenk(8),idenjv(8),idenkv(8),ndenjv(8),ndenkv(8),            11d5s21
     $     ivec(5,8),ivecp(5,8),nvecp(8),ircv(5,8),nrcv(5,8),           2d12s21
     $     itransv(5,8),nh0(*),ih0n(*),i4or(8,8,*),jmatsr(8,8,*),       11d5s21
     $     kmatsr(8,8,*),isopt(*),iifmx(*),iwpb1(4),isy(4),ioxx(2),     2d7s23
     $     iwpk1(4),iwpb2(4),iwpk2(4),mcsf(2),imy(4),igya(4),ipack4(4)  11d7s21
      data loopx/2434540/
      include "common.store"                                            12d12s20
      data icall/0/
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/paddcm/npadddi                                             6d7s22
      save icall
      icall=icall+1
      irori=irorip-1                                                    11d5s21
      if(isopt(3).ne.0)then                                             10d28s21
       ieoro=1-irori                                                    10d28s21
      else                                                              10d28s21
       ieoro=irori                                                      10d28s21
      end if                                                            10d28s21
      if(lprt)                                                          1d25s21
     $    write(6,*)('hello! my name is hcssbk4'),isymbra,isymket,
     $ nec,mdon,mdoobp,mdookp,nsymb,ixw1,ixw2,norb,maxbx,ibcoff
      loop=0
c
c     maxbx is for input wavefcns ...
c     but here bra has ket roots, so maxbx may not be large enough.     6d7s22
c
      maxbxb=0                                                          6d7s22
      do isb=1,nsymb                                                    6d7s22
       isbv=multh(isb,isymbra)                                          6d7s22
       do ii=mdon+1,mdoobp                                              6d30s23
        iarg=ii-mdon                                                    6d7s22
        maxbxb=                                                         6d8s22
     $       max(maxbxb,ncsfb(iarg)*nff1b(ii,isb,1)*nvirt(isbv)*nrootu) 6d8s22
       end do                                                           6d7s22
      end do                                                            6d7s22
      maxbxb=maxbxb+npadddi                                             6d7s22
      iacc=ibcoff                                                       2d12s21
      ibcoff=iacc+mynprocg                                              2d12s21
      do isb=1,nsymb                                                    2d12s21
       nvecp(isb)=0                                                     2d12s21
       do i=1,5                                                         2d12s21
        ivec(i,isb)=ibcoff                                              2d12s21
        ircv(i,isb)=ivec(i,isb)+maxbx                                   2d12s21
        ibcoff=ircv(i,isb)+mynprocg                                     2d12s21
       end do                                                           2d12s21
      end do                                                            2d12s21
      last8(1)=-1                                                       2d8s21
      igg=ibcoff                                                        2d12s21
      ibcoff=igg+maxbxb                                                 6d7s22
      nacc=0                                                            1d30s21
      itransgg=1                                                        2d24s22
      call enough('hcssbk4.  1',bc,ibc)
      do iz=igg,ibcoff-1                                                2d24s22
       bc(iz)=0d0                                                       2d24s22
      end do                                                            2d24s22
      nrootm=nrootu-1                                                   12d15s20
      nsing=0                                                           12d15s20
      loop=0
      norbx=norb+1                                                      12d14s20
      norbxx=norbx+1                                                    12d14s20
      do isb=1,nsymb                                                    12d12s20
       isbv=multh(isb,isymbra)                                          7d9s21
       if(lprt)write(6,*)('for isb = '),isb,isbv,loop
       do jsb=1,nsymb                                                   2d12s21
        nvecp(jsb)=0                                                    2d12s21
       end do                                                           2d12s21
       do ncloip=mdon+1,mdoobp                                          11d5s21
        if(min(nff1b(ncloip,isb,1),nvirt(isbv)).gt.0)then               7d9s21
         if(lprt)write(6,*)('for ncloip = '),ncloip,loop
         ncloi=ncloip-1                                                  12d12s20
         iarg=ncloip-mdon                                                12d12s20
         nsing=nsing+nff1b(ncloip,isb,1)*ncsfb(iarg)*nvirt(isbv)        11d5s21
         ncolti=nff1b(ncloip,isb,1)*ncsfb(iarg)*nvirt(isbv)*nrootu      11d5s21
         ngg=nvirt(isbv)*nrootu                                         12d15s20
         nnc=nff1b(ncloip,isb,1)*ncsfb(iarg)                            11d5s21
         nopeni=nec-2*ncloi                                              12d12s20
         do jsb=1,nsymb                                                   12d12s20
          jsbv=multh(jsb,isymket)                                       7d9s21
          if(lprt)write(6,*)('jsb '),jsb,jsbv,loop
          i28=nvirt(jsbv)*nrootu                                         12d14s20
          call ilimts(nvirt(isbv),nvirt(jsbv),mynprocg,mynowprog,il,ih, 12d14s20
     $         i1s,i1e,i2s,i2e)                                         12d14s20
          nhere=ih+1-il                                                 12d15s20
          nherev=i2e+1-i2s                                              12d14s20
          if(jsbv.eq.isbv)then                                          7d19s21
           if(mynowprog.lt.mynprocg-1)then                              12d15s20
            ip=mynowprog+1                                              12d15s20
            call ilimts(nvirt(isbv),nvirt(jsbv),mynprocg,ip,jl,jh,      12d15s20
     $         j1s,j1e,j2s,j2e)                                         12d15s20
            i2top=j2s-1                                                 12d15s20
           else                                                         12d15s20
            i2top=i2e                                                   12d15s20
           end if                                                       12d15s20
          end if                                                        12d15s20
          do nclojp=max(mdon+1,ncloip-2),min(mdookp,ncloip+2)           11d5s21
           if(min(nff1k(nclojp,jsb,1),nherev).gt.0)then                 7d9s21
            if(lprt)write(6,*)('for nclojp '),nclojp,loop
            ncloj=nclojp-1                                               12d12s20
            jarg=nclojp-mdon                                             12d12s20
            nopenj=nec-2*ncloj                                           12d14s20
            nrow=nherev*nrootu                                          12d14s20
            ncoltj=nff1k(nclojp,jsb,1)*ncsfk(jarg)*nrow                 11d5s21
            ivecbc=ibcoff                                               2d12s21
            i48=nff1k(nclojp,jsb,1)*ncsfk(jarg)                         11d5s21
            ncol=i48                                                    2d12s21
c
c     has vector alread been loaded?
c
            nminc=mdookp                                                11d5s21
            do i=1,nvecp(jsb)                                           2d12s21
             if(nclojp.eq.ivecp(i,jsb))then                             2d12s21
              nuse=i                                                    2d12s21
              go to 22                                                  2d12s21
             end if                                                     2d12s21
             if(ivecp(i,jsb).lt.nminc)then                              2d12s21
              nminc=ivecp(i,jsb)                                        2d12s21
              iminc=i                                                   2d12s21
             end if                                                     2d12s21
            end do                                                      2d12s21
c
c     no ... are we building up library?
c
            if(nvecp(jsb).lt.5)then                                     2d12s21
c     build
             nuse=nvecp(jsb)+1                                          2d12s21
             nvecp(jsb)=nvecp(jsb)+1                                    2d12s21
             ivecp(nvecp(jsb),jsb)=nclojp                               2d12s21
            else
c
c     smallest one out
c
             nuse=iminc                                                 2d12s21
             ivecp(iminc,jsb)=nclojp                                    2d12s21
            end if
            i18=1+nrootu*(i2s-1)                                        12d14s20
            i28=nrootu*i2e                                              12d14s20
            i38=1                                                       12d14s20
            itransv(nuse,jsb)=0                                         2d12s21
            call ddi_iget(bc,ibc,ihsdiagk(nclojp,jsb,1),i18,i28,i38,i48,11d15s22
     $           bc(ivec(nuse,jsb)),ibc(ircv(nuse,jsb)),nrcv(nuse,jsb)) 2d12s21
   22       continue                                                    2d12s21
            nrowi=nrow*ncsfb(iarg)                                      11d5s21
            nrowim=nrowi-1                                               12d14s20
            idenh=ibcoff                                                 12d14s20
            idend=idenh+nrowi                                           12d14s20
            ibcoff=idend+nrowi                                          12d14s20
            ijsbv=multh(isbv,jsbv)                                       12d14s20
            do isc=1,nsymb                                               12d14s20
             isd=multh(isymc,multh(isc,ijsbv))                          11d5s21
             nn=irefo(isc)*irefo(isd)                                   12d14s20
             idenk(isc)=ibcoff                                          11d5s21
             ndenk(isc)=idenk(isc)+nn*nrowi*ntype                       11d5s21
             idenkv(isc)=ndenk(isc)+nn*ntype                            11d5s21
             ndenkv(isc)=idenkv(isc)+nn*nrowi*ntype                     11d5s21
             ibcoff=ndenkv(isc)+nn*ntype                                11d5s21
             idenj(isc)=ibcoff                                          11d5s21
             ndenj(isc)=idenj(isc)+nn*nrowi*ntype                       11d5s21
             idenjv(isc)=ndenj(isc)+nn*ntype                            11d5s21
             ndenjv(isc)=idenjv(isc)+nn*nrowi*ntype                     11d5s21
             ibcoff=ndenjv(isc)+nn*ntype                                11d5s21
            end do                                                       12d14s20
            call enough('hcssbk4.  2',bc,ibc)
            ibctop=ibcoff-1                                              12d14s20
            if1o=nff1b(ncloip,isb,2)                                      12d14s20
            jgg=igg                                                     12d14s20
            jvecer=0                                                    2d12s21
            do if1=1,nff1b(ncloip,isb,1)                                7d9s21
             if(lprt)write(6,*)('for if1 = '),if1,if1o
             ndenh=0                                                    12d14s20
             do i=idenh,ibctop                                           12d14s20
              bc(i)=0d0                                                  12d14s20
             end do                                                      12d14s20
             i1c=iff1b(if1o)                                                11d25s20
             ntest=popcnt(i1c)
             if(ntest.ne.ncloi)then
              write(6,*)('ntest = '),ntest,if1
              call dcbit(i1c,32,'dorb')
              call dws_synca
              call dws_finalize
              stop
             end if
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
             if1o=if1o+1                                                   11d25s20
             i1o=iff1b(if1o)                                                11d25s20
             i1o=ibset(i1o,norbx)                                          11d25s20
             ii=1                                                          11d25s20
             do i=1,norb                                                   11d25s20
              if(btest(i1o,i))then                                         11d25s20
               isorb1(ii)=i                                                11d25s20
               itest(i,1)=1                                                11d25s20
               ii=ii+1                                                     11d25s20
              end if                                                       11d25s20
             end do                                                        11d25s20
             if1o=if1o+1                                                   11d25s20
             jf1o=nff1k(nclojp,jsb,2)                                   7d9s21
             jvec=ivec(nuse,jsb)                                        2d12s21
             do jf1=1,nff1k(nclojp,jsb,1)                               7d9s21
              if(lprt)write(6,*)('for jf1 '),jf1,jf1o,loop
              j1c=iff1k(jf1o)                                           7d9s21
              jf1o=jf1o+1                                                   11d25s20
              j1o=iff1k(jf1o)                                           7d9s21
              j1o=ibset(j1o,norbxx)                                          11d25s20
              jf1o=jf1o+1                                                   11d25s20
              gandcc=ieor(j1c,i1c)                                        2d6s23
              gandco=ieor(j1o,i1o)                                        2d6s23
              gandcb=ior(gandcc,gandco)                                     2d6s23
              ndifb=popcnt(gandcb)                                          2d6s23
c
c     vec is ordered nrootu,nherev. re-order to be nherev,nrootu.
              if(itransv(nuse,jsb).eq.0.and.ndifb.le.4)then             2d8s23
               call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))        2d12s21
               itmp=ibcoff                                                    1d27s21
               ibcoff=itmp+nrow                                            1d27s21
               call enough('hcssbk4.  3',bc,ibc)
               jjvec=ivec(nuse,jsb)                                     2d12s21
               nrowm=nrow-1                                                1d27s21
               do i=0,ncol-1                                               1d27s21
                do iv=0,nherev-1                                        4d15s21
                 do ir=0,nrootm                                               1d27s21
                  ij=jjvec+ir+nrootu*iv                                     1d27s21
                  ji=itmp+iv+nherev*ir                                  4d15s21
                  bc(ji)=bc(ij)                                            1d27s21
                 end do                                                    1d27s21
                end do                                                     1d27s21
                do j=0,nrowm                                               1d27s21
                 bc(jjvec+j)=bc(itmp+j)                                        1d27s21
                end do                                                        1d27s21
                jjvec=jjvec+nrow                                             1d27s21
               end do                                                         1d27s21
               ibcoff=itmp                                                    1d27s21
               itransv(nuse,jsb)=1                                      2d12s21
              end if                                                     2d12s21
              if(ndifb.le.4)then                                        2d7s23
               ivtrans=ibcoff                                           2d19s21
               ibcoff=ivtrans+nrow*ncsfk(jarg)                          11d5s21
               call enough('hcssbk4.  4',bc,ibc)
               do j=0,ncsfk(jarg)-1                                     11d5s21
                do i=0,nrow-1                                           2d19s21
                 ij=jvec+i+nrow*j                                       2d19s21
                 ji=ivtrans+j+ncsfk(jarg)*i                             11d5s21
                 bc(ji)=bc(ij)                                          2d19s21
                end do                                                  2d19s21
               end do                                                   2d19s21
               ndifs=popcnt(gandco)                                     2d6s23
               ndifd=popcnt(gandcc)                                     2d6s23
               if(ndifs.eq.2.and.ndifb.eq.2)then                        2d6s23
                nnot=2                                                  2d7s23
                do i=1,norbxx
                 if(btest(gandco,i))then
                  if((btest(i1c,i).and.btest(j1o,i)).or.                 10d21s22
     $             (btest(i1o,i).and..not.btest(j1c,i)))then            10d21s22
                   nab4(1,1)=i
                  else
                   nab4(2,1)=i
                  end if
                 end if
                end do
                call gandcr(i1c,i1o,j1c,j1o,nopeni,nopenj,norbxx,nnot1,  11d7s21
     $             nab1,icode1,imap1,nx1,irw1,irw2,iwpb1,iwpk1,bc,ibc)  11d14s22
                call gencup(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,nab1,    11d5s21
     $              iwpb1,iwpk1,iouta,imata,ntypea,mcsf,bc(ivtrans),    11d5s21
     $              ncsfk(jarg),nrow,bc,ibc)                            11d14s22
                nok=0                                                         11d13s20
                do i=1,norb                                                   11d25s20
                 itest(i,2)=0                                                 11d25s20
                end do                                                        11d25s20
                ii=1                                                          11d25s20
                do i=1,norb                                                   11d25s20
                 if(btest(j1c,i))then                                         11d25s20
                  itest(i,2)=2                                                11d25s20
                  ii=ii+1                                                     11d25s20
                 end if                                                       11d25s20
                end do                                                        11d25s20
                ii=1                                                          11d25s20
                do i=1,norb                                                   11d25s20
                 if(btest(j1o,i))then                                         11d25s20
                  itest(i,2)=1                                                11d25s20
                  ii=ii+1                                                     11d25s20
                 end if                                                       11d25s20
                end do                                                        11d25s20
                 do i=1,norb                                                   11d13s20
                 ixn=min(itest(i,1),itest(i,2))
                 if(ixn.gt.0)then                                             11d13s20
                  nok=nok+1                                                   11d13s20
                  itest(nok,3)=ixn                                            11d13s20
                  itest(nok,4)=i                                         11d5s21
                 end if                                                       11d13s20
                end do                                                        11d13s20
                do ii=0,ntypea-1                                           10d20s21
                 ipackc=ibc(iouta+ii)                                      9d10s21
                 imat=imata+nrowi*ii                                     11d9s21
                 if(ipackc1(1).gt.0)then                                   9d10s21
                  imy(3)=0                                                 10d12s21
                 else
                  imy(3)=1                                                 10d12s21
                 end if
                 if(ipackc1(2).gt.0)then                                   9d10s21
                  imy(4)=0                                                 10d12s21
                 else
                  imy(4)=1                                                 10d12s21
                 end if
                 phs=1d0                                                   10d20s21
                 if(imy(3).eq.1.and.isopt(4).ne.0)phs=-1d0                 10d20s21
                 ndenh=1                                                 11d7s21
                 do i=0,nrowim                                            11d5s21
                  bc(idenh+i)=bc(idenh+i)+phs*bc(imat+i)                   10d20s21
                 end do                                                     11d26s20
                 do i=1,norb                                               9d10s21
                  if(btest(i1c,i).and.btest(j1c,i).and.l2e.eq.0)then     2d18s22
                   js=ism(i)                                               9d10s21
                   jg=irel(i)-1                                            9d10s21
                   do jpass=0,1                                            9d10s21
c
c                                          a    b    c  d
c     integral is (jg,jg|imy(3)imy(4))=(imy(3)imy(4)|jg,jg)
c
                    ltest=imy(3)+2*(imy(4)+2*(jpass+2*jpass))             11d5s21
                    itestp=ltest+1                                         10d12s21
                    jtest=1-imy(3)+2*(1-imy(4)+2*(1-jpass+2*(1-jpass)))   11d5s21
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
                     iicolj=jg+irefo(js)*(jg+irefo(js)*iuse)              11d5s21
                     jden=idenj(js)+nrowi*iicolj                          11d5s21
                     ibc(ndenj(js)+iicolj)=1                              11d5s21
                     do iii=0,nrowim                                      11d5s21
                      bc(jden+iii)=bc(jden+iii)+bc(imat+iii)*phs           10d20s21
                     end do                                                10d22s21
                    end if                                                 10d20s21
c
c                   a   b     c      d
c     integral is (jg,imy(4)|imy(3),jg)=+/-(jg,imy(3)|imy(4),jg)
c
                    ltest=jpass+2*(imy(3)+2*(imy(4)+2*jpass))            11d7s21
                    itestp=ltest+1                                         10d12s21
                    jtest=1-jpass+2*(1-imy(3)+2*(1-imy(4)+2*(1-jpass)))  11d7s21
                    jtestp=jtest+1                                         10d12s21
                    iuse=-1                                               11d5s21
                    phs=1d0                                               10d20s21
                    if(isopt(2).ne.0)phs=-phs                            11d7s21
                    if(iifmx(itestp).ge.0)then                             10d20s21
                     iuse=iifmx(itestp)                                    10d20s21
                     kchoice=iuse
                    else if(iifmx(jtestp).ge.0)then                        10d20s21
                     iuse=iifmx(jtestp)                                    10d20s21
                     if(isopt(4).ne.0)phs=-phs                             10d20s21
                     kchoice=100+iuse
                    end if                                                 10d20s21
                    if(iuse.ge.0)then                                      10d20s21
                     iicolk=jg+irefo(js)*(jg+irefo(js)*iuse)              11d5s21
                     jden=idenk(js)+nrowi*iicolk                          11d5s21
                     ibc(ndenk(js)+iicolk)=1                              11d5s21
                     do iii=0,nrowim                                      11d5s21
                      bc(jden+iii)=bc(jden+iii)-bc(imat+iii)*phs           10d22s21
                     end do                                                   11d26s20
                    end if                                                 10d20s21
                   end do                                                  10d22s21
                  end if                                                   10d22s21
                 end do                                                    10d22s21
                end do                                                     10d22s21
                ibcoff=iouta                                               10d22s21
                do i=1,nok                                                    11d13s20
                 if(itest(i,3).eq.1.and.l2e.eq.0)then                    2d18s22
                  itestc=i1c                                             11d7s21
                  itesto=i1o                                             11d7s21
                  nopenk=nopeni                                          11d7s21
c
c     anihilate common
c
                  if(btest(itestc,itest(i,4)))then                             11d13s20
                   itestc=ibclr(itestc,itest(i,4))                        11d5s21
                   itesto=ibset(itesto,itest(i,4))                        11d5s21
                   nopenk=nopenk+1                                        11d5s21
                  else                                                         11d13s20
                   itesto=ibclr(itesto,itest(i,4))                        11d5s21
                   nopenk=nopenk-1                                             11d13s20
                  end if                                                       11d13s20
c
c     create ket
c
                  if(btest(itesto,nab4(2,1)))then                         12d14s20
                   itestc=ibset(itestc,nab4(2,1))                         12d14s20
                   itesto=ibclr(itesto,nab4(2,1))                         12d14s20
                   nopenk=nopenk-1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibset(itesto,nab4(2,1))                         12d14s20
                   nopenk=nopenk+1                                             11d13s20
                  end if                                                       11d13s20
                  call gandcr(i1c,i1o,itestc,itesto,nopeni,nopenk,        10d22s21
     $         norbxx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,iwpb1,iwpk1,11d14s22
     $                bc,ibc)                                           11d14s22
                  call gandcr(itestc,itesto,j1c,j1o,nopenk,nopenj,        11d5s21
     $         norbxx,nnot2,nab2,icode1,imap1,nx1,irw1,irw2,iwpb2,iwpk2,11d14s22
     $                bc,ibc)                                           11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                       11d5s21
                   ibail=0
                   if(nrow.gt.ncsfk(jarg).and.
     $                  min(ncsfk(jarg),ncsfb(iarg)).gt.1)then
                    ibail=1
                    iunit=ibcoff
                    ibcoff=iunit+ncsfk(jarg)*ncsfk(jarg)
                    call enough('hcssbk4.unit',bc,ibc)
                    do ix=iunit,ibcoff-1
                     bc(ix)=0d0
                    end do
                    do ix=0,ncsfk(jarg)-1
                     ii=iunit+ix*(ncsfk(jarg)+1)
                     bc(ii)=1d0
                    end do
                    ibcoff=max(ibcoff,iunit+ncsfb(iarg)*nrow)           3d14s23
                   end if
                   if(ibail.eq.0)then
                    call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,     11d5s21
     $               nopenk,ncsfb(iarg),nrow,itype,imatx,ntypex,nab1,   11d7s21
     $               iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(ivtrans),
     $               ncsfk(jarg),ieoro,bc,ibc)                          11d14s22
                   else                                                 3d14s23
                    call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,  3d14s23
     $               nopenk,ncsfb(iarg),ncsfk(jarg),itype,imatx,ntypex, 3d14s23
     $               nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),       3d14s23
     $               ncsfk(jarg),ieoro,bc,ibc)                          11d14s22
                   end if
                   do ini=0,ntypex-1                                      10d24s21
                    ipack=ibc(itype+ini)
                    if(ibail.eq.0)then                                  3d14s23
                     imat=imatx+nrowi*ini                                 11d7s21
                    else                                                3d14s23
                     imatt=imatx+ncsfb(iarg)*ncsfk(jarg)*ini            3d14s23
                     call dgemm('n','n',ncsfb(iarg),nrow,ncsfk(jarg),   3d14s23
     $                    1d0,bc(imatt),ncsfb(iarg),bc(ivtrans),
     $                    ncsfk(jarg),0d0,bc(iunit),ncsfb(iarg),        3d14s23
     $                    'hcssbk4.matt')                               3d14s23
                     imat=iunit                                         3d14s23
                    end if                                              3d14s23
                    do j=1,4                                             11d7s21
                     ipack4(j)=ipack2(j)                                 11d7s21
                     ipack4(j)=iabs(ipack4(j))                           11d7s21
                     if(ipack2(j).gt.0)then                              11d7s21
                      imy(j)=0                                           11d7s21
                     else                                                11d7s21
                      imy(j)=1                                           11d7s21
                     end if                                              11d7s21
                     if(ipack4(j).le.norb)then                           11d7s21
                      isy(j)=ism(ipack4(j))                              11d7s21
                      igya(j)=irel(ipack4(j))-1                          11d7s21
                     end if                                              11d7s21
                    end do                                               11d7s21
                    if(ipack4(3).eq.norbx.and.ipack4(2).eq.norbxx)then   11d7s21
c     v=norbx=isbv, v'=norbxx=jsbv. kmatsr((isc) isbv,jsbv,isd)
c
c     (1v'|v4)=+/-(4v|v'1)
c
                     phs=1d0                                             11d7s21
                     if(isopt(2).ne.0)phs=-phs                           11d7s21
                     ltest=imy(4)+2*(imy(3)+2*(imy(2)+2*imy(1)))         11d7s21
                     jtest=1-imy(4)+2*(1-imy(3)+2*(1-imy(2)              11d7s21
     $                   +2*(1-imy(1))))                                11d7s21
                     itestp=ltest+1                                      11d7s21
                     jtestp=jtest+1                                      11d7s21
                     iuse=-1                                             11d7s21
                     if(iifmx(itestp).ge.0)then                          11d7s21
                      iuse=iifmx(itestp)                                 11d7s21
                     else if(iifmx(jtestp).ge.0)then                     11d7s21
                      iuse=iifmx(jtestp)                                 11d7s21
                      if(isopt(4).ne.0)phs=-phs                          11d7s21
                     end if                                              11d7s21
                     if(iuse.ge.0)then                                   11d7s21
                      icolk=igya(4)+irefo(isy(4))*(igya(1)+irefo(isy(1)) 11d7s21
     $                   *iuse)                                         11d7s21
                      ibc(ndenk(isy(4))+icolk)=1                         11d7s21
                      jden=idenk(isy(4))+nrowi*icolk                     11d7s21
                      do j=0,nrowim                                      11d7s21
                       bc(jden+j)=bc(jden+j)+phs*bc(imat+j)              11d7s21
                      end do                                             11d7s21
                     end if                                              11d7s21
c
c     jmatsr((isbv) jsbv,isc,isd)
c     (14|vv')=(vv'|14)
c
                     phs=-1d0                                            11d7s21
                     ltest=imy(3)+2*(imy(2)+2*(imy(1)+2*imy(4)))         11d7s21
                     jtest=1-imy(3)+2*(1-imy(2)+2*(1-imy(1)              11d7s21
     $                    +2*(1-imy(4))))                                11d7s21
                     itestp=ltest+1                                      11d7s21
                     jtestp=jtest+1                                      11d7s21
                     iuse=-1                                             11d7s21
                     if(iifmx(itestp).ge.0)then                          11d7s21
                      iuse=iifmx(itestp)                                 11d7s21
                     else if(iifmx(jtestp).ge.0)then                          11d7s21
                      iuse=iifmx(jtestp)                                 11d7s21
                      if(isopt(4).ne.0)phs=-phs                          11d7s21
                     end if                                              11d7s21
                     if(iuse.ge.0)then                                   11d7s21
                      icolj=igya(1)+irefo(isy(1))*(igya(4)+irefo(isy(4)) 11d7s21
     $                   *iuse)                                         11d7s21
                      ibc(ndenj(isy(1))+icolj)=1                         11d7s21
                      jden=idenj(isy(1))+nrowi*icolj                     11d7s21
                      do j=0,nrowim                                      11d7s21
                       bc(jden+j)=bc(jden+j)+phs*bc(imat+j)              11d7s21
                      end do                                             11d7s21
                     end if                                              11d7s21
                    else                                                 11d7s21
                     write(6,*)('I don''t know how to handle this set'),
     $                    ipack4
                     call dws_synca
                     call dws_finalize
                     stop 'hcssbk4'
                    end if                                               11d7s21
                   end do                                                 10d22s21
                   ibcoff=itype                                           10d22s21
                  end if                                                   12d14s20
                 end if                                                  11d7s21
                end do                                                    12d14s20
               else if(l2e.eq.0)then                                    2d7s23
                nnot=0                                                  2d7s23
                if(ndifs.eq.4.and.ndifb.eq.4)then                       2d6s23
                 nnot=4                                                 2d6s23
                 ioxx(1)=1                                              2d6s23
                 ioxx(2)=1                                              2d6s23
                 do i=1,norbxx                                          2d6s23
                  if(btest(gandcb,i))then                               2d6s23
                   if((btest(j1c,i).and.btest(i1o,i)).or.               2d6s23
     $                   (btest(j1o,i).and..not.btest(i1c,i)))then      2d6s23
                    nab4(2,ioxx(2))=i                                   2d6s23
                    ioxx(2)=ioxx(2)+1                                   2d6s23
                   else                                                 2d6s23
                    nab4(1,ioxx(1))=i                                   2d6s23
                    ioxx(1)=ioxx(1)+1                                   2d6s23
                   end if                                               2d6s23
                  end if                                                2d6s23
                 end do                                                 2d6s23
                else if(ndifb.eq.3)then                                 2d6s23
                 nnot=3                                                 2d6s23
                 ioxx(1)=1                                              2d6s23
                 ioxx(2)=1                                              2d6s23
                 iswap=0                                                2d6s23
                 do i=1,norbxx                                          2d6s23
                  if(btest(gandcb,i))then                               2d6s23
                   if(btest(gandcc,i).and.                              2d6s23
     $                   ((btest(i1c,i).and..not.btest(j1o,i)).or.      2d6s23
     $                   (btest(j1c,i).and..not.btest(i1o,i))))then     2d6s23
                    if(btest(j1c,i))iswap=1                             2d6s23
                    nab4(1,1)=i                                         2d6s23
                    nab4(1,2)=i                                         2d6s23
                   else                                                 2d6s23
                    nab4(2,ioxx(2))=i                                   2d6s23
                    ioxx(2)=ioxx(2)+1                                   2d6s23
                   end if                                               2d6s23
                  end if                                                2d6s23
                 end do                                                 2d6s23
                 if(iswap.ne.0)then                                     2d6s23
                  icpy=nab4(1,1)                                        2d6s23
                  nab4(1,1)=nab4(2,1)                                   2d6s23
                  nab4(2,1)=icpy                                        2d6s23
                  icpy=nab4(1,2)                                        2d6s23
                  nab4(1,2)=nab4(2,2)                                   2d6s23
                  nab4(2,2)=icpy                                        2d6s23
                  nbt=0                                                 2d6s23
                  if(btest(j1c,nab4(2,2)).and.                          2d6s23
     $                   .not.btest(j1c,nab4(2,1)))nbt=1                2d6s23
                 else                                                   2d6s23
                  nbt=0                                                 2d6s23
                  if(btest(i1c,nab4(1,2)).and.                          2d6s23
     $                 .not.btest(i1c,nab4(1,1)))nbt=1                  2d6s23
                 end if                                                 2d6s23
                 if(nbt.ne.0)then                                       2d6s23
                  nab4(1,1)=nab4(1,2)                                   2d6s23
                  nab4(2,1)=nab4(2,2)                                   2d6s23
                 end if                                                 2d6s23
                else if(ndifs.eq.0.and.ndifd.eq.2)then                  2d6s23
                 nnot=3                                                 2d6s23
                 do i=1,norbxx                                          2d6s23
                  if(btest(gandcb,i))then                               2d6s23
                   if(btest(i1c,i))then                                 2d6s23
                    nab4(1,1)=i                                         2d6s23
                    nab4(1,2)=i                                         2d6s23
                   else                                                 2d6s23
                    nab4(2,1)=i                                         2d6s23
                    nab4(2,2)=i                                         2d6s23
                   end if                                               2d6s23
                  end if                                                2d6s23
                 end do                                                 2d6s23
                end if                                                  2d6s23
                if(nnot.gt.0)then                                       2d7s23
                 iu1=1
                 iu2=1
                 itestc=i1c                                               11d7s21
                 itesto=i1o                                               11d7s21
                 if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                  itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                  itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopeni+1                                         11d7s21
                 else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                  itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopeni-1                                         11d7s21
                 end if                                                        11d13s20
                 if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                  itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                  itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                  nopenk=nopenk-1                                              11d13s20
                 else                                                          11d13s20
                  itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                  nopenk=nopenk+1                                              11d13s20
                 end if                                                        11d13s20
                 call gandcr(i1c,i1o,itestc,itesto,nopeni,nopenk,            11d26s20
     $         norbxx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,iwpb1,iwpk1,11d14s22
     $              bc,ibc)                                             11d14s22
                 call gandcr(itestc,itesto,j1c,j1o,nopenk,nopenj,             11d26s20
     $         norbxx,nnot2,nab2,icode1,imap1,nx1,irw1,irw2,iwpb2,iwpk2,11d14s22
     $              bc,ibc)                                             11d14s22
                 if(nnot1.eq.2.and.nnot2.eq.2)then                              11d13s20
                  ibail=0
                  if(nrow.gt.ncsfk(jarg).and.
     $                min(ncsfk(jarg),ncsfb(iarg)).gt.1)then
                   ibail=1
                   iunit=ibcoff
                   ibcoff=iunit+ncsfk(jarg)*ncsfk(jarg)
                   call enough('hcssbk4.unit',bc,ibc)
                   do ix=iunit,ibcoff-1
                    bc(ix)=0d0
                   end do
                   do ix=0,ncsfk(jarg)-1
                    ii=iunit+ix*(ncsfk(jarg)+1)
                    bc(ii)=1d0
                   end do
                   ibcoff=max(ibcoff,iunit+ncsfb(iarg)*nrow)            3d14s23
                  end if
                  if(ibail.eq.0)then                                    3d14s23
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,        10d22s21
     $               nopenk,ncsfb(iarg),nrow,itype,imatx,ntypex,nab1,   11d7s21
     $             iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(ivtrans),ncsfk(jarg),11d5s21
     $             ieoro,bc,ibc)                                        11d14s22
                  else                                                  3d14s23
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,   3d14s23
     $               nopenk,ncsfb(iarg),ncsfk(jarg),itype,imatx,ntypex, 3d14s23
     $               nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),       3d14s23
     $               ncsfk(jarg),ieoro,bc,ibc)                          3d14s23
                  end if                                                3d14s23
                  do i=0,ntypex-1                                          10d24s21
                   ipack=ibc(itype+i)
                   imat=imatx+nrowi*i                                     11d7s21
                   if(ibail.ne.0)then                                   3d14s23
                    imatt=imatx+ncsfb(iarg)*ncsfk(jarg)*i               3d14s23
                    call dgemm('n','n',ncsfb(iarg),nrow,ncsfk(jarg),    3d14s23
     $                1d0,bc(imatt),ncsfb(iarg),bc(ivtrans),ncsfk(jarg),3d14s23
     $                0d0,bc(iunit),ncsfb(iarg))                        3d14s23
                    imat=iunit                                          3d14s23
                   end if                                               3d14s23
                   do j=1,4                                               11d7s21
                    ipack4(j)=ipack2(j)                                   11d7s21
                    ipack4(j)=iabs(ipack4(j))                             11d7s21
                    if(ipack2(j).gt.0)then                                11d7s21
                     imy(j)=0                                             11d7s21
                    else                                                  11d7s21
                     imy(j)=1                                             11d7s21
                    end if                                                11d7s21
                    if(ipack4(j).le.norb)then                             11d7s21
                     isy(j)=ism(ipack4(j))                                11d7s21
                     igya(j)=irel(ipack4(j))-1                            11d7s21
                    end if                                                11d7s21
                   end do                                                 11d7s21
                   if(ipack4(3).eq.norbx.and.ipack4(4).eq.norbxx)then     11d7s21
                    phs=1d0                                               11d7s21
                    ltest=imy(3)+2*(imy(4)+2*(imy(1)+2*imy(2)))           11d7s21
                    jtest=1-imy(3)+2*(1-imy(4)+2*(1-imy(1)              2d7s23
     $                   +2*(1-imy(2))))                                2d7s23
                    itestp=Ltest+1                                        11d7s21
                    jtestp=jtest+1                                        11d7s21
                    iuse=-1                                               11d7s21
                    if(iifmx(itestp).ge.0)then                             11d7s21
                     iuse=iifmx(itestp)                                    11d7s21
                    else if(iifmx(jtestp).ge.0)then                       11d7s21
                     iuse=iifmx(jtestp)                                   11d7s21
                     if(isopt(4).ne.0)phs=-phs                            11d7s21
                    end if                                                11d7s21
                    if(iuse.ge.0)then                                     11d7s21
                     iicolj=igya(1)+irefo(isy(1))*(igya(2)+irefo(isy(2))  11d7s21
     $                 *iuse)                                           11d7s21
                     ibc(ndenj(isy(1))+iicolj)=1                          11d7s21
                     jden=idenj(isy(1))+nrowi*iicolj                      11d7s21
                     do j=0,nrowim                                        11d7s21
                      bc(jden+j)=bc(jden+j)+phs*bc(imat+j)                11d7s21
                     end do                                               11d7s21
                    end if                                                11d7s21
                    if(nnot.eq.4)then                                     11d7s21
c               v'     v                 a       b     c     d
c     (imy(1),imy(4)|imy(3),imy(2))=+/-(imy(2)imy(3)|imy(4)imy(1))
                     phs=-1d0                                             11d7s21
                     if(isopt(2).ne.0)phs=-phs                            11d7s21
                     ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))           11d7s21
                     jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)              2d7s23
     $                   +2*(1-imy(1))))                                2d7s23
                     itestp=ltest+1                                        11d7s21
                     jtestp=jtest+1                                        11d7s21
                     iuse=-1                                               11d7s21
                     if(iifmx(itestp).ge.0)then                             11d7s21
                      iuse=iifmx(itestp)                                    11d7s21
                     else if(iifmx(jtestp).ge.0)then                       11d7s21
                      iuse=iifmx(jtestp)                                   11d7s21
                      if(isopt(4).ne.0)phs=-phs                            11d7s21
                     end if                                                11d7s21
                     if(iuse.ge.0)then                                     11d7s21
                      iicolk=igya(2)+irefo(isy(2))*(igya(1)             2d7s23
     $                    +irefo(isy(1))*iuse)                          2d7s23
                      ibc(ndenk(isy(2))+iicolk)=1                          11d7s21
                      jden=idenk(isy(2))+nrowi*iicolk                      11d7s21
                      do j=0,nrowim                                        11d7s21
                       bc(jden+j)=bc(jden+j)+phs*bc(imat+j)                11d7s21
                      end do                                              11d7s21
                     end if                                                11d7s21
                    end if                                                11d7s21
                   else                                                   11d7s21
                    write(6,*)('I don''t know how to handle case '),
     $                 ipack2
                    call dws_synca
                    call dws_finalize
                    stop
                   end if                                                 11d7s21
                  end do                                                   11d5s21
                  ibcoff=itype                                             11d5s21
                 end if                                                 2d7s23
                end if                                                   11d5s21
               end if                                                   2d7s23
               ibcoff=ivtrans                                           2d8s23
              end if                                                     12d14s20
              if(jsbv.eq.isbv)then                                      7d19s21
               j1o=ibclr(j1o,norbxx)                                     12d14s20
               j1o=ibset(j1o,norbx)                                      12d14s20
               gandco=ieor(j1o,i1o)                                        2d6s23
               gandcb=ior(gandcc,gandco)                                     2d6s23
               ndifb=popcnt(gandcb)                                          2d6s23
c
c     vec is ordered nrootu,nherev. re-order to be nherev,nrootu.
               if(itransv(nuse,jsb).eq.0.and.ndifb.le.4)then            2d7s23
                call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))       2d12s21
                itmp=ibcoff                                                    1d27s21
                ibcoff=itmp+nrow                                            1d27s21
                call enough('hcssbk4.  5',bc,ibc)
                jjvec=ivec(nuse,jsb)                                    2d12s21
                nrowm=nrow-1                                                1d27s21
                do i=0,ncol-1                                               1d27s21
                 do iv=0,nherev-1                                       4d15s21
                  do ir=0,nrootm                                               1d27s21
                   ij=jjvec+ir+nrootu*iv                                     1d27s21
                   ji=itmp+iv+nherev*ir                                 4d15s21
                   bc(ji)=bc(ij)                                            1d27s21
                  end do                                                    1d27s21
                 end do                                                     1d27s21
                 do j=0,nrowm                                               1d27s21
                  bc(jjvec+j)=bc(itmp+j)                                        1d27s21
                 end do                                                        1d27s21
                 jjvec=jjvec+nrow                                             1d27s21
                end do                                                         1d27s21
                ibcoff=itmp                                                    1d27s21
                itransv(nuse,jsb)=1                                     2d12s21
               end if                                                     2d12s21
               if(ndifb.le.4)then                                       2d7s23
                ivtrans=ibcoff                                           2d19s21
                ibcoff=ivtrans+nrow*ncsfk(jarg)                         11d5s21
                call enough('hcssbk4.  6',bc,ibc)
                do j=0,ncsfk(jarg)-1                                    11d5s21
                 do i=0,nrow-1                                           2d19s21
                  ij=jvec+i+nrow*j                                       2d19s21
                  ji=ivtrans+j+ncsfk(jarg)*i                            11d5s21
                  bc(ji)=bc(ij)                                          2d19s21
                 end do                                                  2d19s21
                end do                                                   2d19s21
                ndifs=popcnt(gandco)                                    2d6s23
                ndifd=popcnt(gandcc)                                    2d6s23
                if(ndifs.eq.0.and.ndifd.eq.0)then                       2d7s23
                 nnot=1                                                 2d8s23
                 idelta=i2smk-i2smb                                      11d5s21
                 idelta=iabs(idelta)                                     11d5s21
                 ill=i2e+1                                               7d15s21
                 ihh=i2s-1                                               7d15s21
                 do i2=i2s,i2e                                           7d15s21
                  icol=i2+nvirt(isbv)*(i2-1)                             7d15s21
                  if(icol.ge.il.and.icol.le.ih)then                      7d15s21
                   ill=min(ill,i2)                                       7d15s21
                   ihh=max(ihh,i2)                                       7d15s21
                  end if                                                 7d15s21
                 end do                                                  7d15s21
                 nnhere=ihh+1-ill                                        7d15s21
                 if(idelta.eq.0)then                                        9d13s21
                  if(i2sb.eq.i2sk)then                                      9d13s21
                   sum=0d0                                                  9d13s21
                   do i=1,norb                                               9d1s21
                    if(btest(i1c,i))then                                    9d1s21
                     is=ism(i)                                               9d1s21
                     ig=irel(i)-1                                            9d1s21
                     ig0=ig                                                 10d8s21
                     do ipass=1,2                                            9d1s21
                      h0=-2d0                                               10d14s21
                      if(ih0n(is).gt.0)then                                 10d14s21
                       iadn=ih0n(is)+(nh0(is)+1)*ig0                        10d14s21
                       h0=bc(iadn)                                           10d14s21
                       if(ipass.eq.2.and.isopt(4).ne.0)h0=-h0            2d24s22
                       sum=sum+h0                                           10d14s21
                      end if                                                10d14s21
                      isp=ig/irefo(is)                                      10d8s21
                      ispm=1-isp                                            10d8s21
                      do j=1,norb                                            9d1s21
                       if(btest(i1c,j).and.l2e.eq.0)then                 2d18s22
                        js=ism(j)                                            9d1s21
                        jg=irel(j)-1                                         9d1s21
                        jg0=jg                                              10d8s21
                        do jpass=1,2                                         9d1s21
                         jsp=jg/irefo(js)                                   10d8s21
                         jspm=1-jsp                                         10d8s21
                         ltest=jsp+2*(jsp+2*(isp+2*isp))                    10d8s21
                         itestp=ltest+1                                     10d8s21
                         jtest=jspm+2*(jspm+2*(ispm+2*ispm))                    10d8s21
                         jtestp=jtest+1                                     10d8s21
                         if(iifmx(itestp).ge.0)then                         10d8s21
                          iad=i4or(js,is,is)+jg0                        2d7s23
     $                        +irefo(js)*(jg0+irefo(js)                 2d7s23
     $                  *(ig0+irefo(is)*(ig0+irefo(is)*iifmx(itestp)))) 10d8s21
                          sum=sum+bc(iad)*0.5d0                             10d8s21
                         else if(iifmx(jtestp).ge.0)then                    10d8s21
                          iad=i4or(js,is,is)+jg0                        2d7s23
     $                         +irefo(js)*(jg0+irefo(js)                2d7s23
     $                  *(ig0+irefo(is)*(ig0+irefo(is)*iifmx(jtestp)))) 10d8s21
                          xint=bc(iad)                                      10d8s21
                          if(isopt(4).ne.0)then                             10d8s21
                           xint=-xint                                       10d8s21
                          end if
                          sum=sum+xint*0.5d0                                10d8s21
                         end if                                             10d8s21
                         ltest=jsp+2*(isp+2*(isp+2*jsp))                    10d8s21
                         itestp=ltest+1                                     10d8s21
                         jtest=jspm+2*(ispm+2*(ispm+2*jspm))                    10d8s21
                         jtestp=jtest+1                                     10d8s21
                         if(iifmx(itestp).ge.0)then                         10d8s21
                          iad=i4or(is,is,js)+jg0                        2d7s23
     $                        +irefo(js)*(ig0+irefo(is)                 2d7s23
     $                  *(ig0+irefo(is)*(jg0+irefo(js)*iifmx(itestp)))) 10d8s21
                          sum=sum-bc(iad)*0.5d0                             10d8s21
                         else if(iifmx(jtestp).ge.0)then                    10d8s21
                          iad=i4or(is,is,js)+jg0                        2d7s23
     $                         +irefo(js)*(ig0+irefo(is)                2d7s23
     $                  *(ig0+irefo(is)*(jg0+irefo(js)*iifmx(jtestp)))) 10d8s21
                          xint=bc(iad)                                      10d8s21
                          if(isopt(4).ne.0)then                             10d8s21
                           xint=-xint                                       10d8s21
                          end if
                          sum=sum-xint*0.5d0                               10d8s21
                         end if                                             10d8s21
                         jg=jg+irefo(js)                                     9d1s21
                        end do                                               9d1s21
                       end if                                                9d1s21
                      end do                                                 9d1s21
                      ig=ig+irefo(is)                                        9d1s21
                     end do                                                  9d1s21
                    end if                                                   9d1s21
                   end do                                                    9d1s21
                   if(nnhere.gt.0)then                                     7d15s21
                    if(itransgg.eq.0)then                                      2d12s21
                     call ddi_done(ibc(iacc),nacc)                             2d12s21
                     itransgg=1                                                2d12s21
                     do i=igg,igg+ncolti-1                                     2d12s21
                      bc(i)=0d0                                                2d12s21
                     end do                                                    2d12s21
                    end if                                                     2d12s21
                    do i=0,ncsfb(iarg)-1                                  11d5s21
                     do ir=0,nrootm                                        7d15s21
                      iadv=jvec-i2s+nherev*(ir+nrootu*i)                   7d15s21
                      iadg=jgg-1+nvirt(isbv)*(ir+nrootu*i)                 7d15s21
                      do ivp=ill,ihh                                       7d15s21
                       bc(iadg+ivp)=bc(iadg+ivp)+sum*bc(iadv+ivp)         11d5s21
                      end do                                               7d15s21
                     end do                                                7d15s21
                    end do                                                 7d15s21
                   end if                                                 11d5s21
                  end if                                                 11d5s21
                 end if                                                  11d5s21
                 ii=0                                                       9d13s21
                 do i=1,norb                                                9d13s21
                  if(btest(i1o,i))then                                     9d13s21
                   ii=ii+1                                                  9d13s21
                   icode1(ii)=3                                             9d13s21
                   imap1(ii)=i
                  end if                                                    9d13s21
                 end do                                                     9d13s21
                 ii=ii+1                                                 11d5s21
                 icode1(ii)=3                                            11d5s21
                 imap1(ii)=norbx                                         11d5s21
                 mxmat=norb*32                                             9d15s21
                 itype=ibcoff
                 imatx=itype+mxmat
                 ibcoff=imatx+mxmat*ncsfb(iarg)*ncsfk(jarg)              11d7s21
                 call enough('hcssbk4.  7',bc,ibc)
                 do iz=imatx,ibcoff-1                                      9d15s21
                  bc(iz)=0d0                                               9d15s21
                 end do                                                    9d15s21
                 if(idelta.le.2)then                                        9d13s21
                  sigc=1d0                                                  9d21s21
                  ntype1a=0                                                 10d14s21
                  ibcb4=ibcoff                                           11d5s21
                  call getcup(irw0,nopeni,i2sb,i2smb,i2sk,i2smk,ntypeg,     10d18s21
     $             ioutg,imatg,mcsf,imap1,bc,ibc)                       11d14s22
                  if(ntypeg.gt.0)then
                   ncsfa=mcsf(1)                                            10d18s21
                   ncsfc=mcsf(2)                                            10d19s21
                   itmp=ibcoff                                           11d5s21
                   idtmp=itmp+nrowi                                      11d9s21
                   ibcoff=idtmp+nnhere                                    7d15s21
                   call enough('hcssbk4.  8',bc,ibc)
                   jdtmp=idtmp-ill                                        7d15s21
                   do iq=0,ntypeg-1
                    do i=0,nnhere-1                                      11d5s21
                     bc(idtmp+i)=0d0                                     11d5s21
                    end do                                                 7d15s21
                    ipackc=ibc(ioutg+iq)                                    10d19s21
                    if(ipackc1(1).gt.0)then                                 10d19s21
                     mmm=ipackc1(1)                                         10d19s21
                     if(mmm.le.norb)then                                 11d5s21
                      lsb=ism(mmm)                                           10d14s21
                      igya(1)=irel(mmm)-1                                    10d14s21
                     else
                      lsb=isbv                                           11d5s21
                      igya(1)=-1                                         11d5s21
                     end if                                              11d5s21
                     imy(1)=0                                               10d14s21
                    else                                                    10d14s21
                     mmm=-ipackc1(1)                                        10d19s21
                     if(mmm.le.norb)then                                 11d5s21
                      lsb=ism(mmm)                                           10d14s21
                      igya(1)=irel(mmm)-1                                    10d14s21
                     else                                                11d5s21
                      lsb=isbv                                           11d5s21
                      igya(1)=-1                                         11d5s21
                     end if                                              11d5s21
                     imy(1)=1                                               10d14s21
                    end if                                                  10d14s21
                    if(ipackc1(2).gt.0)then                                 10d19s21
                     mmm=ipackc1(2)                                         10d19s21
                     if(mmm.le.norb)then                                 11d5s21
                      lsk=ism(mmm)                                           10d14s21
                      igya(2)=irel(mmm)-1                                    10d14s21
                     else                                                11d5s21
                      lsk=jsbv                                           11d5s21
                      igya(2)=-1                                         11d5s21
                     end if                                              11d5s21
                     imy(2)=0                                               10d14s21
                    else                                                    10d14s21
                     mmm=-ipackc1(2)                                        10d19s21
                     if(mmm.le.norb)then                                 11d5s21
                      lsk=ism(mmm)                                           10d14s21
                      igya(2)=irel(mmm)-1                                    10d14s21
                     else                                                11d5s21
                      lsk=jsbv                                           11d5s21
                      igya(2)=-1                                         11d5s21
                     end if                                              11d5s21
                     imy(2)=1                                               10d14s21
                    end if                                                  10d14s21
                    imatq=imatg+mcsf(1)*mcsf(2)*iq
                    sum=0d0                                                  10d14s21
                    h0=-2d0                                                   10d13s21
                    if(ih0n(lsb).gt.0)then                                    10d13s21
                     phs=1d0
                     if(imy(1).ne.0.and.isopt(4).ne.0)phs=-phs           11d5s21
                     if(igya(1).ge.0.and.igya(2).ge.0)then               11d5s21
                      iad=ih0n(lsb)+igya(1)+nh0(lsb)*igya(2)                  10d14s21
                      sum=bc(iad)*phs                                    11d5s21
                     else if(igya(1).ge.0)then                           11d5s21
                      do ivp=ill,ihh                                     11d5s21
                       iv=ivp-1+irefo(lsk)                               11d5s21
                       iad=ih0n(lsb)+igya(1)+nh0(lsb)*iv                  11d5s21
                       bc(jdtmp+ivp)=bc(iad)*phs                         11d5s21
                      end do                                             11d5s21
                     else if(igya(2).ge.0)then                           11d5s21
                      do ivp=ill,ihh                                     11d5s21
                       iv=ivp-1+irefo(lsb)                               11d5s21
                       iad=ih0n(lsb)+iv+nh0(lsb)*igya(2)                  11d5s21
                       bc(jdtmp+ivp)=bc(iad)*phs                         11d5s21
                      end do                                             11d5s21
                     else                                                11d9s21
                      do ivp=ill,ihh                                     11d5s21
                       iv=ivp-1+irefo(lsb)                               11d5s21
                       iad=ih0n(lsb)+iv+nh0(lsb)*iv                      11d9s21
                       bc(jdtmp+ivp)=bc(iad)*phs                         11d5s21
                      end do                                             11d5s21
                     end if                                              11d5s21
                    else if(ih0n(lsk).gt.0)then                          11d5s21
                     phs=1d0                                             11d5s21
                     if(isopt(2).ne.0)phs=-phs                           11d5s21
                     if(imy(2).ne.0.and.isopt(4).ne.0)phs=-phs           11d5s21
                     if(igya(2).ge.0.and.igya(1).ge.0)then               11d5s21
                      iad=ih0n(lsk)+igya(2)+nh0(lsk)*igya(1)                  10d14s21
                      sum=bc(iad)*phs                                    11d5s21
                     else if(igya(2).ge.0)then                           11d5s21
                      do ivp=ill,ihh                                     11d5s21
                       iv=ivp-1+irefo(lsb)                               11d5s21
                       iad=ih0n(lsk)+igya(2)+nh0(lsk)*iv                  11d5s21
                       bc(jdtmp+ivp)=bc(iad)*phs                         11d5s21
                      end do                                             11d5s21
                     else if(igya(1).ge.0)then                           11d5s21
                      do ivp=ill,ihh                                     11d5s21
                       iv=ivp-1+irefo(lsk)                               11d5s21
                       iad=ih0n(lsk)+iv+nh0(lsk)*igya(1)                 11d5s21
                       bc(jdtmp+ivp)=bc(iad)*phs                         11d5s21
                      end do                                             11d5s21
                     else                                                11d9s21
                      do ivp=ill,ihh                                     11d5s21
                       iv=ivp-1+irefo(lsk)                               11d5s21
                       iad=ih0n(lsk)+iv+nh0(lsk)*iv                      11d9s21
                       bc(jdtmp+ivp)=bc(iad)*phs                         11d5s21
                      end do                                             11d5s21
                     end if                                              11d5s21
                    end if                                                    10d13s21
                    do i=1,norb                                             10d14s21
                     if((btest(j1c,i).or.btest(j1o,i)).and.l2e.eq.0)then 2d18s22
                      is=ism(i)                                               9d13s21
                      ig=irel(i)-1                                            9d13s21
                      if(btest(j1c,i))then                               11d5s21
                       ff=-1d0                                              10d14s21
                      else                                                  10d14s21
                       ff=-0.5d0                                            10d14s21
                      end if                                                10d14s21
                      do isp=0,1                                             10d14s21
                       ltest=isp+2*(imy(1)+2*(imy(2)+2*isp))                10d14s21
                       itestp=ltest+1                                        10d14s21
                       jtest=1-isp+2*(1-imy(1)+2*(1-imy(2)+2*(1-isp)))      10d14s21
                       jtestp=jtest+1                                        10d14s21
                       phs=ff                                            11d12s21
                       iuse=-1                                           11d5s21
                       if(iifmx(itestp).ge.0)then                           10d14s21
                        iuse=iifmx(itestp)                               11d5s21
                       else if(iifmx(jtestp).ge.0)then                           10d14s21
                        iuse=iifmx(jtestp)                               11d5s21
                        if(isopt(4).ne.0)phs=-phs                        11d5s21
                       end if                                               10d14s21
                       if(iuse.ge.0)then                                 11d5s21
                        if(igya(1).ge.0.and.igya(2).ge.0)then             11d5s21
                         iad=i4or(lsb,lsb,is)+ig+irefo(is)*(igya(1)       11d5s21
     $                     +irefo(lsb)*(igya(2)+irefo(lsb)*             11d5s21
     $                     (ig+irefo(is)*iuse)))                        11d5s21
                         sum=sum+bc(iad)*phs                             11d12s21
                        else                                             11d12s21
c
c     (gv|vg)
                         do ivp=ill,ihh                                     11d5s21
                          irow=ivp+nvirt(isbv)*(ivp-1)                     11d9s21
                          if(irow.ge.il.and.irow.le.ih)then                 11d5s21
                           irow=irow-il                                     11d5s21
                           iad=kmatsr(isbv,jsbv,is)+irow+nhere           11d12s21
     $                      *(ig+irefo(is)*(ig+irefo(is)*iuse))         11d12s21
                           bc(jdtmp+ivp)=bc(jdtmp+ivp)+bc(iad)*phs         11d9s21
                          end if                                         11d12s21
                         end do                                          11d12s21
                        end if                                            11d5s21
                       end if                                            11d5s21
                      end do                                                10d14s21
                     end if                                                 10d14s21
                    end do                                                  10d14s21
                    if(igya(1).ge.0.and.igya(2).ge.0.and.l2e.eq.0)then   2d18s22
c
c     (v1|2v)=+/-(1v|v2)
c
                     ff=-0.5d0                                            11d5s21
                     if(isopt(2).ne.0)ff=-ff                             11d9s21
                     do isp=0,1                                          11d5s21
                      ltest=imy(1)+2*(isp+2*(isp+2*imy(2)))              11d12s21
                      itestp=ltest+1                                        10d14s21
                      jtest=1-imy(1)+2*(1-isp+2*(1-isp+2*(1-imy(2))))    11d12s21
                      jtestp=jtest+1                                        10d14s21
                      phs=ff                                             11d12s21
                      iuse=-1                                            11d5s21
                      if(iifmx(itestp).ge.0)then                           10d14s21
                       iuse=iifmx(itestp)                                11d5s21
                      else if(iifmx(jtestp).ge.0)then                           10d14s21
                       iuse=iifmx(jtestp)                                11d5s21
                       if(isopt(4).ne.0)phs=-phs                         11d5s21
                      end if                                               10d14s21
                      if(iuse.ge.0)then                                  11d5s21
                       do ivp=ill,ihh                                     11d5s21
                        irow=ivp+nvirt(isbv)*(ivp-1)                     11d9s21
                        if(irow.ge.il.and.irow.le.ih)then                 11d5s21
                         irow=irow-il                                     11d5s21
                         iad=kmatsr(isbv,jsbv,lsk)+irow+nhere            11d9s21
     $                      *(igya(1)+irefo(lsb)*(igya(2)               11d12s21
     $                       +irefo(lsk)*iuse))                         11d12s21
                         bc(jdtmp+ivp)=bc(jdtmp+ivp)+bc(iad)*phs         11d9s21
                        end if                                            11d5s21
                       end do                                             11d5s21
                      end if                                             11d5s21
                     end do                                              11d5s21
                    end if                                               11d5s21
                    do ivp=ill,ihh                                       11d5s21
                     bc(jdtmp+ivp)=bc(jdtmp+ivp)+sum                     11d5s21
                    end do                                               11d5s21
                    call dgemm('n','n',ncsfb(iarg),nrow,ncsfk(jarg),     11d9s21
     $                  1d0,bc(imatq),ncsfb(iarg),bc(ivtrans),          11d9s21
     $                  ncsfk(jarg),0d0,bc(itmp),ncsfb(iarg),           11d9s21
     d' hcssbk4.  1')
                    itmptt=ibcoff                                        2d2s22
                    ibcoff=itmptt+ncsfb(iarg)*nrow                       2d2s22
                    call enough('hcssbk4.  9',bc,ibc)
                    do i=0,ncsfb(iarg)-1                                 2d2s22
                     do iv=0,nrow-1                                      2d2s22
                      iiv=itmp+i+ncsfb(iarg)*iv                          2d2s22
                      ivi=itmptt+iv+nrow*i                               2d2s22
                      bc(ivi)=bc(iiv)                                    2d2s22
                     end do                                              2d2s22
                    end do                                               2d2s22
                    do i=0,ncsfb(iarg)-1                                 11d5s21
                     do ir=0,nrootm                                      11d5s21
                      iadg=jgg-1+nvirt(isbv)*(ir+nrootu*i)                 7d15s21
                      jtmp=itmptt-i2s+nherev*(ir+nrootu*i)               2d22s22
                      do ivp=ill,ihh                                       7d15s21
                       bc(iadg+ivp)=bc(iadg+ivp)                         11d5s21
     $                     +bc(jdtmp+ivp)*bc(jtmp+ivp)                  11d5s21
                      end do                                             11d5s21
                     end do                                              11d5s21
                    end do                                               11d5s21
                    ibcoff=itmptt                                        2d2s22
                   end do                                                   9d20s21
                  end if                                                    9d20s21
                  ibcoff=ibcb4                                           11d5s21
                 end if                                                     9d13s21
                 if(idelta.le.4.and.l2e.eq.0)then                        2d18s22
c
c     intermediate state = bra or ket
c
                  nx1=0                                                     9d15s21
                  do i=1,norb                                               9d15s21
                   if(btest(i1o,i))then                                    9d15s21
                    nx1=nx1+1                                               9d15s21
                    imap1(nx1)=i                                            9d15s21
                    icode1(nx1)=3                                           9d15s21
                   end if                                                   9d15s21
                  end do                                                    9d15s21
                  nx1=nx1+1                                              11d5s21
                  imap1(nx1)=norbx                                       11d5s21
                  icode1(nx1)=3                                          11d5s21
                  ibail=0                                               3d14s23
                  if(nrow.gt.ncsfk(jarg).and.                           3d14s23
     $                 min(ncsfb(iarg),ncsfk(jarg)).gt.1)then           3d14s23
                   ibail=1                                              3d14s23
                  end if                                                3d14s23
                  call spinloop1(i2sb,i2smb,i2sk,i2smk,nopeni,           11d9s21
     $                ncsfb(iarg),imap1,nrow,itypeg,imatg,ntypeg,irw0,  11d5s21
     $                bc(ivtrans),ncsfk(jarg),ieoro,bc,ibc,bc,ibc)      11d14s22
                  do i=0,ntypeg-1
                   ipack=ibc(itypeg+i)
                   imat=imatg+nrowi*i                                    11d5s21
                   do j=1,4                                                 10d20s21
                    if(ipack2(j).gt.0)then                                  10d20s21
                     if(ipack2(j).le.norb)then                           11d5s21
                      isy(j)=ism(ipack2(j))                                  10d20s21
                      igya(j)=irel(ipack2(j))-1                               10d20s21
                     else                                                11d5s21
                      isy(j)=isbv                                        11d5s21
                      igya(j)=-1                                         11d5s21
                     end if                                              11d5s21
                     imy(j)=0                                               10d20s21
                    else                                                    10d20s21
                     if(-ipack2(j).le.norb)then                          11d5s21
                      isy(j)=ism(-ipack2(j))                                 10d20s21
                      igya(j)=irel(-ipack2(j))-1                             10d20s21
                     else                                                11d5s21
                      isy(j)=isbv                                        11d5s21
                      igya(j)=-1                                         11d5s21
                     end if                                              11d5s21
                     imy(j)=1                                               10d20s21
                    end if                                                  10d20s21
                   end do                                                   10d20s21
                   phs=0.5d0                                             11d5s21
                   ispincase=0                                           11d9s21
                   if(igya(1).lt.0.and.igya(2).lt.0)then                 11d9s21
                    if(igya(3).lt.0.and.igya(4).lt.0)then                11d9s21
                     go to 2021                                          11d9s21
                    end if                                               11d9s21
                    ispincase=1                                          11d9s21
                   else if(igya(3).lt.0.and.igya(4).lt.0)then            11d9s21
                    ispincase=-1                                         11d9s21
                   else if                                              2d7s23
     $                  (min(igya(1),igya(2),igya(3),igya(4)).ge.0)then 2d7s23
                   else                                                  11d9s21
                    write(6,*)                                           11d9s21
     $                 ('unknown igya case after spinloop1 in hcssbk4'),11d9s21
     $                  igya                                            11d9s21
                    call dws_synca                                       11d9s21
                    call dws_finalize                                    11d9s21
                    stop                                                 11d9s21
                   end if                                                11d9s21
                   if(ispincase.ge.0)then                                11d9s21
                    ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))              10d20s21
                    jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)              2d7s23
     $                   +2*(1-imy(4))))                                2d7s23
                   else                                                  11d5s21
                    ltest=imy(2)+2*(imy(1)+2*(imy(4)+2*imy(3)))              10d20s21
                    jtest=1-imy(2)+2*(1-imy(1)+2*(1-imy(4)              2d7s23
     $                   +2*(1-imy(3))))                                2d7s23
                    if(isopt(2).ne.0)phs=-phs                            11d5s21
                   end if                                                11d5s21
                   itestp=ltest+1                                           10d20s21
                   jtestp=jtest+1                                           10d20s21
                   iuse=-1                                               11d5s21
                   if(iifmx(itestp).ge.0)then                               10d20s21
                    iuse=iifmx(itestp)                                   11d5s21
                   else if(iifmx(jtestp).ge.0)then                       11d5s21
                    iuse=iifmx(jtestp)                                   11d5s21
                    if(isopt(4).ne.0)phs=-phs                            11d5s21
                   end if                                                11d5s21
                   if(iuse.ge.0)then                                     11d5s21
                    if(ispincase.eq.0)then                               11d9s21
                     iad=i4or(isy(2),isy(3),isy(4))+igya(1)             2d7s23
     $                   +irefo(isy(1))                                 2d7s23
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d20s21
     $                (igya(4)+irefo(isy(4))*iuse)))                    11d5s21
                     xint=bc(iad)*phs                                    11d5s21
                     do j=0,nrowim                                       11d5s21
                      bc(idend+j)=bc(idend+j)+xint*bc(imat+j)            11d5s21
                     end do                                              11d5s21
                    else if(ispincase.eq.1)then                          11d9s21
                     icolj=igya(3)+irefo(isy(3))*(igya(4)                11d9s21
     $                   +irefo(isy(4))*iuse)                           11d9s21
                     ibc(ndenjv(isy(3))+icolj)=1                         11d9s21
                     jden=idenjv(isy(3))+nrowi*icolj                     11d9s21
                     do j=0,nrowim                                       11d5s21
                      bc(jden+j)=bc(jden+j)+bc(imat+j)*phs               11d5s21
                     end do                                              11d5s21
                    else if(ispincase.eq.-1)then                         11d9s21
                     icolj=igya(1)+irefo(isy(1))*(igya(2)                11d9s21
     $                   +irefo(isy(2))*iuse)                           11d9s21
                     ibc(ndenjv(isy(1))+icolj)=1                         11d9s21
                     jden=idenjv(isy(1))+nrowi*icolj                     11d9s21
                     do j=0,nrowim                                       11d5s21
                      bc(jden+j)=bc(jden+j)+bc(imat+j)*phs               11d5s21
                     end do                                              11d5s21
                    end if                                               11d5s21
                   end if                                                11d5s21
 2021              continue                                              11d9s21
                  end do                                                    10d20s21
                  ibcoff=itypeg
                  do i1=1,norb                                           11d5s21
                   if(btest(i1o,i1))then                                   9d15s21
                    do i2=i1+1,norbx                                             9d13s21
                     if(btest(i1o,i2))then                                 9d15s21
                      itestc=i1c                                              5d7s21
                      itesto=i1o                                              5d7s21
                      if(btest(itestc,i1))then                               9d13s21
                       itestc=ibclr(itestc,i1)                               9d13s21
                       itesto=ibset(itesto,i1)                               9d13s21
                       nopenkk=nopeni+1                                           1d22s21
                      else if(btest(itesto,i1))then                          9d13s21
                       itesto=ibclr(itesto,i1)                               9d13s21
                       nopenkk=nopeni-1                                              11d13s20
                      else                                                          11d13s20
                       write(6,*)('bit not set for i1 = '),i1
                       stop 'nab4(1,1)'                                          11d27s20
                      end if                                                        11d13s20
                      if(btest(itesto,i2))then                               9d13s21
                       itestc=ibset(itestc,i2)                               9d13s21
                       itesto=ibclr(itesto,i2)                               9d13s21
                       nopenkk=nopenkk-1                                              11d13s20
                      else if(btest(itestc,i2))then                         9d13s21
                       go to 1848
                      else                                                          11d13s20
                       itesto=ibset(itesto,i2)                               9d13s21
                       nopenkk=nopenkk+1                                              11d13s20
                      end if                                                        11d13s20
                      call gandcr(i1c,i1o,itestc,itesto,nopeni,            9d13s21
     $               nopenkk,norbx,nnot1,nab1,icode1,imap1,nx1,irw1,     9d16s21
     $                 irw2,iwpb1,iwpk1,bc,ibc)                         11d14s22
                      call gandcr(itestc,itesto,j1c,j1o,nopenkk,nopenj,        5d7s21
     $                 norbx,nnot2,nab2,icode2,imap2,nx2,irw1,irw2,      9d16s21
     $                 iwpb2,iwpk2,bc,ibc)                              11d14s22
                      call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,    10d19s21
     $                 nopenkk,ncsfb(iarg),nrow,itype,imatx,ntypeg,nab1,11d5s21
     $                 iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(ivtrans),        11d5s21
     $                 ncsfk(jarg),ieoro,bc,ibc)                        11d14s22
                      do i=0,ntypeg-1                                    11d5s21
                       ipack=ibc(itype+i)
                       imat=imatx+nrowi*i                                 11d5s21
                       do j=1,4                                               9d9s21
                        if(ipack2(j).gt.0)then                                9d9s21
                         if(ipack2(j).le.norb)then                        11d5s21
                          isy(j)=ism(ipack2(j))                                9d9s21
                          igya(j)=irel(ipack2(j))-1                             9d9s21
                         else                                             11d5s21
                          isy(j)=isbv                                     11d5s21
                          igya(j)=-1                                      11d5s21
                         end if                                           11d5s21
                         imy(j)=0                                             10d13s21
                        else                                                  9d9s21
                         if(-ipack2(j).le.norb)then                       11d5s21
                          isy(j)=ism(-ipack2(j))                                9d9s21
                          igya(j)=irel(-ipack2(j))-1                           10d13s21
                         else                                             11d5s21
                          isy(j)=isbv                                     11d5s21
                          igya(j)=-1                                      11d5s21
                         end if                                           11d5s21
                         imy(j)=1                                             10d13s21
                        end if                                                9d9s21
                       end do                                                 9d9s21
                       phs=1d0                                            11d5s21
                       if(imy(2).ge.0)then                                11d5s21
                        ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))       11d5s21
                        jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*         11d5s21
     $                     (1-imy(4))))                                 11d5s21
                       else                                               11d5s21
                        ltest=imy(2)+2*(imy(1)+2*(imy(4)+2*imy(3)))       11d5s21
                        jtest=1-imy(2)+2*(1-imy(1)+2*(1-imy(4)+2*         11d5s21
     $                     (1-imy(3))))                                 11d5s21
                        if(isopt(2).ne.0)phs=-phs                         11d5s21
                       end if                                             11d5s21
                       itestp=ltest+1                                         10d13s21
                       jtestp=jtest+1                                         10d13s21
                       iuse=-1                                            11d5s21
                       if(iifmx(itestp).ge.0)then                             10d13s21
                        iuse=iifmx(itestp)                                11d5s21
                       else if(iifmx(jtestp).ge.0)then                        10d13s21
                        iuse=iifmx(jtestp)                                11d5s21
                        if(isopt(4).ne.0)phs=-phs                         11d5s21
                       end if                                             11d5s21
                       if(iuse.ge.0)then                                  11d5s21
                        if(min(igya(1),igya(2),igya(3),igya(4)).ge.0)   2d7s23
     $                      then                                        2d7s23
                         iad=i4or(isy(2),isy(3),isy(4))+igya(1)           11d5s21
     $                       +irefo(isy(1))                               11d5s21
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d20s21
     $                (igya(4)+irefo(isy(4))*iuse)))                    11d5s21
                         xint=bc(iad)*phs                                    11d5s21
                         do j=0,nrowim                                       11d5s21
                          bc(idend+j)=bc(idend+j)+xint*bc(imat+j)            11d5s21
                         end do                                              11d5s21
                        else if(igya(1).ge.0)then                            11d5s21
                         icolk=igya(1)+irefo(isy(1))*(igya(4)                11d5s21
     $                        +irefo(isy(4))*iuse)                           11d5s21
                         ibc(ndenkv(isy(1))+icolk)=1                         11d5s21
                         jden=idenkv(isy(1))+nrowi*icolk                     11d5s21
                         do j=0,nrowim                                       11d5s21
                          bc(jden+j)=bc(jden+j)+bc(imat+j)*phs               11d5s21
                         end do                                              11d5s21
                        else if(igya(2).ge.0)then                            11d5s21
                         icolk=igya(2)+irefo(isy(2))*(igya(3)                11d5s21
     $                       +irefo(isy(3))*iuse)                           11d5s21
                         ibc(ndenkv(isy(2))+icolk)=1                         11d5s21
                         jden=idenkv(isy(2))+nrowi*icolk                     11d5s21
                         do j=0,nrowim                                       11d5s21
                          bc(jden+j)=bc(jden+j)+bc(imat+j)*phs               11d5s21
                         end do                                              11d5s21
                        end if                                               11d5s21
                       end if                                             11d5s21
                      end do
                      ibcoff=itype                                            9d6s21
 1848                 continue
                     end if
                    end do                                                   9d13s21
                   end if                                                   9d13s21
                  end do                                                    9d13s21
                 end if                                                  11d7s21
                else if(ndifs.eq.2.and.ndifb.eq.2)then
                 nnot=2                                                 2d8s23
                 do i=1,norbx                                           2d6s23
                  if(btest(gandco,i))then
                   if((btest(i1c,i).and.btest(j1o,i)).or.                 10d21s22
     $             (btest(i1o,i).and..not.btest(j1c,i)))then            10d21s22
                    nab4(1,1)=i
                   else
                    nab4(2,1)=i
                   end if
                  end if
                 end do
                 nok=0                                                         11d13s20
                 do i=1,norbx                                                   11d25s20
                  itest(i,2)=0                                                 11d25s20
                 end do                                                        11d25s20
                 ii=1                                                          11d25s20
                 do i=1,norb                                                   11d25s20
                  if(btest(j1c,i))then                                         11d25s20
                   itest(i,2)=2                                                11d25s20
                   ii=ii+1                                                     11d25s20
                  end if                                                       11d25s20
                 end do                                                        11d25s20
                 ii=1                                                          11d25s20
                 do i=1,norb                                             11d7s21
                  if(btest(j1o,i))then                                         11d25s20
                   itest(i,2)=1                                                11d25s20
                   ii=ii+1                                                     11d25s20
                  end if                                                       11d25s20
                 end do                                                        11d25s20
                 do i=1,norb                                             11d7s21
                  ixn=min(itest(i,1),itest(i,2))
                  if(ixn.gt.0)then                                             11d13s20
                   nok=nok+1                                                   11d13s20
                   itest(nok,3)=ixn                                            11d13s20
                   itest(nok,4)=i                                         11d5s21
                  end if                                                       11d13s20
                 end do                                                        11d13s20
                 call gandcr(i1c,i1o,j1c,j1o,nopeni,                     11d7s21
     $          nopenj,norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,     11d5s21
     $            iwpb1,iwpk1,bc,ibc)                                   11d14s22
                 call gencup(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,nab1,
     $               iwpb1,iwpk1,ioutg,imatg,ntypeg,mcsf,bc(ivtrans),   11d7s21
     $               ncsfk(jarg),nrow,bc,ibc)                           11d14s22
                 ntype1a=ntypeg                                             10d12s21
                 iouta=ioutg                                                10d12s21
                 imata=imatg                                                10d12s21
                 do ii=0,ntype1a-1                                          9d10s21
                  ipackc=ibc(iouta+ii)                                      9d10s21
                  if(ipackc1(1).gt.0)then                                   9d10s21
                   lsb=ism(ipackc1(1))                                      10d12s21
                   lgb=irel(ipackc1(1))-1                                   10d12s21
                   igya(3)=lgb                                              10d12s21
                   imy(3)=0                                                 10d12s21
                  else
                   lsb=ism(-ipackc1(1))                                     10d12s21
                   igya(3)=irel(-ipackc1(1))-1                              10d12s21
                   lgb=irel(-ipackc1(1))-1+irefo(lsb)                       10d12s21
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
                   lgk=irel(-ipackc1(2))-1+irefo(lsk)                       10d12s21
                   imy(4)=1                                                 10d12s21
                  end if
                  h0=-2d0                                                   10d13s21
                  sum=0d0                                                   10d14s21
                  if(ih0n(lsb).gt.0)then                                    10d13s21
                   iad=ih0n(lsb)+igya(3)+nh0(lsb)*igya(4)                   10d14s21
                   h0=bc(iad)                                                10d13s21
                   if(imy(3).ne.0.and.isopt(4).ne.0)h0=-h0                   10d13s21
                   sum=h0                                                   10d14s21
                  else if(ih0n(lsk).gt.0)then                               10d13s21
                   iad=ih0n(lsk)+igya(4)+nh0(lsk)*igya(3)                   10d14s21
                   h0=bc(iad)                                                10d13s21
                   if(isopt(2).ne.0)h0=-h0                                  10d13s21
                   if(imy(4).ne.0.and.isopt(4).ne.0)h0=-h0                  10d13s21
                   sum=h0                                                   10d14s21
                  end if                                                    10d13s21
                  do i=1,norb                                               9d10s21
                   if(btest(i1c,i).and.btest(j1c,i).and.l2e.eq.0)then    2d18s22
                    js=ism(i)                                               9d10s21
                    jg=irel(i)-1                                            9d10s21
                    do jpass=0,1                                            9d10s21
                     ltest=jpass+2*(jpass+2*(imy(3)+2*imy(4)))           11d5s21
                     itestp=ltest+1                                         10d12s21
                     jtest=1-jpass+2*(1-jpass+2*(1-imy(3)+2*(1-imy(4))))    10d12s21
                     jtestp=jtest+1                                         10d12s21
                     xintj=-2d0                                             10d12s21
                     if(iifmx(itestp).ge.0)then                             10d12s21
                      iad=i4or(js,lsb,lsk)+jg+irefo(js)*(jg+irefo(js)*   11d7s21
     $                (igya(3)+irefo(lsb)*(igya(4)+irefo(lsk)*          10d12s21
     $                iifmx(itestp))))                                  10d12s21
                      xintj=bc(iad)                                         10d12s21
                      orig=sum
                      sum=sum+xintj                                         10d12s21
                     else if(iifmx(jtestp).ge.0)then                        10d12s21
                      iad=i4or(js,lsb,lsk)+jg+irefo(js)*(jg+irefo(js)*   11d7s21
     $                (igya(3)+irefo(lsb)*(igya(4)+irefo(lsk)*          10d12s21
     $                iifmx(jtestp))))                                  10d12s21
                      xintj=bc(iad)                                         10d12s21
                      if(isopt(4).ne.0)xintj=-xintj                         10d12s21
                      orig=sum
                      sum=sum+xintj                                         10d12s21
                     end if                                                 10d12s21
                     ltest=jpass+2*(imy(4)+2*(imy(3)+2*jpass))           11d5s21
                     itestp=ltest+1                                         10d12s21
                     jtest=1-jpass+2*(1-imy(4)+2*(1-imy(3)+2*(1-jpass)))    10d12s21
                     jtestp=jtest+1                                         10d12s21
                     xintj=-2d0                                             10d12s21
                     if(iifmx(itestp).ge.0)then                             10d12s21
                      iad=i4or(lsk,lsb,js)+jg+irefo(js)*(igya(4)         11d7s21
     $                +irefo(lsk)*(igya(3)+irefo(lsb)*(jg+irefo(js)*    11d7s21
     $                iifmx(itestp))))                                  10d12s21
                      xintk=bc(iad)                                         10d12s21
                      orig=sum
                      sum=sum-xintk                                         10d12s21
                     else if(iifmx(jtestp).ge.0)then                        10d12s21
                      iad=i4or(lsk,lsb,js)+jg+irefo(js)*(igya(4)         11d7s21
     $                 +irefo(lsk)*(igya(3)+irefo(lsb)*(jg+irefo(js)*   11d7s21
     $                iifmx(jtestp))))                                  10d12s21
                      xintk=bc(iad)                                         10d12s21
                      if(isopt(4).ne.0)xintk=-xintk                         10d12s21
                      orig=sum
                      sum=sum-xintk                                         10d12s21
                     end if                                                 10d12s21
                    end do                                                  9d10s21
                   end if                                                   9d10s21
                  end do                                                    9d10s21
                  if(ipackc1(1).gt.0)then                                   10d12s21
                   in=ipackc1(1)                                            10d12s21
                   imy(1)=1                                                 10d12s21
                  else                                                      9d10s21
                   in=-ipackc1(1)                                           10d12s21
                   imy(1)=0                                                 10d12s21
                  end if                                                    9d13s21
                  if(btest(i1c,in).and.l2e.eq.0)then                     2d18s22
                   ltest=imy(1)+2*(imy(1)+2*(imy(3)+2*imy(4)))              10d12s21
                   itestp=ltest+1                                           10d12s21
                   jtest=1-imy(1)+2*(1-imy(1)+2*(1-imy(3)+2*(1-imy(4))))    10d12s21
                   jtestp=jtest+1                                           10d12s21
                   xintj=-2d0                                               10d12s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iad=i4or(lsb,lsb,lsk)+igya(3)+irefo(lsb)*(igya(3)       10d12s21
     $              +irefo(lsb)*(igya(3)+irefo(lsb)*(igya(4)+irefo(lsk) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                    xintj=bc(iad)                                           10d12s21
                    orig=sum
                    sum=sum+xintj                                           10d13s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iad=i4or(lsb,lsb,lsk)+igya(3)+irefo(lsb)*(igya(3)       10d12s21
     $               +irefo(lsb)*(igya(3)+irefo(lsb)*(igya(4)           10d12s21
     $               +irefo(lsk)*iifmx(jtestp))))                       10d12s21
                    xintj=bc(iad)                                           10d12s21
                    if(isopt(4).ne.0)xintj=-xintj                           10d12s21
                    orig=sum
                    sum=sum+xintj                                           10d13s21
                   end if                                                   10d12s21
                   ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(1)))              10d12s21
                   itestp=ltest+1                                           10d12s21
                   jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(1))))    10d12s21
                   jtestp=jtest+1                                           10d12s21
                   xintk=-2d0                                               10d12s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iad=i4or(lsk,lsb,lsb)+igya(3)+irefo(lsb)*(igya(4)       10d12s21
     $              +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(3)+irefo(lsb) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                    xintk=bc(iad)                                           10d12s21
                    orig=sum
                    sum=sum-xintk                                           10d13s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iad=i4or(lsk,lsb,lsb)+igya(3)+irefo(lsb)*(igya(4)       10d12s21
     $               +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(3)           10d12s21
     $               +irefo(lsb)*iifmx(jtestp))))                       10d12s21
                    xintk=bc(iad)                                           10d12s21
                    if(isopt(4).ne.0)xintk=-xintk                           10d12s21
                    orig=sum
                    sum=sum-xintk                                           10d13s21
                   end if                                                   10d12s21
                  end if                                                    9d13s21
                  if(ipackc1(2).gt.0)then                                   10d12s21
                   im=ipackc1(2)                                            10d12s21
                   imy(2)=1                                                 10d13s21
                  else                                                      10d12s21
                   im=-ipackc1(2)                                           10d12s21
                   imy(2)=0                                                 10d13s21
                  end if                                                    10d12s21
                  if(btest(j1c,im).and.l2e.eq.0)then                     2d18s22
                   ltest=imy(2)+2*(imy(2)+2*(imy(3)+2*imy(4)))              10d12s21
                   itestp=ltest+1                                           10d12s21
                   jtest=1-imy(2)+2*(1-imy(2)+2*(1-imy(3)+2*(1-imy(4))))    10d12s21
                   jtestp=jtest+1                                           10d12s21
                   xintj=-2d0                                               10d12s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $              +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)+irefo(lsk) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                    xintj=bc(iad)                                           10d12s21
                    orig=sum
                    sum=sum+xintj                                           10d13s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $               +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)           10d12s21
     $               +irefo(lsk)*iifmx(jtestp))))                       10d12s21
                    xintj=bc(iad)                                           10d12s21
                    if(isopt(4).ne.0)xintj=-xintj                           10d12s21
                    orig=sum
                    sum=sum+xintj                                           10d13s21
                   end if                                                   10d12s21
                   ltest=imy(2)+2*(imy(4)+2*(imy(3)+2*imy(2)))              10d12s21
                   itestp=ltest+1                                           10d12s21
                   jtest=1-imy(2)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(2))))    10d12s21
                   jtestp=jtest+1                                           10d12s21
                   xintk=-2d0                                               10d12s21
                   if(iifmx(itestp).ge.0)then                               10d12s21
                    iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $              +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)+irefo(lsk) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                    xintk=bc(iad)                                           10d12s21
                    orig=sum
                    sum=sum-xintk                                           10d13s21
                   else if(iifmx(jtestp).ge.0)then                          10d12s21
                    iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $               +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)           10d12s21
     $               +irefo(lsk)*iifmx(jtestp))))                       10d12s21
                    xintk=bc(iad)                                           10d12s21
                    if(isopt(4).ne.0)xintk=-xintk                           10d12s21
                    orig=sum
                    sum=sum-xintk                                           10d13s21
                   end if                                                   10d12s21
                  end if                                                    9d13s21
                  imat=imata+nrowi*ii                                    11d5s21
                  do j=0,nrowim                                          11d5s21
                   bc(idend+j)=bc(idend+j)+sum*bc(imat+j)                11d7s21
                  end do                                                 11d5s21
                 end do                                                  11d5s21
                 ibcoff=ioutg                                            11d5s21
                 do i=1,nok                                                    11d13s20
                  if(itest(i,3).eq.1.and.l2e.eq.0)then                   2d18s22
                   itestc=i1c                                            11d7s21
                   itesto=i1o                                            11d7s21
                   nopenk=nopeni                                         11d7s21
c
c     anihilate common
c
                   if(btest(itestc,itest(i,4)))then                      11d5s21
                    itestc=ibclr(itestc,itest(i,4))                      11d5s21
                    itesto=ibset(itesto,itest(i,4))                      11d5s21
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,4))                      11d5s21
                    nopenk=nopenk-1                                             11d13s20
                   end if                                                       11d13s20
c
c     create ket
c
                   if(btest(itesto,nab4(2,1)))then                        12d14s20
                    itestc=ibset(itestc,nab4(2,1))                        12d14s20
                    itesto=ibclr(itesto,nab4(2,1))                        12d14s20
                    nopenk=nopenk-1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibset(itesto,nab4(2,1))                        12d14s20
                    nopenk=nopenk+1                                             11d13s20
                   end if                                                       11d13s20
                   call gandcr(i1c,i1o,itestc,itesto,nopeni,             11d5s21
     $              nopenk,norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2, 11d7s21
     $              iwpb1,iwpk1,bc,ibc)                                 11d14s22
                   call gandcr(itestc,itesto,j1c,j1o,nopenk,nopenj,      11d7s21
     $         norbx,nnot2,nab2,icode2,imap2,nx2,irw1,irw2,iwpb2,iwpk2, 11d14s22
     $                 bc,ibc)                                          11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                        9d10s21
                    call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,      10d19s21
     $               nopenk,ncsfb(iarg),nrow,itype,imatx,ntypeg,nab1,   11d7s21
     $                 iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(ivtrans),        11d5s21
     $                 ncsfk(jarg),ieoro,bc,ibc)                        11d14s22
                    do ini=0,ntypeg-1                                      11d5s21
                     ipack=ibc(itype+ini)                                11d7s21
                     imat=imatx+nrowi*ini                                11d7s21
                     do j=1,4                                               9d9s21
                      if(ipack2(j).gt.0)then                                9d9s21
                       if(ipack2(j).le.norb)then                         11d5s21
                        isy(j)=ism(ipack2(j))                                9d9s21
                        igya(j)=irel(ipack2(j))-1                             9d9s21
                       else                                              11d5s21
                        igya(j)=-1                                       11d5s21
                       end if                                            11d5s21
                       imy(j)=0                                             10d13s21
                      else                                                  9d9s21
                       if(-ipack2(j).le.norb)then                        11d5s21
                        isy(j)=ism(-ipack2(j))                                9d9s21
                        igya(j)=irel(-ipack2(j))-1                           10d13s21
                       else                                              11d5s21
                        igya(j)=-1                                       11d5s21
                       end if                                            11d5s21
                       imy(j)=1                                             10d13s21
                      end if                                                9d9s21
                     end do                                                 9d9s21
                     phs=1d0                                             11d5s21
                     if(igya(2).ge.0)then                                11d5s21
                      ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))            10d13s21
                      jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*          11d5s21
     $                    (1-imy(4))))                                  11d5s21
                     else                                                11d5s21
                      ltest=imy(2)+2*(imy(1)+2*(imy(4)+2*imy(3)))        11d5s21
                      jtest=1-imy(2)+2*(1-imy(1)+2*(1-imy(4)+2*          11d5s21
     $                    (1-imy(3))))                                  11d5s21
                      if(isopt(2).ne.0)phs=-phs                          11d5s21
                     end if                                              11d5s21
                     itestp=ltest+1                                         10d13s21
                     jtestp=jtest+1                                         10d13s21
                     iuse=-1                                             11d5s21
                     if(iifmx(itestp).ge.0)then                             10d13s21
                      iuse=iifmx(itestp)                                 11d5s21
                     else if(iifmx(jtestp).ge.0)then                        10d13s21
                      iuse=iifmx(jtestp)                                 11d5s21
                      if(isopt(4).ne.0)phs=-phs                          11d5s21
                     end if                                                 10d13s21
                     if(iuse.ge.0)then                                   11d7s21
                      if(min(igya(1),igya(2),igya(3),igya(4)).ge.0)then   11d5s21
                       iad=i4or(isy(2),isy(3),isy(4))+igya(1)            11d5s21
     $                    +irefo(isy(1))*(igya(2)+irefo(isy(2))*(igya(3)11d5s21
     $                    +irefo(isy(3))*(igya(4)+irefo(isy(4))*iuse))) 11d5s21
                       xint=bc(iad)*phs                                  11d5s21
                       do j=0,nrowim                                     11d5s21
                        bc(idend+j)=bc(idend+j)+xint*bc(imat+j)          11d5s21
                       end do                                            11d5s21
                      else if(igya(1).lt.0)then                          11d5s21
                       icolk=igya(2)+irefo(isy(2))*(igya(3)              11d5s21
     $                     +irefo(isy(3))*iuse)                         11d5s21
                       ibc(ndenkv(isy(2))+icolk)=1                       11d5s21
                       jden=idenkv(isy(2))+nrowi*icolk                   11d5s21
                       do j=0,nrowim                                     11d5s21
                        bc(jden+j)=bc(jden+j)+phs*bc(imat+j)             11d5s21
                       end do                                            11d5s21
                      else if(igya(2).lt.0)then                          11d5s21
                       icolk=igya(1)+irefo(isy(1))*(igya(4)              11d5s21
     $                     +irefo(isy(4))*iuse)                         11d5s21
                       ibc(ndenkv(isy(1))+icolk)=1                       11d5s21
                       jden=idenkv(isy(1))+nrowi*icolk                   11d5s21
                       do j=0,nrowim                                     11d5s21
                        bc(jden+j)=bc(jden+j)+phs*bc(imat+j)             11d5s21
                       end do                                            11d5s21
                      end if                                             11d5s21
                     end if                                              11d5s21
                     phs=1d0                                             11d5s21
                     if(imy(2).ge.0)then                                 11d5s21
                      ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            10d13s21
                      jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)             11d5s21
     $                    +2*(1-imy(2))))                               11d5s21
                     else                                                11d5s21
                      ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))            10d13s21
                      jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)             11d5s21
     $                    +2*(1-imy(1))))                               11d5s21
                      if(isopt(2).ne.0)phs=-phs                          11d5s21
                     end if                                              11d5s21
                     itestp=ltest+1                                         10d13s21
                     jtestp=jtest+1                                         10d13s21
                     iuse=-1                                             11d5s21
                     if(iifmx(itestp).ge.0)then                             10d13s21
                      iuse=iifmx(itestp)                                 11d5s21
                      xintk=bc(iad)                                         10d13s21
                      ihit=1                                                10d13s21
                      xint=xint-xintk                                       10d13s21
                     else if(iifmx(jtestp).ge.0)then                        10d13s21
                      iuse=iifmx(jtestp)
                      if(isopt(4).ne.0)phs=-phs                          11d5s21
                     end if                                                 10d13s21
                     if(iuse.ge.0)then                                   11d5s21
                      if(min(igya(1),igya(2),igya(3),igya(4)).ge.0)then  11d5s21
                       iad=i4or(isy(4),isy(3),isy(2))+igya(1)            11d5s21
     $                    +irefo(isy(1))                                11d5s21
     $                *(igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(2)+irefo(isy(2))*iuse)))                    11d5s21
                       xint=-bc(iad)*phs                                 11d5s21
                       do j=0,nrowim                                     11d5s21
                        bc(idend+j)=bc(idend+j)+xint*bc(imat+j)          11d5s21
                       end do                                            11d5s21
                      else if(igya(1).ge.0)then                          11d5s21
                       icolj=igya(4)+irefo(isy(4))*(igya(1)              11d5s21
     $                      +irefo(isy(1))*iuse)                         11d5s21
                       ibc(ndenjv(isy(4))+icolj)=1                       11d5s21
                       jden=idenjv(isy(4))+nrowi*icolj                   11d5s21
                       do j=0,nrowim                                     11d5s21
                        bc(jden+j)=bc(jden+j)-bc(imat+j)*phs             11d5s21
                       end do                                            11d5s21
                      else if(igya(2).ge.0)then                          11d5s21
                       icolj=igya(3)+irefo(isy(3))*(igya(2)              11d5s21
     $                     +irefo(isy(2))*iuse)                         11d5s21
                       ibc(ndenjv(isy(3))+icolj)=1                       11d5s21
                       jden=idenjv(isy(3))+nrowi*icolj                   11d5s21
                       do j=0,nrowim                                     11d5s21
                        bc(jden+j)=bc(jden+j)-bc(imat+j)*phs             11d5s21
                       end do                                            11d5s21
                      end if                                             11d5s21
                     end if                                              11d5s21
                    end do                                               11d5s21
                    ibcoff=itype                                         11d5s21
                   end if                                                11d5s21
                  end if                                                   12d14s20
                 end do                                                    12d14s20
                 if(l2e.eq.0)then                                        2d18s22
                  itestc=i1c                                              11d7s21
                  itesto=i1o                                              11d7s21
                  nopenk=nopeni                                           11d7s21
c
c     anihilate common
c
                  if(btest(itestc,norbx))then                             11d13s20
                   itestc=ibclr(itestc,norbx)                             11d13s20
                   itesto=ibset(itesto,norbx)                             11d13s20
                   nopenk=nopenk+1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibclr(itesto,norbx)                             11d13s20
                   nopenk=nopenk-1                                             11d13s20
                  end if                                                       11d13s20
c
c     create ket
c
                  if(btest(itesto,nab4(2,1)))then                         12d14s20
                   itestc=ibset(itestc,nab4(2,1))                         12d14s20
                   itesto=ibclr(itesto,nab4(2,1))                         12d14s20
                   nopenk=nopenk-1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibset(itesto,nab4(2,1))                         12d14s20
                   nopenk=nopenk+1                                             11d13s20
                  end if                                                       11d13s20
                  call gandcr(i1c,i1o,itestc,itesto,nopeni,               11d5s21
     $              nopenk,norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2, 11d7s21
     $              iwpb1,iwpk1,bc,ibc)                                 11d14s22
                  call gandcr(itestc,itesto,j1c,j1o,nopenk,nopenj,        11d7s21
     $         norbx,nnot2,nab2,icode2,imap2,nx2,irw1,irw2,iwpb2,iwpk2, 11d14s22
     $                bc,ibc)                                           11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                        9d10s21
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,      10d19s21
     $               nopenk,ncsfb(iarg),nrow,itype,imatx,ntypeg,nab1,   11d7s21
     $                 iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(ivtrans),        11d5s21
     $                 ncsfk(jarg),ieoro,bc,ibc)                        11d14s22
                   do i=0,ntypeg-1                                        11d5s21
                    ipack=ibc(itype+i)
                    imat=imatx+nrowi*i                                    11d5s21
                    do j=1,4                                               9d9s21
                     if(ipack2(j).gt.0)then                                9d9s21
                      if(ipack2(j).le.norb)then                           11d5s21
                       isy(j)=ism(ipack2(j))                                9d9s21
                       igya(j)=irel(ipack2(j))-1                             9d9s21
                      else                                                11d5s21
                       igya(j)=-1                                         11d5s21
                      end if                                              11d5s21
                      imy(j)=0                                             10d13s21
                     else                                                  9d9s21
                      if(-ipack2(j).le.norb)then                          11d5s21
                       isy(j)=ism(-ipack2(j))                                9d9s21
                       igya(j)=irel(-ipack2(j))-1                            10d13s21
                      else                                                11d5s21
                       igya(j)=-1                                         11d5s21
                      end if                                              11d5s21
                      imy(j)=1                                             10d13s21
                     end if                                                9d9s21
                    end do                                                 9d9s21
                    phs=1d0                                               11d5s21
c
c     (v2|3v)=(3v|v2)=+/-(2v|v3)
                    if(igya(2).ge.0)then                                  11d5s21
                     ltest=imy(2)+2*(imy(1)+2*(imy(4)+2*imy(3)))          11d9s21
                     jtest=1-imy(2)+2*(1-imy(1)+2*(1-imy(4)+2*            11d9s21
     $                    (1-imy(3))))                                  11d9s21
                     if(isopt(2).ne.0)phs=-phs                            11d9s21
                    else                                                  11d5s21
                     ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))          11d9s21
                     jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*            11d9s21
     $                    (1-imy(4))))                                  11d9s21
                    end if                                                11d5s21
                    itestp=ltest+1                                         10d13s21
                    jtestp=jtest+1                                         10d13s21
                    iuse=-1                                               11d5s21
                    if(iifmx(itestp).ge.0)then                             10d13s21
                     iuse=iifmx(itestp)                                   11d5s21
                    else if(iifmx(jtestp).ge.0)then                        10d13s21
                     iuse=iifmx(jtestp)                                   11d5s21
                     if(isopt(4).ne.0)phs=-phs                            11d5s21
                    end if                                                 10d13s21
                    if(iuse.ge.0)then                                     11d7s21
                     if(min(igya(1),igya(2),igya(3),igya(4)).ge.0)then    11d5s21
                      iad=i4or(isy(2),isy(3),isy(4))+igya(1)              11d5s21
     $                    +irefo(isy(1))*(igya(2)+irefo(isy(2))*(igya(3)11d5s21
     $                    +irefo(isy(3))*(igya(4)+irefo(isy(4))*iuse))) 11d5s21
                      xint=bc(iad)*phs                                    11d5s21
                      do j=0,nrowim                                       11d5s21
                       bc(idend+j)=bc(idend+j)+xint*bc(imat+j)            11d5s21
                      end do                                              11d5s21
                     else if(igya(1).lt.0)then                            11d5s21
                      icolk=igya(2)+irefo(isy(2))*(igya(3)                11d5s21
     $                     +irefo(isy(3))*iuse)                         11d5s21
                      ibc(ndenkv(isy(2))+icolk)=1                         11d5s21
                      jden=idenkv(isy(2))+nrowi*icolk                     11d5s21
                      do j=0,nrowim                                       11d5s21
                       bc(jden+j)=bc(jden+j)+phs*bc(imat+j)               11d5s21
                      end do                                              11d5s21
                     else if(igya(2).lt.0)then                            11d5s21
                      icolk=igya(1)+irefo(isy(1))*(igya(4)                11d5s21
     $                     +irefo(isy(4))*iuse)                           11d5s21
                      ibc(ndenkv(isy(1))+icolk)=1                         11d5s21
                      jden=idenkv(isy(1))+nrowi*icolk                     11d5s21
                      do j=0,nrowim                                       11d5s21
                       bc(jden+j)=bc(jden+j)+phs*bc(imat+j)               11d5s21
                      end do                                              11d5s21
                     end if                                               11d5s21
                    end if                                                11d5s21
                    phs=1d0                                               11d5s21
                    if(imy(2).ge.0)then                                   11d5s21
                     ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            10d13s21
                     jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)               11d7s21
     $                    +2*(1-imy(2))))                               11d5s21
                    else                                                  11d7s21
                     ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))            10d13s21
                     jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)               11d7s21
     $                    +2*(1-imy(1))))                               11d5s21
                     if(isopt(2).ne.0)phs=-phs                            11d7s21
                    end if                                                11d7s21
                    itestp=ltest+1                                         10d13s21
                    jtestp=jtest+1                                         10d13s21
                    iuse=-1                                               11d5s21
                    if(iifmx(itestp).ge.0)then                             10d13s21
                     iuse=iifmx(itestp)                                   11d5s21
                    else if(iifmx(jtestp).ge.0)then                        10d13s21
                     iuse=iifmx(jtestp)
                     if(isopt(4).ne.0)phs=-phs                            11d5s21
                    end if                                                 10d13s21
                    if(iuse.ge.0)then                                     11d5s21
                     if(min(igya(1),igya(2),igya(3),igya(4)).ge.0)then    11d5s21
                      iad=i4or(isy(4),isy(3),isy(2))+igya(1)              11d5s21
     $                    +irefo(isy(1))                                11d5s21
     $                  *(igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(2)+irefo(isy(2))*iuse)))                    11d5s21
                      xint=-bc(iad)*phs                                   11d5s21
                      do j=0,nrowim                                       11d5s21
                       bc(idend+j)=bc(idend+j)+xint*bc(imat+j)            11d5s21
                      end do                                              11d5s21
                     else if(igya(1).ge.0)then                            11d5s21
                      icolj=igya(4)+irefo(isy(4))*(igya(1)                11d5s21
     $                     +irefo(isy(1))*iuse)                         11d5s21
                      ibc(ndenjv(isy(4))+icolj)=1                         11d5s21
                      jden=idenjv(isy(4))+nrowi*icolj                     11d5s21
                      do j=0,nrowim                                       11d5s21
                       bc(jden+j)=bc(jden+j)-bc(imat+j)*phs               11d5s21
                      end do                                              11d5s21
                     else if(igya(2).ge.0)then                            11d5s21
                      icolj=igya(3)+irefo(isy(3))*(igya(2)                11d5s21
     $                     +irefo(isy(2))*iuse)                         11d5s21
                      ibc(ndenjv(isy(3))+icolj)=1                         11d5s21
                      jden=idenjv(isy(3))+nrowi*icolj                     11d5s21
                      do j=0,nrowim                                       11d5s21
                       bc(jden+j)=bc(jden+j)-bc(imat+j)*phs               11d5s21
                      end do                                              11d5s21
                     end if                                               11d5s21
                    end if                                                11d5s21
                   end do                                                 11d5s21
                   ibcoff=itype                                           11d5s21
                  end if                                                  11d5s21
                 end if                                                  2d18s22
                else if(l2e.eq.0)then                                   2d7s23
                 nnot=0                                                 2d7s23
                 if(ndifs.eq.4.and.ndifb.eq.4)then                      2d6s23
                  nnot=4                                                2d6s23
                  ioxx(1)=1                                             2d6s23
                  ioxx(2)=1                                             2d6s23
                  do i=1,norbx                                          2d6s23
                   if(btest(gandcb,i))then                              2d6s23
                    if((btest(j1c,i).and.btest(i1o,i)).or.              2d8s23
     $                (btest(j1o,i).and..not.btest(i1c,i)))then         2d8s23
                     nab4(2,ioxx(2))=i                                  2d6s23
                     ioxx(2)=ioxx(2)+1                                  2d6s23
                    else                                                2d6s23
                     nab4(1,ioxx(1))=i                                  2d6s23
                     ioxx(1)=ioxx(1)+1                                  2d6s23
                    end if                                              2d6s23
                   end if                                               2d6s23
                  end do                                                2d6s23
                 else if(ndifb.eq.3)then                                2d6s23
                  nnot=3                                                2d6s23
                  ioxx(1)=1                                             2d6s23
                  ioxx(2)=1                                             2d6s23
                  iswap=0                                               2d6s23
                  do i=1,norbx                                          2d6s23
                   if(btest(gandcb,i))then                              2d6s23
                   if(btest(gandcc,i).and.                              2d6s23
     $                  ((btest(i1c,i).and..not.btest(j1o,i)).or.       2d6s23
     $                  (btest(j1c,i).and..not.btest(i1o,i))))then      2d6s23
                     if(btest(j1c,i))iswap=1                            2d6s23
                     nab4(1,1)=i                                        2d6s23
                     nab4(1,2)=i                                        2d6s23
                    else                                                2d6s23
                     nab4(2,ioxx(2))=i                                  2d6s23
                     ioxx(2)=ioxx(2)+1                                  2d6s23
                    end if                                              2d6s23
                   end if                                               2d6s23
                  end do                                                2d6s23
                  if(iswap.ne.0)then                                    2d6s23
                   icpy=nab4(1,1)                                       2d6s23
                   nab4(1,1)=nab4(2,1)                                  2d6s23
                   nab4(2,1)=icpy                                       2d6s23
                   icpy=nab4(1,2)                                       2d6s23
                   nab4(1,2)=nab4(2,2)                                  2d6s23
                   nab4(2,2)=icpy                                       2d6s23
                   nbt=0                                                2d6s23
                   if(btest(j1c,nab4(2,2)).and.                         2d8s23
     $                  .not.btest(j1c,nab4(2,1)))nbt=1                 2d8s23
                  else                                                  2d6s23
                   nbt=0                                                2d6s23
                   if(btest(i1c,nab4(1,2)).and.                         2d8s23
     $                  .not.btest(i1c,nab4(1,1)))nbt=1                 2d8s23
                  end if                                                2d6s23
                  if(nbt.ne.0)then                                      2d6s23
                   nab4(1,1)=nab4(1,2)                                  2d6s23
                   nab4(2,1)=nab4(2,2)                                  2d6s23
                  end if                                                2d6s23
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                 2d6s23
                  nnot=3                                                2d6s23
                  do i=1,norbx                                          2d6s23
                   if(btest(gandcb,i))then                              2d6s23
                    if(btest(i1c,i))then                                2d8s23
                     nab4(1,1)=i                                        2d6s23
                     nab4(1,2)=i                                        2d6s23
                    else                                                2d6s23
                     nab4(2,1)=i                                        2d6s23
                     nab4(2,2)=i                                        2d6s23
                    end if                                              2d6s23
                   end if                                               2d6s23
                  end do                                                2d6s23
                 end if                                                 2d6s23
                 if(nnot.ne.0)then                                      2d7s23
                  iu1=1                                                     1d22s21
                  iu2=1                                                     1d22s21
                  itestc=i1c                                              11d7s21
                  itesto=i1o                                              11d7s21
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopeni+1                                        11d7s21
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopeni-1                                        11d7s21
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(i1c,i1o,itestc,itesto,nopeni,               11d5s21
     $              nopenk,norbx,nnot1,nab1,icode1,imap1,nx1,irw1,irw2, 11d7s21
     $             iwpb1,iwpk1,bc,ibc)                                  11d14s22
                  call gandcr(itestc,itesto,j1c,j1o,nopenk,nopenj,        11d7s21
     $         norbx,nnot2,nab2,icode2,imap2,nx2,irw1,irw2,iwpb2,iwpk2, 11d14s22
     $               bc,ibc)                                            11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                         1d22s21
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopeni,nopenj,      10d19s21
     $               nopenk,ncsfb(iarg),nrow,itype,imatx,ntypeg,nab1,   11d7s21
     $               iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(ivtrans),          11d5s21
     $               ncsfk(jarg),ieoro,bc,ibc)                          11d14s22
                   do i=0,ntypeg-1                                        11d5s21
                    ipack=ibc(itype+i)
                    ltest=0                                                10d12s21
                    jtest=0                                                10d12s21
                    mtest=1                                                10d12s21
                    do j=1,4                                               9d9s21
                     if(ipack2(j).gt.0)then                                9d9s21
                      isy(j)=ism(ipack2(j))                                9d9s21
                      igya(j)=irel(ipack2(j))-1                             9d9s21
                      jtest=jtest+mtest                                    10d12s21
                      imy(j)=0                                             10d12s21
                     else                                                  9d9s21
                      isy(j)=ism(-ipack2(j))                                9d9s21
                      igya(j)=irel(-ipack2(j))-1                           10d12s21
                      ltest=ltest+mtest                                    10d12s21
                      imy(j)=1                                             10d12s21
                     end if                                                9d9s21
                     mtest=mtest*2                                         10d12s21
                    end do                                                 9d9s21
                    itestp=ltest+1                                         10d12s21
                    jtestp=jtest+1                                         10d12s21
                    ihit=0                                                11d5s21
                    xint=0d0                                              11d5s21
                    if(iifmx(itestp).ge.0)then                             10d12s21
                     iad=i4or(isy(2),isy(3),isy(4))+igya(1)              2d7s23
     $                  +irefo(isy(1))*                                 2d7s23
     $          (igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(4)+irefo(isy(4))*iifmx(itestp))))            10d12s21
                     xintj=bc(iad)                                          10d12s21
                     ihit=1                                                10d12s21
                     xint=xintj                                            10d12s21
                    else if(iifmx(jtestp).ge.0)then                        10d12s21
                     iad=i4or(isy(2),isy(3),isy(4))+igya(1)              2d7s23
     $                    +irefo(isy(1))*                                2d7s23
     $          (igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(4)+irefo(isy(4))*iifmx(jtestp))))            10d12s21
                     xintj=bc(iad)                                          10d12s21
                     if(isopt(4).ne.0)then                                 10d12s21
                      xintj=-xintj                                           10d12s21
                     end if
                     ihit=1                                                10d12s21
                     xint=xintj                                            10d12s21
                    end if                                                 10d12s21
                    if(nnot.eq.4)then                                      10d12s21
                     ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            10d12s21
                     itestp=ltest+1                                         10d12s21
                     jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)              2d7s23
     $                    +2*(1-imy(2))))                                2d7s23
                     jtestp=jtest+1                                         10d12s21
                     if(iifmx(itestp).ge.0)then                             10d12s21
                      iad=i4or(isy(4),isy(3),isy(2))+igya(1)              11d5s21
     $                  +irefo(isy(1))*                                 11d5s21
     $          (igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(2)+irefo(isy(2))*iifmx(itestp))))            10d12s21
                      xintk=bc(iad)                                          10d12s21
                      xint=xint-xintk                                       10d12s21
                      ihit=1                                              11d9s21
                     else if(iifmx(jtestp).ge.0)then                        10d12s21
                      iad=i4or(isy(4),isy(3),isy(2))+igya(1)
     $                   +irefo(isy(1))*                                11d5s21
     $          (igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(2)+irefo(isy(2))*iifmx(jtestp))))            10d12s21
                      xintk=bc(iad)                                          10d12s21
                      if(isopt(4).ne.0)then                                 10d12s21
                       xintk=-xintk                                           10d12s21
                      end if
                      xint=xint-xintk                                       10d12s21
                      ihit=1                                              11d9s21
                     end if                                                 10d12s21
                    end if                                                 10d12s21
                    if(ihit.ne.0)then                                      9d13s21
                     xint0=xint
                     if(nab4(1,1).eq.nab4(1,2).and.                      2d7s23
     $                   nab4(2,1).eq.nab4(2,2))                        2d7s23
     $                 xint=xint*0.5d0                                  9d9s21
                     itmpxq=imatx+nrowi*i                                 11d5s21
                     do j=0,nrowim                                        11d5s21
                      bc(idend+j)=bc(idend+j)+xint*bc(itmpxq+j)           11d5s21
                     end do                                               11d5s21
                    end if                                                 9d9s21
                   end do
                   ibcoff=itype                                            9d6s21
                  end if                                                    12d14s20
                 end if                                                 2d7s23
                end if                                                  2d7s23
               end if                                                   11d5s21
              end if                                                     12d14s20
              jvec=jvec+ncsfk(jarg)*nrow                                11d5s21
             end do                                                      12d14s20
             if(itransgg.eq.0)then                                      2d12s21
              call ddi_done(ibc(iacc),nacc)                             2d12s21
              itransgg=1                                                2d12s21
              do i=igg,igg+ncolti-1                                     2d12s21
               bc(i)=0d0                                                2d12s21
              end do                                                    2d12s21
             end if                                                     2d12s21
             if(isbv.eq.jsbv)then                                       7d19s21
              ltmp1=ibcoff                                              2d22s21
              ibcoff=ltmp1+nrow*ncsfb(iarg)                             11d5s21
              call enough('hcssbk4. 10',bc,ibc)
              do j=0,nrow-1                                             2d22s21
               do ii=0,ncsfb(iarg)-1                                    11d5s21
                ij=idend+ii+ncsfb(iarg)*j                               11d5s21
                ji=ltmp1+j+nrow*ii                                      2d22s21
                bc(ji)=bc(ij)                                           2d22s21
               end do                                                   2d22s21
              end do                                                    2d22s21
              do i=0,ncsfb(iarg)-1                                      11d5s21
               do ir=0,nrootm                                           1d27s21
                iadg=jgg-1+nvirt(isbv)*(ir+nrootu*i)                    1d27s21
                iadd=ltmp1-i2s+nherev*(ir+nrootu*i)                     2d22s21
                do iv=i2s,i2top                                           12d15s20
                 iadgv=iadg+iv                                          7d19s21
                 iaddv=iadd+iv                                          7d19s21
                 bc(iadgv)=bc(iadgv)+bc(iaddv)                          7d19s21
                end do                                                  12d15s20
               end do                                                   12d15s20
              end do                                                    12d15s20
              ibcoff=ltmp1                                              2d22s21
             end if                                                     12d15s20
             ltmp1=ibcoff                                               2d22s21
             ibcoff=ltmp1+nrowi                                         2d22s21
             call enough('hcssbk4. 11',bc,ibc)
             do j=0,nrow-1                                              2d22s21
              do ii=0,ncsfb(iarg)-1                                      2d22s21
               ij=idenh+ii+ncsfb(iarg)*j                                 2d22s21
               ji=ltmp1+j+nrow*ii                                       2d22s21
               bc(ji)=bc(ij)                                            2d22s21
              end do                                                    2d22s21
             end do                                                     2d22s21
             do j=0,nrowim                                              2d22s21
              bc(idenh+j)=bc(ltmp1+j)                                   2d22s21
             end do                                                     2d22s21
             do isc=1,nsymb                                             2d22s21
              isd=multh(isymc,multh(isc,ijsbv))                         11d5s21
              nn=irefo(isc)*irefo(isd)*ntype                            11d5s21
              do i=0,nn-1                                               2d22s21
               if(ibc(ndenk(isc)+i).ne.0.and.l2e.eq.0)then              2d18s22
                iden=idenk(isc)+nrowi*i                                 11d5s21
                do j=0,nrow-1                                              2d22s21
                 do ii=0,ncsfb(iarg)-1                                      2d22s21
                  ij=iden+ii+ncsfb(iarg)*j                                 2d22s21
                  ji=ltmp1+j+nrow*ii                                       2d22s21
                  bc(ji)=bc(ij)                                            2d22s21
                 end do                                                    2d22s21
                end do                                                     2d22s21
                do j=0,nrowim                                              2d22s21
                 bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                end do                                                     2d22s21
               end if                                                   2d22s21
               if(ibc(ndenkv(isc)+i).ne.0.and.l2e.eq.0)then             2d18s22
                iden=idenkv(isc)+nrowi*i                                 2d22s21
                do j=0,nrow-1                                              2d22s21
                 do ii=0,ncsfb(iarg)-1                                  11d5s21
                  ij=iden+ii+ncsfb(iarg)*j                              11d5s21
                  ji=ltmp1+j+nrow*ii                                       2d22s21
                  bc(ji)=bc(ij)                                            2d22s21
                 end do                                                    2d22s21
                end do                                                     2d22s21
                do j=0,nrowim                                              2d22s21
                 bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                end do                                                     2d22s21
               end if                                                   2d22s21
              end do                                                    2d22s21
              do i=0,nn-1                                               2d22s21
               if(ibc(ndenj(isc)+i).ne.0.and.l2e.eq.0)then              2d18s22
                iden=idenj(isc)+nrowi*i                                 11d5s21
                do j=0,nrow-1                                              2d22s21
                 do ii=0,ncsfb(iarg)-1                                      2d22s21
                  ij=iden+ii+ncsfb(iarg)*j                                 2d22s21
                  ji=ltmp1+j+nrow*ii                                       2d22s21
                  bc(ji)=bc(ij)                                            2d22s21
                 end do                                                    2d22s21
                end do                                                     2d22s21
                do j=0,nrowim                                              2d22s21
                 bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                end do                                                     2d22s21
               end if                                                   2d22s21
               if(ibc(ndenjv(isc)+i).ne.0.and.l2e.eq.0)then             2d18s22
                iden=idenjv(isc)+nrowi*i                                 2d22s21
                do j=0,nrow-1                                              2d22s21
                 do ii=0,ncsfb(iarg)-1                                      2d22s21
                  ij=iden+ii+ncsfb(iarg)*j                                 2d22s21
                  ji=ltmp1+j+nrow*ii                                       2d22s21
                  bc(ji)=bc(ij)                                            2d22s21
                 end do                                                    2d22s21
                end do                                                     2d22s21
                do j=0,nrowim                                              2d22s21
                 bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                end do                                                     2d22s21
               end if                                                   2d22s21
              end do                                                    2d22s21
             end do                                                     2d22s21
             ibcoff=ltmp1                                               2d22s21
             if(i2s.eq.i2e)then                                         12d18s20
              n1here=i1e+1-i1s                                          12d18s20
              ii10=i1s                                                  12d18s20
             else                                                       12d18s20
              n1here=nvirt(isbv)                                        12d18s20
              ii10=1                                                    12d18s20
             end if                                                     12d18s20
             mcol=nrootu*ncsfb(iarg)                                     12d17s20
             itmpp=ibcoff                                               12d17s20
             ibcoff=itmpp+n1here*mcol                                   12d18s20
             call enough('hcssbk4. 12',bc,ibc)
             do i=itmpp,ibcoff-1                                        12d17s20
              bc(i)=0d0                                                 12d17s20
             end do                                                     12d17s20
             ihit=0
             if(ndenh.ne.0)then                                         12d15s20
c     g r,v,i = sum v' denh r,v',i h0v v'
              itmp=ibcoff                                               12d17s20
              ibcoff=itmp+n1here*nherev                                 12d18s20
              call enough('hcssbk4. 13',bc,ibc)
              do i=itmp,ibcoff-1                                        12d17s20
               bc(i)=0d0                                                12d17s20
              end do                                                    12d17s20
              i10=i1s                                                   12d15s20
              i1n=nvirt(isbv)                                           12d15s20
              do i2=i2s,i2e                                             12d15s20
               i2m=i2-1                                                 12d17s20
               i2p=i2+1                                                 12d17s20
               if(i2.eq.i2e)i1n=i1e                                     12d15s20
               jtmp=itmp-ii10+n1here*(i2-i2s)                           12d18s20
               if(isymc.eq.1)then                                       11d7s21
                ih0=ih0n(isbv)+irefo(isbv)-1                             11d5s21
     $              +nh0(isbv)*(i2m+irefo(jsbv))                        11d5s21
                do i1=i10,min(i1n,i2m)                                   12d17s20
                 bc(jtmp+i1)=bc(ih0+i1)                                 11d5s21
                end do                                                   12d17s20
                do i1=max(i2p,i10),i1n                                   12d17s20
                 bc(jtmp+i1)=bc(ih0+i1)                                 11d5s21
                end do                                                   12d17s20
               else                                                     7d19s21
                if(ih0n(isbv).gt.0)then                                  11d7s21
                 ih0=ih0n(isbv)+irefo(isbv)-1                             11d5s21
     $              +nh0(isbv)*(i2m+irefo(jsbv))                        11d5s21
                 do i1=i10,i1n                                           7d19s21
                  bc(jtmp+i1)=bc(ih0+i1)                                 11d5s21
                 end do                                                   12d17s20
                else if(ih0n(jsbv).gt.0)then                            11d7s21
                 ih0=ih0n(jsbv)+i2m+irefo(jsbv)                         11d7s21
     $              +nh0(jsbv)*(irefo(isbv)-1)                          11d7s21
                 do i1=i10,i1n                                           7d19s21
                  bc(jtmp+i1)=bc(ih0+i1*nh0(jsbv))                      11d7s21
                 end do                                                   12d17s20
                end if                                                  11d7s21
               end if                                                   7d19s21
               i10=1                                                    12d15s20
              end do                                                    12d15s20
              itmpd=ibcoff                                              12d17s20
              ibcoff=itmpd+nrowi                                        12d17s20
              call enough('hcssbk4. 14',bc,ibc)
              jtmp=idenh                                                12d17s20
              do i=0,ncsfb(iarg)-1                                      11d5s21
               do ir=0,nrootm                                           1d28s21
                iad=itmpd+nherev*(ir+nrootu*i)                          1d28s21
                do iv=0,nherev-1                                        1d28s21
                 bc(iad+iv)=bc(jtmp+iv)                                 1d28s21
                end do                                                  1d28s21
                jtmp=jtmp+nherev                                        1d28s21
               end do                                                   12d17s20
              end do                                                    12d17s20
              ihit=1
              if(ih0n(isbv).gt.0)then                                   11d9s21
               phasf=1d0                                                11d9s21
              else                                                      11d9s21
               phasf=1d0                                                11d9s21
               if(isopt(2).ne.0)phasf=-phasf                            11d9s21
               if(isopt(3).ne.0.and.isopt(4).ne.0)phasf=-phasf          11d9s21
              end if                                                    11d9s21
              call dgemm('n','n',n1here,mcol,nherev,phasf,              11d7s21
     $             bc(itmp),n1here,bc(itmpd),nherev,1d0,                12d18s20
     $             bc(itmpp),n1here,                                    12d18s20
     d' hcssbk4.  2')
              ibcoff=itmp
             end if                                                     12d15s20
             do isc=1,nsymb                                             12d15s20
              isd=multh(isymc,multh(isc,ijsbv))                         11d5s21
              nnp=irefo(isc)*irefo(isd)                                 11d5s21
              nn=nnp*ntype                                              11d5s21
              ikeep=ibcoff                                              11d5s21
              ibcoff=ikeep+nn                                           11d5s21
              call enough('hcssbk4. 15',bc,ibc)
              nokk=0                                                    11d5s21
              nokkv=0                                                   11d5s21
              do i=0,nn-1                                               11d5s21
               ibc(ikeep+i)=0                                           11d5s21
               if(ibc(ndenkv(isc)+i).ne.0.and.l2e.eq.0)then             2d18s22
                iden=idenkv(isc)+nrowi*i                                 11d5s21
                rms=1d0                                                  11d5s21
                if(ibc(ndenk(isc)+i).ne.0)then                           11d5s21
                 jden=idenk(isc)+nrowi*i                                 11d5s21
                 rms=0d0                                                 11d5s21
                 do k=0,nrowim                                           11d5s21
                  rms=rms+(bc(jden+k)-bc(iden+k))**2                     11d5s21
                 end do                                                  11d5s21
                 rms=sqrt(rms/dfloat(nrowi))                             11d5s21
                end if                                                   11d5s21
                if(rms.gt.1d-10)then                                     11d5s21
                 jden=idenkv(isc)+nrowi*nokkv                            11d5s21
                 do k=0,nrowim                                           11d5s21
                  bc(jden+k)=bc(iden+k)                                  11d5s21
                 end do                                                  11d5s21
                 ibc(ndenkv(isc)+nokkv)=i                                11d5s21
                 nokkv=nokkv+1                                           11d5s21
                else                                                     11d5s21
                 ibc(ikeep+i)=1                                          11d5s21
                end if                                                   11d5s21
               end if                                                    11d5s21
               if(ibc(ndenk(isc)+i).ne.0.and.l2e.eq.0)then              2d18s22
                iden=idenk(isc)+nrowi*i                                  11d5s21
                jden=idenk(isc)+nrowi*nokk                               11d5s21
                do k=0,nrowim                                            11d5s21
                 bc(jden+k)=bc(iden+k)                                   11d5s21
                end do                                                   11d5s21
                ibc(ndenk(isc)+nokk)=i                                   11d5s21
                nokk=nokk+1                                              11d5s21
               end if                                                    11d5s21
              end do                                                    11d7s21
              if(nokk.gt.0)then                                         11d5s21
               itmp=ibcoff                                               12d17s20
               ibcoff=itmp+nokk*n1here*nherev                           11d5s21
               call enough('hcssbk4. 16',bc,ibc)
               ihit=1                                                   12d17s20
               do i=itmp,ibcoff-1                                       12d17s20
                bc(i)=0d0                                               12d17s20
               end do                                                   12d17s20
               do i=0,nokk-1                                             12d17s20
                icol=ibc(ndenk(isc)+i)                                  12d17s20
                icol0=icol                                              11d5s21
                iq=icol/nnp                                             11d5s21
                icol=icol-nnp*iq                                        11d5s21
                id=icol/irefo(isc)                                      7d12s21
                ic=icol-irefo(isc)*id                                   7d12s21
                i10=i1s                                                   12d15s20
                i1n=nvirt(isbv)                                           12d15s20
                if(ibc(ikeep+icol0).eq.1.or.isbv.ne.jsbv)then           11d9s21
                 kint=kmatsr(isbv,jsbv,isd)+nhere*icol0                 11d5s21
                 do i2=i2s,i2e                                             12d15s20
                  jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)               11d7s21
                  if(i2.eq.i2e)i1n=i1e                                     12d15s20
                  do i1=i10,i1n                                         12d17s20
                   bc(jtmp+i1)=bc(jtmp+i1)+bc(kint)                     11d5s21
                   kint=kint+1                                          11d5s21
                  end do                                                12d17s20
                  i10=1                                                 12d17s20
                 end do                                                 12d17s20
                else                                                    12d17s20
                 kint=kmatsr(isbv,jsbv,isd)+nhere*icol0                 11d5s21
                 do i2=i2s,i2e                                             12d15s20
                  i2m=i2-1                                              12d17s20
                  i2p=i2+1                                              12d17s20
                  jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)               11d7s21
                  if(i2.eq.i2e)i1n=i1e                                     12d15s20
                  kkint=kint-i10                                        11d5s21
                  do i1=i10,min(i1n,i2m)                                12d18s20
                   bc(jtmp+i1)=bc(jtmp+i1)+bc(kkint+i1)                 11d5s21
                  end do                                                12d18s20
                  do i1=max(i10,i2p),i1n                                12d18s20
                   bc(jtmp+i1)=bc(jtmp+i1)+bc(kkint+i1)                 11d5s21
                  end do                                                12d18s20
                  kint=kint+i1n+1-i10                                   11d5s21
                  i10=1                                                 12d17s20
                 end do                                                 12d17s20
                end if                                                  12d17s20
               end do                                                   12d17s20
               itmpd=ibcoff                                             12d17s20
               ibcoff=itmpd+nherev*nokk*nrootu*ncsfb(iarg)              11d5s21
               call enough('hcssbk4. 17',bc,ibc)
               jden=idenk(isc)                                          11d5s21
               do i=0,nokk-1                                            11d5s21
                do m=0,ncsfb(iarg)-1                                    11d5s21
                 do ir=0,nrootm                                         1d28s21
                  iad=itmpd+nherev*(i+nokk*(ir+nrootu*m))               11d5s21
                  do iv=0,nherev-1                                      1d28s21
                   bc(iad+iv)=bc(jden+iv)                               1d28s21
                  end do                                                12d17s20
                  jden=jden+nherev                                      1d28s21
                 end do                                                 12d17s20
                end do                                                  12d17s20
               end do                                                   12d17s20
               mm=nherev*nokk                                           11d5s21
               call dgemm('n','n',n1here,mcol,mm,1d0,                   12d18s20
     $              bc(itmp),n1here,bc(itmpd),mm,1d0,bc(itmpp),         12d18s20
     $              n1here,                                             12d18s20
     d' hcssbk4.  3')
               ibcoff=itmp                                              12d17s20
              end if                                                    12d17s20
              if(nokkv.gt.0)then                                        11d12s21
               itmp=ibcoff                                              12d17s20
               ibcoff=itmp+nherev*nokkv                                 11d5s21
               call enough('hcssbk4. 18',bc,ibc)
               do i=itmp,ibcoff-1                                       12d17s20
                bc(i)=0d0                                               12d17s20
               end do                                                   12d17s20
               do i=0,nokkv-1                                           11d5s21
                icol=ibc(ndenkv(isc)+i)                                 11d5s21
                icol0=icol                                              11d5s21
                iq=icol/nnp                                             11d5s21
                icol=icol-nnp*iq                                        11d5s21
                id=icol/irefo(isc)                                      7d12s21
                ic=icol-irefo(isc)*id                                   7d12s21
                i10=i1s                                                   12d15s20
                i1n=nvirt(isbv)                                           12d15s20
                kint=kmatsr(isbv,jsbv,isd)+nhere*icol0                  11d12s21
                do i2=i2s,i2e                                           12d17s20
                 if(i2.eq.i2e)i1n=i1e                                   12d17s20
                 jtmp=itmp+i2-i2s+nherev*i                              12d17s20
                 if(i2.ge.i10.and.i2.le.i1n)then                        12d18s20
                  kkint=kint-i10
                  bc(jtmp)=bc(jtmp)+bc(kkint+i2)                        11d5s21
                 end if                                                 12d18s20
                 kint=kint+i1n+1-i10                                    11d5s21
                 i10=1                                                  12d17s20
                end do                                                  12d17s20
               end do                                                   12d17s20
               jtmp=itmp-i2s                                            12d17s20
               jtmpd=idenkv(isc)                                        11d5s21
               do i=0,nokkv-1                                           11d5s21
                do m=0,ncsfb(iarg)-1                                     12d17s20
                 do ir=0,nrootm                                         12d17s20
                  jjgg=jgg-1+nvirt(isbv)*ir+ngg*m                       1d27s21
                  ktmpd=jtmpd-i2s+nherev*ir                             1d27s21
                  do iv=i2s,i2e                                          12d17s20
                   bc(jjgg+iv)=bc(jjgg+iv)+bc(jtmp+iv)*bc(ktmpd+iv)     1d27s21
                  end do                                                12d17s20
                 end do                                                 12d17s20
                 jtmpd=jtmpd+nrootu*nherev                              1d27s21
                end do                                                  12d17s20
                jtmp=jtmp+nherev                                        12d17s20
               end do                                                   12d17s20
               ibcoff=itmp                                              12d17s20
              end if                                                    12d17s20
              ikeep=ibcoff                                              11d5s21
              ibcoff=ikeep+nn                                           11d5s21
              call enough('hcssbk4. 19',bc,ibc)
              nok=0                                                     11d5s21
              nokv=0                                                    11d5s21
              do i=0,nn-1                                               12d17s20
               ibc(ikeep+i)=0                                           11d5s21
               if(ibc(ndenjv(isc)+i).ne.0.and.l2e.eq.0)then             2d18s22
                iden=idenjv(isc)+nrowi*i                                11d5s21
                rms=1d0                                                 12d17s20
                if(ibc(ndenj(isc)+i).ne.0)then                          11d5s21
                 jden=idenj(isc)+nrowi*i                                11d5s21
                 rms=0d0                                                12d17s20
                 do k=0,nrowim                                          12d17s20
                  rms=rms+(bc(jden+k)-bc(iden+k))**2                    12d17s20
                 end do                                                 12d17s20
                 rms=sqrt(rms/dfloat(nrowi))                            12d17s20
                end if                                                  12d17s20
                if(rms.gt.1d-10)then
                 jden=idenjv(isc)+nrowi*nokv                            11d5s21
                 do k=0,nrowim                                          12d17s20
                  bc(jden+k)=bc(iden+k)                                 12d17s20
                 end do                                                 12d17s20
                 ibc(ndenjv(isc)+nokv)=i                                11d5s21
                 nokv=nokv+1                                            11d5s21
                else                                                    12d17s20
                 ibc(ikeep+i)=1                                         11d5s21
                end if                                                  12d17s20
               end if                                                   12d17s20
               if(ibc(ndenj(isc)+i).ne.0.and.l2e.eq.0)then              2d18s22
                iden=idenj(isc)+nrowi*i                                 11d5s21
                jden=idenj(isc)+nrowi*nok                               11d5s21
                do k=0,nrowim                                           12d17s20
                 bc(jden+k)=bc(iden+k)                                  12d17s20
                end do                                                  12d17s20
                ibc(ndenj(isc)+nok)=i                                   11d5s21
                nok=nok+1                                               11d5s21
               end if                                                   12d17s20
              end do                                                    12d17s20
              if(nok.gt.0)then                                          12d17s20
               itmp=ibcoff                                               12d17s20
               ibcoff=itmp+nok*n1here*nherev                            12d18s20
               call enough('hcssbk4. 20',bc,ibc)
               ihit=1                                                   12d17s20
               do i=itmp,ibcoff-1                                       12d17s20
                bc(i)=0d0                                               12d17s20
               end do                                                   12d17s20
               do i=0,nok-1                                             12d17s20
                icol=ibc(ndenj(isc)+i)                                  12d17s20
                icol0=icol                                              11d5s21
                iq=icol/nnp                                             11d5s21
                icol=icol-iq*nnp                                        11d5s21
                id=icol/irefo(isc)                                      7d12s21
                ic=icol-id*irefo(isc)                                   7d12s21
                i10=i1s                                                   12d15s20
                i1n=nvirt(isbv)                                           12d15s20
                jint=jmatsr(jsbv,isc,isd)+nhere*icol0                   11d7s21
                if(ibc(ikeep+icol0).eq.1.or.isbv.ne.jsbv)then           11d9s21
                 do i2=i2s,i2e                                             12d15s20
                  if(i2.eq.i2e)i1n=i1e                                     12d15s20
                  jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)               12d18s20
                  do i1=i10,i1n                                         12d17s20
                   bc(jtmp+i1)=bc(jtmp+i1)+bc(jint)                     11d5s21
                   jint=jint+1                                          11d5s21
                  end do                                                12d17s20
                  i10=1                                                 12d17s20
                 end do                                                 12d17s20
                else                                                    12d17s20
                 do i2=i2s,i2e                                             12d15s20
                  i2m=i2-1                                              12d17s20
                  i2p=i2+1                                              12d17s20
                  if(i2.eq.i2e)i1n=i1e                                     12d15s20
                  jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)               11d5s21
                  jjint=jint-i10                                        11d5s21
                  do i1=i10,min(i1n,i2m)                                12d18s20
                   bc(jtmp+i1)=bc(jtmp+i1)+bc(jjint+i1)                 11d5s21
                  end do                                                12d18s20
                  do i1=max(i10,i2p),i1n                                12d18s20
                   bc(jtmp+i1)=bc(jtmp+i1)+bc(jjint+i1)                 11d5s21
                  end do                                                12d18s20
                  jint=jint+i1n+1-i10                                   11d5s21
                  i10=1                                                 12d17s20
                 end do                                                 12d17s20
                end if                                                  12d17s20
               end do                                                   11d5s21
               itmpd=ibcoff                                             12d17s20
               ibcoff=itmpd+nherev*nok*nrootu*ncsfb(iarg)                12d17s20
               call enough('hcssbk4. 21',bc,ibc)
               jden=idenj(isc)                                          12d17s20
               sz=0d0
               do i=0,nok-1                                             12d17s20
                do m=0,ncsfb(iarg)-1                                     12d17s20
                 do ir=0,nrootm                                         12d17s20
                  iad=itmpd+nherev*(i+nok*(ir+nrootu*m))                1d28s21
                  do iv=0,nherev-1                                       12d17s20
                   bc(iad+iv)=bc(jden+iv)                               1d28s21
                   sz=sz+bc(jden+iv)**2
                  end do                                                12d17s20
                  jden=jden+nherev                                      1d28s21
                 end do                                                 12d17s20
                end do                                                  12d17s20
               end do                                                   12d17s20
               mm=nherev*nok                                            12d17s20
               call dgemm('n','n',n1here,mcol,mm,1d0,                   12d18s20
     $              bc(itmp),n1here,bc(itmpd),mm,1d0,bc(itmpp),n1here,  12d18s20
     d' hcssbk4.  4')
               ibcoff=itmp                                              12d17s20
               if(nokv.gt.0)then                                         12d17s20
                itmp=ibcoff                                              12d17s20
                ibcoff=itmp+nherev*nokv                                  12d17s20
                call enough('hcssbk4. 22',bc,ibc)
                do i=itmp,ibcoff-1                                       12d17s20
                 bc(i)=0d0                                               12d17s20
                end do                                                   12d17s20
                do i=0,nokv-1                                            12d17s20
                 icol=ibc(ndenjv(isc)+i)                                 12d17s20
                 i10=i1s                                                   12d15s20
                 i1n=nvirt(isbv)                                           12d15s20
                 jint=jmatsr(jsbv,isc,isd)+nhere*icol                   11d7s21
                 do i2=i2s,i2e                                           12d17s20
                  if(i2.eq.i2e)i1n=i1e                                   12d17s20
                  if(i2.ge.i10.and.i2.le.i1n)then                        12d18s20
                   jtmp=itmp+i2-i2s+nherev*i                              12d17s20
                   jjint=jint-i10                                       11d5s21
                   bc(jtmp)=bc(jtmp)+bc(jjint+i2)                       11d5s21
                  end if                                                 12d18s20
                  jint=jint+i1n+1-i10                                   11d5s21
                  i10=1                                                  12d17s20
                 end do                                                  12d17s20
                end do                                                   12d17s20
                jtmp=itmp-i2s                                            12d17s20
                jtmpd=idenjv(isc)
                do i=0,nokv-1                                            12d17s20
                 do m=0,ncsfb(iarg)-1                                     12d17s20
                  do ir=0,nrootm                                        12d17s20
                   jjgg=jgg-1+nvirt(isbv)*ir+ngg*m                      1d27s21
                   ktmpd=jtmpd-i2s+nherev*ir                            1d27s21
                   do iv=i2s,i2e                                          12d17s20
                    bc(jjgg+iv)=bc(jjgg+iv)+bc(jtmp+iv)*bc(ktmpd+iv)    1d27s21
                   end do                                                12d17s20
                  end do                                                 12d17s20
                  jtmpd=jtmpd+nrootu*nherev                             1d27s21
                 end do                                                  12d17s20
                 jtmp=jtmp+nherev                                        12d17s20
                end do                                                   12d17s20
                ibcoff=itmp                                              12d17s20
               end if                                                    12d17s20
              end if                                                    12d15s20
             end do                                                     12d15s20
             if(ihit.ne.0)then                                          12d17s20
              jjgg=jgg                                                  12d17s20
              do i=0,ncsfb(iarg)-1                                       12d17s20
               do ir=0,nrootm                                           1d28s21
                jjjgg=jjgg+ii10-1+nvirt(isbv)*ir                         1d27s21
                iad=itmpp+n1here*(ir+nrootu*i)                          1d28s21
                do iv=0,n1here-1                                         12d18s20
                 bc(jjjgg+iv)=bc(jjjgg+iv)+bc(iad+iv)                   1d27s21
                end do                                                  12d17s20
               end do                                                   12d17s20
               jjgg=jjgg+nrootu*nvirt(isbv)                             12d18s20
              end do                                                    12d17s20
             end if                                                     12d17s20
             ibcoff=itmpp                                               12d17s20
             jgg=jgg+ngg*ncsfb(iarg)                                     12d16s20
            end do                                                       12d14s20
            ibcoff=ivecbc                                               2d12s21
           end if                                                       12d14s20
          end do                                                          12d12s20
         end do                                                           12d12s20
         i18=1                                                          12d15s20
         i28=ngg
         i38=1                                                          12d15s20
         i48=nff1b(ncloip,isb,1)*ncsfb(iarg)                             7d12s21
         if(itransgg.ne.0)then                                          2d12s21
c
c     igg is virt,root,... so transpose to root,virt,...
          itmp=ibcoff                                                    1d27s21
          ibcoff=itmp+ngg                                                1d27s21
          call enough('hcssbk4. 23',bc,ibc)
          jgg=igg                                                        1d27s21
          do icol=1,nnc                                                  1d28s21
           do i=0,nrootm                                                 1d27s21
            do j=0,nvirt(isbv)-1                                         1d27s21
             ji=jgg+j+nvirt(isbv)*i                                      1d27s21
             ij=itmp+i+nrootu*j                                          1d27s21
             bc(ij)=bc(ji)                                               1d27s21
            end do                                                       1d27s21
           end do                                                        1d27s21
           do i=0,ngg-1                                                  1d28s21
            bc(jgg+i)=bc(itmp+i)                                         1d27s21
           end do                                                        1d27s21
           jgg=jgg+ngg                                                   1d28s21
          end do                                                         1d27s21
          ibcoff=itmp                                                    1d27s21
          call ddi_iacc(bc,ibc,ihsdiagb(ncloip,isb),i18,i28,i38,i48,    11d15s22
     $         bc(igg),ibc(iacc),nacc)                                  11d15s22
          nwiacc=nwiacc+i28*(i48+1-i38)                                 1d26s23
          last8(1)=ihsdiagb(ncloip,isb)                                 8d21s21
          last8(2)=i48                                                  2d17s21
          itransgg=0                                                     2d12s21
         else                                                           2d12s21
          nacc=0                                                        2d12s21
         end if                                                         2d12s21
        end if                                                           12d12s20
       end do                                                           12d14s20
      end do                                                            12d12s20
      ibcoff=iacc                                                       2d12s21
       call dws_synca                                                   2d22s21
      return                                                            12d12s20
      end                                                               12d12s20
