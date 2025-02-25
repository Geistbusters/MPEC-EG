c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcssbk(ihsdiagb,nff1b,iff1b,isymbra,ihsdiagk,nff1k,
     $     iff1k,isymket,                                               7d9s21
     $     ncsf,nec,mdon,mdoop,nsymb,multh,                             7d9s21
     $     ixw1,ixw2,lxmt,ixmt,phase,phase2,nrootu,ism,irel,             7d9s21
     $     irefo,idoubo,nvirt,nbasdws,norb,lprt,maxbx,isymop,n2e,i2eop, 7d9s21
     $     ixmtf,bc,ibc,nwiacc)                                         11d21s22
      implicit real*8 (a-h,o-z)
      external second                                                   2d18s21
      integer*8 ihsdiagb(mdoop,nsymb),i18,i28,i38,i48,i1c,i1o,j1c,j1o,  8d21s21
     $     itestc,itesto,last8(2),ihsdiagk(mdoop,nsymb,2),gandcc,gandco,2d6s23
     $     gandcb                                                       2d6s23
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d19s21
      logical lkeep,lprt,lchoice,lkey                                        3d17s21
      dimension nff1b(mdoop,nsymb,2),ncsf(*),multh(8,8),iff1k(*),       7d9s21
     $     nff1k(mdoop,nsymb,2),nvirt(*),ism(*),irel(*),idoubo(*),      7d9s21
     $     itest(32,3),isorb1(32),idorb1(32),jsorb1(32),jdorb1(32),
     $     nab4(2,3),idenj(8,3),idenk(8,3),irefo(*),iff1b(*),ndenj(8,3),7d15s21
     $     ndenk(8,3),idenjv(8,3),idenkv(8,3),ndenjv(8,3),ndenkv(8,3),  7d15s21
     $     ivec(5,8),ivecp(5,8),nvecp(8),ircv(5,8),nrcv(5,8),           2d12s21
     $     itransv(5,8),nbasdws(*),ixmt(8,*),isymop(*),i2eop(2,3),      7d12s21
     $     ixmtf(*),ikeep(3),nokk(3),nokkv(3),ioxx(2)                   2d6s23
      data loopx/2434540/
      include "common.store"                                            12d12s20
      data icall/0/
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/paddcm/npadddi                                             6d7s22
      save icall
      icall=icall+1
      if(lprt)                                                          1d25s21
     $    write(6,*)('hello! my name is hcssbk'),isymbra,isymket,
     $ nec,mdon,mdoop,nsymb,ixw1,ixw2,norb,maxbx,ibcoff                                     7d9s21
          lkeep=icall.eq.-7
      loop=0
      mdoo=mdoop-1                                                      7d12s21
c
c     maxbx is for input wavefcns ...
c     but here bra has ket roots, so maxbx may not be large enough.     6d7s22
c
      maxbxb=0                                                          6d7s22
      do isb=1,nsymb                                                    6d7s22
       isbv=multh(isb,isymbra)                                          6d7s22
       do ii=mdon+1,mdoop                                               6d30s23
        iarg=ii-mdon                                                    6d7s22
        maxbxb=max(maxbxb,ncsf(iarg)*nff1b(ii,isb,1)*nvirt(isbv)*nrootu)6d8s22
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
      itransgg=0                                                        2d12s21
      call enough('hcssbk.  1',bc,ibc)
      nrootm=nrootu-1                                                   12d15s20
      nsing=0                                                           12d15s20
      loop=0
      norbx=norb+1                                                      12d14s20
      norbxx=norbx+1                                                    12d14s20
      do isb=1,nsymb                                                    12d12s20
       isbv=multh(isb,isymbra)                                          7d9s21
       do jsb=1,nsymb                                                   2d12s21
        nvecp(jsb)=0                                                    2d12s21
       end do                                                           2d12s21
       do ncloip=mdon+1,mdoop                                           7d12s21
        if(min(nff1b(ncloip,isb,1),nvirt(isbv)).gt.0)then               7d9s21
         ncloi=ncloip-1                                                  12d12s20
         iarg=ncloip-mdon                                                12d12s20
         nsing=nsing+nff1b(ncloip,isb,1)*ncsf(iarg)*nvirt(isbv)         7d9s21
         ncolti=nff1b(ncloip,isb,1)*ncsf(iarg)*nvirt(isbv)*nrootu       7d9s21
         ngg=nvirt(isbv)*nrootu                                         12d15s20
         nnc=nff1b(ncloip,isb,1)*ncsf(iarg)                             7d9s21
         nopeni=nec-2*ncloi                                              12d12s20
         do jsb=1,nsymb                                                   12d12s20
          jsbv=multh(jsb,isymket)                                       7d9s21
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
          do nclojp=max(mdon+1,ncloip-2),min(mdoop,ncloip+2)            7d12s21
           if(min(nff1k(nclojp,jsb,1),nherev).gt.0)then                 7d9s21
            ncloj=nclojp-1                                               12d12s20
            jarg=nclojp-mdon                                             12d12s20
            nopenj=nec-2*ncloj                                           12d14s20
            nrow=nherev*nrootu                                          12d14s20
            ncoltj=nff1k(nclojp,jsb,1)*ncsf(jarg)*nrow                  7d9s21
            ivecbc=ibcoff                                               2d12s21
            i48=nff1k(nclojp,jsb,1)*ncsf(jarg)                          7d9s21
            ncol=i48                                                    2d12s21
c
c     has vector alread been loaded?
c
            nminc=mdoop                                                 7d12s21
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
            nrowi=nrow*ncsf(iarg)                                        12d14s20
            nrowim=nrowi-1                                               12d14s20
            idenh=ibcoff                                                 12d14s20
            idend=idenh+nrowi                                           12d14s20
            ibcoff=idend+nrowi                                          12d14s20
            ijsbv=multh(isbv,jsbv)                                       12d14s20
            do isc=1,nsymb                                               12d14s20
             isd=multh(isymop(lxmt),multh(isc,ijsbv))                   7d14s21
             nn=irefo(isc)*irefo(isd)                                   12d14s20
             do ii2e=1,n2e                                              7d15s21
              idenk(isc,ii2e)=ibcoff                                    7d15s21
              ndenk(isc,ii2e)=idenk(isc,ii2e)+nn*nrowi                  7d15s21
              idenkv(isc,ii2e)=ndenk(isc,ii2e)+nn                       7d15s21
              ndenkv(isc,ii2e)=idenkv(isc,ii2e)+nn*nrowi                7d15s21
              ibcoff=ndenkv(isc,ii2e)+nn                                7d15s21
              idenj(isc,ii2e)=ibcoff                                    7d15s21
              ndenj(isc,ii2e)=idenj(isc,ii2e)+nn*nrowi                  7d15s21
              idenjv(isc,ii2e)=ndenj(isc,ii2e)+nn                       7d15s21
              ndenjv(isc,ii2e)=idenjv(isc,ii2e)+nn*nrowi                7d15s21
              ibcoff=ndenjv(isc,ii2e)+nn                                7d15s21
             end do                                                     7d15s21
            end do                                                       12d14s20
            call enough('hcssbk.  2',bc,ibc)
            ibctop=ibcoff-1                                              12d14s20
            if1o=nff1b(ncloip,isb,2)                                      12d14s20
            jgg=igg                                                     12d14s20
            jvecer=0                                                    2d12s21
            do if1=1,nff1b(ncloip,isb,1)                                7d9s21
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
              if(itransv(nuse,jsb).eq.0.and.ndifb.le.4)then             2d6s23
               call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))        2d12s21
               itmp=ibcoff                                                    1d27s21
               ibcoff=itmp+nrow                                            1d27s21
               call enough('hcssbk.  3',bc,ibc)
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
              if(ndifb.le.4)then                                        2d6s23
               ivtrans=ibcoff                                           2d19s21
               ibcoff=ivtrans+nrow*ncsf(jarg)                           2d19s21
               call enough('hcssbk.  4',bc,ibc)
               do j=0,ncsf(jarg)-1                                      2d19s21
                do i=0,nrow-1                                           2d19s21
                 ij=jvec+i+nrow*j                                       2d19s21
                 ji=ivtrans+j+ncsf(jarg)*i                              2d19s21
                 bc(ji)=bc(ij)                                          2d19s21
                end do                                                  2d19s21
               end do                                                   2d19s21
               ndifs=popcnt(gandco)                                     2d6s23
               ndifd=popcnt(gandcc)                                     2d6s23
               if(ndifs.eq.2.and.ndifb.eq.2)then                        2d6s23
                do i=1,norbxx
                 if(btest(gandco,i))then
                  if((btest(j1c,i).and.btest(i1o,i)).or.                 10d21s22
     $             (btest(j1o,i).and..not.btest(i1c,i)))then            10d21s22
                   nab4(1,1)=i
                  else
                   nab4(2,1)=i
                  end if
                 end if
                end do
                call gandc(i1c,i1o,j1c,j1o,nopeni,nopenj,iarg,jarg,ncsf, 2d19s21
     $            norbxx,ixw1,ixw2,nnot,nab1,iwpb1,iwpk1,ncsfmid1,bc,   11d14s22
     $             ibc)                                                 11d14s22
                itmp=ibcoff                                               12d14s20
                ibcoff=itmp+nrowi                                         12d14s20
                call enough('hcssbk.  5',bc,ibc)
                call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid1,nrow,        2d19s21
     $              iwpb1,iwpk1,bc(ivtrans),ncsf(jarg),bc(itmp),        2d19s21
     $              ncsf(iarg),1d0,0d0,bc,ibc)                          11d10s22
                if(isymop(lxmt).eq.multh(jsbv,isbv))then                 7d19s21
                 ndenh=1                                                  12d14s20
c
c     g r,v,i = sum v' denh r,v',i h0v v'
                 do i=0,nrowim                                             12d14s20
                  bc(idenh+i)=bc(idenh+i)+bc(itmp+i)                      2d22s21
                 end do                                                    12d14s20
                end if                                                   7d15s21
                nok=0                                                         11d13s20
                if(n2e.gt.0)then                                         8d23s21
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
                   itest(nok,2)=i                                              11d13s20
                  end if                                                       11d13s20
                 end do                                                        11d13s20
                end if                                                   8d23s21
                do i=1,nok                                                    11d13s20
                 js=ism(itest(i,2))                                           11d13s20
                 jg=irel(itest(i,2))-1                                        11d13s20
                 do ii2e=1,n2e                                           7d15s21
                  if(isymop(i2eop(1,ii2e)).eq.1.and.isbv.eq.jsbv)then    7d15s21
                   icolj=jg+irefo(js)*jg                                   7d14s21
c     (vu|gg)
                   jdenj=idenj(js,ii2e)+nrowi*icolj                      7d15s21
                   ibc(ndenj(js,ii2e)+icolj)=1                           7d15s21
                   fct=1d0                                                7d14s21
                   if(itest(i,3).eq.2)fct=2d0                             7d14s21
                   do j=0,nrowim                                           12d14s20
                    bc(jdenj+j)=bc(jdenj+j)+bc(itmp+j)                    2d22s21
                   end do                                                  12d14s20
                  end if                                                   12d14s20
                 end do                                                  7d15s21
                 if(itest(i,3).eq.2)then                                  12d14s20
                  do ii2e=1,n2e                                          7d15s21
                   if(multh(js,jsbv).eq.isymop(i2eop(1,ii2e)))then       7d15s21
                    icolk=jg+irefo(js)*jg                                  12d14s20
                    jdenk=idenk(js,ii2e)+nrowi*icolk                     7d15s21
                    ibc(ndenk(js,ii2e)+icolk)=1                          7d15s21
                    do j=0,nrowim                                           12d14s20
                     bc(jdenk+j)=bc(jdenk+j)-bc(itmp+j)                    2d22s21
                    end do                                                  12d14s20
                   end if                                                 7d14s21
                  end do                                                 7d15s21
                 end if                                                  7d14s21
                end do                                                    12d14s20
                ibcoff=itmp                                              2d19s21
                do i=1,nok                                                    11d13s20
                 if(itest(i,3).eq.1)then                                  12d14s20
                  itestc=j1c                                              12d14s20
                  itesto=j1o                                              12d14s20
                  nopenk=nopenj                                                11d13s20
c
c     anihilate common
c
                  if(btest(itestc,itest(i,2)))then                             11d13s20
                   itestc=ibclr(itestc,itest(i,2))                             11d13s20
                   itesto=ibset(itesto,itest(i,2))                             11d13s20
                   karg=jarg-1                                                11d13s20
                   nopenk=nopenk+1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibclr(itesto,itest(i,2))                             11d13s20
                   karg=jarg                                                  11d13s20
                   nopenk=nopenk-1                                             11d13s20
                  end if                                                       11d13s20
c
c     create ket
c
                  if(btest(itesto,nab4(2,1)))then                        12d14s20
                   itestc=ibset(itestc,nab4(2,1))                        12d14s20
                   itesto=ibclr(itesto,nab4(2,1))                        12d14s20
                   karg=karg+1                                                11d13s20
                   nopenk=nopenk-1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibset(itesto,nab4(2,1))                        12d14s20
                   nopenk=nopenk+1                                             11d13s20
                  end if                                                       11d13s20
                  nqq=karg+mdon-1
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                   call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,       2d19s21
     $            iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,  2d19s21
     $                iwpk1b,ncsfmid1b,bc,ibc)                          11d14s22
                   call gandc(itestc,itesto,j1c,j1o,nopenk,nopenj,       2d19s21
     $            karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2b,nab2b,iwpb2b,     12d14s20
     $               iwpk2b,ncsfmid2b,bc,ibc)                           11d14s22
                   if(nnot1b.eq.2.and.nnot2b.eq.2)then                         11d13s20
                    itmp1=ibcoff                                          2d19s21
                    itmp2=itmp1+ncsf(iarg)*nrow                           2d19s21
                    ibcoff=itmp2+ncsf(iarg)*nrow                          2d19s21
                    call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),
     $                ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b,  2d19s21
     $                bc(ivtrans),ncsf(jarg),nrow,bc(itmp2),0d0,bc,ibc) 11d10s22
c     (vb|au)
                    lsa=ism(nab2b(2))                                      12d14s20
                    lga=irel(nab2b(2))-1                                   12d14s20
                    lsb=ism(nab1b(1))                                      12d14s20
                    lgb=irel(nab1b(1))-1                                   12d14s20
                    icol=lga+irefo(lsa)*lgb                               12d14s20
                    do ii2e=1,n2e                                        7d15s21
                     if(multh(lsa,jsbv).eq.isymop(i2eop(1,ii2e)))then    7d15s21
                      kden=idenk(lsa,ii2e)+nrowi*icol                    7d15s21
                      ibc(ndenk(lsa,ii2e)+icol)=1                        7d15s21
                      do j=0,nrowim                                         2d19s21
                       bc(kden+j)=bc(kden+j)+bc(itmp2+j)                    2d22s21
                      end do                                                2d19s21
                     end if                                               7d14s21
                    end do                                               7d15s21
                    ibcoff=itmp1                                          2d19s21
                   end if                                                  12d14s20
                  end if                                                 4d14s21
                 end if                                                   12d14s20
                end do                                                    12d14s20
               else if(n2e.gt.0)then                                    2d6s23
                nnot=0                                                  2d6s23
                ipssx=0                                                 2d6s23
                if(ndifs.eq.4.and.ndifb.eq.4)then                       2d6s23
                 nnot=4                                                 2d6s23
                 ioxx(1)=1                                              2d6s23
                 ioxx(2)=1                                              2d6s23
                 do i=1,norbxx                                          2d6s23
                  if(btest(gandcb,i))then                               2d6s23
                   if((btest(i1c,i).and.btest(j1o,i)).or.               2d6s23
     $                   (btest(i1o,i).and..not.btest(j1c,i)))then      2d6s23
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
     $                   ((btest(j1c,i).and..not.btest(i1o,i)).or.      2d6s23
     $                   (btest(i1c,i).and..not.btest(j1o,i))))then     2d6s23
                    if(btest(i1c,i))iswap=1                             2d6s23
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
                  if(btest(j1c,nab4(1,2)).and.                          2d6s23
     $                 .not.btest(j1c,nab4(1,1)))nbt=1                  2d6s23
                 else                                                   2d6s23
                  nbt=0                                                 2d6s23
                  if(btest(i1c,nab4(2,2)).and.                          2d6s23
     $                   .not.btest(i1c,nab4(2,1)))nbt=1                2d6s23
                 end if                                                 2d6s23
                 if(nbt.ne.0)then                                       2d6s23
                  nab4(1,1)=nab4(1,2)                                   2d6s23
                  nab4(2,1)=nab4(2,2)                                   2d6s23
                 end if                                                 2d6s23
                else if(ndifs.eq.0.and.ndifd.eq.2)then                  2d6s23
                 nnot=3                                                 2d6s23
                 do i=1,norbxx                                          2d6s23
                  if(btest(gandcb,i))then                               2d6s23
                   if(btest(j1c,i))then                                 2d6s23
                    nab4(1,1)=i                                         2d6s23
                    nab4(1,2)=i                                         2d6s23
                   else                                                 2d6s23
                    nab4(2,1)=i                                         2d6s23
                    nab4(2,2)=i                                         2d6s23
                   end if                                               2d6s23
                  end if                                                2d6s23
                 end do                                                 2d6s23
                end if                                                  2d6s23
                if(nnot.eq.3)then                                          12d8s20
                 ipssx=1                                                   12d8s20
                else if(nnot.eq.4)then                                  2d6s23
                 ipssx=2                                                   12d8s20
                end if                                                     12d8s20
                do ipss=1,ipssx                                            12d8s20
                 if(ipss.eq.1)then
                  iu1=1
                  iu2=1
                 else if(ipss.eq.2)then
                  iu1=1
                  iu2=2
                 else
                  iu1=2
                  iu2=1
                 end if                                                    12d8s20
                 itestc=j1c                                              12d14s20
                 itesto=j1o                                              12d14s20
                 if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                  itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                  itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopenj+1                                              11d13s20
                  karg=jarg-1                                             12d8s20
                 else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                  itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopenj-1                                              11d13s20
                  karg=jarg                                              12d14s20
                 end if                                                        11d13s20
                 if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                  itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                  itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                  nopenk=nopenk-1                                              11d13s20
                  karg=karg+1                                                  11d13s20
                 else                                                          11d13s20
                  itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                  nopenk=nopenk+1                                              11d13s20
                 end if                                                        11d13s20
                 nqq=karg+mdon-1                                         4d14s21
                 if(nqq.ge.mdon.and.nqq.le.mdoo)then                     4d14s21
                  call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,         2d19s21
     $         iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,
     $               iwpk1b,ncsfmid1b,bc,ibc)                           11d14s22
                  call gandc(itestc,itesto,j1c,j1o,nopenk,nopenj,         2d19s21
     $           karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2b,nab2b,iwpb2b,
     $               iwpk2b,ncsfmid2b,bc,ibc)                           11d14s22
                  if(nnot1b.eq.2.and.nnot2b.eq.2)then                              11d13s20
                   nab1(1)=nab2b(2)                                       2d19s21
                   nab1(2)=nab2b(1)                                       2d19s21
                   nab2(1)=nab1b(2)                                       2d19s21
                   nab2(2)=nab1b(1)                                       2d19s21
                   igotit=0                                              7d14s21
                   if(max(nab1(1),nab1(2)).le.norb)then                    12d14s20
                    lsa=ism(nab1(1))                                       12d14s20
                    lga=irel(nab1(1))-1                                    12d14s20
                    lsb=ism(nab1(2))                                       12d14s20
                    lgb=irel(nab1(2))-1                                    12d14s20
                    icol=lgb+irefo(lsb)*lga                              7d16s21
                    do ii2e=1,n2e                                        7d15s21
                     if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)).and.     7d15s21
     $                    multh(isbv,jsbv).eq.isymop(i2eop(2,ii2e)))then 7d15s21
                      iden=idenj(lsb,ii2e)+nrowi*icol                    7d15s21
                      ibc(ndenj(lsb,ii2e)+icol)=1                        7d15s21
                      itmp1=ibcoff                                           2d19s21
                      itmp2=itmp1+ncsf(iarg)*nrow                            2d19s21
                      ibcoff=itmp2+ncsf(iarg)*nrow                           2d19s21
                      call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),
     $                ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b,  2d19s21
     $                bc(ivtrans),ncsf(jarg),nrow,bc(itmp2),0d0,bc,ibc) 11d10s22
                      do j=0,nrowim                                          2d19s21
                       bc(iden+j)=bc(iden+j)+bc(itmp2+j)                     2d22s21
                      end do                                                 2d19s21
                      ibcoff=itmp1                                           2d19s21
                     end if                                               7d14s21
                    end do                                               7d15s21
                   else                                                    12d14s20
                    lsb=ism(nab2(2))                                       12d14s20
                    lsa=ism(nab1(1))                                       12d14s20
                    lga=irel(nab1(1))-1                                    12d14s20
                    lgb=irel(nab2(2))-1                                    12d14s20
                    icol=lgb+irefo(lsb)*lga                                12d14s20
                    do ii2e=1,n2e                                        7d15s21
                     if(multh(lsb,jsbv).eq.isymop(i2eop(1,ii2e)))then    7d15s21
                      iden=idenk(lsb,ii2e)+nrowi*icol                             12d14s20
                      ibc(ndenk(lsb,ii2e)+icol)=1                        7d15s21
                      itmp1=ibcoff                                           2d19s21
                      itmp2=itmp1+ncsf(iarg)*nrow                            2d19s21
                      ibcoff=itmp2+ncsf(iarg)*nrow                           2d19s21
                      call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),
     $                ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b,  2d19s21
     $                bc(ivtrans),ncsf(jarg),nrow,bc(itmp2),0d0,bc,ibc) 11d10s22
                      do j=0,nrowim                                          2d19s21
                       bc(iden+j)=bc(iden+j)+bc(itmp2+j)                     2d22s21
                      end do                                                 2d19s21
                      ibcoff=itmp1                                           2d19s21
                     end if                                               7d14s21
                    end do                                               7d15s21
                   end if                                                  12d14s20
                  end if                                                  12d14s20
                 end if                                                  4d14s21
                end do                                                   12d14s20
               end if                                                     12d14s20
               ibcoff=ivtrans                                           2d6s23
              end if                                                    2d6s23
              if(jsbv.eq.isbv)then                                      7d19s21
               j1o=ibclr(j1o,norbxx)                                     12d14s20
               j1o=ibset(j1o,norbx)                                      12d14s20
               gandco=ieor(j1o,i1o)                                        2d6s23
               gandcb=ior(gandcc,gandco)                                     2d6s23
               ndifb=popcnt(gandcb)                                          2d6s23
c
c     vec is ordered nrootu,nherev. re-order to be nherev,nrootu.
               if(itransv(nuse,jsb).eq.0.and.ndifb.le.4)then            2d6s23
                call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))       2d12s21
                itmp=ibcoff                                                    1d27s21
                ibcoff=itmp+nrow                                            1d27s21
                call enough('hcssbk.  6',bc,ibc)
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
               if(ndifb.le.4)then                                       2d6s23
                ivtrans=ibcoff                                           2d19s21
                ibcoff=ivtrans+nrow*ncsf(jarg)                           2d19s21
                call enough('hcssbk.  7',bc,ibc)
                do j=0,ncsf(jarg)-1                                      2d19s21
                 do i=0,nrow-1                                           2d19s21
                  ij=jvec+i+nrow*j                                       2d19s21
                  ji=ivtrans+j+ncsf(jarg)*i                              2d19s21
                  bc(ji)=bc(ij)                                          2d19s21
                 end do                                                  2d19s21
                end do                                                   2d19s21
                ndifs=popcnt(gandco)                                    2d6s23
                ndifd=popcnt(gandcc)                                    2d6s23
                if(ndifs.eq.0.and.ndifd.eq.0)then                       2d6s23
                 sum=0d0                                                 7d14s21
                 if(isymop(lxmt).eq.1)then                               7d14s21
                  do i=1,norb                                            7d14s21
                   js=ism(i)                                             7d14s21
                   jg=irel(i)-1+idoubo(js)                               7d14s21
                   if(btest(i1c,i))then                                    7d14s21
                    iad=ixmtf(js)+jg*(nbasdws(js)+1)                     7d14s21
                    sum=sum+2d0*bc(iad)*phase                            7d15s21
                   else if(btest(i1o,i))then                               7d14s21
                    iad=ixmtf(js)+jg*(nbasdws(js)+1)                     7d14s21
                    sum=sum+bc(iad)*phase                                7d15s21
                   end if                                                7d14s21
                  end do                                                 7d14s21
                 end if                                                  7d14s21
                 if(n2e.gt.0)then                                        8d23s21
                  do i=1,norb                                             7d14s21
                   is=ism(i)                                              7d14s21
                   ig=irel(i)-1+idoubo(is)                                7d14s21
                   if(btest(i1c,i))then                                   7d14s21
                    do j=1,norb                                           7d14s21
                     js=ism(j)                                            7d14s21
                     jg=irel(j)-1+idoubo(js)                              7d14s21
                     if(btest(i1c,j))then                                 7d14s21
                      do ii2e=1,n2e                                       7d14s21
                       if(isymop(i2eop(1,ii2e)).eq.1)then                 7d14s21
                        iii=ixmt(is,i2eop(1,ii2e))+ig*(nbasdws(is)+1)     7d14s21
                        jjj=ixmt(js,i2eop(2,ii2e))+jg*(nbasdws(js)+1)     7d14s21
                        sum=sum+2d0*phase2*bc(iii)*bc(jjj)                7d14s21
                       end if                                             7d14s21
                       if(isymop(i2eop(1,ii2e)).eq.multh(is,js))then      7d14s21
                        ij=ixmt(is,i2eop(1,ii2e))+ig+nbasdws(is)*jg       7d14s21
                        ji=ixmt(js,i2eop(2,ii2e))+jg+nbasdws(js)*ig       7d14s21
                        sum=sum-phase2*bc(ij)*bc(ji)                      7d14s21
                       end if                                             7d14s21
                      end do                                              7d14s21
                     else if(btest(i1o,j))then                            7d14s21
                      do ii2e=1,n2e                                       7d14s21
                       if(isymop(i2eop(1,ii2e)).eq.1)then                 7d14s21
                        iii=ixmt(is,i2eop(1,ii2e))+ig*(nbasdws(is)+1)     7d14s21
                        jjj=ixmt(js,i2eop(2,ii2e))+jg*(nbasdws(js)+1)     7d14s21
                        sum=sum+2d0*phase2*bc(iii)*bc(jjj)                7d14s21
                       end if                                             7d14s21
                       if(isymop(i2eop(1,ii2e)).eq.multh(is,js))then      7d14s21
                        ij=ixmt(is,i2eop(1,ii2e))+ig+nbasdws(is)*jg       7d14s21
                        ji=ixmt(js,i2eop(2,ii2e))+jg+nbasdws(js)*ig       7d14s21
                        term=bc(ij)*bc(ji)
                        sum=sum-phase2*bc(ij)*bc(ji)                      7d14s21
                       end if                                             7d14s21
                      end do                                              7d14s21
                     end if                                               7d14s21
                    end do                                                7d14s21
                   end if                                                 7d14s21
                   if(btest(i1o,i))then                                   7d14s21
                    do j=1,norb                                           7d14s21
                     if(j.ne.i.and.btest(i1o,j))then                      7d14s21
                      js=ism(j)                                           7d14s21
                      jg=irel(j)-1+idoubo(js)                             7d14s21
                      if(isymop(i2eop(1,1)).eq.1)then                     7d14s21
                       do ii2e=1,n2e                                      7d14s21
                        iii=ixmt(is,i2eop(1,ii2e))+ig*(nbasdws(is)+1)     7d14s21
                        jjj=ixmt(js,i2eop(2,ii2e))+jg*(nbasdws(js)+1)     7d14s21
                        sum=sum+0.5d0*phase2*bc(iii)*bc(jjj)              7d14s21
                       end do                                             7d14s21
                      end if                                              7d14s21
                     end if                                               7d14s21
                    end do                                                7d14s21
                   end if                                                 7d14s21
                  end do                                                  7d14s21
                 end if                                                  8d23s21
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
                 if(nnhere.gt.0)then                                     7d15s21
                  if(itransgg.eq.0)then                                      2d12s21
                   call ddi_done(ibc(iacc),nacc)                             2d12s21
                   itransgg=1                                                2d12s21
                   do i=igg,igg+ncolti-1                                     2d12s21
                    bc(i)=0d0                                                2d12s21
                   end do                                                    2d12s21
                  end if                                                     2d12s21
                  idtmp=ibcoff                                           7d15s21
                  ibcoff=idtmp+nnhere                                    7d15s21
                  call enough('hcssbk.  8',bc,ibc)
                  do i=idtmp,ibcoff-1                                    7d15s21
                   bc(i)=sum                                             7d15s21
                  end do                                                 7d15s21
                  jdtmp=idtmp-ill                                        7d15s21
                  if(isymop(lxmt).eq.1)then                              7d15s21
                   do ivp=ill,ihh                                         7d15s21
                    iv=ivp-1+idoubo(isbv)+irefo(isbv)                    7d15s21
                    iad=ixmtf(isbv)+iv*(nbasdws(isbv)+1)                 7d15s21
                    bc(jdtmp+ivp)=bc(jdtmp+ivp)+phase*bc(iad)            7d15s21
                   end do                                                 7d15s21
                  end if                                                 7d15s21
                  if(n2e.gt.0)then                                       8d23s21
                   if(isymop(i2eop(1,1)).eq.1)then                        7d15s21
                    do i=1,norb                                           7d15s21
                     if(btest(i1c,i).or.btest(i1o,i))then                 7d15s21
                      is=ism(i)                                           7d15s21
                      ig=irel(i)-1+idoubo(is)                             7d15s21
                      ff=phase2                                           7d15s21
                      if(btest(i1c,i))ff=ff*2d0                           7d15s21
                      do ii2e=1,n2e                                       7d15s21
                       iaa=ixmt(is,i2eop(1,ii2e))+ig*(nbasdws(is)+1)      7d15s21
                       faa=ff*bc(iaa)                                     7d15s21
                       do ivp=ill,ihh
                        iv=ivp-1+idoubo(isbv)+irefo(isbv)                 7d15s21
                        ivv=ixmt(isbv,i2eop(2,ii2e))                    2d6s23
     $                       +iv*(nbasdws(isbv)+1)                      2d6s23
                        bc(jdtmp+ivp)=bc(jdtmp+ivp)+faa*bc(ivv)           7d15s21
                       end do                                             7d15s21
                      end do                                              7d15s21
                     end if                                               7d15s21
                    end do                                                7d15s21
                   end if                                                 7d15s21
                   do i=1,norb                                            7d15s21
                    if(btest(i1c,i))then                                  7d15s21
                     is=ism(i)                                            7d15s21
                     ig=irel(i)-1+idoubo(is)                              7d15s21
                     do ii2e=1,n2e                                        7d15s21
                      if(multh(is,isbv).eq.isymop(i2eop(1,ii2e)))then     7d15s21
                       do ivp=ill,ihh                                     7d15s21
                        iv=ivp-1+idoubo(isbv)+irefo(isbv)                 7d15s21
                        ivg=ixmt(isbv,i2eop(1,ii2e))+iv+nbasdws(isbv)*ig  7d15s21
                        igv=ixmt(is,i2eop(2,ii2e))+ig+nbasdws(is)*iv      7d15s21
                        term=bc(ivg)*bc(igv)
                        bc(jdtmp+ivp)=bc(jdtmp+ivp)                      8d23s21
     $                      -phase2*bc(ivg)*bc(igv)                     8d23s21
                       end do                                             7d15s21
                      end if                                              7d15s21
                     end do                                               7d15s21
                    end if                                                7d15s21
                   end do                                                 7d15s21
                  end if                                                 8d23s21
                  do i=0,ncsf(iarg)-1                                    7d15s21
                   do ir=0,nrootm                                        7d15s21
                    iadv=jvec-i2s+nherev*(ir+nrootu*i)                   7d15s21
                    iadg=jgg-1+nvirt(isbv)*(ir+nrootu*i)                 7d15s21
                    do ivp=ill,ihh                                       7d15s21
                     bc(iadg+ivp)=bc(iadg+ivp)                          2d6s23
     $                    +bc(jdtmp+ivp)*bc(iadv+ivp)                   2d6s23
                    end do                                               7d15s21
                   end do                                                7d15s21
                  end do                                                 7d15s21
                  ibcoff=idtmp                                           7d15s21
                 end if                                                  7d15s21
                 if(n2e.gt.0)then                                        8d23s21
                  do i1=1,nopeni-1                                        12d17s20
                   do i2=i1+1,nopeni                                      12d17s20
                    if(i2.eq.nopeni)then                                  12d17s20
                     nab1(1)=norbx                                        12d17s20
                    else                                                  12d17s20
                     nab1(1)=isorb1(i2)                                      12d6s20
                    end if                                                12d17s20
                    nab1(2)=isorb1(i1)                                      12d6s20
                    nopenk=nopeni-2                                              11d13s20
                    karg=iarg+1                                                 11d13s20
                    nab2(1)=nab1(2)                                             11d13s20
                    nab2(2)=nab1(1)                                             11d13s20
                    lsb=ism(nab1(2))                                              11d13s20
                    lgb=irel(nab1(2))-1                                           11d13s20
                    lsc=ism(nab2(1))                                              11d13s20
                    lgc=irel(nab2(1))-1                                           11d13s20
                    nqq=karg+mdon-1
                    igotit=0                                              7d14s21
                    if(i2.ne.nopeni)then                                  12d17s20
                     lsa=ism(nab1(1))                                              11d13s20
                     lga=irel(nab1(1))-1                                           11d13s20
                     lsd=ism(nab2(2))                                              11d13s20
                     lgd=irel(nab2(2))-1                                           11d13s20
                     xint=0d0                                             7d9s21
                     do ii2e=1,n2e                                        7d9s21
                      if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)))then         7d14s21
                       igotit=1                                            7d14s21
                       iab=ixmt(lsa,i2eop(1,ii2e))+lga+idoubo(lsa)         7d9s21
     $                   +nbasdws(lsa)*(lgb+idoubo(lsb))                7d9s21
                       icd=ixmt(lsc,i2eop(2,ii2e))+lgc+idoubo(lsc)         7d9s21
     $                   +nbasdws(lsc)*(lgd+idoubo(lsd))                7d9s21
                       xint=xint+bc(iab)*bc(icd)*phase2                    7d9s21
                      end if                                              7d16s21
                     end do                                               7d9s21
                     if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                      itestc=i1c                                            12d14s20
                      itesto=i1o                                            12d14s20
                      if(i2.ne.nopeni)then                                 12d17s20
                       itesto=ibclr(itesto,isorb1(i2))                                     11d13s20
                      else                                                 12d17s20
                       itesto=ibclr(itesto,norbx)                          12d17s20
                      end if                                               12d17s20
                      itesto=ibclr(itesto,isorb1(i1))                                     11d13s20
                      itestc=ibset(itestc,isorb1(i1))                           11d13s20
                      call gandc(itestc,itesto,i1c,i1o,nopenk,nopeni,      2d19s21
     $         karg,iarg,ncsf,norbx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d14s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                      call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,       12d14s20
     $         iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   12d14s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                      ltmp1=ibcoff                                          2d19s21
                      ltmp2=ltmp1+ncsf(iarg)*nrow                           2d19s21
                      ibcoff=ltmp2+ncsf(iarg)*nrow                          2d19s21
                      call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),        2d19s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        2d19s21
     $                bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc) 11d10s22
                      do j=0,nrowim                                         2d19s21
                       bc(ltmp2+j)=bc(ltmp2+j)-bc(ivtrans+j)               2d22s21
                      end do                                                2d19s21
                      do j=0,nrowim                                       2d19s21
                       bc(idend+j)=bc(idend+j)+xint*bc(ltmp2+j)           2d22s21
                      end do                                              2d19s21
                     else
                      do i=0,nrowim                                         12d14s20
                       bc(idend+i)=bc(idend+i)-xint*bc(ivtrans+i)         2d22s21
                      end do                                                12d14s20
                     end if                                               7d15s21
                    else                                                  12d17s20
                     do ii2e=1,n2e                                        7d15s21
                      if(multh(lsb,jsbv).eq.isymop(i2eop(1,ii2e)))then    7d15s21
                       igotit=1                                            7d14s21
                       icol=lgb+irefo(lsb)*lgc                              12d17s20
                       iden=idenkv(lsb,ii2e)+nrowi*icol                   7d15s21
                       ibc(ndenkv(lsb,ii2e)+icol)=1                       7d15s21
                       if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                        itestc=i1c                                            12d14s20
                        itesto=i1o                                            12d14s20
                        if(i2.ne.nopeni)then                                 12d17s20
                         itesto=ibclr(itesto,isorb1(i2))                                     11d13s2
                        else                                                 12d17s20
                         itesto=ibclr(itesto,norbx)                          12d17s20
                        end if                                               12d17s20
                        itesto=ibclr(itesto,isorb1(i1))                                     11d13s20
                        itestc=ibset(itestc,isorb1(i1))                           11d13s20
                        call gandc(itestc,itesto,i1c,i1o,nopenk,nopeni,      2d19s21
     $         karg,iarg,ncsf,norbx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d14s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                        call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,       12d14s20
     $         iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   12d14s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                        ltmp1=ibcoff                                          2d19s21
                        ltmp2=ltmp1+ncsf(iarg)*nrow                           2d19s21
                        ibcoff=ltmp2+ncsf(iarg)*nrow                          2d19s21
                        call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),        2d19s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        2d19s21
     $                bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc) 11d10s22
                        do j=0,nrowim                                         2d19s21
                         bc(ltmp2+j)=bc(ltmp2+j)-bc(ivtrans+j)               2d22s21
                        end do                                                2d19s21
                        do j=0,nrowim                                       2d19s21
                         bc(iden+j)=bc(iden+j)+bc(ltmp2+j)                  2d22s21
                        end do                                              2d19s21
                        ibcoff=ltmp1                                         2d19s21
                       else
                        do i=0,nrowim                                         12d14s20
                         bc(iden+i)=bc(iden+i)-bc(ivtrans+i)                2d22s21
                        end do                                                12d14s20
                       end if                                             7d15s21
                      end if                                               7d14s21
                     end do                                               7d15s21
                    end if                                                12d17s20
                   end do                                                  12d14s20
                  end do                                                  12d14s20
                 end if                                                  8d23s21
                else if(ndifs.eq.2.and.ndifb.eq.2)then
                 do i=1,norbx                                           2d6s23
                  if(btest(gandco,i))then
                   if((btest(j1c,i).and.btest(i1o,i)).or.                 10d21s22
     $             (btest(j1o,i).and..not.btest(i1c,i)))then            10d21s22
                    nab4(1,1)=i
                   else
                    nab4(2,1)=i
                   end if
                  end if
                 end do
                 call gandc(i1c,i1o,j1c,j1o,nopeni,nopenj,iarg,jarg,    2d6s23
     $         ncsf,norbx,ixw1,ixw2,nnotb,nab1b,iwpb1b,iwpk1b,ncsfmid1b,2d6s23
     $               bc,ibc)                                            11d14s22
c     i is bra, j is ket
                 nab1(1)=nab1b(2)                                        2d22s21
                 nab1(2)=nab1b(1)                                        2d22s21
                 lsa=ism(nab1(1))                                        12d15s20
                 lsb=ism(nab1(2))                                        7d9s21
                 lga=irel(nab1(1))-1                                     12d15s20
                 lgb=irel(nab1(2))-1                                     12d15s20
                 if(multh(lsa,lsb).eq.isymop(lxmt))then                  7d15s21
                  ih0=ixmtf(lsb)+lgb+idoubo(lsb)                         7d19s21
     $                +nbasdws(lsb)*(lga+idoubo(lsa))                    7d19s21
                  sum=bc(ih0)*phase                                       7d9s21
                 else                                                    7d15s21
                  sum=0d0                                                7d15s21
                 end if                                                  7d15s21
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
                 nok=0                                                         11d13s20
                 if(n2e.gt.0)then                                        8d23s21
                  do i=1,norb                                                   11d13s20
                   ixn=min(itest(i,1),itest(i,2))
                   if(ixn.gt.0)then                                             11d13s20
                    nok=nok+1                                                   11d13s20
                    itest(nok,3)=ixn                                            11d13s20
                    itest(nok,2)=i                                              11d13s20
                   end if                                                       11d13s20
                  end do                                                        11d13s20
                 end if                                                  8d23s21
                 do ii2e=1,n2e                                           7d9s21
                  do i=1,nok                                                    11d13s20
                   js=ism(itest(i,2))                                           11d13s20
                   jg=irel(itest(i,2))-1                                        11d13s20
                   if(isymop(i2eop(1,ii2e)).eq.1)then                    7d14s21
                    jab=ixmt(js,i2eop(1,ii2e))+jg+idoubo(js)              7d9s21
     $                 +nbasdws(js)*(jg+idoubo(js))                     7d9s21
                    jcd=ixmt(lsa,i2eop(2,ii2e))+lga+idoubo(lsa)           7d9s21
     $                 +nbasdws(lsa)*(lgb+idoubo(lsa))                  7d9s21
                    if(itest(i,3).eq.2)then                              7d14s21
                     sum=sum+phase2*2d0*bc(jab)*bc(jcd)                  7d14s21
                    else                                                 7d14s21
                     sum=sum+phase2*bc(jab)*bc(jcd)                      7d14s21
                    end if                                               7d14s21
                   end if                                                7d14s21
                   if(itest(i,3).eq.2)then                                  12d14s20
                    if(multh(lsa,js).eq.isymop(i2eop(1,ii2e)))then       7d14s21
                     kba=ixmt(lsa,i2eop(1,ii2e))+lga+idoubo(lsa)         7d16s21
     $                   +nbasdws(lsa)*(jg+idoubo(js))                   7d16s21
                     kcd=ixmt(js,i2eop(2,ii2e))+jg+idoubo(js)             7d9s21
     $                  +nbasdws(js)*(lgb+idoubo(lsb))                  7d9s21
                     sum=sum-phase2*bc(kba)*bc(kcd)                      7d16s21
                    end if                                               7d14s21
                   end if                                                   12d14s20
                  end do                                                    12d14s20
                 end do                                                  7d9s21
                 ltmp1=ibcoff                                            2d19s21
                 ltmp2=ltmp1+ncsf(iarg)*nrow                             2d19s21
                 ibcoff=ltmp2+ncsf(iarg)*nrow                             2d19s21
                 call enough('hcssbk.  9',bc,ibc)
                 call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid1b,nrow,      2d19s21
     $              iwpb1b,iwpk1b,bc(ivtrans),ncsf(jarg),bc(ltmp2),     2d19s21
     $              ncsf(iarg),1d0,0d0,bc,ibc)                          11d10s22
                 do ii=0,nrow*ncsf(iarg)-1                               2d22s21
                  bc(idend+ii)=bc(idend+ii)+sum*bc(ltmp2+ii)             2d22s21
                 end do                                                  2d22s21
                 do ii2e=1,n2e                                           7d15s21
                  if(isymop(i2eop(1,ii2e)).eq.1)then                     7d15s21
                   icol=lga+irefo(lsa)*lga                                7d14s21
                   iden=idenjv(lsa,ii2e)+nrowi*icol                      7d15s21
                   ibc(ndenjv(lsa,ii2e)+icol)=1                          7d15s21
                   do ii=0,nrow*ncsf(iarg)-1                               2d22s21
                    bc(iden+ii)=bc(iden+ii)+bc(ltmp2+ii)                   2d22s21
                   end do                                                  2d22s21
                  end if                                                  7d14s21
                 end do                                                  7d15s21
                 ibcoff=ltmp1                                            2d22s21
                 do i=1,nok                                                    11d13s20
                  if(itest(i,3).eq.1)then                                  12d14s20
                   itestc=j1c                                              12d14s20
                   itesto=j1o                                              12d14s20
                   nopenk=nopenj                                                11d13s20
c
c     anihilate common
c
                   if(btest(itestc,itest(i,2)))then                             11d13s20
                    itestc=ibclr(itestc,itest(i,2))                             11d13s20
                    itesto=ibset(itesto,itest(i,2))                             11d13s20
                    karg=jarg-1                                                11d13s20
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,2))                             11d13s20
                    karg=jarg                                                  11d13s20
                    nopenk=nopenk-1                                             11d13s20
                   end if                                                       11d13s20
c
c     create ket
c
                   if(btest(itesto,nab4(2,1)))then                        12d14s20
                    itestc=ibset(itestc,nab4(2,1))                        12d14s20
                    itesto=ibclr(itesto,nab4(2,1))                        12d14s20
                    karg=karg+1                                                11d13s20
                    nopenk=nopenk-1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibset(itesto,nab4(2,1))                        12d14s20
                    nopenk=nopenk+1                                             11d13s20
                   end if                                                       11d13s20
                   nqq=karg+mdon-1
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                  2d6s23
                    call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,      2d22s21
     $               iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,2d22s21
     $                iwpk1b,ncsfmid1b,bc,ibc)                          11d14s22
                    call gandc(itestc,itesto,j1c,j1o,nopenk,nopenj,      2d22s21
     $            karg,jarg,ncsf,norbx,ixw1,ixw2,nnot2b,nab2b,iwpb2b,   2d22s21
     $               iwpk2b,ncsfmid2b,bc,ibc)                           11d14s22
                    if(nnot1b.eq.2.and.nnot2b.eq.2)then                   2d22s21
c     (dc|ba)
                     lsa=ism(nab2b(2))                                      12d14s20
                     lga=irel(nab2b(2))-1                                   12d14s20
                     lsb=ism(nab2b(1))                                      12d14s20
                     lgb=irel(nab2b(1))-1                                   12d14s20
                     lsc=ism(nab1b(2))                                      12d14s20
                     lgc=irel(nab1b(2))-1                                   12d14s20
                     lsd=ism(nab1b(1))                                      12d14s20
                     lgd=irel(nab1b(1))-1                                   12d14s20
                     xint=0d0                                            7d9s21
                     do ii2e=1,n2e                                       7d9s21
                      if(multh(lsa,lsb).eq.isymop(i2eop(2,ii2e)))then    7d16s21
                       iba=ixmt(lsb,i2eop(2,ii2e))+lgb+idoubo(lsb)        7d9s21
     $                    +nbasdws(lsb)*(lga+idoubo(lsa))               7d9s21
                       idc=ixmt(lsd,i2eop(1,ii2e))+lgd+idoubo(lsd)        7d12s21
     $                    +nbasdws(lsd)*(lgc+idoubo(lsc))               7d12s21
                       xint=xint+bc(iba)*bc(idc)                          7d12s21
                      end if                                             7d16s21
                     end do                                              7d9s21
                     xint=xint*phase2                                    7d12s21
                     ltmp1=ibcoff                                         2d22s21
                     ltmp2=ltmp1+nrow*ncsf(iarg)                          2d22s21
                     ibcoff=ltmp2+nrow*ncsf(iarg)                         2d22s21
                     call enough('hcssbk. 10',bc,ibc)
                     call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),       2d22s21
     $                 ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b,
     $                 bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc)11d10s22
                     do ii=0,nrow*ncsf(iarg)-1                            2d22s21
                      bc(idend+ii)=bc(idend+ii)+xint*bc(ltmp2+ii)         2d22s21
                     end do                                               2d22s21
                     ibcoff=ltmp1                                         2d22s21
                    end if                                                  12d14s20
                   end if                                                4d14s21
                  end if                                                   12d14s20
                 end do                                                    12d14s20
                 if(n2e.gt.0)then                                        8d23s21
                  itestc=j1c                                              12d14s20
                  itesto=j1o                                              12d14s20
                  nopenk=nopenj                                                11d13s20
c
c     anihilate common
c
                  if(btest(itestc,norbx))then                             11d13s20
                   itestc=ibclr(itestc,norbx)                             11d13s20
                   itesto=ibset(itesto,norbx)                             11d13s20
                   karg=jarg-1                                                11d13s20
                   nopenk=nopenk+1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibclr(itesto,norbx)                             11d13s20
                   karg=jarg                                                  11d13s20
                   nopenk=nopenk-1                                             11d13s20
                  end if                                                       11d13s20
c
c     create ket
c
                  if(btest(itesto,nab4(2,1)))then                         12d14s20
                   itestc=ibset(itestc,nab4(2,1))                         12d14s20
                   itesto=ibclr(itesto,nab4(2,1))                         12d14s20
                   karg=karg+1                                                11d13s20
                   nopenk=nopenk-1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibset(itesto,nab4(2,1))                         12d14s20
                   nopenk=nopenk+1                                             11d13s20
                  end if                                                       11d13s20
                  nqq=karg+mdon-1
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                   call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,        2d22s21
     $               iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,2d22s21
     $                iwpk1b,ncsfmid1b,bc,ibc)                          11d14s22
                   call gandc(itestc,itesto,j1c,j1o,nopenk,nopenj,        2d22s21
     $            karg,jarg,ncsf,norbx,ixw1,ixw2,nnot2b,nab2b,iwpb2b,   2d22s21
     $               iwpk2b,ncsfmid2b,bc,ibc)                           11d14s22
                   if(nnot1b.eq.2.and.nnot2b.eq.2)then                     2d22s21
                    ltmp1=ibcoff                                           2d22s21
                    ltmp2=ltmp1+nrow*ncsf(iarg)                            2d22s21
                    ibcoff=ltmp2+nrow*ncsf(iarg)                           2d22s21
                    call enough('hcssbk. 11',bc,ibc)
                    call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),         2d22s21
     $                 ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b, 2d22s21
     $                bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc) 11d10s22
                    lsb=ism(nab2b(1))                                       12d14s20
                    lgb=irel(nab2b(1))-1                                    12d14s20
                    lsc=ism(nab1b(2))                                       12d14s20
                    lgc=irel(nab1b(2))-1                                   2d22s21
                    icol=lgc+irefo(lsc)*lgb                                12d17s20
                    do ii2e=1,n2e                                         7d15s21
                     if(multh(lsb,jsbv).eq.isymop(i2eop(1,ii2e)))then         7d14s21
                      iden=idenkv(lsc,ii2e)+nrowi*icol                    7d15s21
                      ibc(ndenkv(lsc,ii2e)+icol)=1                        7d15s21
                      do ii=0,nrow*ncsf(iarg)-1                              2d22s21
                       bc(iden+ii)=bc(iden+ii)+bc(ltmp2+ii)                  2d22s21
                      end do                                                 2d22s21
                     end if                                                7d14s21
                    end do                                                7d15s21
                    ibcoff=ltmp1                                           2d22s21
                   end if                                                  12d17s20
                  end if                                                  4d14s21
                 end if                                                  8d23s21
                else if(n2e.gt.0)then                                   2d6s23
                 nnot=0                                                 2d6s23
                 ipssx=0                                                2d6s23
                 if(ndifs.eq.4.and.ndifb.eq.4)then                      2d6s23
                  nnot=4                                                2d6s23
                  ioxx(1)=1                                             2d6s23
                  ioxx(2)=1                                             2d6s23
                  do i=1,norbx                                          2d6s23
                   if(btest(gandcb,i))then                              2d6s23
                    if((btest(i1c,i).and.btest(j1o,i)).or.              2d6s23
     $                (btest(i1o,i).and..not.btest(j1c,i)))then         2d6s23
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
     $                  ((btest(j1c,i).and..not.btest(i1o,i)).or.       2d6s23
     $                  (btest(i1c,i).and..not.btest(j1o,i))))then      2d6s23
                     if(btest(i1c,i))iswap=1                            2d6s23
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
                   if(btest(j1c,nab4(1,2)).and.                         2d6s23
     $                  .not.btest(j1c,nab4(1,1)))nbt=1                 2d6s23
                  else                                                  2d6s23
                   nbt=0                                                2d6s23
                   if(btest(i1c,nab4(2,2)).and.                         2d6s23
     $                  .not.btest(i1c,nab4(2,1)))nbt=1                 2d6s23
                  end if                                                2d6s23
                  if(nbt.ne.0)then                                      2d6s23
                   nab4(1,1)=nab4(1,2)                                  2d6s23
                   nab4(2,1)=nab4(2,2)                                  2d6s23
                  end if                                                2d6s23
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                 2d6s23
                  nnot=3                                                2d6s23
                  do i=1,norbx                                          2d6s23
                   if(btest(gandcb,i))then                              2d6s23
                    if(btest(j1c,i))then                                2d6s23
                     nab4(1,1)=i                                        2d6s23
                     nab4(1,2)=i                                        2d6s23
                    else                                                2d6s23
                     nab4(2,1)=i                                        2d6s23
                     nab4(2,2)=i                                        2d6s23
                    end if                                              2d6s23
                   end if                                               2d6s23
                  end do                                                2d6s23
                 end if                                                 2d6s23
                 if(nnot.eq.3)then                                      2d6s23
                  ipssx=1                                                   12d8s20
                 else if(nnot.eq.4)then                                 2d6s23
                  ipssx=3                                                12d16s20
                 end if                                                     12d8s20
                 do ipss=1,ipssx                                            12d8s20
                  if(ipss.eq.1)then
                   iu1=1
                   iu2=1
                  else if(ipss.eq.2)then
                   iu1=1
                   iu2=2
                  else
                   iu1=2
                   iu2=1
                  end if                                                    12d8s20
                  itestc=j1c                                              12d14s20
                  itesto=j1o                                              12d14s20
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenj+1                                              11d13s20
                   karg=jarg-1                                             12d8s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenj-1                                              11d13s20
                   karg=jarg                                              12d14s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                   karg=karg+1                                                  11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  nqq=karg+mdon-1
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                   call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,        2d22s21
     $        iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,2d22s21
     $         ncsfmid1b,bc,ibc)                                        11d14s22
                   call gandc(itestc,itesto,j1c,j1o,nopenk,nopenj,        2d22s21
     $        karg,jarg,ncsf,norbx,ixw1,ixw2,nnot2b,nab2b,iwpb2b,iwpk2b,2d22s21
     $               ncsfmid2b,bc,ibc)                                  11d14s22
                   if(nnot1b.eq.2.and.nnot2b.eq.2)then                    2d22s21
                    nab1(1)=nab2b(2)                                       2d22s21
                    nab1(2)=nab2b(1)                                       2d22s21
                    nab2(1)=nab1b(2)                                       2d22s21
                    nab2(2)=nab1b(1)                                       2d22s21
                    lsa=ism(nab1(1))                                       12d14s20
                    lga=irel(nab1(1))-1                                    12d14s20
                    lsb=ism(nab1(2))                                       12d14s20
                    lgb=irel(nab1(2))-1                                    12d14s20
                    lsc=ism(nab2(1))                                       12d14s20
                    lgc=irel(nab2(1))-1                                    12d14s20
                    lsd=ism(nab2(2))                                       12d14s20
                    lgd=irel(nab2(2))-1                                    12d14s20
                    if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2)).and.  12d16s20
     $               max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))then  12d16s20
                     fact=0.5d0                                           12d16s20
                    else                                                  12d16s20
                     fact=1d0                                             12d16s20
                    end if                                                12d16s20
                    xint=0d0                                             7d12s21
                    do ii2e=1,n2e                                        7d12s21
                     if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)))then     7d16s21
                      iab=ixmt(lsa,i2eop(1,ii2e))+lga+idoubo(lsa)         7d12s21
     $                   +nbasdws(lsa)*(lgb+idoubo(lsb))                7d12s21
                      icd=ixmt(lsc,i2eop(2,ii2e))+lgc+idoubo(lsc)         7d12s21
     $                   +nbasdws(lsc)*(lgd+idoubo(lsd))                7d12s21
                      trm=bc(iab)*bc(icd)
                      orig=xint
                      xint=xint+bc(iab)*bc(icd)                           7d12s21
                     end if                                              7d16s21
                    end do                                               7d12s21
                    xint=xint*phase2                                     7d12s21
                    xuse=xint*fact                                        12d16s20
                    ltmp1=ibcoff                                          2d22s21
                    ltmp2=ltmp1+nrow*ncsf(iarg)                           2d22s21
                    ibcoff=ltmp2+nrow*ncsf(iarg)                           2d22s21
                    call enough('hcssbk. 12',bc,ibc)
                    call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),        2d22s21
     $                 ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b, 2d22s21
     $                 bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc)11d10s22
                    do ii=0,nrow*ncsf(iarg)-1                             2d22s21
                     bc(idend+ii)=bc(idend+ii)+xuse*bc(ltmp2+ii)          2d22s21
                    end do                                                2d22s21
                    ibcoff=ltmp1                                          2d22s21
                    if(ipss.eq.2)go to 2001                               12d16s20
                   end if                                                  12d14s20
                  end if                                                 4d14s21
                 end do                                                   12d14s20
 2001            continue                                                12d16s20
                end if                                                    12d14s20
               end if                                                   2d6s23
              end if                                                     12d14s20
              jvec=jvec+ncsf(jarg)*nrow                                  12d14s20
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
              ibcoff=ltmp1+nrow*ncsf(iarg)                              2d22s21
              call enough('hcssbk. 13',bc,ibc)
              do j=0,nrow-1                                             2d22s21
               do ii=0,ncsf(iarg)-1                                     2d22s21
                ij=idend+ii+ncsf(iarg)*j                                2d22s21
                ji=ltmp1+j+nrow*ii                                      2d22s21
                bc(ji)=bc(ij)                                           2d22s21
               end do                                                   2d22s21
              end do                                                    2d22s21
              do i=0,ncsf(iarg)-1                                       12d15s20
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
             call enough('hcssbk. 14',bc,ibc)
             do j=0,nrow-1                                              2d22s21
              do ii=0,ncsf(iarg)-1                                      2d22s21
               ij=idenh+ii+ncsf(iarg)*j                                 2d22s21
               ji=ltmp1+j+nrow*ii                                       2d22s21
               bc(ji)=bc(ij)                                            2d22s21
              end do                                                    2d22s21
             end do                                                     2d22s21
             do j=0,nrowim                                              2d22s21
              bc(idenh+j)=bc(ltmp1+j)                                   2d22s21
             end do                                                     2d22s21
             do isc=1,nsymb                                             2d22s21
              isd=multh(isymop(lxmt),multh(isc,ijsbv))                  7d14s21
              nn=irefo(isc)*irefo(isd)                                  2d22s21o
              do i=0,nn-1                                               2d22s21
               do ii2e=1,n2e                                            7d15s21
                if(ibc(ndenk(isc,ii2e)+i).ne.0)then                     7d15s21
                 iden=idenk(isc,ii2e)+nrowi*i                           7d15s21
                 do j=0,nrow-1                                              2d22s21
                  do ii=0,ncsf(iarg)-1                                      2d22s21
                   ij=iden+ii+ncsf(iarg)*j                                 2d22s21
                   ji=ltmp1+j+nrow*ii                                       2d22s21
                   bc(ji)=bc(ij)                                            2d22s21
                  end do                                                    2d22s21
                 end do                                                     2d22s21
                 do j=0,nrowim                                              2d22s21
                  bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                 end do                                                     2d22s21
                end if                                                   2d22s21
                if(ibc(ndenkv(isc,ii2e)+i).ne.0)then                           2d22s21
                 iden=idenkv(isc,ii2e)+nrowi*i                                 2d22s21
                 do j=0,nrow-1                                              2d22s21
                  do ii=0,ncsf(iarg)-1                                      2d22s21
                   ij=iden+ii+ncsf(iarg)*j                                 2d22s21
                   ji=ltmp1+j+nrow*ii                                       2d22s21
                   bc(ji)=bc(ij)                                            2d22s21
                  end do                                                    2d22s21
                 end do                                                     2d22s21
                 do j=0,nrowim                                              2d22s21
                  bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                 end do                                                     2d22s21
                end if                                                   2d22s21
               end do                                                   7d15s21
              end do                                                    2d22s21
              do i=0,nn-1                                               2d22s21
               do ii2e=1,n2e                                            7d15s21
                if(ibc(ndenj(isc,ii2e)+i).ne.0)then                     7d15s21
                 iden=idenj(isc,ii2e)+nrowi*i                           7d15s21
                 do j=0,nrow-1                                              2d22s21
                  do ii=0,ncsf(iarg)-1                                      2d22s21
                   ij=iden+ii+ncsf(iarg)*j                                 2d22s21
                   ji=ltmp1+j+nrow*ii                                       2d22s21
                   bc(ji)=bc(ij)                                            2d22s21
                  end do                                                    2d22s21
                 end do                                                     2d22s21
                 do j=0,nrowim                                              2d22s21
                  bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                 end do                                                     2d22s21
                end if                                                   2d22s21
                if(ibc(ndenjv(isc,ii2e)+i).ne.0)then                           2d22s21
                 iden=idenjv(isc,ii2e)+nrowi*i                                 2d22s21
                 do j=0,nrow-1                                              2d22s21
                  do ii=0,ncsf(iarg)-1                                      2d22s21
                   ij=iden+ii+ncsf(iarg)*j                                 2d22s21
                   ji=ltmp1+j+nrow*ii                                       2d22s21
                   bc(ji)=bc(ij)                                            2d22s21
                  end do                                                    2d22s21
                 end do                                                     2d22s21
                 do j=0,nrowim                                              2d22s21
                  bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                 end do                                                     2d22s21
                end if                                                   2d22s21
               end do                                                   7d15s21
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
             mcol=nrootu*ncsf(iarg)                                     12d17s20
             itmpp=ibcoff                                               12d17s20
             ibcoff=itmpp+n1here*mcol                                   12d18s20
             call enough('hcssbk. 15',bc,ibc)
             do i=itmpp,ibcoff-1                                        12d17s20
              bc(i)=0d0                                                 12d17s20
             end do                                                     12d17s20
             ihit=0
             if(ndenh.ne.0)then                                         12d15s20
c     g r,v,i = sum v' denh r,v',i h0v v'
              itmp=ibcoff                                               12d17s20
              ibcoff=itmp+n1here*nherev                                 12d18s20
              call enough('hcssbk. 16',bc,ibc)
              do i=itmp,ibcoff-1                                        12d17s20
               bc(i)=0d0                                                12d17s20
              end do                                                    12d17s20
              i10=i1s                                                   12d15s20
              i1n=nvirt(isbv)                                           12d15s20
              do i2=i2s,i2e                                             12d15s20
               i2m=i2-1                                                 12d17s20
               i2p=i2+1                                                 12d17s20
               ih0=ixmtf(isbv)+irefo(isbv)+idoubo(isbv)-1               7d9s21
     $              +nbasdws(isbv)*(i2m+irefo(jsbv)+idoubo(jsbv))       7d9s21
               if(i2.eq.i2e)i1n=i1e                                     12d15s20
               jtmp=itmp-ii10+n1here*(i2-i2s)                           12d18s20
               if(isymop(lxmt).eq.1)then                                7d19s21
                do i1=i10,min(i1n,i2m)                                   12d17s20
                 bc(jtmp+i1)=bc(ih0+i1)*phase                            7d9s21
                end do                                                   12d17s20
                do i1=max(i2p,i10),i1n                                   12d17s20
                 bc(jtmp+i1)=bc(ih0+i1)*phase                            7d9s21
                end do                                                   12d17s20
               else                                                     7d19s21
                do i1=i10,i1n                                           7d19s21
                 bc(jtmp+i1)=bc(ih0+i1)*phase                            7d9s21
                end do                                                   12d17s20
               end if                                                   7d19s21
               i10=1                                                    12d15s20
              end do                                                    12d15s20
              itmpd=ibcoff                                              12d17s20
              ibcoff=itmpd+nrowi                                        12d17s20
              call enough('hcssbk. 17',bc,ibc)
              jtmp=idenh                                                12d17s20
              do i=0,ncsf(iarg)-1                                       12d17s20
               do ir=0,nrootm                                           1d28s21
                iad=itmpd+nherev*(ir+nrootu*i)                          1d28s21
                do iv=0,nherev-1                                        1d28s21
                 bc(iad+iv)=bc(jtmp+iv)                                 1d28s21
                end do                                                  1d28s21
                jtmp=jtmp+nherev                                        1d28s21
               end do                                                   12d17s20
              end do                                                    12d17s20
              ihit=1
              call dgemm('n','n',n1here,mcol,nherev,1d0,                12d18s20
     $             bc(itmp),n1here,bc(itmpd),nherev,1d0,                12d18s20
     $             bc(itmpp),n1here,                                    12d18s20
     d' hcssbk.  1')
              ibcoff=itmp
             end if                                                     12d15s20
             do isc=1,nsymb                                             12d15s20
              isd=multh(isymop(lxmt),multh(isc,ijsbv))                  7d14s21
              nn=irefo(isc)*irefo(isd)                                  12d15s20
              do ii2e=1,n2e                                             7d15s21
               ikeep(ii2e)=ibcoff                                       7d15s21
               ibcoff=ikeep(ii2e)+nn                                    7d15s21
               call enough('hcssbk. 18',bc,ibc)
               nokk(ii2e)=0                                             7d15s21
               nokkv(ii2e)=0                                            7d15s21
               do i=0,nn-1                                               12d17s20
                ibc(ikeep(ii2e)+i)=0                                    7d15s21
                if(ibc(ndenkv(isc,ii2e)+i).ne.0)then                    7d15s21
                 iden=idenkv(isc,ii2e)+nrowi*i                          7d15s21
                 rms=1d0                                                 12d17s20
                 if(ibc(ndenk(isc,ii2e)+i).ne.0)then                    7d15s21
                  jden=idenk(isc,ii2e)+nrowi*i                          7d15s21
                  rms=0d0                                                12d17s20
                  do k=0,nrowim                                          12d17s20
                   rms=rms+(bc(jden+k)-bc(iden+k))**2                    12d17s20
                  end do                                                 12d17s20
                  rms=sqrt(rms/dfloat(nrowi))                            12d17s20
                 end if                                                  12d17s20
                 if(rms.gt.1d-10)then                                   7d16s21
                  jden=idenkv(isc,ii2e)+nrowi*nokkv(ii2e)               7d16s21
                  do k=0,nrowim                                          12d17s20
                   bc(jden+k)=bc(iden+k)                                 12d17s20
                  end do                                                 12d17s20
                  ibc(ndenkv(isc,ii2e)+nokkv(ii2e))=i                   7d15s21
                  nokkv(ii2e)=nokkv(ii2e)+1                             7d15s21
                 else                                                    12d17s20
                  ibc(ikeep(ii2e)+i)=1                                  7d15s21
                 end if                                                  12d17s20
                end if                                                   12d17s20
                if(ibc(ndenk(isc,ii2e)+i).ne.0)then                     7d15s21
                 iden=idenk(isc,ii2e)+nrowi*i                           7d15s21
                 jden=idenk(isc,ii2e)+nrowi*nokk(ii2e)                  7d15s21
                 do k=0,nrowim                                           12d17s20
                  bc(jden+k)=bc(iden+k)                                  12d17s20
                 end do                                                  12d17s20
                 ibc(ndenk(isc,ii2e)+nokk(ii2e))=i                      7d15s21
                 nokk(ii2e)=nokk(ii2e)+1                                7d15s21
                end if                                                   12d17s20
               end do                                                    12d17s20
              end do                                                    7d15s21
              do ii2e=1,n2e                                             7d15s21
               if(nokk(ii2e).gt.0)then                                  7d15s21
                itmp=ibcoff                                               12d17s20
                ibcoff=itmp+nokk(ii2e)*n1here*nherev                    7d16s21
                call enough('hcssbk. 19',bc,ibc)
                ihit=1                                                   12d17s20
                do i=itmp,ibcoff-1                                       12d17s20
                 bc(i)=0d0                                               12d17s20
                end do                                                   12d17s20
                do i=0,nokk(ii2e)-1                                             12d17s20
                 icol=ibc(ndenk(isc,ii2e)+i)                                  12d17s20
                 id=icol/irefo(isc)                                      7d12s21
                 ic=icol-irefo(isc)*id                                   7d12s21
                 i10=i1s                                                   12d15s20
                 i1n=nvirt(isbv)                                           12d15s20
                 if(ibc(ikeep(ii2e)+icol).eq.1.or.isbv.ne.jsbv)then            12d17s20
                  do i2=i2s,i2e                                             12d15s20
                   i2p=i2+idoubo(jsbv)+irefo(jsbv)-1                     7d12s21
                   if(i2.eq.i2e)i1n=i1e                                     12d15s20
                   jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)               12d18s20
                   icv=ixmt(isc,i2eop(2,ii2e))+ic+idoubo(isc)            7d12s21
     $                 +nbasdws(isc)*i2p                                7d12s21
                   factk=phase2*bc(icv)                                 7d12s21
                   i1p=idoubo(isbv)+irefo(isbv)-1                       7d14s21
                   ivd=ixmt(isbv,i2eop(1,ii2e))+i1p                     7d14s21
     $                  +nbasdws(isbv)*(id+idoubo(isd))                 7d14s21
                   do i1=i10,i1n                                         12d17s20
                    bc(jtmp+i1)=bc(jtmp+i1)+factk*bc(ivd+i1)            7d14s21
                   end do                                                12d17s20
                   i10=1                                                 12d17s20
                  end do                                                 12d17s20
                 else                                                    12d17s20
                  do i2=i2s,i2e                                             12d15s20
                   i2m=i2-1                                              12d17s20
                   i2p=i2+1                                              12d17s20
                   i2pp=i2+idoubo(jsbv)+irefo(jsbv)-1                    7d14s21
                   if(i2.eq.i2e)i1n=i1e                                     12d15s20
                   jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)               12d18s20
                   icv=ixmt(isc,i2eop(2,ii2e))+ic+idoubo(isc)            7d12s21
     $                 +nbasdws(isc)*i2pp                               7d14s21
                   factk=phase2*bc(icv)                                 7d12s21
                   i1p=idoubo(isbv)+irefo(isbv)-1                       7d14s21
                   ivd=ixmt(isbv,i2eop(1,ii2e))
     $                   +i1p+nbasdws(isbv)*(id+idoubo(isd))            7d14s21
                   do i1=i10,min(i1n,i2m)                                12d18s20
                    bc(jtmp+i1)=bc(jtmp+i1)+factk*bc(ivd+i1)            7d14s21
                   end do                                                12d18s20
                   do i1=max(i10,i2p),i1n                                12d18s20
                    bc(jtmp+i1)=bc(jtmp+i1)+factk*bc(ivd+i1)            7d14s21
                   end do                                                12d18s20
                   i10=1                                                 12d17s20
                  end do                                                 12d17s20
                 end if                                                  12d17s20
                end do                                                   12d17s20
                itmpd=ibcoff                                             12d17s20
                ibcoff=itmpd+nherev*nokk(ii2e)*nrootu*ncsf(iarg)        7d15s21
                call enough('hcssbk. 20',bc,ibc)
                jden=idenk(isc,ii2e)                                    7d15s21
                do i=0,nokk(ii2e)-1                                     7d15s21
                 do m=0,ncsf(iarg)-1                                     12d17s20
                  do ir=0,nrootm                                         1d28s21
                   iad=itmpd+nherev*(i+nokk(ii2e)*(ir+nrootu*m))        7d15s21
                   do iv=0,nherev-1                                      1d28s21
                    bc(iad+iv)=bc(jden+iv)                               1d28s21
                   end do                                                12d17s20
                   jden=jden+nherev                                      1d28s21
                  end do                                                 12d17s20
                 end do                                                  12d17s20
                end do                                                   12d17s20
                mm=nherev*nokk(ii2e)                                    7d15s21
                call dgemm('n','n',n1here,mcol,mm,1d0,                   12d18s20
     $              bc(itmp),n1here,bc(itmpd),mm,1d0,bc(itmpp),         12d18s20
     $              n1here,                                             12d18s20
     d' hcssbk.  2')
                ibcoff=itmp                                              12d17s20
               end if                                                    12d17s20
               if(nokkv(ii2e).gt.0.and.isc.eq.isd)then                  7d15s21
                itmp=ibcoff                                              12d17s20
                ibcoff=itmp+nherev*nokkv(ii2e)                          7d15s21
                call enough('hcssbk. 21',bc,ibc)
                do i=itmp,ibcoff-1                                       12d17s20
                 bc(i)=0d0                                               12d17s20
                end do                                                   12d17s20
                do i=0,nokkv(ii2e)-1                                    7d15s21
                 icol=ibc(ndenkv(isc,ii2e)+i)                           7d15s21
                 id=icol/irefo(isc)                                      7d12s21
                 ic=icol-irefo(isc)*id                                   7d12s21
                 i10=i1s                                                   12d15s20
                 i1n=nvirt(isbv)                                           12d15s20
                 do i2=i2s,i2e                                           12d17s20
                  i2p=i2+idoubo(jsbv)+irefo(jsbv)-1                      7d12s21
                  if(i2.eq.i2e)i1n=i1e                                   12d17s20
                  jtmp=itmp+i2-i2s+nherev*i                              12d17s20
                  if(i2.ge.i10.and.i2.le.i1n)then                        12d18s20
                   ivc=ixmt(jsbv,i2eop(1,ii2e))+i2p                     7d16s21
     $                 +nbasdws(jsbv)*(ic+idoubo(isc))                  7d16s21
                   idv=ixmt(isd,i2eop(2,ii2e))+id+idoubo(isd)           7d12s21
     $                  +nbasdws(isd)*i2p                               7d12s21
                   bc(jtmp)=bc(jtmp)+phase2*bc(ivc)*bc(idv)             7d16s21
                  end if                                                 12d18s20
                  i10=1                                                  12d17s20
                 end do                                                  12d17s20
                end do                                                   12d17s20
                jtmp=itmp-i2s                                            12d17s20
                jtmpd=idenkv(isc,ii2e)                                  7d15s21
                do i=0,nokkv(ii2e)-1                                    7d15s21
                 do m=0,ncsf(iarg)-1                                     12d17s20
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
              end do                                                    7d15s21
              do ii2e=1,n2e                                             7d15s21
               ikeep(ii2e)=ibcoff                                       7d15s21
               ibcoff=ikeep(ii2e)+nn                                    7d15s21
               call enough('hcssbk. 22',bc,ibc)
               nok=0                                                    7d15s21
               nokv=0                                                   7d15s21
               do i=0,nn-1                                               12d17s20
                ibc(ikeep(ii2e)+i)=0                                    7d15s21
                if(ibc(ndenjv(isc,ii2e)+i).ne.0)then                    7d15s21
                 iden=idenjv(isc,ii2e)+nrowi*i                          7d15s21
                 rms=1d0                                                 12d17s20
                 if(ibc(ndenj(isc,ii2e)+i).ne.0)then                    7d15s21
                  jden=idenj(isc,ii2e)+nrowi*i                          7d15s21
                  rms=0d0                                                12d17s20
                  do k=0,nrowim                                          12d17s20
                   rms=rms+(bc(jden+k)-bc(iden+k))**2                    12d17s20
                  end do                                                 12d17s20
                  rms=sqrt(rms/dfloat(nrowi))                            12d17s20
                 end if                                                  12d17s20
                 if(rms.gt.1d-10)then
                  jden=idenjv(isc,ii2e)+nrowi*nokv                      7d15s21
                  do k=0,nrowim                                          12d17s20
                   bc(jden+k)=bc(iden+k)                                 12d17s20
                  end do                                                 12d17s20
                  ibc(ndenjv(isc,ii2e)+nokv)=i                          7d15s21
                  nokv=nokv+1                                           7d15s21
                 else                                                    12d17s20
                  ibc(ikeep(ii2e)+i)=1                                  7d15s21
                 end if                                                  12d17s20
                end if                                                   12d17s20
                if(ibc(ndenj(isc,ii2e)+i).ne.0)then                     7d15s21
                 iden=idenj(isc,ii2e)+nrowi*i                           7d15s21
                 jden=idenj(isc,ii2e)+nrowi*nok                         7d15s21
                 do k=0,nrowim                                           12d17s20
                  bc(jden+k)=bc(iden+k)                                  12d17s20
                 end do                                                  12d17s20
                 ibc(ndenj(isc,ii2e)+nok)=i                             7d15s21
                 nok=nok+1                                              7d15s21
                end if                                                   12d17s20
               end do                                                    12d17s20
               if(nok.gt.0)then                                          12d17s20
                itmp=ibcoff                                               12d17s20
                ibcoff=itmp+nok*n1here*nherev                            12d18s20
                call enough('hcssbk. 23',bc,ibc)
                ihit=1                                                   12d17s20
                do i=itmp,ibcoff-1                                       12d17s20
                 bc(i)=0d0                                               12d17s20
                end do                                                   12d17s20
                do i=0,nok-1                                             12d17s20
                 icol=ibc(ndenj(isc,ii2e)+i)                                  12d17s20
                 id=icol/irefo(isc)                                      7d12s21
                 ic=icol-id*irefo(isc)                                   7d12s21
                 icd=ixmt(isc,i2eop(1,ii2e))+ic+idoubo(isc)             7d12s21
     $                +nbasdws(isc)*(id+idoubo(isd))                    7d12s21
                 factj=bc(icd)*phase2                                   7d12s21
                 i10=i1s                                                   12d15s20
                 i1n=nvirt(isbv)                                           12d15s20
                 if(ibc(ikeep(ii2e)+icol).eq.1.or.isbv.ne.jsbv)then     7d15s21
                  do i2=i2s,i2e                                             12d15s20
                   i2p=i2+idoubo(jsbv)+irefo(jsbv)-1                    7d12s21
                   if(i2.eq.i2e)i1n=i1e                                     12d15s20
                   jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)              12d18s20
                   ivu=ixmt(isbv,i2eop(2,ii2e))-1+idoubo(isbv)          7d12s21
     $                   +irefo(isbv)                                   7d12s21
     $                    +nbasdws(isbv)*i2p                            7d12s21
                   do i1=i10,i1n                                         12d17s20
                    bc(jtmp+i1)=bc(jtmp+i1)+factj*bc(ivu+i1)            7d12s21
                   end do                                                12d17s20
                   i10=1                                                 12d17s20
                  end do                                                 12d17s20
                 else                                                    12d17s20
                  do i2=i2s,i2e                                             12d15s20
                   i2m=i2-1                                              12d17s20
                   i2p=i2+1                                              12d17s20
                   if(i2.eq.i2e)i1n=i1e                                     12d15s20
                   jtmp=itmp-ii10+n1here*(i2-i2s+nherev*i)              12d18s20
                   ivu=ixmt(isbv,i2eop(2,ii2e))-1+idoubo(isbv)
     $                  +irefo(isbv)                                    7d16s21
     $                    +nbasdws(isbv)*(i2-1+irefo(jsbv)+idoubo(jsbv))7d16s21
                   do i1=i10,min(i1n,i2m)                                12d18s20
                    bc(jtmp+i1)=bc(jtmp+i1)+factj*bc(ivu+i1)            7d12s21
                   end do                                                12d18s20
                   do i1=max(i10,i2p),i1n                                12d18s20
                    bc(jtmp+i1)=bc(jtmp+i1)+factj*bc(ivu+i1)            7d12s21
                   end do                                                12d18s20
                   i10=1                                                 12d17s20
                  end do                                                 12d17s20
                 end if                                                  12d17s20
                end do                                                  7d12s21
                itmpd=ibcoff                                             12d17s20
                ibcoff=itmpd+nherev*nok*nrootu*ncsf(iarg)                12d17s20
                call enough('hcssbk. 24',bc,ibc)
                jden=idenj(isc,ii2e)                                          12d17s20
                do i=0,nok-1                                             12d17s20
                 do m=0,ncsf(iarg)-1                                     12d17s20
                  do ir=0,nrootm                                         12d17s20
                   iad=itmpd+nherev*(i+nok*(ir+nrootu*m))                1d28s21
                   do iv=0,nherev-1                                       12d17s20
                    bc(iad+iv)=bc(jden+iv)                               1d28s21
                   end do                                                12d17s20
                   jden=jden+nherev                                      1d28s21
                  end do                                                 12d17s20
                 end do                                                  12d17s20
                end do                                                   12d17s20
                mm=nherev*nok                                            12d17s20
                call dgemm('n','n',n1here,mcol,mm,1d0,                   12d18s20
     $              bc(itmp),n1here,bc(itmpd),mm,1d0,bc(itmpp),n1here,  12d18s20
     d' hcssbk.  3')
                ibcoff=itmp                                              12d17s20
                if(nokv.gt.0)then                                         12d17s20
                 itmp=ibcoff                                              12d17s20
                 ibcoff=itmp+nherev*nokv                                  12d17s20
                 call enough('hcssbk. 25',bc,ibc)
                 do i=itmp,ibcoff-1                                       12d17s20
                  bc(i)=0d0                                               12d17s20
                 end do                                                   12d17s20
                 do i=0,nokv-1                                            12d17s20
                  icol=ibc(ndenjv(isc,ii2e)+i)                                 12d17s20
                  id=icol/irefo(isc)                                     7d12s21
                  ic=icol-id*irefo(isc)                                  7d12s21
                  icd=ixmt(isc,i2eop(1,ii2e))+ic+idoubo(isc)             7d12s21
     $                +nbasdws(isc)*(id+idoubo(isd))                    7d12s21
                  factj=bc(icd)*phase2                                  7d12s21
                  i10=i1s                                                   12d15s20
                  i1n=nvirt(isbv)                                           12d15s20
                  do i2=i2s,i2e                                           12d17s20
                   i2p=i2+idoubo(isbv)+irefo(isbv)-1                    7d12s21
                   if(i2.eq.i2e)i1n=i1e                                   12d17s20
                   jtmp=itmp+i2-i2s+nherev*i                              12d17s20
                   if(i2.ge.i10.and.i2.le.i1n)then                        12d18s20
                    ivu=ixmt(isbv,i2eop(2,ii2e))+i2p*(nbasdws(isbv)+1)  7d12s21
                    bc(jtmp)=bc(jtmp)+factj*bc(ivu)                     7d12s21
                   end if                                                 12d18s20
                   i10=1                                                  12d17s20
                  end do                                                  12d17s20
                 end do                                                   12d17s20
                 jtmp=itmp-i2s                                            12d17s20
                 jtmpd=idenjv(isc,ii2e)                                 7d15s21
                 do i=0,nokv-1                                            12d17s20
                  do m=0,ncsf(iarg)-1                                     12d17s20
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
              end do                                                    7d15s21
             end do                                                     12d15s20
             if(ihit.ne.0)then                                          12d17s20
              jjgg=jgg                                                  12d17s20
              do i=0,ncsf(iarg)-1                                       12d17s20
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
             jgg=jgg+ngg*ncsf(iarg)                                     12d16s20
            end do                                                       12d14s20
            ibcoff=ivecbc                                               2d12s21
           end if                                                       12d14s20
          end do                                                          12d12s20
         end do                                                           12d12s20
         i18=1                                                          12d15s20
         i28=ngg
         i38=1                                                          12d15s20
         i48=nff1b(ncloip,isb,1)*ncsf(iarg)                             7d12s21
         if(itransgg.ne.0)then                                          2d12s21
c
c     igg is virt,root,... so transpose to root,virt,...
          itmp=ibcoff                                                    1d27s21
          ibcoff=itmp+ngg                                                1d27s21
          call enough('hcssbk. 26',bc,ibc)
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
          nwiacc=nwiacc+i28*i48                                         11d21s22
          last8(1)=ihsdiagb(ncloip,isb)                                 8d21s21
          last8(2)=i48                                                  2d17s21
          itransgg=0                                                     2d12s21
         else                                                           2d12s21
          nacc=0                                                        2d12s21
         end if                                                         2d12s21
        end if                                                           12d12s20
       end do                                                           12d14s20
      end do                                                            12d12s20
      if(nacc.ne.0)call ddi_done(ibc(iacc),nacc)                        11d21s22
      ibcoff=iacc                                                       2d12s21
      return                                                            12d12s20
      end                                                               12d12s20
