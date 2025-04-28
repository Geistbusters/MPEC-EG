c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcssd12(ihsdiag,nff1,iff1,ncsf,nec,mdon,mdoo,nsymb,    5d8s23
     $     multh,ixw1,ixw2,ih0av,nh0av,ioooo,jmats,kmats,nvirt,nrootu,  5d8s23
     $     ism,irel,irefo,isymmrci,norb,lprt,maxbx,tdenss,tovr,         5d12s23
     $     bc,ibc,igoal)                                                      11d10s22
      implicit real*8 (a-h,o-z)
      external second                                                   2d18s21
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48,i1c,i1o,j1c,j1o,12d14s20
     $     itestc,itesto,last8(2),gandcc,gandco,gandcb                  11d1s22
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d19s21
      logical lkeep,lprt,lchoice                                        3d17s21
      dimension nff1(mdoo+1,nsymb,2),ncsf(*),multh(8,8),ih0av(*),       12d12s20
     $     nh0av(*),jmats(*),kmats(*),nvirt(*),ism(*),irel(*),
     $     itest(32,3),isorb1(32),idorb1(32),jsorb1(32),jdorb1(32),
     $     nab4(2,3),idenj(8),idenk(8),irefo(*),iff1(*),ndenj(8),       12d14s20
     $     ndenk(8),ioooo(*),idenjv(8),idenkv(8),ndenjv(8),ndenkv(8),   2d12s21
     $     ivec(5,8),ivecp(5,8),nvecp(8),ircv(5,8),nrcv(5,8),           2d12s21
     $     itransv(5,8),ioxx(2)                                         5d12s23
      include "common.store"                                            12d12s20
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      loop=0
      loopx=18800000
      if(lprt)                                                          1d25s21
     $    write(6,*)('hello my name is hcssd12'),isymmrci,nec,mdon,mdoo,
     $     nsymb,ixw1,ixw2,norb
c
c     can we hold all of singles vector that we need in memory?
c
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
      ntri=(nrootu*(nrootu+1))/2                                        5d12s23
      ivv=ibcoff                                                        5d12s23
      ibcoff=ivv+ntri                                                   5d12s23
      call enough('hcssd12.ivv',bc,ibc)                                 5d12s23
      last8(1)=-1                                                       2d8s21
      igg=ibcoff                                                        2d12s21
      ibcoff=igg+maxbx                                                  1d30s21
      nacc=0                                                            1d30s21
      itransgg=0                                                        2d12s21
      call enough('hcssd12.1',bc,ibc)                                    5d12s23
      nrootm=nrootu-1                                                   12d15s20
      nsing=0                                                           12d15s20
      loop=0
      norbx=norb+1                                                      12d14s20
      norbxx=norbx+1                                                    12d14s20
      do isb=1,nsymb                                                    12d12s20
       isbv=multh(isb,isymmrci)                                         12d12s20
       do jsb=1,nsymb                                                   2d12s21
        nvecp(jsb)=0                                                    2d12s21
       end do                                                           2d12s21
       do ncloip=mdon+1,mdoo+1                                           12d12s20
        if(min(nff1(ncloip,isb,1),nvirt(isbv)).gt.0)then                 12d14s20
         ncloi=ncloip-1                                                  12d12s20
         iarg=ncloip-mdon                                                12d12s20
         nsing=nsing+nff1(ncloip,isb,1)*ncsf(iarg)*nvirt(isbv)          12d15s20
         ncolti=nff1(ncloip,isb,1)*ncsf(iarg)*nvirt(isbv)*nrootu         12d12s20
         i18=1                                                          5d12s23
         i28=nrootu*nvirt(isbv)                                         5d12s23
         i48=nff1(ncloip,isb,1)*ncsf(iarg)                              5d12s23
         if(ncolti.gt.maxbx)then                                        5d12s23
          call dws_synca                                                5d12s23
          call dws_finalize                                             5d12s23
          stop 'hcssd12:maxbx'                                          5d12s23
         end if                                                         5d12s23
         igcpy=ibcoff                                                   5d12s23
         ibcoff=igcpy+ncolti                                            5d12s23
         call enough('hcssd12.gcpy',bc,ibc)                             5d12s23
         ngg=nvirt(isbv)*nrootu                                         12d15s20
         nnc=nff1(ncloip,isb,1)*ncsf(iarg)                              1d28s21
         call ddi_get(bc,ibc,ihsdiag(ncloip,isb,1),i18,i28,i18,i48,     5d12s23
     $        bc(igcpy))                                                5d12s23
         do icol=0,nnc-1                                                5d18s23
          do iv=0,nvirt(isbv)-1                                         5d12s23
           iad1=igcpy+nrootu*(iv+nvirt(isbv)*icol)                      5d12s23
           iad2=igg+iv+ngg*icol                                         5d12s23
           do ir=0,nrootm                                               5d12s23
            bc(iad2)=bc(iad1+ir)                                        5d12s23
            iad2=iad2+nvirt(isbv)                                       5d12s23
           end do                                                       5d12s23
          end do                                                        5d12s23
         end do                                                         5d12s23
         ibcoff=igcpy                                                   5d12s23
         nopeni=nec-2*ncloi                                              12d12s20
         do jsb=1,nsymb                                                   12d12s20
          jsbv=multh(jsb,isymmrci)                                        12d12s20
          i28=nvirt(jsbv)*nrootu                                         12d14s20
          call ilimts(nvirt(isbv),nvirt(jsbv),mynprocg,mynowprog,il,ih, 12d14s20
     $         i1s,i1e,i2s,i2e)                                         12d14s20
          nhere=ih+1-il                                                 12d15s20
          nherev=i2e+1-i2s                                              12d14s20
          if(jsb.eq.isb)then                                            12d15s20
           if(mynowprog.lt.mynprocg-1)then                              12d15s20
            ip=mynowprog+1                                              12d15s20
            call ilimts(nvirt(isbv),nvirt(jsbv),mynprocg,ip,jl,jh,      12d15s20
     $         j1s,j1e,j2s,j2e)                                         12d15s20
            i2top=j2s-1                                                 12d15s20
            nherevu=i2top+1-i2s                                         5d12s23
           else                                                         12d15s20
            i2top=i2e                                                   12d15s20
            nherevu=nherev                                              5d12s23
           end if                                                       12d15s20
          end if                                                        12d15s20
          do nclojp=max(mdon+1,ncloip-2),min(mdoo+1,ncloip+2)            12d12s20
           if(min(nff1(nclojp,jsb,1),nherev).gt.0)then                  12d14s20
            ncloj=nclojp-1                                               12d12s20
            jarg=nclojp-mdon                                             12d12s20
            nopenj=nec-2*ncloj                                           12d14s20
            nrow=nherev*nrootu                                          12d14s20
            ncoltj=nff1(nclojp,jsb,1)*ncsf(jarg)*nrow                    12d14s20
            ivecbc=ibcoff                                               2d12s21
            i48=nff1(nclojp,jsb,1)*ncsf(jarg)                           12d14s20
            ncol=i48                                                    2d12s21
c
c     has vector alread been loaded?
c
            nminc=mdoo+1                                                2d12s21
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
            call ddi_iget(bc,ibc,ihsdiag(nclojp,jsb,1),i18,i28,i38,i48, 11d15s22
     $           bc(ivec(nuse,jsb)),ibc(ircv(nuse,jsb)),nrcv(nuse,jsb)) 2d12s21
   22       continue                                                    2d12s21
            nrowi=nrow*ncsf(iarg)                                        12d14s20
            nrowim=nrowi-1                                               12d14s20
            idenh=ibcoff                                                 12d14s20
            idend=idenh+nrowi                                           12d14s20
            ibcoff=idend+nrowi                                          12d14s20
            ijsbv=multh(isbv,jsbv)                                       12d14s20
            do isc=1,nsymb                                               12d14s20
             isd=multh(isc,ijsbv)                                        12d14s20
             idenk(isc)=ibcoff                                          12d14s20
             nn=irefo(isc)*irefo(isd)                                   12d14s20
             ndenk(isc)=idenk(isc)+nn*nrowi                             12d14s20
             idenkv(isc)=ndenk(isc)+nn                                  12d17s20
             ndenkv(isc)=idenkv(isc)+nn*nrowi                           12d17s20
             ibcoff=ndenkv(isc)+nn                                      12d17s20
             if(isd.ge.isc)then                                         12d14s20
              if(isd.eq.isc)nn=(irefo(isc)*(irefo(isc)+1))/2            12d14s20
              idenj(isc)=ibcoff                                           12d14s20
              ndenj(isc)=idenj(isc)+nn*nrowi                            12d14s20
              idenjv(isc)=ndenj(isc)+nn                                 12d17s20
              ndenjv(isc)=idenjv(isc)+nn*nrowi                          12d17s20
              ibcoff=ndenjv(isc)+nn                                     12d17s20
             end if                                                     12d14s20
            end do                                                       12d14s20
            call enough('hcssd12.2',bc,ibc)                              5d12s23
            ibctop=ibcoff-1                                              12d14s20
            if1o=nff1(ncloip,isb,2)                                      12d14s20
            jgg=igg                                                     12d14s20
            jvecer=0                                                    2d12s21
            do if1=1,nff1(ncloip,isb,1)                                  12d14s20
             ndenh=0                                                    12d14s20
             do i=idenh,ibctop                                           12d14s20
              bc(i)=0d0                                                  12d14s20
             end do                                                      12d14s20
             i1c=iff1(if1o)                                                11d25s20
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
             i1o=iff1(if1o)                                                11d25s20
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
             jf1o=nff1(nclojp,jsb,2)                                      12d14s20
             jvec=ivec(nuse,jsb)                                        2d12s21
             do jf1=1,nff1(nclojp,jsb,1)                                  12d14s20
              j1c=iff1(jf1o)                                                11d25s20
              jf1o=jf1o+1                                                   11d25s20
              j1o=iff1(jf1o)                                                11d25s20
              j1o=ibset(j1o,norbxx)                                          11d25s20
              jf1o=jf1o+1                                                   11d25s20
              gandcc=ieor(j1c,i1c)                                      10d13s22
              gandco=ieor(j1o,i1o)                                      10d13s22
              gandcb=ior(gandcc,gandco)                                 10d20s22
              ndifb=popcnt(gandcb)                                      10d20s22
              if(ndifb.le.4)then                                        10d20s22
               ndifs=popcnt(gandco)                                      10d13s22
               ndifd=popcnt(gandcc)                                      10d13s22
c
c     vec is ordered nrootu,nherev. re-order to be nherev,nrootu.
              if(itransv(nuse,jsb).eq.0)then                            11d1s22
               call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))        2d12s21
               itmp=ibcoff                                                    1d27s21
               ibcoff=itmp+nrow                                            1d27s21
               call enough('hcssd12.3',bc,ibc)                           5d12s23
               jjvec=ivec(nuse,jsb)                                     2d12s21
               nrowm=nrow-1                                                1d27s21
               szj=0d0
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
              ivtrans=ibcoff                                            2d19s21
              ibcoff=ivtrans+nrow*ncsf(jarg)                            2d19s21
              call enough('hcssd12.4',bc,ibc)                            5d12s23
              do j=0,ncsf(jarg)-1                                       2d19s21
               do i=0,nrow-1                                            2d19s21
                ij=jvec+i+nrow*j                                        2d19s21
                ji=ivtrans+j+ncsf(jarg)*i                               2d19s21
                bc(ji)=bc(ij)                                           2d19s21
               end do                                                   2d19s21
              end do                                                    2d19s21
              if(ndifs.eq.2.and.ndifb.eq.2)then
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
     $              ibc)                                                11d14s22
               itmp=ibcoff                                               12d14s20
               ibcoff=itmp+nrowi                                         12d14s20
               call enough('hcssd12.5',bc,ibc)                           5d12s23
               call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid1,nrow,        2d19s21
     $              iwpb1,iwpk1,bc(ivtrans),ncsf(jarg),bc(itmp),        2d19s21
     $              ncsf(iarg),1d0,0d0,bc,ibc)                          11d10s22
               ndenh=1                                                  12d14s20
c
c     g r,v,i = sum v' denh r,v',i h0v v'
               do i=0,nrowim                                             12d14s20
                bc(idenh+i)=bc(idenh+i)+bc(itmp+i)                      2d22s21
               end do                                                    12d14s20
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
               do i=1,norb                                                   11d13s20
                ixn=min(itest(i,1),itest(i,2))
                if(ixn.gt.0)then                                             11d13s20
                 nok=nok+1                                                   11d13s20
                 itest(nok,3)=ixn                                            11d13s20
                 itest(nok,2)=i                                              11d13s20
                end if                                                       11d13s20
               end do                                                        11d13s20
               do i=1,nok                                                    11d13s20
                js=ism(itest(i,2))                                           11d13s20
                jg=irel(itest(i,2))-1                                        11d13s20
                icolj=((jg*(jg+1))/2)+jg                                12d14s20
                jdenj=idenj(js)+nrowi*icolj                             12d14s20
                ibc(ndenj(js)+icolj)=1                                  12d14s20
                if(itest(i,3).eq.2)then                                  12d14s20
                 icolk=jg+irefo(js)*jg                                  12d14s20
                 jdenk=idenk(js)+nrowi*icolk                            12d14s20
                 ibc(ndenk(js)+icolk)=1                                 12d14s20
                 do j=0,nrowim                                           12d14s20
                  bc(jdenj+j)=bc(jdenj+j)+2d0*bc(itmp+j)                2d22s21
                  bc(jdenk+j)=bc(jdenk+j)-bc(itmp+j)                    2d22s21
                 end do                                                  12d14s20
                else                                                     12d14s20
                 do j=0,nrowim                                           12d14s20
                  bc(jdenj+j)=bc(jdenj+j)+bc(itmp+j)                    2d22s21
                 end do                                                  12d14s20
                end if                                                   12d14s20
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
                   lsa=ism(nab2b(2))                                      12d14s20
                   lga=irel(nab2b(2))-1                                   12d14s20
                   lsb=ism(nab1b(1))                                      12d14s20
                   lgb=irel(nab1b(1))-1                                   12d14s20
                   icol=lga+irefo(lsa)*lgb                               12d14s20
                   kden=idenk(lsa)+nrowi*icol                            12d14s20
                   ibc(ndenk(lsa)+icol)=1                                12d14s20
                   do j=0,nrowim                                         2d19s21
                    bc(kden+j)=bc(kden+j)+bc(itmp2+j)                    2d22s21
                   end do                                                2d19s21
                   ibcoff=itmp1                                          2d19s21
                  end if                                                  12d14s20
                 end if                                                 4d14s21
                end if                                                   12d14s20
               end do                                                    12d14s20
              else                                                      10d21s22
               nnot=0                                                   10d21s22
               ipssx=0                                                  10d21s22
               if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s22
                nnot=4                                                          10d14s22
                ioxx(1)=1                                                     10d17s22
                ioxx(2)=1                                                     10d17s22
                do i=1,norbxx                                                     10d17s22
                 if(btest(gandcb,i))then                                         10d14s22
                  if((btest(i1c,i).and.btest(j1o,i)).or.                         10d17s22
     $                (btest(i1o,i).and..not.btest(j1c,i)))then                   10d14s22
                   nab4(2,ioxx(2))=i                                       10d17s22
                   ioxx(2)=ioxx(2)+1                                       10d17s22
                  else                                                           10d14s22
                   nab4(1,ioxx(1))=i                                       10d17s22
                   ioxx(1)=ioxx(1)+1                                       10d17s22
                  end if                                                         10d14s22
                 end if
                end do
               else if(ndifb.eq.3)then                                           10d14s22
                nnot=3
                ioxx(1)=1                                                        10d14s22
                ioxx(2)=1                                                        10d14s22
                iswap=0                                                          10d17s22
                do i=1,norbxx
                 if(btest(gandcb,i))then                                         10d14s22
                  if(btest(gandcc,i).and.                                        10d14s22
     $        ((btest(j1c,i).and..not.btest(i1o,i)).or.                 10d14s22
     $         (btest(i1c,i).and..not.btest(j1o,i))))then                     10d14s22
                   if(btest(i1c,i))iswap=1                                        10d17s22
                   nab4(1,1)=i                                                10d17s22
                   nab4(1,2)=i                                                10d17s22
                  else                                                           10d14s22
                   nab4(2,ioxx(2))=i                                             10d14s22
                   ioxx(2)=ioxx(2)+1
                  end if
                 end if                                                          10d14s22
                end do
                if(iswap.ne.0)then                                               10d17s22
                 icpy=nab4(1,1)                                               10d17s22
                 nab4(1,1)=nab4(2,1)                                       10d17s22
                 nab4(2,1)=icpy                                               10d17s22
                 icpy=nab4(1,2)                                               10d17s22
                 nab4(1,2)=nab4(2,2)                                       10d17s22
                 nab4(2,2)=icpy                                               10d17s22
                 nbt=0                                                           10d17s22
                 if(btest(j1c,nab4(1,2)).and..not.btest(j1c,nab4(1,1)))    10d17s22
     $       nbt=1                                                      10d17s22
                else                                                             10d17s22
                 nbt=0                                                           10d17s22
                 if(btest(i1c,nab4(2,2)).and..not.btest(i1c,nab4(2,1)))
     $              nbt=1                                                      10d17s22
                end if                                                           10d17s22
                if(nbt.ne.0)then                                                 10d17s22
                 nab4(1,1)=nab4(1,2)                                             10d17s22
                 nab4(2,1)=nab4(2,2)                                             10d17s22
                end if                                                           10d17s22
               else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                nnot=3
                do i=1,norbxx
                 if(btest(gandcb,i))then                                         10d14s22
                  if(btest(j1c,i))then
                   nab4(1,1)=i
                   nab4(1,2)=i
                  else                                                           10d14s22
                   nab4(2,1)=i
                   nab4(2,2)=i
                  end if
                 end if                                                          10d14s22
                end do
               end if                                                            10d14s22
               if(nnot.eq.3)then                                          12d8s20
                ipssx=1                                                   12d8s20
               else if(nnot.eq.4)then                                   11d1s22
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
                  if(max(nab1(1),nab1(2)).le.norb)then                    12d14s20
                   if(nab1(2).lt.nab1(1))then                             12d14s20
                    icpy=nab1(2)                                          12d14s20
                    nab1(2)=nab1(1)                                       12d14s20
                    nab1(1)=icpy                                          12d14s20
                   end if                                                 12d14s20
                   lsa=ism(nab1(1))                                       12d14s20
                   lga=irel(nab1(1))-1                                    12d14s20
                   lsb=ism(nab1(2))                                       12d14s20
                   lgb=irel(nab1(2))-1                                    12d14s20
                   if(lsa.eq.lsb)then                                     12d14s20
                    icol=((lgb*(lgb+1))/2)+lga                             12d14s20
                   else                                                   12d14s20
                    icol=lga+irefo(lsa)*lgb                               12d14s20
                   end if                                                 12d14s20
                   iden=idenj(lsa)+nrowi*icol                             12d14s20
                   ibc(ndenj(lsa)+icol)=1                                 12d14s20
                  else                                                    12d14s20
                   lsa=ism(nab1(1))                                       12d14s20
                   lga=irel(nab1(1))-1                                    12d14s20
                   lsb=ism(nab2(2))                                       12d14s20
                   lgb=irel(nab2(2))-1                                    12d14s20
                   icol=lgb+irefo(lsb)*lga                                12d14s20
                   iden=idenk(lsb)+nrowi*icol                             12d14s20
                   ibc(ndenk(lsb)+icol)=1                                 12d14s20
                  end if                                                  12d14s20
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
                 end if                                                  12d14s20
                end if                                                  4d14s21
               end do                                                   12d14s20
              end if                                                     12d14s20
              if(nnot.ne.0)ibcoff=ivtrans                               2d19s21
              end if                                                    10d13s22
              if(jsb.eq.isb)then                                         12d14s20
               j1o=ibclr(j1o,norbxx)                                     12d14s20
               j1o=ibset(j1o,norbx)                                      12d14s20
               gandcc=ieor(j1c,i1c)                                     10d13s22
               gandco=ieor(j1o,i1o)                                     10d13s22
               gandcb=ior(gandcc,gandco)                                10d21s22
               ndifb=popcnt(gandcb)                                     10d21s22
               if(ndifb.le.4)then                                       10d21s22
                ndifs=popcnt(gandco)                                     10d13s22
                ndifd=popcnt(gandcc)                                     10d13s22
c
c     vec is ordered nrootu,nherev. re-order to be nherev,nrootu.
               if(itransv(nuse,jsb).eq.0)then                           11d1s22
                call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))       2d12s21
                itmp=ibcoff                                                    1d27s21
                ibcoff=itmp+nrow                                            1d27s21
                call enough('hcssd12.7',bc,ibc)                          5d12s23
                jjvec=ivec(nuse,jsb)                                    2d12s21
                nrowm=nrow-1                                                1d27s21
               szj=0d0
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
               ivtrans=ibcoff                                           2d19s21
               ibcoff=ivtrans+nrow*ncsf(jarg)                           2d19s21
               call enough('hcssd12.8',bc,ibc)                           5d12s23
              sztz=0d0
               do j=0,ncsf(jarg)-1                                      2d19s21
                do i=0,nrow-1                                           2d19s21
                 ij=jvec+i+nrow*j                                       2d19s21
                 ji=ivtrans+j+ncsf(jarg)*i                              2d19s21
                 bc(ji)=bc(ij)                                          2d19s21
                end do                                                  2d19s21
               end do                                                   2d19s21
               if(ndifd.eq.0.and.ndifs.eq.0)then                        10d21s22
                jvv=ivv                                                 5d12s23
                nnad=nh0av(isbv)*nh0av(isbv)                            5d25s23
                jh0av=ih0av(isbv)                                       5d25s23
                ivadd=irefo(isbv)+i2s-1                                 5d25s23
                do ir1=0,nrootm                                         5d12s23
                 do ir2=0,ir1                                           5d12s23
                  bc(jvv)=0d0                                           5d12s23
                  do j=0,ncsf(jarg)-1                                   5d12s23
                   iad1=jvec+nherev*(ir2+nrootu*j)                      5d12s23
                   iad2=jgg+i2s-1+nvirt(isbv)*(ir1+nrootu*j)            5d18s23
                   do iv=0,nherevu-1                                    5d12s23
                    term=bc(iad1+iv)*bc(iad2+iv)                        5d12s23
                    bc(jvv)=bc(jvv)+term                                5d12s23
                    ivp=iv+ivadd                                        5d12s23
                    iad=jh0av+ivp*(nh0av(isbv)+1)                       5d25s23
                    bc(iad)=bc(iad)+term                                5d12s23
                   end do                                               5d12s23
                  end do                                                5d12s23
                  jvv=jvv+1                                             5d12s23
                  jh0av=jh0av+nnad                                      5d12s23
                 end do                                                 5d12s23
                end do                                                  5d12s23
                do i=1,norb                                             5d12s23
                 ig=irel(i)-1                                           5d12s23
                 is=ism(i)                                              5d12s23
                 if(btest(i1c,i))then                                   5d12s23
                  i2eu=inv(1,is,is,isbv)                                5d12s23
                  icol=((ig*(ig+1))/2)+ig
                  jj=jmats(i2eu)+nhere*icol                             5d12s23
                  nnadj=((irefo(is)*(irefo(is)+1))/2)*nhere             5d12s23
                  i2eu=invk1(1,is,is,isbv,1)                            5d12s23
                  icol=ig*(irefo(is)+1)                                 5d31s23
                  kk=kmats(i2eu)+nhere*icol                             5d12s23
                  nnad=irefo(is)*irefo(is)*nhere                        5d12s23
                  do ir1=0,nrootm                                       5d12s23
                   do ir2=0,ir1                                         5d12s23
                    do iv=0,nherev-1                                    6d5s23
                     ivp=iv+i2s                                         5d12s23
                     icol=ivp+nvirt(isbv)*(ivp-1)                       5d12s23
                     if(icol.ge.il.and.icol.le.ih)then                  5d12s23
                      kkk=kk+icol-il                                    5d12s23
                      jjj=jj+icol-il                                    5d12s23
                      sum=0d0                                           5d12s23
                      do j=0,ncsf(jarg)-1                               5d12s23
                       iad1=jvec+iv+nherev*(ir2+nrootu*j)               5d12s23
                       iad2=jgg+ivp-1+nvirt(isbv)*(ir1+nrootu*j)        5d12s23
                       sum=sum+bc(iad1)*bc(iad2)                        5d12s23
                      end do                                             5d12s23
                      bc(jjj)=bc(jjj)+2d0*sum                           5d12s23
                      bc(kkk)=bc(kkk)-sum                               5d12s23
                     end if                                             5d12s23
                    end do                                              5d12s23
                    kk=kk+nnad                                          5d12s23
                    jj=jj+nnadj                                         5d12s23
                   end do                                               5d12s23
                  end do                                                5d12s23
                 else if(btest(i1o,i))then                              5d12s23
                  i2eu=inv(1,is,is,isbv)                                5d12s23
                  icol=((ig*(ig+1))/2)+ig
                  jj=jmats(i2eu)+nhere*icol                             5d12s23
                  nnadj=((irefo(is)*(irefo(is)+1))/2)*nhere             5d12s23
                  do ir1=0,nrootm                                       5d12s23
                   do ir2=0,ir1                                         5d12s23
                    do iv=0,nherev-1                                    6d5s23
                     ivp=iv+i2s                                         5d12s23
                     icol=ivp+nvirt(isbv)*(ivp-1)                       5d12s23
                     if(icol.ge.il.and.icol.le.ih)then                  5d12s23
                      jjj=jj+icol-il                                    5d12s23
                      sum=0d0                                           5d12s23
                      do j=0,ncsf(jarg)-1                               5d12s23
                       iad1=jvec+iv+nherev*(ir2+nrootu*j)               5d12s23
                       iad2=jgg+ivp-1+nvirt(isbv)*(ir1+nrootu*j)        5d12s23
                       sum=sum+bc(iad1)*bc(iad2)                        5d12s23
                      end do                                             5d12s23
                      bc(jjj)=bc(jjj)+sum                               5d25s23
                     end if                                             5d12s23
                    end do                                              5d12s23
                    jj=jj+nnadj                                         5d12s23
                   end do                                               5d12s23
                  end do                                                5d12s23
                 end if                                                 5d12s23
                end do                                                  5d12s23
c     ic jc 2J-K
c     ic jo 2J-k
c     io jc
c     io jo i ne j J/2
                do i=1,norb                                             5d12s23
                 ff=0d0                                                 5d12s23
                 ig=irel(i)-1                                           5d12s23
                 is=ism(i)                                              5d12s23
                 if(btest(i1c,i))then                                   5d12s23
                  ff=2d0                                                5d12s23
                  do j=1,norb                                           5d12s23
                   if(btest(i1c,j).or.btest(i1o,j))then                 5d12s23
                    jg=irel(j)-1                                        5d12s23
                    js=ism(j)                                           5d12s23
                    iadint=igetint(ioooo,is,is,js,js,ig,ig,jg,jg,       5d12s23
     $                  nnad,.false.,bc,ibc)                            5d12s23
                    do itri=0,ntri-1                                    5d12s23
                     bc(iadint)=bc(iadint)+2d0*bc(ivv+itri)             5d12s23
                     iadint=iadint+nnad                                 5d12s23
                    end do                                              5d12s23
                    iadint=igetint(ioooo,is,js,js,is,ig,jg,jg,ig,       5d12s23
     $                  nnad,.false.,bc,ibc)                            5d12s23
                    do itri=0,ntri-1                                    5d12s23
                     bc(iadint)=bc(iadint)-bc(ivv+itri)                 5d12s23
                     iadint=iadint+nnad                                 5d12s23
                    end do                                              5d12s23
                   end if
                  end do                                                5d12s23
                 else if(btest(i1o,i))then                              5d12s23
                  ff=1d0                                                5d12s23
                  do j=1,norb                                           5d12s23
                   if(btest(i1o,j).and.i.ne.j)then                      5d12s23
                    jg=irel(j)-1                                        5d12s23
                    js=ism(j)                                           5d12s23
                    iadint=igetint(ioooo,is,is,js,js,ig,ig,jg,jg,       5d12s23
     $                  nnad,.false.,bc,ibc)                            5d12s23
                    do itri=0,ntri-1                                    5d12s23
                     bc(iadint)=bc(iadint)+0.5d0*bc(ivv+itri)           5d12s23
                     iadint=iadint+nnad                                 5d12s23
                    end do                                              5d12s23
                   end if                                               5d12s23
                  end do                                                5d12s23
                 end if                                                 5d12s23
                 if(ff.gt.1d-4)then                                     5d12s23
                  iad=ih0av(is)+ig*(1+nh0av(is))                        5d12s23
                  nnad=nh0av(is)*nh0av(is)                              5d12s23
                  do itri=0,ntri-1                                      5d12s23
                   bc(iad)=bc(iad)+ff*bc(ivv+itri)                      5d12s23
                   iad=iad+nnad                                         5d12s23
                  end do                                                5d12s23
                 end if                                                 5d12s23
                end do                                                  5d12s23
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
                  if(i2.ne.nopeni)then                                  12d17s20
                   lsa=ism(nab1(1))                                              11d13s20
                   lga=irel(nab1(1))-1                                           11d13s20
                   lsd=ism(nab2(2))                                              11d13s20
                   lgd=irel(nab2(2))-1                                           11d13s20
                   iadint=igetint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd, 5d12s23
     $                  nnad,.false.,bc,ibc)                            5d12s23
                  else                                                  12d17s20
                   icol=lgb+irefo(lsb)*lgc                              12d17s20
                   iden=idenkv(lsb)+nrowi*icol                          12d17s20
                   ibc(ndenkv(lsb)+icol)=1                              12d17s20
                  end if                                                12d17s20
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
                   if(i2.ne.nopeni)then                                 12d17s20
                    jadint=iadint                                       5d12s23
                    do ir1=0,nrootm                                     5d12s23
                     do ir2=0,ir1                                       5d12s23
                      sumvv=0d0                                         5d12s23
                      do j=0,ncsf(jarg)-1                               5d12s23
                       iad1=ltmp2+j+ncsf(jarg)*nherev*ir2               5d22s23
                       iad2=jgg+i2s-1+nvirt(isbv)*(ir1+nrootu*j)        5d12s23
                       do iv=0,nherevu-1                                5d12s23
                        sumvv=sumvv+bc(iad2+iv)*bc(iad1)                5d22s23
                        iad1=iad1+ncsf(jarg)                            5d22s23
                       end do                                           5d12s23
                      end do                                            5d12s23
                      bc(jadint)=bc(jadint)+sumvv                       5d12s23
                      jadint=jadint+nnad                                5d12s23
                     end do                                             5d12s23
                    end do                                              5d12s23
                   else                                                 12d17s20
                    do j=0,nrowim                                       2d19s21
                     bc(iden+j)=bc(iden+j)+bc(ltmp2+j)                  2d22s21
                    end do                                              2d19s21
                   end if                                               12d17s20
                   ibcoff=ltmp1                                         2d19s21
                  else                                                   12d14s20
                   if(i2.ne.nopeni)then                                 12d17s20
                    jadint=iadint                                       5d12s23
                    do ir1=0,nrootm                                     5d12s23
                     do ir2=0,ir1                                       5d12s23
                      sumvv=0d0                                         5d12s23
                      do j=0,ncsf(jarg)-1                               5d12s23
                       iad1=jgg+i2s-1+nvirt(isbv)*(ir2+nrootu*j)        5d22s23
                       iad2=jgg+i2s-1+nvirt(isbv)*(ir1+nrootu*j)        5d12s23
                       do iv=0,nherevu-1                                5d12s23
                        sumvv=sumvv+bc(iad2+iv)*bc(iad1+iv)             5d22s23
                       end do                                           5d12s23
                      end do                                            5d12s23
                      bc(jadint)=bc(jadint)+sumvv                       5d12s23
                      jadint=jadint+nnad                                5d12s23
                     end do                                             5d12s23
                    end do                                              5d12s23
                   else                                                 12d17s20
                    do i=0,nrowim                                         12d14s20
                     bc(iden+i)=bc(iden+i)-bc(ivtrans+i)                2d22s21
                    end do                                                12d14s20
                   end if                                               12d17s20
                  end if                                                 12d14s20
                 end do                                                  12d14s20
                end do                                                  12d14s20
               else if(ndifs.eq.2.and.ndifb.eq.2)then                   10d21s22
                do i=1,norbx
                 if(btest(gandco,i))then                                10d21s22
                  if((btest(j1o,i).and..not.btest(i1c,i)).or.
     $                (btest(j1c,i).and.btest(i1o,i)))then                           10d14s22
                   nab4(1,1)=i                                          10d21s22
                  else                                                           10d14s22
                   nab4(2,1)=i                                          10d21s22
                  end if
                 end if                                                          10d14s22
                end do                                                           10d14s22
                call gandc(i1c,i1o,j1c,j1o,nopeni,nopenj,iarg,jarg,ncsf,2d19s21
     $            norbx,ixw1,ixw2,nnotb,nab1b,iwpb1b,iwpk1b,ncsfmid1b,  11d14s22
     $               bc,ibc)                                            11d14s22
                nab1(1)=nab1b(2)                                        2d22s21
                nab1(2)=nab1b(1)                                        2d22s21
                if(nab1(1).gt.nab1(2))then                              12d17s20
                 icpy=nab1(1)                                           12d17s20
                 nab1(1)=nab1(2)                                        12d17s20
                 nab1(2)=icpy                                           12d17s20
                end if                                                  12d17s20
                lsa=ism(nab1(1))                                        12d15s20
                lga=irel(nab1(1))-1                                     12d15s20
                lgb=irel(nab1(2))-1                                     12d15s20
                nnad=nh0av(isb)*nh0av(isb)                              5d12s23
                jh0av=ih0av(isb)                                        5d12s23
                ivadd=irefo(isb)+i2s-1                                  5d12s23
                ih0=ih0av(lsa)+lga+nh0av(lsa)*lgb                       12d15s20
                ltmp1=ibcoff                                            2d19s21
                ltmp2=ltmp1+ncsf(iarg)*nrow                             2d19s21
                ibcoff=ltmp2+ncsf(iarg)*nrow                             2d19s21
                call enough('hcss.  9',bc,ibc)                          5d12s23
                call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid1b,nrow,      2d19s21
     $              iwpb1b,iwpk1b,bc(ivtrans),ncsf(jarg),bc(ltmp2),     2d19s21
     $              ncsf(iarg),1d0,0d0,bc,ibc)                          11d10s22
                jvv=ivv                                                 5d12s23
                do ir1=0,nrootm                                         5d12s23
                 do ir2=0,ir1                                           5d12s23
                  bc(jvv)=0d0                                           5d12s23
                  do j=0,ncsf(iarg)-1                                   5d12s23
                   iad1=ltmp2+j+ncsf(iarg)*nherev*ir2                   5d18s23
                   iad2=jgg+i2s-1+nvirt(isbv)*(ir1+nrootu*j)            5d12s23
                   do iv=0,nherevu-1                                    5d12s23
                    term=bc(iad1)*bc(iad2+iv)                           5d18s23
                    iad1=iad1+ncsf(iarg)                                5d18s23
                    bc(jvv)=bc(jvv)+term                                5d12s23
                   end do                                               5d12s23
                  end do                                                5d12s23
                  bc(ih0)=bc(ih0)+bc(jvv)                               5d12s23
                  ih0=ih0+nnad                                          5d12s23
                  jvv=jvv+1                                             5d12s23
                 end do                                                 5d12s23
                end do                                                  5d12s23
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
                do i=1,norb                                                   11d13s20
                 ixn=min(itest(i,1),itest(i,2))
                 if(ixn.gt.0)then                                             11d13s20
                  nok=nok+1                                                   11d13s20
                  itest(nok,3)=ixn                                            11d13s20
                  itest(nok,2)=i                                              11d13s20
                 end if                                                       11d13s20
                end do                                                        11d13s20
                do i=1,nok                                                    11d13s20
                 js=ism(itest(i,2))                                           11d13s20
                 jg=irel(itest(i,2))-1                                        11d13s20
                 iadj=igetint(ioooo,js,js,lsa,lsa,jg,jg,lga,lgb,nnadj,  5d12s23
     $                .false.,bc,ibc)                                   5d12s23
                 if(itest(i,3).eq.2)then                                5d12s23
                  iadk=igetint(ioooo,js,lsa,js,lsa,jg,lga,jg,lgb,nnadk, 5d12s23
     $                .false.,bc,ibc)                                   5d12s23
                  do j=0,ntri-1                                         5d12s23
                   bc(iadj)=bc(iadj)+2d0*bc(ivv+j)                      5d12s23
                   iadj=iadj+nnadj                                      5d12s23
                   bc(iadk)=bc(iadk)-bc(ivv+j)                          5d12s23
                   iadk=iadk+nnadk                                      5d12s23
                  end do                                                5d12s23
                 else                                                     12d14s20
                  do j=0,ntri-1                                         5d12s23
                   bc(iadj)=bc(iadj)+bc(ivv+j)                          5d12s23
                   iadj=iadj+nnadj                                      5d12s23
                  end do                                                5d12s23
                 end if                                                   12d14s20
                end do                                                    12d14s20
                icol=((lgb*(lgb+1))/2)+lga                              12d17s20
                iden=idenjv(lsa)+nrowi*icol                             12d17s20
                ibc(ndenjv(lsa)+icol)=1                                 12d17s20
                do ii=0,nrow*ncsf(iarg)-1                               2d22s21
                 bc(iden+ii)=bc(iden+ii)+bc(ltmp2+ii)                   2d22s21
                end do                                                  2d22s21
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
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                   call gandc(i1c,i1o,itestc,itesto,nopeni,nopenk,      2d22s21
     $               iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,2d22s21
     $                iwpk1b,ncsfmid1b,bc,ibc)                          11d14s22
                   call gandc(itestc,itesto,j1c,j1o,nopenk,nopenj,      2d22s21
     $            karg,jarg,ncsf,norbx,ixw1,ixw2,nnot2b,nab2b,iwpb2b,   2d22s21
     $               iwpk2b,ncsfmid2b,bc,ibc)                           11d14s22
                   if(nnot1b.eq.2.and.nnot2b.eq.2)then                   2d22s21
                    lsa=ism(nab2b(2))                                      12d14s20
                    lga=irel(nab2b(2))-1                                   12d14s20
                    lsb=ism(nab2b(1))                                      12d14s20
                    lgb=irel(nab2b(1))-1                                   12d14s20
                    lsc=ism(nab1b(2))                                      12d14s20
                    lgc=irel(nab1b(2))-1                                   12d14s20
                    lsd=ism(nab1b(1))                                      12d14s20
                    lgd=irel(nab1b(1))-1                                   12d14s20
                    iadint=igetint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,   5d12s23
     $                   lgd,nnad,.false.,bc,ibc)                       5d12s23
                    ltmp1=ibcoff                                         2d22s21
                    ltmp2=ltmp1+nrow*ncsf(iarg)                          2d22s21
                    ibcoff=ltmp2+nrow*ncsf(iarg)                         2d22s21
                    call enough('hcssd12.10',bc,ibc)                      5d12s23
                    call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),       2d22s21
     $                 ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b,
     $                 bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc)11d10s22
                    do ir1=0,nrootm                                     5d12s23
                     do ir2=0,ir1                                       5d12s23
                      sum=0d0                                           5d12s23
                      do j=0,ncsf(iarg)-1                               5d12s23
                       iad1=ltmp2+j+ncsf(iarg)*nherev*ir2               5d25s23
                       iad2=jgg+i2s-1+nvirt(isbv)*(ir1+nrootu*j)        5d12s23
                       do iv=0,nherevu-1                                5d12s23
                        sum=sum+bc(iad1)*bc(iad2+iv)                    5d25s23
                        iad1=iad1+ncsf(iarg)                            5d25s23
                       end do                                           5d12s23
                      end do                                            5d12s23
                      bc(iadint)=bc(iadint)+sum                         5d12s23
                      iadint=iadint+nnad                                5d12s23
                     end do                                             5d12s23
                    end do                                              5d12s23
                    ibcoff=ltmp1                                         2d22s21
                   end if                                                  12d14s20
                  end if                                                4d14s21
                 end if                                                   12d14s20
                end do                                                    12d14s20
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
                  call enough('hcssd12.11',bc,ibc)                        5d12s23
                  call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),         2d22s21
     $                 ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b, 2d22s21
     $                bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc) 11d10s22
                  lsb=ism(nab2b(1))                                       12d14s20
                  lgb=irel(nab2b(1))-1                                    12d14s20
                  lsc=ism(nab1b(2))                                       12d14s20
                  lgc=irel(nab1b(2))-1                                   2d22s21
                  icol=lgc+irefo(lsc)*lgb                                12d17s20
                  iden=idenkv(lsc)+nrowi*icol                            12d17s20
                  ibc(ndenkv(lsc)+icol)=1                                12d17s20
                  sic=0d0
                  do ii=0,nrow*ncsf(iarg)-1                              2d22s21
                   bc(iden+ii)=bc(iden+ii)+bc(ltmp2+ii)                  2d22s21
                  end do                                                 2d22s21
                  ibcoff=ltmp1                                           2d22s21
                 end if                                                  12d17s20
                end if                                                  4d14s21
               else                                                     11d1s22
                nnot=0
                if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s22
                 nnot=4                                                          10d14s22
                 ioxx(1)=1                                                     10d17s22
                 ioxx(2)=1                                                     10d17s22
                 do i=1,norbx                                                     10d17s22
                  if(btest(gandcb,i))then                                         10d14s22
                   if((btest(i1c,i).and.btest(j1o,i)).or.                         10d17s22
     $                 (btest(i1o,i).and..not.btest(j1c,i)))then                   10d14s22
                    nab4(2,ioxx(2))=i                                       10d17s22
                    ioxx(2)=ioxx(2)+1                                       10d17s22
                   else                                                           10d14s22
                    nab4(1,ioxx(1))=i                                       10d17s22
                    ioxx(1)=ioxx(1)+1                                       10d17s22
                   end if                                                         10d14s22
                  end if
                 end do
                else if(ndifb.eq.3)then                                           10d14s22
                 nnot=3
                 ioxx(1)=1                                                        10d14s22
                 ioxx(2)=1                                                        10d14s22
                 iswap=0                                                          10d17s22
                 do i=1,norbx
                  if(btest(gandcb,i))then                                         10d14s22
                   if(btest(gandcc,i).and.                                        10d14s22
     $                 ((btest(j1c,i).and..not.btest(i1o,i)).or.                 10d14s22
     $                 (btest(i1c,i).and..not.btest(j1o,i))))then                     10d14s22
                   if(btest(i1c,i))iswap=1                                        10d17s22
                    nab4(1,1)=i                                                10d17s22
                    nab4(1,2)=i                                                10d17s22
                   else                                                           10d14s22
                    nab4(2,ioxx(2))=i                                             10d14s22
                    ioxx(2)=ioxx(2)+1
                   end if
                  end if                                                          10d14s22
                 end do
                 if(iswap.ne.0)then                                               10d17s22
                  icpy=nab4(1,1)                                               10d17s22
                  nab4(1,1)=nab4(2,1)                                       10d17s22
                  nab4(2,1)=icpy                                               10d17s22
                  icpy=nab4(1,2)                                               10d17s22
                  nab4(1,2)=nab4(2,2)                                       10d17s22
                  nab4(2,2)=icpy                                               10d17s22
                  nbt=0                                                           10d17s22
                  if(btest(j1c,nab4(1,2)).and..not.btest(j1c,nab4(1,1)))    10d17s22
     $                 nbt=1                                                      10d17s22
                 else                                                             10d17s22
                  nbt=0                                                           10d17s22
                  if(btest(i1c,nab4(2,2)).and..not.btest(i1c,nab4(2,1)))
     $                 nbt=1                                                      10d17s22
                 end if                                                           10d17s22
                 if(nbt.ne.0)then                                                 10d17s22
                  nab4(1,1)=nab4(1,2)                                             10d17s22
                  nab4(2,1)=nab4(2,2)                                             10d17s22
                 end if                                                           10d17s22
                else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                 nnot=3
                 do i=1,norbx
                  if(btest(gandcb,i))then                                         10d14s22
                   if(btest(j1c,i))then
                    nab4(1,1)=i
                    nab4(1,2)=i
                   else                                                           10d14s22
                    nab4(2,1)=i
                    nab4(2,2)=i
                   end if
                  end if                                                          10d14s22
                 end do
                end if                                                            10d14s22
                ipssx=0
                if(nnot.eq.3)then                                          12d8s20
                 ipssx=1                                                   12d8s20
                else if(nnot.eq.4)then                                  11d1s22
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
                   iadint=igetint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,5d12s23
     $                  nnad,.false.,bc,ibc)                            5d12s23
                   ltmp1=ibcoff                                          2d22s21
                   ltmp2=ltmp1+nrow*ncsf(iarg)                           2d22s21
                   ibcoff=ltmp2+nrow*ncsf(iarg)                           2d22s21
                   call enough('hcssd12.12',bc,ibc)                       5d12s23
                   call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),        2d22s21
     $                 ncsfmid1b,iwpb1b,iwpk1b,ncsfmid2b,iwpb2b,iwpk2b, 2d22s21
     $                 bc(ivtrans),ncsf(jarg),nrow,bc(ltmp2),0d0,bc,ibc)11d10s22
                   do ir1=0,nrootm                                      5d12s23
                    do ir2=0,ir1                                        5d12s23
                     sum=0d0                                            5d12s23
                     do j=0,ncsf(iarg)-1                                5d12s23
                      iad1=ltmp2+j+ncsf(iarg)*nherev*ir2                5d25s23
                      iad2=jgg+i2s-1+nvirt(isbv)*(ir1+nrootu*j)         5d12s23
                      do iv=0,nherevu-1                                 5d12s23
                       sum=sum+bc(iad1)*bc(iad2+iv)                     5d25s23
                       iad1=iad1+ncsf(iarg)
                      end do                                            5d12s23
                     end do                                             5d12s23
                     bc(iadint)=bc(iadint)+sum*fact                     5d12s23
                     iadint=iadint+nnad                                 5d12s23
                    end do                                              5d12s23
                   end do                                               5d12s23
                   ibcoff=ltmp1                                          2d22s21
                   if(ipss.eq.2)go to 2001                               12d16s20
                  end if                                                  12d14s20
                 end if                                                 4d14s21
                end do                                                   12d14s20
 2001           continue                                                12d16s20
               end if                                                    12d14s20
               end if                                                   10d13s22
              end if                                                     12d14s20
              jvec=jvec+ncsf(jarg)*nrow                                  12d14s20
             end do                                                      12d14s20
             ltmp1=ibcoff                                               2d22s21
             ibcoff=ltmp1+nrowi                                         2d22s21
             call enough('hcssd12.14',bc,ibc)                             5d12s23
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
             igotj=0
             igotjv=0
             do isc=1,nsymb                                             2d22s21
              isd=multh(isc,ijsbv)                                      2d22s21
              nn=irefo(isc)*irefo(isd)                                  2d22s21o
              do i=0,nn-1                                               2d22s21
               if(ibc(ndenk(isc)+i).ne.0)then                           2d22s21
                iden=idenk(isc)+nrowi*i                                 2d22s21
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
               if(ibc(ndenkv(isc)+i).ne.0)then                           2d22s21
                iden=idenkv(isc)+nrowi*i                                 2d22s21
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
              end do                                                    2d22s21
              if(isd.ge.isc)then                                        2d22s21
               if(isd.eq.isc)nn=(irefo(isc)*(irefo(isc)+1))/2           2d22s21
               do i=0,nn-1                                               2d22s21
                if(ibc(ndenj(isc)+i).ne.0)then                           2d22s21
                 iden=idenj(isc)+nrowi*i                                 2d22s21
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
                if(ibc(ndenjv(isc)+i).ne.0)then                           2d22s21
                 iden=idenjv(isc)+nrowi*i                                 2d22s21
                 szjv=0d0
                 do j=0,nrow-1                                              2d22s21
                  do ii=0,ncsf(iarg)-1                                      2d22s21
                   ij=iden+ii+ncsf(iarg)*j                                 2d22s21
                   ji=ltmp1+j+nrow*ii                                       2d22s21
                   bc(ji)=bc(ij)                                            2d22s21
                   szjv=szjv+bc(ij)**2
                  end do                                                    2d22s21
                 end do                                                     2d22s21
                 if(szjv.gt.1d-4)igotjv=isc
                 do j=0,nrowim                                              2d22s21
                  bc(iden+j)=bc(ltmp1+j)                                 2d22s21
                 end do                                                     2d22s21
                end if                                                   2d22s21
               end do                                                    2d22s21
              end if                                                    2d22s21
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
             call enough('hcssd12. 15',bc,ibc)                             5d12s23
             do i=itmpp,ibcoff-1                                        12d17s20
              bc(i)=0d0                                                 12d17s20
             end do                                                     12d17s20
             ihit=0
             if(ndenh.ne.0)then                                         12d15s20
c     g r,v,i = sum v' denh r,v',i h0v v'
              itmp=ibcoff                                               12d17s20
              irr=0                                                     5d12s23
              do ir1=0,nrootm                                           5d12s23
               do ir2=0,ir1                                             5d12s23
                do i=0,ncsf(iarg)-1                                     5d12s23
                 i10=i1s                                                5d12s23
                 i1n=nvirt(isbv)
                 do i2=i2s,i2e                                          5d12s23
                  if(i2.eq.i2e)i1n=i1e                                  5d12s23
                  i2m=i2-1                                              5d25s23
                  i2p=i2+1                                              5d25s23
                  iadd=idenh+i2-i2s+nherev*(ir2+nrootu*i)               5d12s23
                  iadg=jgg-1+nvirt(isbv)*(ir1+nrootu*i)                 5d12s23
                  iadh=ih0av(isbv)+irefo(isbv)-1+nh0av(isbv)            5d12s23
     $                 *(i2+irefo(isbv)-1+nh0av(isbv)*irr)              5d12s23
                  do i1=i10,min(i1n,i2m)                                5d25s23
                   bc(iadh+i1)=bc(iadh+i1)+bc(iadd)*bc(iadg+i1)         5d18s23
                  end do                                                5d12s23
                  do i1=max(i2p,i10),i1n                                5d25s23
                   bc(iadh+i1)=bc(iadh+i1)+bc(iadd)*bc(iadg+i1)         5d18s23
                  end do                                                5d12s23
                  i10=1                                                 5d12s23
                 end do                                                 5d12s23
                end do                                                  5d12s23
                irr=irr+1                                               5d12s23
               end do                                                   5d12s23
              end do                                                    5d12s23
             end if                                                     12d15s20
             do isc=1,nsymb                                             12d15s20
              isd=multh(isc,ijsbv)                                      12d15s20
              nn=irefo(isc)*irefo(isd)                                  12d15s20
              i2eu=invk1(1,isc,isd,isbv,1)                              12d15s20
              ikeep=ibcoff                                              12d17s20
              ibcoff=ikeep+nn                                           12d17s20
              call enough('hcssd12.18',bc,ibc)                            5d12s23
              nok=0                                                     12d17s20
              nokv=0                                                    12d17s20
              szk=0d0
              nszkk=0
              do i=0,nn-1                                               12d17s20
               ibc(ikeep+i)=0                                           12d17s20
               if(ibc(ndenkv(isc)+i).ne.0)then                          12d17s20
                iden=idenkv(isc)+nrowi*i
                rms=1d0                                                 12d17s20
                if(ibc(ndenk(isc)+i).ne.0)then                          12d17s20
                 jden=idenk(isc)+nrowi*i                                12d17s20
                 rms=0d0                                                12d17s20
                 szkk=0d0
                 do k=0,nrowim                                          12d17s20
                  rms=rms+(bc(jden+k)-bc(iden+k))**2                    12d17s20
                  szkk=szkk+bc(iden+k)**2
                 end do                                                 12d17s20
                 rms=sqrt(rms/dfloat(nrowi))                            12d17s20
                end if                                                  12d17s20
                if(rms.gt.1d-10)then
                 jden=idenkv(isc)+nrowi*nokv                            12d17s20
                 do k=0,nrowim                                          12d17s20
                  bc(jden+k)=bc(iden+k)                                 12d17s20
                 end do                                                 12d17s20
                 ibc(ndenkv(isc)+nokv)=i                                12d17s20
                 nokv=nokv+1                                            12d17s20
                else                                                    12d17s20
                 ibc(ikeep+i)=1                                         12d17s20
                end if                                                  12d17s20
               end if                                                   12d17s20
               if(ibc(ndenk(isc)+i).ne.0)then                           12d17s20
                iden=idenk(isc)+nrowi*i                                 12d17s20
                jden=idenk(isc)+nrowi*nok                               12d17s20
                do k=0,nrowim                                           12d17s20
                 bc(jden+k)=bc(iden+k)                                  12d17s20
                end do                                                  12d17s20
                ibc(ndenk(isc)+nok)=i                                   12d17s20
                nok=nok+1                                               12d17s20
               end if                                                   12d17s20
              end do                                                    12d17s20
              if(nok.gt.0)then                                          12d17s20
               irr=0                                                    5d12s23
               do ir1=0,nrootm                                          5d12s23
                do ir2=0,ir1                                            5d12s23
                 do ii=0,nok-1                                          5d12s23
                  icol=ibc(ndenk(isc)+ii)                                  12d17s20
                  kint=kmats(i2eu)+nhere*(icol+nn*irr)                  5d12s23
                  do i=0,ncsf(iarg)-1                                    5d12s23
                   i10=i1s                                               5d12s23
                   i1n=nvirt(isbv)                                      5d22s23
                   kkint=kint                                           5d12s23
                   iadg=jgg-1+nvirt(isbv)*(ir1+nrootu*i)                 5d12s23
                   if(ibc(ikeep+icol).eq.1.or.isbv.ne.jsbv)then         5d24s23
                    do i2=i2s,i2e                                         5d12s23
                     if(i2.eq.i2e)i1n=i1e                                 5d12s23
                     iadk=idenk(isc)+i2-i2s                              5d12s23
     $                   +nherev*(ir2+nrootu*(i+ncsf(iarg)*ii))         5d31s23
                     do i1=i10,i1n                                        5d12s23
                      bc(kkint)=bc(kkint)+bc(iadg+i1)*bc(iadk)           5d12s23
                      kkint=kkint+1                                      5d12s23
                     end do                                               5d12s23
                     i10=1                                                5d12s23
                    end do                                               5d12s23
                   else                                                 5d24s23
                    do i2=i2s,i2e                                       5d24s23
                     i2m=i2-1                                           5d24s23
                     i2p=i2+1                                           5d24s23
                     if(i2.eq.i2e)i1n=i1e                               5d24s23
                     iadk=idenk(isc)+i2-i2s                              5d12s23
     $                   +nherev*(ir2+nrootu*(i+ncsf(iarg)*ii))         5d31s23
                     kintu=kkint-i10                                        12d18s20
                     do i1=i10,min(i1n,i2m)                                12d18s20
                      bc(kintu+i1)=bc(kintu+i1)+bc(iadg+i1)*bc(iadk)    5d24s23
                     end do                                                12d18s20
                     do i1=max(i10,i2p),i1n                                12d18s20
                      bc(kintu+i1)=bc(kintu+i1)+bc(iadg+i1)*bc(iadk)    5d24s23
                     end do
                     kkint=kkint+i1n+1-i10                                   12d18s20
                     i10=1                                                 12d17s20
                    end do                                              5d24s23
                   end if                                               5d24s23
                  end do                                                5d12s23
                 end do                                                 5d12s23
                 irr=irr+1                                              5d12s23
                end do                                                  5d12s23
               end do                                                   5d12s23
              end if                                                    12d17s20
              if(nokv.gt.0)then                                         12d17s20
c
               do i=0,nokv-1                                            12d17s20
                icol=ibc(ndenkv(isc)+i)                                 12d17s20
                do m=0,ncsf(iarg)-1                                     12d17s20
                 irr=0                                                   5d12s23
                 do ir1=0,nrootm
                  do ir2=0,ir1                                           5d12s23
                   kint=kmats(i2eu)+nhere*(icol+nn*irr)                  5d12s23
                   jjgg=jgg-1+nvirt(isbv)*ir1+ngg*m                       1d27s21
                   ktmpd=idenkv(isc)-i2s+nherev*(ir2                    5d12s23
     $                  +nrootu*(m+ncsf(iarg)*i))                       5d12s23
                   do i2=i2s,i2e                                        6d5s23
                    if(i2.eq.i2e)i1n=i1e                                   12d17s20
                    irow=i2+nvirt(isbv)*(i2-1)                          6d5s23
                    if(irow.ge.il.and.irow.le.ih)then                   6d5s23
                     irow=irow-il                                       6d5s23
                     bc(kint+irow)=bc(kint+irow)                        5d12s23
     $                    +bc(jjgg+i2)*bc(ktmpd+i2)                     6d5s23
                    end if
                   end do                                                12d17s20
                   irr=irr+1
                  end do                                                 12d17s20
                 end do
                end do                                                  12d17s20
               end do                                                   12d17s20
c
              end if                                                    12d17s20
              if(isd.ge.isc)then                                        12d15s20
               if(isd.eq.isc)nn=(irefo(isc)*(irefo(isc)+1))/2           12d15s20
               i2eu=inv(1,isc,isd,isbv)                                 12d15s20
               icase=inv(2,isc,isd,isbv)                                 12d15s20
               ikeep=ibcoff                                              12d17s20
               ibcoff=ikeep+nn                                           12d17s20
               call enough('hcssd12.22',bc,ibc)                           5d12s23
               nok=0                                                     12d17s20
               nokv=0                                                    12d17s20
               do i=0,nn-1                                               12d17s20
                ibc(ikeep+i)=0                                           12d17s20
                if(ibc(ndenjv(isc)+i).ne.0)then                          12d17s20
                 iden=idenjv(isc)+nrowi*i
                 rms=1d0                                                 12d17s20
                 if(ibc(ndenj(isc)+i).ne.0)then                          12d17s20
                  jden=idenj(isc)+nrowi*i                                12d17s20
                  rms=0d0                                                12d17s20
                  do k=0,nrowim                                          12d17s20
                   rms=rms+(bc(jden+k)-bc(iden+k))**2                    12d17s20
                  end do                                                 12d17s20
                  rms=sqrt(rms/dfloat(nrowi))                            12d17s20
                 end if                                                  12d17s20
                 if(rms.gt.1d-10)then
                  jden=idenjv(isc)+nrowi*nokv                            12d17s20
                  do k=0,nrowim                                          12d17s20
                   bc(jden+k)=bc(iden+k)                                 12d17s20
                  end do                                                 12d17s20
                  ibc(ndenjv(isc)+nokv)=i                                12d17s20
                  nokv=nokv+1                                            12d17s20
                 else                                                    12d17s20
                  ibc(ikeep+i)=1                                         12d17s20
                 end if                                                  12d17s20
                end if                                                   12d17s20
                if(ibc(ndenj(isc)+i).ne.0)then                           12d17s20
                 iden=idenj(isc)+nrowi*i                                 12d17s20
                 jden=idenj(isc)+nrowi*nok                               12d17s20
                 do k=0,nrowim                                           12d17s20
                  bc(jden+k)=bc(iden+k)                                  12d17s20
                 end do                                                  12d17s20
                 ibc(ndenj(isc)+nok)=i                                   12d17s20
                 nok=nok+1                                               12d17s20
                end if                                                   12d17s20
               end do                                                    12d17s20
               if(nok.gt.0)then                                          12d17s20
                irr=0                                                   5d12s23
                do ir1=0,nrootm                                         5d12s23
                 do ir2=0,ir1                                           5d12s23
                  do i=0,nok-1                                             12d17s20
                   icol=ibc(ndenj(isc)+i)                                  12d17s20
                   if(icase.eq.1)then                                     12d18s20
                    jint=jmats(i2eu)+nhere*(icol+nn*irr)                5d12s23
                   else                                                   12d17s20
                    id=icol/irefo(isc)                                    12d18s20
                    ic=icol-id*irefo(isc)                                 12d18s20
                    it=id+irefo(isd)*ic                                   12d17s20
                    jint=jmats(i2eu)+nhere*(it+nn*irr)                  5d12s23
                   end if                                                 12d17s20
                   do m=0,ncsf(iarg)-1                                  5d12s23
                    i10=i1s                                                   12d15s20
                    i1n=nvirt(isbv)                                           12d15s20
                    jjint=jint                                          5d24s23
                    if(ibc(ikeep+icol).eq.1.or.isbv.ne.jsbv)then            12d17s20
                     do i2=i2s,i2e                                             12d15s20
                      if(i2.eq.i2e)i1n=i1e                                     12d15s20
                      iad=idenj(isc)+i2-i2s+nherev*(ir2+nrootu*(m       5d25s23
     $                     +ncsf(iarg)*i))                              5d25s23
                      iadg=jgg-1+nvirt(isbv)*(ir1+nrootu*m)                 5d12s23
                      do i1=i10,i1n                                         12d17s20
                       bc(jjint)=bc(jjint)+bc(iadg+i1)*bc(iad)          5d24s23
                       jjint=jjint+1                                    5d24s23
                      end do                                                12d17s20
                      i10=1                                                 12d17s20
                     end do                                                 12d17s20
                    else                                                    12d17s20
                     do i2=i2s,i2e                                             12d15s20
                      i2m=i2-1                                              12d17s20
                      i2p=i2+1                                              12d17s20
                      if(i2.eq.i2e)i1n=i1e                                     12d15s20
                      iad=idenj(isc)+i2-i2s+nherev*(ir2+nrootu*(m       5d25s23
     $                     +ncsf(iarg)*i))                              5d25s23
                      jintu=jjint-i10                                   5d24s23
                      iadg=jgg-1+nvirt(isbv)*(ir1+nrootu*m)                 5d12s23
                      do i1=i10,min(i1n,i2m)                                12d18s20
                       bc(jintu+i1)=bc(jintu+i1)+bc(iadg+i1)*bc(iad)    5d12s23
                      end do                                                12d18s20
                      do i1=max(i10,i2p),i1n                                12d18s20
                       bc(jintu+i1)=bc(jintu+i1)+bc(iadg+i1)*bc(iad)    5d12s23
                      end do                                                12d18s20
                      jjint=jjint+i1n+1-i10                             5d24s23
                      i10=1                                                 12d17s20
                     end do                                                 12d17s20
                    end if                                                  12d17s20
                   end do                                               5d12s23
                  end do                                                5d12s23
                  irr=irr+1                                             5d12s23
                 end do                                                 5d12s23
                end do                                                  5d12s23
               end if                                                    12d17s20
               if(nokv.gt.0)then                                         12d17s20
c
                do i=0,nokv-1                                            12d17s20
                 icol=ibc(ndenjv(isc)+i)                                 12d17s20
                 do m=0,ncsf(iarg)-1                                     12d17s20
                  irr=0                                                 5d12s23
                  do ir1=0,nrootm                                        12d17s20
                   do ir2=0,ir1                                         5d12s23
                    jjgg=jgg-1+nvirt(isbv)*ir1+ngg*m                      1d27s21
                    ktmpd=idenjv(isc)-i2s+nherev*(ir2                   5d12s23
     $                   +nrootu*(m+ncsf(iarg)*i))                      5d12s23
                    jint=jmats(i2eu)+nhere*(icol+nn*irr)                5d12s23
                    do i2=i2s,i2e                                       6d5s23
                     irow=i2+nvirt(isbv)*(i2-1)                         6d5s23
                     if(irow.ge.il.and.irow.le.ih)then                  6d5s23
                      irow=irow-il                                      6d5s23
                      bc(jint+irow)=bc(jint+irow)+bc(jjgg+i2)           5d12s23
     $                     *bc(ktmpd+i2)                                5d12s23
                     end if
                    end do                                                12d17s20
                    irr=irr+1                                           5d12s23
                   end do                                               5d12s23
                  end do                                                 12d17s20
                 end do                                                  12d17s20
                end do                                                   12d17s20
c
               end if                                                    12d17s20
              end if                                                    12d15s20
             end do                                                     12d15s20
             jgg=jgg+ngg*ncsf(iarg)                                     12d16s20
            end do                                                       12d14s20
            ibcoff=ivecbc                                               2d12s21
           end if                                                       12d14s20
          end do                                                          12d12s20
         end do                                                           12d12s20
        end if                                                           12d12s20
       end do                                                           12d14s20
      end do                                                            12d12s20
      ibcoff=iacc                                                       2d12s21
      return                                                            12d12s20
      end                                                               12d12s20
