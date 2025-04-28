c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcssd1(ihsdiag,nff1,iff1,ncsf,mdon,mdoo,nsymb,multh,   3d4s21
     $    ixw1,ixw2,iden,nh0av,nrootu,ism,irel,                         3d4s21
     $     irefo,nvirt,isymmrci,norb,maxbx,bc,ibc,iroff)                1d27s23
      implicit real*8 (a-h,o-z)
c
c     1-e density
c
      external second                                                   2d18s21
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48,i1c,i1o,j1c,j1o,12d14s20
     $     itestc,itesto,last8(2),gandcc,gandco,gandcb                  2d6s23
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d19s21
      logical lkeep,lprt                                                1d25s21
      dimension nff1(mdoo+1,nsymb,2),ncsf(*),multh(8,8),iden(*),        13d4s21
     $     nh0av(*),nvirt(*),ism(*),irel(*),                            3d4s21
     $     nab4(2,3),irefo(*),iff1(*),ivec(5,8),ivecp(5,8),nvecp(8),    3d4s21
     $     ircv(5,8),nrcv(5,8),itransv(5,8),idend(8)                    3d4s21
      include "common.store"                                            12d12s20
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      iacc=ibcoff                                                       2d12s21
      ibcoff=iacc+mynprocg                                              2d12s21
      idoit=0                                                           3d4s21
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
      ibcoff=igg+maxbx                                                  1d30s21
      nacc=0                                                            1d30s21
      itransgg=0                                                        2d12s21
      call enough('hcssd1.  1',bc,ibc)
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
         ngg=nvirt(isbv)*nrootu                                         12d15s20
         nnc=nff1(ncloip,isb,1)*ncsf(iarg)                              1d28s21
         i18=1                                                          12d15s20
         i28=ngg
         i38=1                                                          12d15s20
         i48=nnc                                                        4d21s21
         call ddi_iget(bc,ibc,ihsdiag(ncloip,isb,1),i18,i28,i38,i48,    11d15s22
     $        bc(igg),ibc(iacc),nacc)                                   11d15s22
         itransgg=0                                                     2d12s21
         do jsb=isb,isb                                                 3d4s21
          jsbv=multh(jsb,isymmrci)                                        12d12s20
          i28=nvirt(jsbv)*nrootu                                         12d14s20
          do nclojp=max(mdon+1,ncloip-1),min(mdoo+1,ncloip+1)           3d4s21
           if(min(nff1(nclojp,jsb,1),nvirt(jsbv)).gt.0)then             3d4s21
            ncloj=nclojp-1                                               12d12s20
            jarg=nclojp-mdon                                             12d12s20
            nrow=nvirt(jsbv)*nrootu                                     3d4s21
            ncoltj=nff1(nclojp,jsb,1)*ncsf(jarg)*nrow                    12d14s20
            ivecbc=ibcoff                                               2d12s21
            i48=nff1(nclojp,jsb,1)*ncsf(jarg)                           12d14s20
            ncol=i48                                                    2d12s21
c
c     has vector already been loaded?
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
            i18=1                                                       3d4s21
            i28=nrootu*nvirt(jsbv)                                      3d4s21
            i38=1                                                       12d14s20
            itransv(nuse,jsb)=0                                         2d12s21
            call ddi_iget(bc,ibc,ihsdiag(nclojp,jsb,1),i18,i28,i38,i48, 11d15s22
     $           bc(ivec(nuse,jsb)),ibc(ircv(nuse,jsb)),nrcv(nuse,jsb)) 2d12s21
   22       continue                                                    2d12s21
            nrowi=nrow*ncsf(iarg)                                        12d14s20
            nrowim=nrowi-1                                               12d14s20
            idenh=ibcoff                                                 12d14s20
            idenhd=idenh+nrowi                                          3d4s21
            ibcoff=idenhd+nrowi                                         3d4s21
            do isc=1,nsymb                                              3d4s21
             idend(isc)=ibcoff                                          3d4s21
             nvvn=(irefo(isc)*(irefo(isc)+1))/2                         3d4s21
             ibcoff=idend(isc)+nvvn*nrowi                               3d4s21
            end do                                                      3d4s21
            ijsbv=multh(isbv,jsbv)                                       12d14s20
            call enough('hcssd1.  2',bc,ibc)
            ibctop=ibcoff-1                                              12d14s20
            if1o=nff1(ncloip,isb,2)                                      12d14s20
            jgg=igg                                                     12d14s20
            jvecer=0                                                    2d12s21
            do if1=1,nff1(ncloip,isb,1)                                  12d14s20
             call second(timea)                                         2d18s21
             ndenh=0                                                    12d14s20
             do i=idenh,ibctop                                           12d14s20
              bc(i)=0d0                                                  12d14s20
             end do                                                      12d14s20
             ndenhd=0                                                   3d5s21
             i1c=iff1(if1o)                                                11d25s20
             ntest=popcnt(i1c)
             if(ntest.ne.ncloi)then
              write(6,*)('ntest = '),ntest,if1
              call dcbit(i1c,32,'dorb')
              call dws_synca
              call dws_finalize
              stop
             end if
             if1o=if1o+1                                                   11d25s20
             i1o=iff1(if1o)                                                11d25s20
             i1o=ibset(i1o,norbx)                                          11d25s20
             nopeni=popcnt(i1o)                                         3d4s21
             if1o=if1o+1                                                   11d25s20
             jf1o=nff1(nclojp,jsb,2)                                      12d14s20
             jvec=ivec(nuse,jsb)                                        2d12s21
             do jf1=1,nff1(nclojp,jsb,1)                                  12d14s20
              j1c=iff1(jf1o)                                                11d25s20
              jf1o=jf1o+1                                                   11d25s20
              j1o=iff1(jf1o)                                                11d25s20
              j1o=ibset(j1o,norbxx)                                          11d25s20
              nopenj=popcnt(j1o)                                        3d4s21
              jf1o=jf1o+1                                                   11d25s20
              if(mod(idoit,mynprocg).eq.mynowprog)then                  3d4s21
               gandcc=ieor(j1c,i1c)                                     2d6s23
               gandco=ieor(j1o,i1o)                                     2d6s23
               gandcb=ior(gandcc,gandco)                                2d6s23
               ndifb=popcnt(gandcb)                                     2d6s23
c
c     vec is ordered nrootu,nvirt. re-order to be nvirt,nrootu.
               if(itransv(nuse,jsb).eq.0.and.ndifb.le.2)then              2d12s21
                call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))        2d12s21
                itmp=ibcoff                                                    1d27s21
                ibcoff=itmp+nrow                                            1d27s21
                call enough('hcssd1.  3',bc,ibc)
                jjvec=ivec(nuse,jsb)                                     2d12s21
                nrowm=nrow-1                                                1d27s21
                do i=0,ncol-1                                               1d27s21
                 do iv=0,nvirt(jsbv)-1                                         1d27s21
                  do ir=0,nrootm                                               1d27s21
                   ij=jjvec+ir+nrootu*iv                                     1d27s21
                   ji=itmp+iv+nvirt(jsbv)*ir                                1d27s21
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
               ndifs=popcnt(gandco)                                     10d14s22
               if(ndifs.eq.2.and.ndifb.eq.2)then                        2d6s23
                do i=1,norbxx                                           2d6s23
                 if(btest(gandco,i))then                                2d6s23
                  if((btest(i1o,i).and..not.btest(j1c,i)).or.           2d6s23
     $                (btest(i1c,i).and.btest(j1o,i)))then              2d6s23
                   nab4(1,1)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,1)=i                                          2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
                call gandc(i1c,i1o,j1c,j1o,nopeni,nopenj,iarg,jarg,ncsf, 2d19s21
     $            norbxx,ixw1,ixw2,nnot,nab1,iwpb1,iwpk1,ncsfmid1,bc,   11d14s22
     $              ibc)                                                11d14s22
                itmp=ibcoff                                               12d14s20
                ibcoff=itmp+nrowi                                         12d14s20
                call enough('hcssd1.  4',bc,ibc)
                ivtrans=ibcoff                                           2d19s21
                ibcoff=ivtrans+nrow*ncsf(jarg)                           2d19s21
                call enough('hcssd1.  5',bc,ibc)
                do j=0,ncsf(jarg)-1                                      2d19s21
                 do i=0,nrow-1                                           2d19s21
                  ij=jvec+i+nrow*j                                       2d19s21
                  ji=ivtrans+j+ncsf(jarg)*i                              2d19s21
                  bc(ji)=bc(ij)                                          2d19s21
                 end do                                                  2d19s21
                end do                                                   2d19s21
                call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid1,nrow,        2d19s21
     $              iwpb1,iwpk1,bc(ivtrans),ncsf(jarg),bc(itmp),        2d19s21
     $              ncsf(iarg),1d0,0d0,bc,ibc)                          11d10s22
                ndenh=1                                                  12d14s20
                do i=0,nrowim                                             12d14s20
                 bc(idenh+i)=bc(idenh+i)+bc(itmp+i)                      2d22s21
                end do                                                    12d14s20
                ibcoff=itmp                                             3d5s21
               end if                                                     12d14s20
               j1o=ibclr(j1o,norbxx)                                     12d14s20
               j1o=ibset(j1o,norbx)                                      12d14s20
               gandcc=ieor(j1c,i1c)                                     2d6s23
               gandco=ieor(j1o,i1o)                                     2d6s23
               gandcb=ior(gandcc,gandco)                                2d6s23
               ndifb=popcnt(gandcb)                                     2d6s23
               ndifs=popcnt(gandco)                                     10d14s22
c
c     vec is ordered nrootu,nherev. re-order to be nherev,nrootu.
               if(itransv(nuse,jsb).eq.0.and.ndifb.le.2)then            2d6s23
                call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))       2d12s21
                itmp=ibcoff                                                    1d27s21
                ibcoff=itmp+nrow                                            1d27s21
                call enough('hcssd1.  6',bc,ibc)
                jjvec=ivec(nuse,jsb)                                    2d12s21
                nrowm=nrow-1                                                1d27s21
                do i=0,ncol-1                                               1d27s21
                 do iv=0,nvirt(jsbv)-1                                         1d27s21
                  do ir=0,nrootm                                               1d27s21
                   ij=jjvec+ir+nrootu*iv                                     1d27s21
                   ji=itmp+iv+nvirt(jsbv)*ir                                1d27s21
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
               nnot=0                                                   2d9s23
               if(ndifb.le.2)then                                       2d6s23
                ivtrans=ibcoff                                           2d19s21
                ibcoff=ivtrans+nrow*ncsf(jarg)                           2d19s21
                call enough('hcssd1.  7',bc,ibc)
                do j=0,ncsf(jarg)-1                                      2d19s21
                 do i=0,nrow-1                                           2d19s21
                  ij=jvec+i+nrow*j                                       2d19s21
                  ji=ivtrans+j+ncsf(jarg)*i                              2d19s21
                  bc(ji)=bc(ij)                                          2d19s21
                 end do                                                  2d19s21
                end do                                                   2d19s21
               end if                                                    2d19s21
               if(ndifb.eq.0)then                                       2d6s23
                nnot=1                                                  2d9s23
                ndenhd=1                                                3d5s21
                do i=1,norb                                             3d4s21
                 fact=0d0                                               3d4s21
                 if(btest(i1c,i))then                                   3d4s21
                  fact=2d0                                              3d4s21
                 else if(btest(i1o,i))then                              3d4s21
                  fact=1d0                                              3d4s21
                 end if                                                 3d4s21
                 if(fact.ne.0d0)then                                    3d4s21
                  jsg=ism(i)                                            3d4s21
                  jbg=irel(i)-1                                         3d4s21
                  ixx=((jbg*(jbg+1))/2)+jbg                             3d4s21
                  jdend=idend(jsg)+nrowi*ixx                            3d4s21
                  do j=0,nrowim                                         3d4s21
                   bc(jdend+j)=bc(jdend+j)+fact*bc(ivtrans+j)           3d5s21
                  end do                                                3d4s21
                 end if                                                 3d4s21
                end do                                                  3d4s21
               else if(ndifs.eq.2.and.ndifb.eq.2)then                   2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandco,i))then                                2d6s23
                  if((btest(i1o,i).and..not.btest(j1c,i)).or.           2d6s23
     $                (btest(i1c,i).and.btest(j1o,i)))then              2d6s23
                   nab4(1,1)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,1)=i                                          2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
                nnot=2                                                  2d9s23
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
                ixx=((lgb*(lgb+1))/2)+lga                               3d4s21
                ndenhd=1                                                3d5s21
                jdend=idend(lsa)+nrowi*ixx                              3d4s21
                ltmp1=ibcoff                                            2d19s21
                ltmp2=ltmp1+ncsf(iarg)*nrow                             2d19s21
                ibcoff=ltmp2+ncsf(iarg)*nrow                             2d19s21
                call enough('hcssd1.  8',bc,ibc)
                call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid1b,nrow,      2d19s21
     $              iwpb1b,iwpk1b,bc(ivtrans),ncsf(jarg),bc(ltmp2),     2d19s21
     $              ncsf(iarg),1d0,0d0,bc,ibc)                          11d10s22
                do ii=0,nrowim                                          4d21s21
                 bc(jdend+ii)=bc(jdend+ii)+bc(ltmp2+ii)                 3d4s21
                end do                                                  2d22s21
                ibcoff=ltmp1                                            3d5s21
               end if                                                    12d14s20
               if(nnot.gt.0.and.nnot.le.2)then
                ibcoff=ivtrans
               end if
              end if                                                     12d14s20
              idoit=idoit+1                                             3d4s21
              jvec=jvec+ncsf(jarg)*nrow                                  12d14s20
             end do                                                      12d14s20
             if(itransgg.eq.0)then                                      2d12s21
              call ddi_done(ibc(iacc),nacc)                             2d12s21
              itransgg=1                                                2d12s21
c
c     transpose igg to virt,root                                     3d4s21
              itmp=ibcoff                                                    1d27s21
              ibcoff=itmp+ngg                                                1d27s21
              call enough('hcssd1.  9',bc,ibc)
              jggg=igg                                                        1d27s21
              do icol=1,nnc                                                  1d28s21
               do i=0,nrootm                                                 1d27s21
                do j=0,nvirt(isbv)-1                                         1d27s21
                 ji=jggg+i+nrootu*j                                     3d4s21
                 ij=itmp+j+nvirt(isbv)*i                                          1d27s21
                 bc(ij)=bc(ji)                                               1d27s21
                end do                                                       1d27s21
               end do                                                        1d27s21
               do i=0,ngg-1                                                  1d28s21
                bc(jggg+i)=bc(itmp+i)                                         1d27s21
               end do                                                        1d27s21
               jggg=jggg+ngg                                            3d4s21
              end do                                                         1d27s21
              ibcoff=itmp                                                    1d27s21
             end if                                                     2d12s21
             if(ndenhd.ne.0)then                                        3d5s21
              do isc=1,nsymb                                             3d4s21
               jden=iden(isc)
               do ir=0,nrootm                                            3d4s21
                do j=0,irefo(isc)-1                                      3d4s21
                 do jj=0,j                                               3d4s21
                  jjden2=jden+j+nh0av(isc)*jj                           3d5s21
                  sum=0d0                                                3d4s21
                  ixx=((j*(j+1))/2)+jj                                     3d4s21
                  do i=0,ncsf(iarg)-1                                      3d4s21
                   jggg=jgg+nvirt(isbv)*(ir+nrootu*i)                    3d4s21
                   do iv=0,nvirt(isbv)-1                                  3d4s21
                    iad=idend(isc)+i                                      3d4s21
     $                 +ncsf(iarg)*(iv+nvirt(isbv)*(ir+nrootu*ixx))       3d4s21
                    sum=sum+bc(jggg+iv)*bc(iad)                          3d4s21
                   end do                                                 3d4s21
                  end do                                                  3d4s21
                  if(jj.ne.j)sum=sum*0.5d0                              3d19s21
                  bc(jjden2)=bc(jjden2)+sum                              3d4s21
                 end do                                                  3d4s21
                end do                                                   3d4s21
                jden=jden+iroff*nh0av(isc)*nh0av(isc)                   1d27s23
               end do                                                    3d4s21
              end do                                                     3d4s21
             end if
             if(ndenh.ne.0)then                                         3d5s21
              jden=iden(isbv)+irefo(isbv)+nh0av(isbv)*irefo(isbv)        3d4s21
              do ir=0,nrootm                                             3d4s21
               do i=0,ncsf(iarg)-1                                       3d4s21
                jggg=jgg+nvirt(isbv)*(ir+nrootu*i)                       3d4s21
                do iv=0,nvirt(isbv)-1                                    3d4s21
                 jjden=jden+nh0av(isbv)*iv                               3d4s21
                 iad=idenh+i+ncsf(iarg)*(iv+nvirt(isbv)*ir)              3d4s21
                 do jv=iv,nvirt(isbv)-1                                 3d5s21
                  bc(jjden+jv)=bc(jjden+jv)+bc(iad)*bc(jggg+jv)          3d4s21
                 end do                                                  3d4s21
                end do                                                   3d4s21
               end do                                                    3d4s21
               jden=jden+iroff*nh0av(isbv)*nh0av(isbv)                  1d27s23
              end do                                                     3d4s21
             end if                                                     3d5s21
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
