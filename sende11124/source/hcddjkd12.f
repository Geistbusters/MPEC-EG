c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcddjko(nff22,nfdat,gd,vd,nsymb,mdon,mdoo,nec,multh,    1d8s21
     $     isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,ismult,ixw1,ixw2,   1d8s21
     $     norb,nrootu,ih0av,nh0av,ioooo,jmats,kmats,ndoub,mdoub,shift, 2d18s21
     $     tdendd,tovr,sr2,srh,bc,ibc)                                  9d28s23
      implicit real*8 (a-h,o-z)                                         1d8s21
c                                                                       1d8s21
c     0 and 2 virt contribution to gd=hdd*vd                            1d8s21
c                                                                       1d8s21
      logical lprt,lchoice                                              3d19s21
      include "common.store"                                            1d8s21
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d22s21
      integer*8 ipack8,itc,ito,jtc,jto,itestc,itesto,itesta,itestb,     1d17s21
     $     jtesta,jtestb,gandcc,gandco,gandcb                           11d1s22
      external second                                                   2d18s21
      equivalence (ipack8,ipack4)                                       1d8s21
      dimension nff22(mdoo+1,2,nsymb),nfdat(5,4,*),vd(*),gd(*),         1d8s21
     $     multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),irel(*),ism(*),       1d8s21
     $     irefo(*),ih0av(*),nh0av(*),ioooo(*),jmats(*),kmats(*),       1d8s21
     $     ipack4(2),nprei(4),mprei(4),nli(4),nprej(4),                 1d12s21
     $     mprej(4),nlj(4),nab4(2,3),idenj(4,8),ndenj(4,8),             1d9s21
     $     mdenj(4,8),idenk(4,8),ndenk(4,8),mdenk(4,8),loope(6),                 1d9s21
     $     idenhvv(4),mdenhvv(4),itest(32,3),iden1e(4),mden1e(4),       1d8s21
     $     ldenj(4,8),ldenk(4,8),iden1ef(4),idenhvvf(4),ivects(4),      11d3s22
     $     idenjf(4,8),idenkf(4,4,8),ndenjf(4,8),ndenkf(4,4,8),
     $     nokj(4,8),nokk(4,4,8),isorb(32),idorb(32),idenjn(4,8),       2d23s21
     $     idenkn(4,4,8),iden1en(4),idenhvvn(4),tdendd(*),ioxx(2)
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/fnd2cm/inv(2,8,8,8)                                        9d2s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/1300000/
      write(6,*)('Hi, my name is hcddjko'),ibcoff
      igoal=2133967
      loop=0
      nnes=0
      ndelta2=0
      nsamex=0
      ntotcalc=0
      nnzznn=0
      mplanb=0
      ibcoffo=ibcoff
      norbx=norb+1                                                      1d8s21
      norbxx=norbx+1                                                    1d8s21
      norbxxx=norbxx+1                                                  1d8s21
      ioffdnon=1                                                        1d8s21
      nrootm=nrootu-1                                                   1d11s21
      ioffd=1                                                           1d11s21
      sumjj=0d0
      do isb=1,nsymb                                                    1d8s21
       isbv12=multh(isb,isymmrci)                                       1d8s21
       joffdnon=1                                                       1d8s21
       do jsb=1,isb                                                     1d8s21
        ijsb=multh(jsb,isb)                                             1d9s21
        jsbv12=multh(jsb,isymmrci)                                      1d8s21
        ibc0=ibcoff                                                     1d8s21
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    2d24s21
          do ll=1,4                                                     1d8s21
           if(nfdat(2,ll,jsb).gt.0)then                                 2d24s21
            nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                          2d24s21
            if(l.eq.ll)then                                             1d11s21
             iden1en(l)=ibcoff                                          2d24s21
             idenhvvn(l)=iden1en(l)+nll                                 2d24s21
             ibcoff=idenhvvn(l)+nll                                     2d24s21
            end if                                                      1d11s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(lsb.ge.lsa.and.l.eq.ll)then                             1d11s21
              if(lsb.eq.lsa)then                                        1d9s21
               nn=(irefo(lsa)*(irefo(lsa)+1))/2                         1d9s21
              else                                                      1d9s21
               nn=irefo(lsa)*irefo(lsb)                                 1d9s21
              end if                                                    1d9s21
              ndenjf(l,lsa)=ibcoff                                      2d25s21
              ibcoff=ndenjf(l,lsa)+nn                                    1d11s21
              idenjn(l,lsa)=ibcoff                                      1d11s21
              ibcoff=idenjn(l,lsa)+nn*nll                               2d24s21
             end if                                                     1d9s21
             nn=irefo(lsa)*irefo(lsb)                                   1d9s21
             ndenkf(l,ll,lsa)=ibcoff                                    2d25s21
             ibcoff=ndenkf(l,ll,lsa)+nn                                     1d11s21
             idenkn(l,ll,lsa)=ibcoff                                    1d9s21
             ibcoff=idenkn(l,ll,lsa)+nn*nll                             2d24s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call enough('hcddjk.  1',bc,ibc)
        nwds=ibcoff-ibc0                                                1d8s21
        do i=ibc0,ibcoff-1                                              1d8s21
         bc(i)=0d0                                                      1d8s21
        end do                                                          1d8s21
        ibclast=ibcoff                                                  1d17s21
        idoit=0                                                         2d23s21
        do ncloi=mdon,mdoo                                              1d8s21
         ncloip=ncloi+1                                                 1d8s21
         if(nff22(ncloip,1,isb).gt.0)then                               1d8s21
          iarg=ncloip-mdon                                              1d8s21
          nopeni=nec-2*ncloip                                           1d8s21
          nopenip=nopeni+2                                              1d8s21
          ibcst=ibcoff                                                  1d8s21
          ibcnd=ibcoff-1                                                1d8s21
          ivcv=nfdat(5,1,isb)+nff22(ncloip,2,isb)                       1d8s21
          do if=1,nff22(ncloip,1,isb)                                   1d8s21
           ipack8=ibc(ivcv)                                             1d8s21
           itc=ipack4(1)                                                1d8s21
           ito=ipack4(2)                                                1d8s21
           ito=ibset(ito,norbx)                                         1d8s21
           ito=ibset(ito,norbxx)                                        1d8s21
           nspacei=ibc(ivcv+1)                                          3d19s21
           ibcvects=ibcoff                                              11d3s22
           do l=1,4                                                     3d19s21
            nli(l)=ibc(ivcv+1+l)                                        3d19s21
            if(nli(l).gt.0)then                                         11d3s22
             iad1=ivcv+ibc(ivcv+5+l)                                    11d3s22
             iad2=iad1+nli(l)                                           11d3s22
             ivects(l)=ibcoff                                           11d3s22
             ibcoff=ivects(l)+nli(l)*ncsf2(l,iarg)                      11d3s22
             call enough('hcddjk.1.1',bc,ibc)
             do i=0,nli(l)-1                                            11d3s22
              do j=0,ncsf2(l,iarg)                                      11d3s22
               ji=iad2+j+ncsf2(l,iarg)*i                                11d3s22
               ij=ivects(l)+i+nli(l)*j                                  11d3s22
               bc(ij)=bc(ji)                                            11d3s22
              end do                                                    11d3s22
             end do                                                     11d3s22
            end if                                                      11d3s22
           end do                                                       3d19s21
           do i=1,norb                                                  1d8s21
            itest(i,1)=0                                                1d8s21
           end do                                                       1d8s21
           do i=1,norb                                                  1d8s21
            if(btest(itc,i))itest(i,1)=2                                1d8s21
            if(btest(ito,i))itest(i,1)=1                                1d8s21
           end do                                                       1d8s21
           itest(norbx,1)=1                                             1d8s21
           itest(norbxx,1)=1                                            1d8s21
           do ncloj=max(ncloi-2,mdon),min(ncloi+2,mdoo)                 1d8s21
            nclojp=ncloj+1                                                1d8s21
            if(nff22(nclojp,1,jsb).gt.0)then                            1d14s21
             if(mod(idoit,mynprocg).eq.mynowprog)then                   3d1s21
              jarg=nclojp-mdon                                             1d8s21
              nopenj=nec-2*nclojp                                          1d8s21
              nopenjp=nopenj+2                                           1d8s21
              jvcv=nfdat(5,1,jsb)+nff22(nclojp,2,jsb)                      1d8s21
              do jf=1,nff22(nclojp,1,jsb)                                  1d8s21
               ipack8=ibc(jvcv)                                            1d8s21
               jtc=ipack4(1)                                               1d8s21
               jto=ipack4(2)                                               1d8s21
               jto=ibset(jto,norbx)                                        1d8s21
               jto=ibset(jto,norbxx)                                       1d8s21
               nspacej=ibc(jvcv+1)                                       3d19s21
               do l=1,4                                                  3d19s21
                nlj(l)=ibc(jvcv+1+l)                                     3d19s21
               end do                                                    3d19s21
               if(isb.eq.jsb)then                                        1d17s21
c                                                                       1d8s21
                gandcc=ieor(itc,jtc)                                     10d14s22
                gandco=ieor(ito,jto)                                     10d14s22
                gandcb=ior(gandcc,gandco)                                10d21s22
                ndifb=popcnt(gandcb)                                     10d21s22
                if(ndifb.le.4)then                                       10d21s22
                 ndifs=popcnt(gandco)                                     10d14s22
                 ndifd=popcnt(gandcc)                                     10d14s22
                 if(ndifd.eq.0.and.ndifs.eq.0)then                       10d21s22
                  ii=1                                                           11d21s20
                  do i=1,norb                                                    11d21s20
                   if(btest(ipack4(1),i))then                                    11d21s20
                    idorb(ii)=i                                                  11d21s20
                    ii=ii+1                                                      11d21s20
                   end if                                                        11d21s20
                  end do                                                         11d21s20
                  ii=1                                                           11d21s20
                  do i=1,norb                                                    11d21s20
                   if(btest(ipack4(2),i))then                                    11d21s20
                    isorb(ii)=i                                                  11d21s20
                    ii=ii+1                                                      11d21s20
                   end if                                                        11d21s20
                  end do                                                         11d21s20
                  sumx=0d0                                                 2d24s21
                  sumc=0d0                                                         9d23s19
                  do i=1,ncloi                                                   11d23s20
                   ig=irel(idorb(i))-1                                           11d23s20
                   is=ism(idorb(i))                                              11d23s20
                   iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
                   h0=bc(iad)
                   do j=1,ncloi                                                  11d23s20
                    jg=irel(idorb(j))-1                                          11d23s20
                    js=ism(idorb(j))                                             11d23s20
                    xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)      9d28s23
                    xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)      9d28s23
                    sumc=sumc+2d0*xj-xk                                            7d15s19
                   end do                                                          7d15s19
                   sumc=sumc+h0*2d0                                                7d15s19
                  end do                                                           7d15s19
                  sumo=0d0                                                         7d15s19
                  do i=1,nopeni                                                  11d23s20
                   ig=irel(isorb(i))-1                                           11d23s20
                   is=ism(isorb(i))                                              11d23s20
                   iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
                   h0=bc(iad)
                   sumo=sumo+h0                                                    7d15s19
                   do j=1,nopeni                                                 11d23s20
                    if(i.ne.j)then                                                 7d15s19
                     jg=irel(isorb(j))-1                                         11d23s20
                     js=ism(isorb(j))                                            11d23s20
                     xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)     9d28s23
                     sumo=sumo+xj*0.5d0                                            7d15s19
                    end if                                                         7d15s19
                   end do                                                          7d15s19
                  end do                                                           7d15s19
                  sum=sumc+sumo
                  do i=1,ncloi                                                   11d23s20
                   is=ism(idorb(i))                                              11d23s20
                   ig=irel(idorb(i))-1                                           11d23s20
                   do j=1,nopeni                                                 11d23s20
                    js=ism(isorb(j))                                             11d23s20
                    jg=irel(isorb(j))-1                                          11d23s20
                    xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)      9d28s23
                    xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)      9d28s23
                    sum=sum+xj*2d0-xk                                              7d15s19
                   end do                                                          7d15s19
                  end do                                                           7d15s19
                  sumx=sum                                               2d24s21
                  do i1=1,nopeni                                                 11d23s20
                   jsa=ism(isorb(i1))                                            11d22s20
                   jga=irel(isorb(i1))-1                                         11d22s20
                   do i2=i1+1,nopenip                                            11d19s20
                    if(i2.le.nopeni)then                                         11d23s20
                     jtesta=itc                                               11d19s20
                     jtestb=ito                                               11d19s20
                     nopenk=nopenip-2                                            11d19s20
                     karg=iarg+1                                                11d19s20
                     nqq=karg+mdon-1
                     xint=0d0
                     ksb=ism(isorb(i2))                                          11d22s20
                     kgb=irel(isorb(i2))-1                                       11d22s20
                     xint=getint(ioooo,jsa,ksb,ksb,jsa,jga,kgb,kgb,jga,  9d28s23
     $                   bc,ibc)                                        9d28s23
                     sumx=sumx-xint                                      2d24s21
                     if(nqq.ge.mdon.and.nqq.le.mdoo)then                         11d19s20
                      jtestb=ibclr(jtestb,isorb(i2))                      2d10s21
                      jtestb=ibclr(jtestb,isorb(i1))                      2d10s21
                      jtesta=ibset(jtesta,isorb(i1))                      2d10s21
                      call gandc(itc,ito,jtesta,jtestb,nopenip,nopenk,      11d19s20
     $             iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot,nab,iwpb,iwpk,       1d17s21
     $                   ncsfmid,bc,ibc)                                9d28s23
                      call gandc(jtesta,jtestb,itc,ito,nopenk,nopenip,   2d25s21
     $             karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,    2d25s21
     $                     iwpk1,ncsfmid1,bc,ibc)                       9d28s23
                      iargo=1                                            2d24s21
                      do l=1,4                                           2d24s21
                       if(nli(l).gt.0)then                               2d24s21
                        iad1=ivcv+ibc(ivcv+5+l)                            3d19s21
                        iad2=iad1+nli(l)                                    3d19s21
                        itmp1=ibcoff                                     2d24s21
                        itmp2=itmp1+nli(l)*ncsf(karg)                    2d24s21
                        itmp3=itmp2+nli(l)*ncsf(karg)                    2d24s21
                        ibcoff=itmp3+nli(l)*nli(l)                       2d24s21
                        call enough('hcddjk.  2',bc,ibc)
                        call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1,    2d25s21
     $                       nli(l),iwpb1,iwpk1,bc(iad2),ncsf2(l,iarg), 2d25s21
     $                       bc(itmp1),ncsf(karg),1d0,0d0,iargo,        2d24s21
     $                       ncsf2(l,iarg),bc,ibc)                      9d28s23
                        do iii=0,nli(l)-1                                2d24s21
                         do k=0,ncsf(karg)-1                             2d24s21
                          ki=itmp1+k+ncsf(karg)*iii                      2d24s21
                          ik=itmp2+iii+nli(l)*k                          2d24s21
                          bc(ik)=bc(ki)                                  2d24s21
                         end do                                          2d24s21
                        end do                                           2d24s21
                        call dgemm('n','n',nli(l),nli(l),ncsf(karg),     2d24s21
     $                       xint,bc(itmp2),nli(l),bc(itmp1),ncsf(karg),2d24s21
     $                       0d0,bc(itmp3),nli(l),                      2d24s21
     d' hcddjk.  1')
                        jtmp3=itmp3                                      2d24s21
                        do j=0,nli(l)-1                                   2d23s21
                         jj=ibc(iad1+j)-1                                  1d8s21
                         jden=iden1en(l)+nfdat(2,l,isb)*jj-1             2d24s21
                         do k=0,nli(l)-1                                 2d24s21
                          kk=ibc(iad1+k)                                 2d24s21
                          bc(jden+kk)=bc(jden+kk)+bc(jtmp3+k)            2d24s21
                         end do                                          2d24s21
                         jtmp3=jtmp3+nli(l)                              2d24s21
                        end do                                           2d24s21
                        ibcoff=itmp1                                     2d24s21
                       end if                                            2d24s21
                       iargo=iargo+ncsf2(l,iarg)                         2d24s21
                      end do                                             2d24s21
                     end if                                                      11d19s20
                    end if                                                       11d23s20
                   end do                                                1d18s21
                  end do                                                 1d18s21
                  iargo=1                                                2d24s21
                  do l=1,4                                               2d24s21
                   if(nli(l).gt.0)then                                   2d24s21
                    itmpt=ibcoff                                         2d24s21
                    itmpm=itmpt+ncsf2(l,iarg)*nli(l)                     2d24s21
                    ibcoff=itmpm+nli(l)*nli(l)                           2d24s21
                    call enough('hcddjk.  3',bc,ibc)
                    iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                    iad2=iad1+nli(l)                                     3d19s21
                    do iii=0,nli(l)-1                                    2d24s21
                     do j=0,ncsf2(l,iarg)-1                              2d24s21
                      ji=iad2+j+ncsf2(l,iarg)*iii                        2d24s21
                      ij=itmpt+iii+nli(l)*j                              2d24s21
                      bc(ij)=bc(ji)                                      2d24s21
                     end do                                              2d24s21
                    end do                                               2d24s21
                    call dgemm('n','n',nli(l),nli(l),ncsf2(l,iarg),sumx, 2d24s21
     $                   bc(itmpt),nli(l),bc(iad2),ncsf2(l,iarg),0d0,   2d24s21
     $                   bc(itmpm),nli(l),                              2d24s21
     d' hcddjk.  2')
                    jtmpm=itmpm                                          2d24s21
                    do j=0,nli(l)-1                                      2d24s21
                     jj=ibc(iad1+j)-1                                    2d24s21
                     jden=iden1en(l)+nfdat(2,l,isb)*jj-1                 2d24s21
                     do k=0,nli(l)-1                                     2d24s21
                      kk=ibc(iad1+k)                                     2d24s21
                      bc(jden+kk)=bc(jden+kk)+bc(jtmpm+k)
                     end do                                              2d24s21
                     jtmpm=jtmpm+nli(l)                                  2d24s21
                    end do                                               2d24s21
                    ibcoff=itmpt                                         2d24s21
                   end if                                                2d24s21
                   iargo=iargo+ncsf2(l,iarg)                             2d24s21
                  end do                                                 2d24s21
                 else if(ndifs.eq.2.and.ndifb.eq.2)then                  10d21s22
                  do i=1,norbxx
                   if(btest(gandco,i))then                                          10d14s22
                    if((btest(ito,i).and..not.btest(jtc,i)).or.
     $                  (btest(itc,i).and.btest(jto,i)))then                           10d14s22
                     nab4(1,1)=i
                    else                                                           10d14s22
                     nab4(2,1)=i
                    end if
                   end if                                                          10d14s22
                  end do                                                           10d14s22
                  call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,  1d8s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $               ncsfmid1,bc,ibc)                                   9d28s23
                  call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,  2d24s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,  2d24s21
     $               ncsfmid1b,bc,ibc)                                  9d28s23
                  ksa=ism(nab1(1))                                       1d8s21
                  kga=irel(nab1(1))-1                                    1d8s21
                  ksb=ism(nab1(2))                                       1d8s21
                  kgb=irel(nab1(2))-1                                    1d8s21
                  ih0=ih0av(ksa)+kga+nh0av(ksa)*kgb                      1d8s21
                  sum=bc(ih0)                                            1d8s21
                  do i=1,norb                                            12d18s20
                   itest(i,2)=0                                          12d18s20
                  end do                                                 12d18s20
                  do i=1,norb                                            12d18s20
                   if(btest(jtc,i))then                                  12d18s20
                    itest(i,2)=2                                         12d18s20
                   end if                                                12d18s20
                   if(btest(jto,i))then                                  12d18s20
                    itest(i,2)=1                                         12d18s20
                   end if                                                12d18s20
                  end do                                                 12d18s20
                  nok=0                                                         11d13s20
                  do i=1,norb                                            1d8s21
                   ixn=min(itest(i,1),itest(i,2))
                   if(ixn.gt.0)then                                             11d13s20
                    nok=nok+1                                                   11d13s20
                    itest(nok,3)=ixn                                            11d13s20
                    itest(nok,2)=i                                              11d13s20
                   end if                                                       11d13s20
                  end do                                                        11d13s20
                  do i=1,nok                                             1d8s21
                   lsa=ism(itest(i,2))                                   1d8s21
                   lga=irel(itest(i,2))-1                                1d8s21
                   xj=getint(ioooo,lsa,lsa,ksa,ksb,lga,lga,kga,kgb,bc,   9d28s23
     $                 ibc)                                             9d28s23
                   sum=sum+xj                                            1d8s21
                   if(itest(i,3).eq.2)then                               1d8s21
                    xk=getint(ioooo,lsa,ksa,lsa,ksb,lga,kga,lga,kgb,bc,  9d28s23
     $                 ibc)                                             9d28s23
                    sum=sum+xj-xk                                        1d8s21
                   end if                                                1d8s21
                  end do                                                 1d8s21
                  jargo=1                                                2d23s21
                  iargo=1                                                2d24s21
                  do l=1,4                                                1d8s21
                   if(min(nlj(l),nli(l)).gt.0)then                       2d25s21
                    jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                    jad2=jad1+nlj(l)                                     3d19s21
                    iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                    iad2=iad1+nli(l)                                     3d19s21
                    itmp1=ibcoff                                         2d23s21
                    itmp1b=itmp1+ncsfmid1*nlj(l)                         2d23s21
                    ibcoff=itmp1b+ncsfmid1*nli(l)                        2d23s21
                    call enough('hcddjk.  4',bc,ibc)
                    if(iwpk1.lt.0)then                                   2d23s21
                     nusedi=ibc(-iwpk1)/2                                2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1          2d23s21
                     icmp1=-iwpk1+1                                      2d23s21
                     icmp2=icmp1+ncsf(jarg)                              2d23s21
                     icmp3=icmp2+nusedi                                  2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   bc(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),   2d23s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                    else                                                              12d5s20
                     jwpk1=iwpk1+ncsfmid1*(jargo-1)                      2d23s21
                     call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                   1d0,bc(jwpk1),ncsfmid1,bc(jad2),ncsf2(l,jarg),
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjk.  3')
                    end if                                                            12d5s20
                    if(iwpk1b.lt.0)then                                   2d23s21
                     nusedi=ibc(-iwpk1b)/2                                2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                     icmp1=-iwpk1b+1                                      2d23s21
                     icmp2=icmp1+ncsf(iarg)                              2d23s21
                     icmp3=icmp2+nusedi                                  2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   bc(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),   2d23s21
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                    else                                                              12d5s20
                     jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                     call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                   1d0,bc(jwpk1b),ncsfmid1,bc(iad2),ncsf2(l,iarg),
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjk.  4')
                    end if                                                            12d5s20
                    itmpt=ibcoff                                         2d23s21
                    itmpm=itmpt+nli(l)*ncsfmid1                          2d24s21
                    ibcoff=itmpm+nli(l)*nlj(l)                           2d24s21
                    call enough('hcddjk.  5',bc,ibc)
                    do iii=0,nli(l)-1                                    2d23s21
                     do j=0,ncsfmid1-1                                   2d23s21
                      ji=itmp1b+j+ncsfmid1*iii                           2d23s21
                      ij=itmpt+iii+nli(l)*j                              2d23s21
                      bc(ij)=bc(ji)                                      2d23s21
                     end do                                              2d23s21
                    end do                                               2d23s21
                    call dgemm('n','n',nli(l),nlj(l),ncsfmid1,sum,       2d24s21
     $                   bc(itmpt),nli(l),bc(itmp1),ncsfmid1,0d0,        2d23s21
     $                  bc(itmpm),nli(l),                               2d23s21
     d' hcddjk.  5')
                    jtmpm=itmpm                                          2d24s21
                    do iii=0,nlj(l)-1                                    2d24s21
                     jj=ibc(jad1+iii)-1                                  2d24s21
                     jjden=iden1en(l)+nfdat(2,l,isb)*jj-1                2d24s21
                     do j=0,nli(l)-1                                     2d24s21
                      kk=ibc(iad1+j)                                     2d24s21
                      bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)              2d24s21
                     end do                                              2d24s21
                     jtmpm=jtmpm+nli(l)                                  2d24s21
                    end do                                               2d24s21
                    ibcoff=itmp1                                         2d24s21
                   end if                                                 1d8s21
                   jargo=jargo+ncsf2(l,jarg)                             2d23s21
                   iargo=iargo+ncsf2(l,iarg)                             2d24s21
                  end do                                                  1d8s21
                  do i=1,nok                                             1d8s21
                   if(itest(i,3).eq.1)then                                  12d14s20
                    itestc=itc                                              12d14s20
                    itesto=ito                                              12d14s20
                    nopenk=nopenip                                                11d13s20
c
c     anihilate common
c
                    if(btest(itestc,itest(i,2)))then                             11d13s20
                     itestc=ibclr(itestc,itest(i,2))                             11d13s20
                     itesto=ibset(itesto,itest(i,2))                             11d13s20
                     karg=iarg-1                                                11d13s20
                     nopenk=nopenk+1                                             11d13s20
                    else                                                         11d13s20
                     itesto=ibclr(itesto,itest(i,2))                             11d13s20
                     karg=iarg                                                  11d13s20
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
                    nnot1=0                                              8d4s22
                    nnot2=0                                              8d4s22
                    if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                     call gandc(itestc,itesto,itc,ito,nopenk,nopenip,    2d24s21
     $                karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,    2d24s21
     $                iwpb1b,iwpk1b,ncsfmid1b,bc,ibc)                   9d28s23
                     nnot1=nnot1b                                        2d25s21
                     nab1(1)=nab1b(2)                                    2d25s21
                     nab1(2)=nab1b(1)                                    2d25s21
                     call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,        12d14s20
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            9d28s23
                    end if
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     lsa=ism(nab1(1))                                    1d8s21
                     lga=irel(nab1(1))-1                                 1d11s21
                     lsb=ism(nab1(2))                                    1d8s21
                     lgb=irel(nab1(2))-1                                 1d11s21
                     lsc=ism(nab2(1))                                    1d8s21
                     lgc=irel(nab2(1))-1                                 1d11s21
                     lsd=ism(nab2(2))                                    1d8s21
                     lgd=irel(nab2(2))-1                                 1d11s21
                     xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  9d28s23
     $                   bc,ibc)                                        9d28s23
                     jargo=1                                             2d23s21
                     iargo=1                                             2d24s21
                     do l=1,4                                             1d8s21
                      if(min(ncsf(karg),ncsf(jarg),ncsfmid2,             4d12s22
     $                     nli(l),nlj(l)).gt.0)then                      4d12s22
                       jad1=jvcv+ibc(jvcv+5+l)                           3d19s21
                       jad2=jad1+nlj(l)                                  3d19s21
                       iad1=ivcv+ibc(ivcv+5+l)                           3d19s21
                       iad2=iad1+nli(l)                                  3d19s21
                       itmpj=ibcoff                                      2d24s21
                       ibcoff=itmpj+ncsf(karg)*nlj(l)                    2d24s21
                       call enough('hcddjk.  6',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,bc(jad2),ncsf2(l,jarg),     2d24s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                       itmpi=ibcoff                                      2d24s21
                       itmpt=itmpi+ncsf(karg)*nli(l)                     2d24s21
                       ibcoff=itmpt+ncsf(karg)*nli(l)                    2d24s21
                       call enough('hcddjk.  7',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,bc(iad2),ncsf2(l,iarg),   2d24s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                       do iii=0,nli(l)-1                                 2d24s21
                        do k=0,ncsf(karg)-1                                2d23s21
                         ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                         ik=itmpi+iii+nli(l)*k                            2d23s21
                         bc(ik)=bc(ki)                                     2d23s21
                        end do                                             2d23s21
                       end do                                              2d23s21
                       itmpm=itmpt                                       2d24s21
                       ibcoff=itmpm+nli(l)*nlj(l)                        2d24s21
                       call enough('hcddjk.  8',bc,ibc)
                       call dgemm('n','n',nli(l),nlj(l),ncsf(karg),      2d25s21
     $                     xint,                                        2d25s21
     $                      bc(itmpi),nli(l),bc(itmpj),ncsf(karg),0d0,  2d24s21
     $                      bc(itmpm),nli(l),                           2d24s21
     d' hcddjk.  6')
                       jtmpm=itmpm                                       2d24s21
                       do iii=0,nlj(l)-1                                   2d24s21
                        jj=ibc(jad1+iii)-1                                 2d24s21
                        jjden=iden1en(l)+nfdat(2,l,isb)*jj-1               2d24s21
                        do j=0,nli(l)-1                                    2d24s21
                         kk=ibc(iad1+j)                                    2d24s21
                         bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)             2d24s21
                        end do                                             2d24s21
                        jtmpm=jtmpm+nli(l)                                 2d24s21
                       end do                                              2d24s21
                       ibcoff=itmpj                                      2d24s21
                      end if                                              1d8s21
                      jargo=jargo+ncsf2(l,jarg)                          2d23s21
                      iargo=iargo+ncsf2(l,iarg)                          2d24s21
                     end do                                               1d8s21
                    end if                                                  12d14s20
                   end if                                                   12d14s20
                  end do                                                 1d8s21
                 else                                                    1d8s21
                  nnot=0                                                 10d21s22
                  if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s
                   nnot=4                                                          10d14s22
                   ioxx(1)=1                                                     10d17s22
                   ioxx(2)=1                                                     10d17s22
                   do i=1,norbxx                                                     10d17s22
                    if(btest(gandcb,i))then                                         10d14s22
                     if((btest(jtc,i).and.btest(ito,i)).or.                         10d17s22
     $                   (btest(jto,i).and..not.btest(itc,i)))then                   10d14s22
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
     $                   ((btest(itc,i).and..not.btest(jto,i)).or.                 10d14s22
     $                   (btest(jtc,i).and..not.btest(ito,i))))then                     10d14s22
                      if(btest(jtc,i))iswap=1                                        10d17s22
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
                    if(btest(itc,nab4(1,2)).and.
     $                   .not.btest(itc,nab4(1,1)))nbt=1                10d21s22
                   else                                                             10d17s22
                    nbt=0                                                           10d17s22
                    if(btest(jtc,nab4(2,2)).and.
     $                   .not.btest(jtc,nab4(2,1)))nbt=1                10d21s22
                   end if                                                           10d17s22
                   if(nbt.ne.0)then                                                 10d17s22
                    nab4(1,1)=nab4(1,2)                                             10d17s22
                    nab4(2,1)=nab4(2,2)                                             10d17s22
                   end if                                                           10d17s22
                  else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                   nnot=3
                   do i=1,norbxx
                    if(btest(gandcb,i))then                                         10d14s22
                     if(btest(itc,i))then
                      nab4(1,1)=i
                      nab4(1,2)=i
                     else                                                           10d14s22
                      nab4(2,1)=i
                      nab4(2,2)=i
                     end if
                    end if                                                          10d14s22
                   end do
                  end if                                                            10d14s22
                  ipssx=0                                                10d21s22
                  if(nnot.eq.3)then                                          12d8s20
                   ipssx=1                                                   12d8s20
                  else if(nnot.eq.4)then                                 11d1s22
                   ipssx=3                                                 12d18s20
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
                   itestc=itc                                              12d8s20
                   itesto=ito                                              12d8s20
                   if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                    itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                    itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                    nopenk=nopenip+1                                              11d13s20
                    karg=iarg-1                                             12d8s20
                   else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                    itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                    nopenk=nopenip-1                                              11d13s20
                    karg=iarg
                   else                                                          11d13s20
                    write(6,*)('bit not set for nab4(1,1) ='),           3d1s21
     $                   nab4(1,iu1)                                    3d1s21
                    stop 'nab4(1,1)'                                           11d27s20
                   end if                                                        11d13s20
                   if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                    itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                    itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                    nopenk=nopenk-1                                              11d13s20
                    karg=karg+1                                                  11d13s20
                   else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                    write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                 nab4(2,iu2)                                       12d18s20
                    stop 'nab4(2,1)'                                           11d27s20
                   else                                                          11d13s20
                    itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                    nopenk=nopenk+1                                              11d13s20
                   end if                                                        11d13s20
                   nqq=karg+mdon-1                                       4d20s21
                   nnot1=0                                               8d4s22
                   nnot2=0                                               8d4s22
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                   4d20s21
c
c     plan a: both contraction vectors untransposed, do
c     iwpb2*iwpk2*vecj, iwpb1b*iwpk1b*veci, then multiply the 2.
c     plan b: veci is transposed, form iwpb1bt*iwpb2 outside of l
c     loop, then multiply 5 matrices together via genmatn3
c     which is faster will depend on large ncsf(karg) is (a larger
c     value favoring plan b), and how small nli and nlj are
c     (small values favor plan a).
c
                    call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,       1d8s21
     $         karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d18s20
     $         ncsfmid2,bc,ibc)                                         9d28s23
                    rtb=dfloat(ncsf(karg))/dfloat(ncsfmid2)              11d3s22
                    rta=0d0                                              11d3s22
                    nls=0                                                11d3s22
                    do l=1,4                                             11d3s22
                     if(min(nli(l),nlj(l)).gt.0.and.
     $                   min(ncsf2(l,iarg),ncsf2(l,jarg)).gt.0)then          11d3s22
                      nls=nls+1                                          11d3s22
                      rta=max(rta,dfloat(ncsf2(l,jarg))/dfloat(nlj(l)))  11d3s22
                     end if                                              11d3s22
                    end do                                               11d3s22
                    nplanb=0
                    if(rtb.gt.-rta)then                                   11d3s22
                     call gandc(itc,ito,itestc,itesto,nopenip,nopenk,    11d3s22
     $         iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,     2d25s21
     $                  iwpk1b,ncsfmid1b,bc,ibc)                        9d28s23
                     nab1(1)=nab1b(1)                                    11d3s22
                     nab1(2)=nab1b(2)                                    11d3s22
                     lsa=ism(nab1(1))                                     1d8s21
                     lga=irel(nab1(1))-1                                  1d8s21
                     lsb=ism(nab1(2))                                     1d8s21
                     lgb=irel(nab1(2))-1                                  1d8s21
                     lsc=ism(nab2(1))                                     1d8s21
                     lgc=irel(nab2(1))-1                                  1d8s21
                     lsd=ism(nab2(2))                                     1d8s21
                     lgd=irel(nab2(2))-1                                  1d8s21
                     xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  9d28s23
     $                   bc,ibc)                                        9d28s23
                     if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2))     4d20s21
     $                   .and.                                          4d20s21
     $                max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))     1d12s21
     $                  xint=xint*0.5d0                                 1d12s21
                     if(abs(xint).gt.1d-14)then                          11d3s22
                      nplanb=1
                      nsz=ncsfmid1b*ncsfmid2                             11d3s22
                      iprod=ibcoff                                        11d3s22
                      ibcoff=iprod+nsz                                   11d3s22
                      call enough('hcddjk.8.1',bc,ibc)
                      call prodn(iwpk1b,iwpb2,ncsfmid1b,ncsfmid2,         11d3s22
     $                  ncsf(karg),bc(iprod),bc,ibc,1d0,0d0)            9d28s23
                      itrans=ibcoff
                      ibcoff=itrans+ncsfmid1b*ncsfmid2
                      call enough('hcddjk.8.3',bc,ibc)
                      do i=0,ncsfmid2-1
                       do j=0,ncsfmid1b-1
                        ji=iprod+j+ncsfmid1b*i
                        ij=itrans+i+ncsfmid2*j
                        bc(ij)=bc(ji)
                       end do
                      end do
                      iwprod=ibcoff                                      11d3s22
                      icmp1=iwprod+1                                     11d3s22
                      icmp2=icmp1+ncsfmid2                               11d3s22
                      icmp3=icmp2+nsz                                    11d3s22
                      ibcoff=icmp3+nsz                                   11d3s22
                      call enough('hcddjk.8.2',bc,ibc)
                      call cmpvcsf(ncsfmid2,ncsfmid1b,ibc(icmp1),        11d3s22
     $                    ibc(icmp2),bc(icmp3),bc(itrans),nkeep)         11d3s22
                      nkeeph=nkeep/2                                     11d3s22
                      if(2*nkeeph.ne.nkeep)nkeeph=nkeeph+1               11d3s22
                      if(ncsfmid1b+nkeeph+nkeep.lt.nsz)then              11d3s22
                       icmp3n=icmp2+nkeeph                               11d3s22
                       do i=0,nkeep-1                                    11d3s22
                        bc(icmp3n+i)=bc(icmp3+i)                         11d3s22
                       end do                                            11d3s22
                       ibc(iwprod)=nkeep                                 11d3s22
                       iwprod=-iwprod                                    11d3s22
                      else                                               11d3s22
                       iwprod=iprod                                      11d3s22
                      end if                                             11d3s22
                      icsf0=1                                             11d3s22
                      jcsf0=1                                             11d3s22
                      do l=1,4                                            11d3s22
                       if(min(nlj(l),nli(l)).gt.0)then                    11d3s22
                        jad1=jvcv+ibc(jvcv+5+l)                           3d19s21
                        jad2=jad1+nlj(l)                                  3d19s21
                        iad1=ivcv+ibc(ivcv+5+l)                           3d19s21
                        iad2=iad1+nli(l)                                  3d19s21
                        itmp=ibcoff                                       11d3s22
                        ibcoff=itmp+nli(l)*nlj(l)                         11d3s22
                        call enough('hcddjk.8.2',bc,ibc)
c     nli,ncsf(iarg),ncsfmid1b,ncsf(karg),ncsfmid2,ncsf(jarg),nlj
c       ivects     iwpb1b    iprod     iwpk2    jad2
c      nli,ncsf(iarg),ncsfmid1b,ncsfmid2,ncsf(jarg),nlj
c       1    2           3         4         5       6
                        call genmatn3(nli(l),ncsfmid1b,ncsf(jarg),        11d3s22
     $                     ncsf(iarg),ivects(l),iwpb1b,ncsfmid2,        11d3s22
     $                     iwprod,iwpk2,bc(jad2),ncsf2(l,jarg),nlj(l),  11d3s22
     $                     bc(itmp),0d0,jcsf0,ncsf2(l,jarg),            11d3s22
     $                     icsf0,ncsf2(l,iarg),xint,bc,ibc)             9d28s23
                        do j=1,nlj(l)
                         jm=j-1                                          11d4s22
                         jj=ibc(jad1+jm)-1                               11d4s22
                         jjden=iden1en(l)+nfdat(2,l,isb)*jj-1               2d24s21
                         do iii=1,nli(l)                                 11d4s22
                          iiim=iii-1                                     11d4s22
                          iad=itmp+iiim+nli(l)*jm                        11d4s22
                          kk=ibc(iad1+iiim)                              11d4s22
                          bc(jjden+kk)=bc(jjden+kk)+bc(iad)              11d4s22
                         end do
                        end do
                        ibcoff=itmp                                       11d3s22
                       end if                                             11d3s22
                       icsf0=icsf0+ncsf2(l,iarg)                          11d3s22
                       jcsf0=jcsf0+ncsf2(l,jarg)                          11d3s22
                      end do                                              11d3s22
                      ibcoff=iprod                                       11d3s22
                      if(ipss.eq.2)go to 3                                   12d18s20
                     end if                                              11d3s22
                    else                                                 11d4s22
                     call gandc(itestc,itesto,itc,ito,nopenk,nopenip,      2d25s21
     $         karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,     2d25s21
     $                  iwpk1b,ncsfmid1b,bc,ibc)                        9d28s23
                     nnot1=nnot1b                                          2d25s21
                     nab1(1)=nab1b(2)                                      2d25s21
                     nab1(2)=nab1b(1)                                      2d25s21
                     if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                      lsa=ism(nab1(1))                                     1d8s21
                      lga=irel(nab1(1))-1                                  1d8s21
                      lsb=ism(nab1(2))                                     1d8s21
                      lgb=irel(nab1(2))-1                                  1d8s21
                      lsc=ism(nab2(1))                                     1d8s21
                      lgc=irel(nab2(1))-1                                  1d8s21
                      lsd=ism(nab2(2))                                     1d8s21
                      lgd=irel(nab2(2))-1                                  1d8s21
                      xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  9d28s23
     $                   bc,ibc)                                        9d28s23
                      if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2))     4d20s21
     $                   .and.                                          4d20s21
     $                max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))     1d12s21
     $                  xint=xint*0.5d0                                 1d12s21
                      jargo=1                                              2d23s21
                      iargo=1                                              2d24s21
                      do l=1,4                                               12d19s20
                       if(min(nli(l),nlj(l)).gt.0)then                    3d1s21
                        jad1=jvcv+ibc(jvcv+5+l)                           3d19s21
                        jad2=jad1+nlj(l)                                  3d19s21
                        iad1=ivcv+ibc(ivcv+5+l)                           3d19s21
                        iad2=iad1+nli(l)                                  3d19s21
                        itmpj=ibcoff                                      2d24s21
                        ibcoff=itmpj+ncsf(karg)*nlj(l)                    2d24s21
                        call enough('hcddjk.  9',bc,ibc)
                        call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,bc(jad2),ncsf2(l,jarg),     2d24s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                        itmpi=ibcoff                                      2d24s21
                        itmpt=itmpi+ncsf(karg)*nli(l)                     2d24s21
                        ibcoff=itmpt+ncsf(karg)*nli(l)                    2d24s21
                        call enough('hcddjk. 10',bc,ibc)
                        call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,bc(iad2),ncsf2(l,iarg),   2d24s21
     $                   bc(itmpt),ncsf(karg),xint,0d0,iargo,           11d2s22
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                        do iii=0,nlj(l)-1                                 11d2s22
                         iadj=itmpj+ncsf(karg)*iii                        11d2s22
                         jj=ibc(jad1+iii)-1                                 2d24s21
                         jjden=iden1en(l)+nfdat(2,l,isb)*jj-1               2d24s21
                         do j=0,nli(l)-1                                  11d2s22
                          iadi=itmpt+ncsf(karg)*j                         11d2s22
                          sum=0d0                                         11d2s22
                          do k=0,ncsf(karg)-1                             11d2s22
                           sum=sum+bc(iadj+k)*bc(iadi+k)                  11d2s22
                          end do                                          11d2s22
                          kk=ibc(iad1+j)                                    2d24s21
                          bc(jjden+kk)=bc(jjden+kk)+sum                   11d2s22
                         end do                                           11d2s22
                        end do                                            11d2s22
                        ibcoff=itmpj                                      2d24s21
                       end if                                             2d25s21
                       jargo=jargo+ncsf2(l,jarg)                           2d23s21
                       iargo=iargo+ncsf2(l,iarg)                           2d24s21
                      end do                                                 12d19s20
                      if(ipss.eq.2)go to 3                                   12d18s20
                     end if                                              11d4s22
                    end if                                                  12d18s20
                   end if                                                4d20s21
                  end do                                                   12d18s20
    3             continue                                                 12d18s20
                 end if                                                  1d8s21
                end if                                                   10d14s22
               end if                                                    1d8s21
               jto=ibclr(jto,norbx)                                      1d8s21
               jto=ibset(jto,norbxxx)                                    1d8s21
c                                                                       1d8s21
               gandcc=ieor(itc,jtc)                                      10d14s22
               gandco=ieor(ito,jto)                                      10d14s22
               gandcb=ior(gandcc,gandco)                                 10d21s22
               ndifb=popcnt(gandcb)                                      10d21s22
               if(ndifb.le.4)then                                        10d21s22
                ndifd=popcnt(gandcc)                                      10d14s22
                ndifs=popcnt(gandco)                                      10d14s22
                if(ndifs.eq.2.and.ndifb.eq.2)then                        10d21s22
                 do i=1,norbxxx                                           10d21s22
                  if(btest(gandco,i))then                                          10d14s22
                   if((btest(ito,i).and..not.btest(jtc,i)).or.
     $                 (btest(itc,i).and.btest(jto,i)))then                           10d14s22
                    nab4(1,1)=i
                   else                                                           10d14s22
                    nab4(2,1)=i
                   end if
                  end if                                                          10d14s22
                 end do                                                           10d14s22
                 if(isb.ne.jsb)then                                      1d9s21
                  write(6,*)('hey, what''s up with isb ne jsb ? ')       1d9s21
                  stop 'hcddjk:g'                                                   1d9s21
                 end if                                                  1d9s21
                 call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,   1d8s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $              ncsfmid1,bc,ibc)                                    9d28s23
                call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,   2d23s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,      1d8s21
     $              ncsfmid1b,bc,ibc)                                   9d28s23
                itmpn=ibcoff                                            2d23s21
                nlij=0                                                  2d23s21
                do l=1,4                                                2d23s21
                 nlij=nlij+nli(l)*nlj(l)                                2d23s21
                end do                                                  2d23s21
                itmpn=ibcoff                                            2d23s21
                ibcoff=itmpn+nlij                                       2d23s21
                call enough('hcddjk. 12',bc,ibc)
                jtmpn=itmpn                                             2d23s21
                jargo=1                                                  2d22s21
                iargo=1                                                 2d23s21
                do l=1,4                                                 1d8s21
                 if(min(nlj(l),nli(l)).gt.0)then                         2d23s21
                  iad1=ivcv+ibc(ivcv+5+l)                               3d19s21
                  iad2=iad1+nli(l)                                      3d19s21
                  jad1=jvcv+ibc(jvcv+5+l)                               3d19s21
                  jad2=jad1+nlj(l)                                      3d19s21
                  itmp1=ibcoff                                          2d23s21
                  itmp1b=itmp1+ncsfmid1*nlj(l)                          2d23s21
                  ibcoff=itmp1b+ncsfmid1*nli(l)                         2d23s21
                  call enough('hcddjk. 13',bc,ibc)
                  if(iwpk1.lt.0)then                                    2d23s21
                   nusedi=ibc(-iwpk1)/2                                 2d23s21
                   if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1           2d23s21
                   icmp1=-iwpk1+1                                       2d23s21
                   icmp2=icmp1+ncsf(jarg)                               2d23s21
                   icmp3=icmp2+nusedi                                   2d23s21
                   call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),    2d23s21
     $                   bc(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),   2d23s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                  else                                                              12d5s20
                   jwpk1=iwpk1+ncsfmid1*(jargo-1)                       2d23s21
                   call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                   1d0,bc(jwpk1),ncsfmid1,bc(jad2),ncsf2(l,jarg),
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjk.  8')
                  end if                                                            12d5s20
                  if(iwpk1b.lt.0)then                                   2d23s21
                   nusedi=ibc(-iwpk1b)/2                                2d23s21
                   if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                   icmp1=-iwpk1b+1                                      2d23s21
                   icmp2=icmp1+ncsf(iarg)                               2d23s21
                   icmp3=icmp2+nusedi                                   2d23s21
                   call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),    2d23s21
     $                   bc(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),   2d23s21
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                  else                                                              12d5s20
                   jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                   call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                   1d0,bc(jwpk1b),ncsfmid1,bc(iad2),ncsf2(l,iarg),
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjk.  9')
                  end if                                                            12d5s20
                  itmpt=ibcoff                                          2d23s21
                  ibcoff=itmpt+nli(l)*ncsfmid1                          2d23s21
                  call enough('hcddjk. 14',bc,ibc)
                  do iii=0,nli(l)-1                                     2d23s21
                   do j=0,ncsfmid1-1                                    2d23s21
                    ji=itmp1b+j+ncsfmid1*iii                            2d23s21
                    ij=itmpt+iii+nli(l)*j                               2d23s21
                    bc(ij)=bc(ji)                                       2d23s21
                   end do                                               2d23s21
                  end do                                                2d23s21
                  call dgemm('n','n',nli(l),nlj(l),ncsfmid1,1d0,        2d23s21
     $                  bc(itmpt),nli(l),bc(itmp1),ncsfmid1,0d0,         2d23s21
     $                  bc(jtmpn),nli(l),                                2d23s21
     d' hcddjk. 10')
                  ibcoff=itmp1                                          2d23s21
                  ktmpn=jtmpn                                           2d24s21
                  do iii=0,nlj(l)-1                                     2d24s21
                   jj=ibc(jad1+iii)-1                                   2d24s21
                   jjden=idenhvvn(l)+nfdat(2,l,isb)*jj-1                 2d24s21
                   do j=0,nli(l)-1                                      2d24s21
                    kk=ibc(iad1+j)                                      2d24s21
                    bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)               2d24s21
                   end do                                               2d24s21
                   ktmpn=ktmpn+nli(l)                                   2d24s21
                  end do                                                2d24s21
                  jtmpn=jtmpn+nli(l)*nlj(l)                             2d23s21
                 end if                                                 2d23s21
                 jargo=jargo+ncsf2(l,jarg)                              2d22s21
                 iargo=iargo+ncsf2(l,iarg)                              2d24s21
                end do                                                  1d8s21
                nok=0                                                    1d8s21
                do i=1,norb                                             1d8s21
                 if(itest(i,1).gt.0)then                                1d8s21
                  nok=nok+1                                             1d8s21
                  itest(nok,3)=itest(i,1)                               1d8s21
                  itest(nok,2)=i                                        1d8s21
                 end if                                                 1d8s21
                end do                                                  1d8s21
                do i=1,nok                                              1d8s21
                 lsa=ism(itest(i,2))                                    1d8s21
                 lga=irel(itest(i,2))-1                                 1d8s21
                 icolj=((lga*(lga+1))/2)+lga                            1d8s21
                 icolk=lga+irefo(lsa)*lga                               1d8s21
                 if(itest(i,3).eq.2)then                                1d8s21
                  jtmpn=itmpn                                           2d23s21
                  do l=1,4                                              1d8s21
                   if(min(nli(l),nlj(l)).gt.0)then                      2d23s21
                    jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                    iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                    jden=idenjn(l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)    2d23s21
     $                    *icolj                                        2d23s21
                    kden=idenkn(l,l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)  2d24s21
     $                    *icolk                                        2d24s21
                    do iii=0,nlj(l)-1                                   2d23s21
                     jj=ibc(jad1+iii)-1                                 1d8s21
                     jjden=jden+nfdat(2,l,isb)*jj-1                     2d23s21
                     kkden=kden+nfdat(2,l,isb)*jj-1                     2d24s21
                     do j=0,nli(l)-1                                    2d23s21
                      kk=ibc(iad1+j)                                    2d23s21
                      bc(jjden+kk)=bc(jjden+kk)+bc(jtmpn+j)*2d0         2d23s21
                      bc(kkden+kk)=bc(kkden+kk)-bc(jtmpn+j)             2d24s21
                     end do                                             2d23s21
                     jtmpn=jtmpn+nli(l)                                 2d23s21
                    end do                                              2d23s21
                   end if                                               2d23s21
                  end do                                                1d8s21
                 else                                                   1d8s21
                  jtmpn=itmpn                                           2d23s21
                  do l=1,4                                              1d8s21
                   if(nlj(l).gt.0)then                                  1d8s21
                    jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                    iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                    jden=idenjn(l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)    2d23s21
     $                    *icolj                                        2d23s21
                    do iii=0,nlj(l)-1                                   2d23s21
                     jj=ibc(jad1+iii)-1                                 1d8s21
                     jjden=jden+nfdat(2,l,isb)*jj-1                     2d23s21
                     do j=0,nli(l)-1                                    2d23s21
                      kk=ibc(iad1+j)                                    2d23s21
                      bc(jjden+kk)=bc(jjden+kk)+bc(jtmpn+j)             2d23s21
                     end do                                             2d23s21
                     jtmpn=jtmpn+nli(l)                                 2d23s21
                    end do                                              2d23s21
                   end if                                               1d8s21
                  end do                                                1d8s21
                 end if                                                 1d8s21
                end do                                                  1d8s21
                ibcoff=itmpn                                            4d28s21
                do i=1,nok                                              1d8s21
                 if(itest(i,3).eq.1)then                                  12d14s20
                  itestc=itc                                              12d14s20
                  itesto=ito                                              12d14s20
                  nopenk=nopenip                                                11d13s20
c
c     i is bra, j is ket.
c
c         I   J    K     thus  IK  KJ
c     b   n+1 n    n+1        g   g   (ck|bc)
c     k   m   m+1  m+1         ck  bc
c     c   p   p    p-1
c
c
c     anihilate common
c
                  if(btest(itestc,itest(i,2)))then                             11d13s20
                   itestc=ibclr(itestc,itest(i,2))                             11d13s20
                   itesto=ibset(itesto,itest(i,2))                             11d13s20
                   karg=iarg-1                                                11d13s20
                   nopenk=nopenk+1                                             11d13s20
                  else                                                         11d13s20
                   itesto=ibclr(itesto,itest(i,2))                             11d13s20
                   karg=iarg                                                  11d13s20
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
                  nnot1=0                                               8d4s22
                  nnot2=0                                               8d4s22
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                   call gandc(itestc,itesto,itc,ito,nopenk,nopenip,     2d24s21
     $            karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b, 2d24s21
     $               iwpk1b,ncsfmid1b,bc,ibc)                           9d28s23
                   nnot1=nnot1b                                         2d25s21
                   nab1(1)=nab1b(2)                                     2d25s21
                   nab1(2)=nab1b(1)                                     2d25s21
                   call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,        12d14s20
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            9d28s23
                  end if
                  if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                   if(nab1(2).gt.nab2(1))then                           1d8s21
                    ia=nab1(1)
                    ib=nab2(2)
                   else
                    ia=nab2(2)
                    ib=nab1(1)
                   end if                                               1d8s21
                   lsa=ism(ia)                                          1d8s21
                   lga=irel(ia)-1                                       1d8s21
                   lsb=ism(ib)                                          1d8s21
                   lgb=irel(ib)-1                                       1d8s21
                   lsab=multh(lsa,lsb)                                  1d9s21
                   icol=lga+irefo(lsa)*lgb                              1d8s21
                   if(lsab.ne.ijsb)then                                 1d9s21
                    write(6,*)('hey, lsab = '),lsab,(' ne ijsb = '),
     $                  ijsb
                    stop 'hcddjk:i'
                   end if                                               1d9s21
                   jargo=1                                               2d23s21
                   itmpj=ibcoff                                          2d23s21
                   jtmpj=itmpj                                           2d23s21
                   do lj=1,4                                             2d23s21
                    if(nlj(lj).gt.0)then                                 2d23s21
                     jad1=jvcv+ibc(jvcv+5+lj)                              3d19s21
                     jad2=jad1+nlj(lj)                                     3d19s21
                     ibcoff=jtmpj+ncsf(karg)*nlj(lj)                     2d23s21
                     call enough('hcddjk. 15',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                    nlj(lj),iwpb2,iwpk2,bc(jad2),ncsf2(lj,jarg),   2d23s21
     $                   bc(jtmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(lj,jarg),bc,ibc)                         9d28s23
                     jtmpj=jtmpj+ncsf(karg)*nlj(lj)                      2d23s21
                    end if                                               2d23s21
                    jargo=jargo+ncsf2(lj,jarg)                           2d23s21
                   end do                                                2d23s21
                   iargo=1                                               2d23s21
                   itmpi=ibcoff                                          2d23s21
                   jtmpi=itmpi                                           2d23s21
                   do li=1,4                                             2d23s21
                    if(nli(li).gt.0)then                                 2d23s21
                     iad1=ivcv+ibc(ivcv+5+li)                              3d19s21
                     iad2=iad1+nli(li)                                     3d19s21
                     ibcoff=jtmpi+ncsf(karg)*nli(li)                     2d23s21
                     itmpt=ibcoff                                        2d23s21
                     ibcoff=itmpt+ncsf(karg)*nli(li)                     2d23s21
                     call enough('hcddjk. 16',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                    nli(li),iwpb1b,iwpk1b,bc(iad2),ncsf2(li,iarg), 2d23s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(li,iarg),bc,ibc)                         9d28s23
                     do iii=0,nli(li)-1                                  2d23s21
                      do k=0,ncsf(karg)-1                                2d23s21
                       ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                       ik=jtmpi+iii+nli(li)*k                            2d23s21
                       bc(ik)=bc(ki)                                     2d23s21
                      end do                                             2d23s21
                     end do                                              2d23s21
                     ibcoff=itmpt                                        2d23s21
                     jtmpi=jtmpi+ncsf(karg)*nli(li)                      2d23s21
                    end if                                               2d23s21
                    iargo=iargo+ncsf2(li,iarg)                            2d23s21
                   end do                                                2d23s21
                   jtmpj=itmpj                                          2d23s21
                   do lj=1,4                                             2d23s21
                    if(nlj(lj).gt.0)then                                2d24s21
                     jad1=jvcv+ibc(jvcv+5+lj)                              3d19s21
                     jtmpi=itmpi                                          2d23s21
                     do li=1,4                                          2d24s21
                      if(nli(li).gt.0)then                              2d24s21
                       iad1=ivcv+ibc(ivcv+5+li)                         3d19s21
                       jden=idenkn(li,lj,lsa)+nfdat(2,li,isb)           2d24s21
     $                      *nfdat(2,lj,jsb)*icol                       2d24s21
                       itmpm=ibcoff                                       2d23s21
                       ibcoff=itmpm+nlj(lj)*nli(li)                     2d24s21
                       call enough('hcddjk. 17',bc,ibc)
                       call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),   2d24s21
     $                       1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),2d24s21
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' hcddjk. 11')
                       jtmpm=itmpm                                        2d23s21
                       do iii=0,nlj(lj)-1                               2d24s21
                        jj=ibc(jad1+iii)-1                                 1d8s21
                        jjden=jden+nfdat(2,li,isb)*jj-1                 2d24s21
                        do j=0,nli(li)-1                                2d24s21
                         kk=ibc(iad1+j)                                   2d23s21
                         bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)            2d23s21
                        end do                                            2d23s21
                        jtmpm=jtmpm+nli(li)                                2d23s21
                       end do                                             2d23s21
                       ibcoff=itmpm                                       2d23s21
                      end if                                              2d23s21
                      jtmpi=jtmpi+nli(li)*ncsf(karg)                    2d24s21
                     end do                                             2d24s21
                    end if                                              2d24s21
                    jtmpj=jtmpj+nlj(lj)*ncsf(karg)                      2d24s21
                   end do                                               2d23s21
                  end if                                                  12d14s20
                 end if                                                   12d14s20
                end do                                                  12d18s20
               else                                                     1d8s21
                ipssx=0                                                 10d21s22
                nnot=0                                                  10d21s22
                if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s22
                 nnot=4                                                          10d14s22
                 ioxx(1)=1                                                     10d17s22
                 ioxx(2)=1                                                     10d17s22
                 do i=1,norbxxx                                                     10d17s22
                  if(btest(gandcb,i))then                                         10d14s22
                   if((btest(jtc,i).and.btest(ito,i)).or.                         10d17s22
     $                  (btest(jto,i).and..not.btest(itc,i)))then                   10d14s22
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
                 do i=1,norbxxx
                  if(btest(gandcb,i))then                                         10d14s22
                   if(btest(gandcc,i).and.                                        10d14s22
     $        ((btest(itc,i).and..not.btest(jto,i)).or.                 10d14s22
     $         (btest(jtc,i).and..not.btest(ito,i))))then                     10d14s22
                    if(btest(jtc,i))iswap=1                                        10d17s22
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
                  if(btest(itc,nab4(1,2)).and.
     $                  .not.btest(itc,nab4(1,1)))nbt=1                 10d21s22
                 else                                                             10d17s22
                  nbt=0                                                           10d17s22
                  if(btest(jtc,nab4(2,2)).and.
     $                  .not.btest(jtc,nab4(2,1)))nbt=1                 10d21s22
                 end if                                                           10d17s22
                 if(nbt.ne.0)then                                                 10d17s22
                  nab4(1,1)=nab4(1,2)                                             10d17s22
                  nab4(2,1)=nab4(2,2)                                             10d17s22
                 end if                                                           10d17s22
                else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                 nnot=3
                 do i=1,norbxxx
                  if(btest(gandcb,i))then                                         10d14s22
                   if(btest(itc,i))then
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
                else if(nnot.eq.4)then                                  11d1s22
                 ipssx=3                                                 12d18s20
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
                 itestc=itc                                              12d8s20
                 itesto=ito                                              12d8s20
                 if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                  itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                  itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopenip+1                                              11d13s20
                  karg=iarg-1                                             12d8s20
                 else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                  itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopenip-1                                              11d13s20
                  karg=iarg
                 else                                                          11d13s20
                  write(6,*)('bit not set for nab4(1,1) ='),nab4(1,iu1)       11d27s20
                  stop 'nab4(1,1)'                                           11d27s20
                 end if                                                        11d13s20
                 if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                  itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                  itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                  nopenk=nopenk-1                                              11d13s20
                  karg=karg+1                                                  11d13s20
                 else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                  write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                 nab4(2,iu2)                                       12d18s20
                  stop 'nab4(2,1)'                                           11d27s20
                 else                                                          11d13s20
                  itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                  nopenk=nopenk+1                                              11d13s20
                 end if                                                        11d13s20
                 nqq=karg+mdon-1                                        4d20s21
                 if(nqq.ge.mdon.and.nqq.le.mdoo)then                    4d20s21
                  call gandc(itestc,itesto,itc,ito,nopenk,nopenip,       2d23s21
     $         karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,    2d23s21
     $                iwpk1b,ncsfmid1b,bc,ibc)                          9d28s23
                  nnot1=nnot1b                                           2d25s21
                  nab1(1)=nab1b(2)                                       2d25s21
                  nab1(2)=nab1b(1)                                       2d25s21
                  call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,       1d8s21
     $         karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d18s20
     $         ncsfmid2,bc,ibc)                                         9d28s23
                  if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                   if(nab1(2).gt.norb)then                               1d8s21
                    if(nab2(1).lt.nab2(1))then                           1d8s21
                     ia=nab2(2)
                     ib=nab1(1)                                          1d8s21
                    else                                                 1d8s21
                     ia=nab1(1)                                          1d8s21
                     ib=nab2(2)                                          1d8s21
                    end if                                               1d8s21
                    lsa=ism(ia)                                           1d8s21
                    lsb=ism(ib)                                           1d8s21
                    lga=irel(ia)-1                                        1d8s21
                    lgb=irel(ib)-1                                        1d8s21
                    icol=lga+irefo(lsa)*lgb                              1d8s21
                   else                                                  1d8s21
                    ia=min(nab1(1),nab1(2))                              1d8s21
                    ib=max(nab1(1),nab1(2))                              1d8s21
                    lsa=ism(ia)                                           1d8s21
                    lsb=ism(ib)                                           1d8s21
                    lga=irel(ia)-1                                        1d8s21
                    lgb=irel(ib)-1                                        1d8s21
                    if(lsa.eq.lsb)then                                   1d8s21
                     icol=((lgb*(lgb+1))/2)+lga                          1d10s21
                    else                                                 1d8s21
                     icol=lga+irefo(lsa)*lgb                             1d9s21
                    end if                                               1d8s21
                   end if                                                1d8s21
                   lsab=multh(lsa,lsb)                                   1d9s21
                   jargo=1                                               2d23s21
                   itmpj=ibcoff                                          2d23s21
                   jtmpj=itmpj                                           2d23s21
                   do lj=1,4                                             2d23s21
                    if(nlj(lj).gt.0)then                                 2d23s21
                     jad1=jvcv+ibc(jvcv+5+lj)                           3d19s21
                     jad2=jad1+nlj(lj)                                     3d19s21
                     ibcoff=jtmpj+ncsf(karg)*nlj(lj)                     2d23s21
                     call enough('hcddjk. 18',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(lj),iwpb2,iwpk2,bc(jad2),ncsf2(lj,jarg),   2d23s21
     $                   bc(jtmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(lj,jarg),bc,ibc)                         9d28s23
                     jtmpj=jtmpj+ncsf(karg)*nlj(lj)                      2d23s21
                    end if                                               2d23s21
                    jargo=jargo+ncsf2(lj,jarg)                           2d23s21
                   end do                                                2d23s21
                   iargo=1                                               2d23s21
                   itmpi=ibcoff                                          2d23s21
                   jtmpi=itmpi                                           2d23s21
                   do li=1,4                                             2d23s21
                    if(nli(li).gt.0)then                                 2d23s21
                     iad1=ivcv+ibc(ivcv+5+li)                              3d19s21
                     iad2=iad1+nli(li)                                     3d19s21
                     ibcoff=jtmpi+ncsf(karg)*nli(li)                     2d23s21
                     itmpt=ibcoff                                        2d23s21
                     ibcoff=itmpt+ncsf(karg)*nli(li)                     2d23s21
                     call enough('hcddjk. 19',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(li),iwpb1b,iwpk1b,bc(iad2),ncsf2(li,iarg), 2d23s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(li,iarg),bc,ibc)                         9d28s23
                     do iii=0,nli(li)-1                                  2d23s21
                      do k=0,ncsf(karg)-1                                2d23s21
                       ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                       ik=jtmpi+iii+nli(li)*k                            2d23s21
                       bc(ik)=bc(ki)                                     2d23s21
                      end do                                             2d23s21
                     end do                                              2d23s21
                     ibcoff=itmpt                                        2d23s21
                     jtmpi=jtmpi+ncsf(karg)*nli(li)                      2d23s21
                    end if                                               2d23s21
                    iargo=iargo+ncsf2(li,iarg)                            2d23s21
                   end do                                                2d23s21
                   if(nab1(2).gt.norb)then                               2d23s21
                    jtmpj=itmpj                                          2d23s21
                    do lj=1,4                                           2d24s21
                     if(nlj(lj).gt.0)then                               2d24s21
                      jad1=jvcv+ibc(jvcv+5+lj)                              3d19s21
                      jtmpi=itmpi                                       2d24s21
                      do li=1,4                                         2d24s21
                       if(nli(li).gt.0)then                             2d24s21
                        iad1=ivcv+ibc(ivcv+5+li)                        3d19s21
                        jden=idenkn(li,lj,lsa)+nfdat(2,li,isb)          2d24s21
     $                       *nfdat(2,lj,jsb)*icol                      2d24s21
                        itmpm=ibcoff                                       2d23s21
                        ibcoff=itmpm+nlj(lj)*nli(li)                    2d24s21
                        call enough('hcddjk. 20',bc,ibc)
                        call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),  2d24s21
     $                       1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),2d24s21
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' hcddjk. 12')
                        jtmpm=itmpm                                        2d23s21
                        do iii=0,nlj(lj)-1                                  2d23s21
                         jj=ibc(jad1+iii)-1                                 1d8s21
                         jjden=jden+nfdat(2,li,isb)*jj-1                2d24s21
                         do j=0,nli(li)-1                               2d24s21
                          kk=ibc(iad1+j)                                   2d23s21
                          bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)            2d23s21
                         end do                                            2d23s21
                         jtmpm=jtmpm+nli(li)                            2d24s21
                        end do                                             2d23s21
                        ibcoff=itmpm                                       2d23s21
                       end if                                              2d23s21
                       jtmpi=jtmpi+nli(li)*ncsf(karg)                   2d24s21
                      end do                                            2d24s21
                     end if                                             2d24s21
                     jtmpj=jtmpj+nlj(lj)*ncsf(karg)                       2d23s21
                    end do                                               2d23s21
                   else                                                  2d23s21
                    jtmpj=itmpj                                          2d23s21
                    jtmpi=itmpi                                          2d23s21
                    do l=1,4                                             2d23s21
                     if(min(nli(l),nlj(l)).gt.0)then                     2d23s21
                      jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                      iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                      jden=idenjn(l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)   2d23s21
     $                    *icol                                         2d23s21
                      itmpm=ibcoff                                       2d23s21
                      ibcoff=itmpm+nlj(l)*nli(l)                         2d23s21
                      call enough('hcddjk. 21',bc,ibc)
                      call dgemm('n','n',nli(l),nlj(l),ncsf(karg),1d0,   2d23s21
     $                    bc(jtmpi),nli(l),bc(jtmpj),ncsf(karg),0d0,    2d23s21
     $                    bc(itmpm),nli(l),                             2d23s21
     d' hcddjk. 13')
                      jtmpm=itmpm                                        2d23s21
                      do iii=0,nlj(l)-1                                  2d23s21
                       jj=ibc(jad1+iii)-1                                 1d8s21
                       jjden=jden+nfdat(2,l,isb)*jj-1                    2d23s21
                       do j=0,nli(l)-1                                   2d23s21
                        kk=ibc(iad1+j)                                   2d23s21
                        bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)            2d23s21
                       end do                                            2d23s21
                       jtmpm=jtmpm+nli(l)                                2d23s21
                      end do                                             2d23s21
                      ibcoff=itmpm                                       2d23s21
                     end if                                              2d23s21
                     jtmpj=jtmpj+nlj(l)*ncsf(karg)                       2d23s21
                     jtmpi=jtmpi+nli(l)*ncsf(karg)                       2d23s21
                    end do                                               2d23s21
                   end if                                                2d23s21
                   if(ipss.eq.2)go to 4                                   12d18s20
                  end if                                                  12d18s20
                 end if                                                 4d20s21
                end do                                                   12d18s20
    4           continue                                                 12d18s20
               end if                                                   1d8s21
              end if                                                    10d14s22
              jvcv=jvcv+nspacej                                         3d19s21
             end do                                                       1d8s21
             end if                                                     3d1s21
             idoit=idoit+1                                               3d1s21
            end if                                                      1d8s21
           end do                                                       1d8s21
           ibcoff=ibcvects                                              11d3s22
           ivcv=ivcv+nspacei                                            3d19s21
          end do                                                        1d8s21
          ibcoff=ibcst                                                  1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call dws_gsumf(bc(ibc0),nwds)                                   2d23s21
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    1d8s21
          do ll=1,4                                                     1d8s21
           if(nfdat(2,ll,jsb).gt.0)then                                 1d8s21
            nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                          1d8s21
            if(l.eq.ll.and.nll.gt.0)then                                2d24s21
             iden1ef(l)=iden1en(l)                                      2d25s21
             idenhvvf(l)=idenhvvn(l)                                    2d25s21
            end if                                                      1d11s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(lsb.ge.lsa.and.l.eq.ll.and.                             2d23s21
     $            min(irefo(lsa),irefo(lsb),nll).gt.0)then              2d23s21
              idenjf(l,lsa)=idenjn(l,lsa)                               2d25s21
             end if                                                     1d9s21
             nn=irefo(lsa)*irefo(lsb)                                   1d9s21
             if(min(nn,nll).gt.0)then                                   2d24s21
              idenkf(l,ll,lsa)=idenkn(l,ll,lsa)                         2d25s21
             end if                                                     2d24s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        nall=0                                                          5d10s21
        do ll=1,4                                                       1d8s21
         if(nfdat(2,ll,jsb).gt.0)then                                   1d8s21
          do l=1,4                                                      1d8s21
           if(nfdat(2,l,isb).gt.0)then                                  1d8s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             nokk(l,ll,lsa)=0                                           1d15s21
             if(min(irefo(lsa),irefo(lsb)).gt.0)then                    1d13s21
              nnn=nfdat(2,l,isb)*nfdat(2,ll,jsb)                        1d9s21
              if(lsb.ge.lsa.and.l.eq.ll)then                            1d13s21
               if(lsa.eq.lsb)then                                        1d8s21
                nn=(irefo(lsa)*(irefo(lsa)+1))/2                         1d8s21
               else                                                      1d8s21
                nn=irefo(lsa)*irefo(lsb)                                 1d8s21
               end if                                                    1d8s21
               jden=idenjf(l,lsa)                                        1d11s21
               jdent=idenjf(l,lsa)                                       1d11s21
               nokj(l,lsa)=0                                             1d12s21
               kden=idenjf(l,lsa)                                        1d12s21
               do icol=0,nn-1                                            1d9s21
                sz=0d0                                                    1d8s21
                do i=0,nnn-1                                             1d9s21
                 sz=sz+bc(jden+i)**2                                    1d13s21
                end do                                                    1d8s21
                sz=sqrt(sz/dfloat(nnn))                                  1d9s21
                if(sz.gt.1d-10)then                                       1d8s21
                 ibc(ndenjf(l,lsa)+nokj(l,lsa))=icol                     1d12s21
                 nokj(l,lsa)=nokj(l,lsa)+1                               1d12s21
                 do i=0,nnn-1                                            1d12s21
                  bc(kden+i)=bc(jden+i)                                  1d12s21
                 end do                                                  1d12s21
                 kden=kden+nnn                                           1d12s21
                 nall=nall+1                                            5d10s21
                end if                                                   1d13s21
                jden=jden+nnn                                            1d9s21
               end do                                                    1d9s21
              end if                                                    1d8s21
              nn=irefo(lsa)*irefo(lsb)                                  1d8s21
              jden=idenkf(l,ll,lsa)                                     1d9s21
              jdent=idenkf(ll,l,lsb)                                    1d10s21
              kden=idenkf(l,ll,lsa)                                     1d12s21
              do icol=0,nn-1                                            1d9s21
               ib=icol/irefo(lsa)
               ia=icol-irefo(lsa)*ib
               icolt=ib+irefo(lsb)*ia                                   1d10s21
               jdent=idenkf(ll,l,lsb)+nnn*icolt                         1d10s21
               sz=0d0                                                    1d8s21
               do i=0,nnn-1                                             1d9s21
                sz=sz+bc(jden+i)**2                                     1d9s21
               end do                                                    1d8s21
               sz=sqrt(sz/dfloat(nnn))                                  1d9s21
               if(sz.gt.1d-10)then                                       1d8s21
                ibc(ndenkf(l,ll,lsa)+nokk(l,ll,lsa))=icol               1d12s21
                nokk(l,ll,lsa)=nokk(l,ll,lsa)+1                         1d12s21
                do i=0,nnn-1                                            1d12s21
                 bc(kden+i)=bc(jden+i)                                  1d12s21
                end do                                                  1d12s21
                kden=kden+nnn                                           1d12s21
                nall=nall+1                                             5d10s21
               end if                                                    1d8s21
               jden=jden+nnn                                            1d9s21
              end do                                                    1d9s21
             end if                                                     1d8s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        if(nall.gt.0)then                                               5d10s21
         ioff=ioffdnon                                                   1d12s21
         do isbv1=1,nsymb                                                1d12s21
          isbv2=multh(isbv1,isbv12)                                      1d12s21
          if(isbv2.ge.isbv1)then                                         1d12s21
           if(isbv12.eq.1)then                                           1d12s21
            isw=0                                                        1d12s21
            mvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        1d12s21
            ioffp=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)
           else                                                          1d12s21
            isw=1                                                        1d12s21
            mvv=nvirt(isbv1)*nvirt(isbv2)                                1d12s21
            ioffp=ioff
           end if                                                        1d12s21
           mmvv=mvv*nrootu                                               1d12s21
           joff=joffdnon                                                 1d12s21
           do jsbv1=1,nsymb                                              1d12s21
            jsbv2=multh(jsbv1,jsbv12)                                    1d12s21
            if(jsbv2.ge.jsbv1)then                                       1d12s21
             if(jsbv12.eq.1)then                                         1d12s21
              jsw=0                                                      1d12s21
              nvv=(nvirt(jsbv1)*(nvirt(jsbv2)-1))/2                      1d12s21
              joffp=joff+nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)              1d13s21
             else                                                        1d12s21
              jsw=1                                                      1d12s21
              nvv=nvirt(jsbv1)*nvirt(jsbv2)                              1d12s21
              joffp=joff                                                 1d13s21
             end if                                                      1d12s21
             nnvv=nvv*nrootu                                             1d12s21
c     Gvvrj =djiab (ab|vv") Vvv"ri
c     Gvv'rj=djiab (ab|vv") Vv'v"ri
c     Gvv'rj=djiab (ab|vv") Vv"v'ri
c     Gv'vrj=djiab (ab|vv") Vv'v"ri
c     Gv'vrj=djiab (ab|vv") Vv"v'ri
             if(jsbv2.eq.isbv1.or.jsbv2.eq.isbv2.or.jsbv1.eq.isbv1.or.   1d13s21
     $         jsbv1.eq.isbv2)then                                      1d13s21
              do imatch=1,4                                              1d13s21
               kase=0                                                    1d14s21
               itf=0                                                    8d7s24
               jtf=0                                                    8d7s24
               if(jsbv2.eq.isbv1.and.imatch.eq.1)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv2                                                1d13s21
                itf=0                                                     1d13s21
                jtf=0                                                     1d13s21
                ibl=1                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=1
               end if                                                    1d13s21
               if(jsbv2.eq.isbv2.and.imatch.eq.2)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=0                                                    1d15s21
                ibl=1                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=2
               end if                                                    1d13s21
               if(jsbv1.eq.isbv1.and.imatch.eq.3)then                    1d14s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv2                                                1d13s21
                jtf=1                                                     1d13s21
                itf=0                                                    1d15s21
                ibl=2                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=3
               end if                                                    1d13s21
               if(jsbv1.eq.isbv2.and.imatch.eq.4)then                    1d13s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=1                                                     1d13s21
                ibl=2                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=4                                                    1d13s21
               end if                                                     1d13s21
               iioffp=ioffp                                               1d12s21
               tfi=1d0                                                   1d14s21
               do li=1,4                                                 1d14s21
                if(min(kase,nfdat(2,li,isb)).gt.0)then                   1d14s21
                 jjoffp=joffp                                            1d14s21
                 tfj=1d0                                                 1d14s21
                 do lj=1,4                                               1d14s21
                  if(nfdat(2,lj,jsb).gt.0)then                           1d14s21
                   tf=tfi*tfj
                   call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0.and.isb.ne.jsb)then   3d19s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdat(2,li,isb)*nfdat(2,lj,jsb)                   1d14s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    call enough('hcddjk. 22',bc,ibc)
                    fact=0d0                                                1d13s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then             1d14s21
                      i2eu=invk1(1,lsb,lsa,isbl,1)                              1d13s21
                      if(nokk(li,lj,lsa).gt.0)then                              1d12s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjk. 23',bc,ibc)
                       kint=kmats(i2eu)                                     1d12s21
                       do j=0,nokk(li,lj,lsa)-1                                1d12s21
                        jj=ibc(ndenkf(li,lj,lsa)+j)                            1d12s21
                        jb=jj/irefo(lsa)                                   1d12s21
                        ja=jj-irefo(lsa)*jb                                1d12s21
                        kk=jb+irefo(lsb)*ja                                1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=kint+i+nhere*kk                                1d12s21
                         ji=itmpi+j+nokk(li,lj,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                              1d12s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0,   1d14s21
     $               bc(idenkf(li,lj,lsa)),mn,bc(itmpi),nokk(li,lj,lsa),1d14s21
     $                    fact,bc(intden),mn,                           1d14s21
     d' hcddjk. 14')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
                     nrow=nfdat(2,lj,jsb)*nvirt(isbl)                        1d13s21
                     nmul=nfdat(2,li,isb)*nhere2
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjk. 24',bc,ibc)
                     do i=itmp1,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     i10=i1s                                                1d12s21
                     i1n=nvirt(isbl)                                        1d13s21
                     jntden=intden                                          1d12s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s                                            1d12s21
                      if(i2.eq.i2e)i1n=i1e                                  1d12s21
                      do i1=i10,i1n                                         1d12s21
                       i1m=i1-1                                             1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                              1d12s21
                        do i=0,nfdat(2,li,isb)-1                             1d12s21
                         iad=itmp1+j+nfdat(2,lj,jsb)*(i1m+nvirt(isbl)*(  1d14s21
     $                       i2n+nhere2*i))                             1d14s21
                         bc(iad)=bc(jntden+i)                               1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdat(2,li,isb)                        1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjk. 25',bc,ibc)
                     nx=nhere2*nfdat(2,li,isb)*nrootu                        1d12s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibr.eq.1)then                                       1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do i=0,nfdat(2,li,isb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                         do iv1=0,ntop                                       1d12s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv1=i2-1                                              1d12s21
                       nbot=iv1+1                                            1d12s21
                       nbot=nbot+isw*(0-nbot)                               1d13s21
                       do i=0,nfdat(2,li,isb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                         do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(li.eq.1.and.isbv12.eq.1)then                     1d15s21
                      iioff=iioffp-nvirt(isbv1)*nrootu*nfdat(2,1,isb)    1d15s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do i=0,nfdat(2,1,isb)-1                           1d15s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,1,isb)*(ir    1d15s21
     $                       +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjk. 26',bc,ibc)
                     tff=tf*phs                                          1d14s21
                     call dgemm('n','n',nrow,ncol,nmul,tff,              1d14s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 15')
                     if(ibl.eq.1)then                                       1d13s21
                      do iv2=0,nvirt(isbvc)-1                                1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+jsw*(nvirt(jsbv1)-1-ntop)                   1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv1=0,ntop                                        1d12s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,lj,jsb)*(iv1+nvirt(jsbv1)*( 1d14s21
     $                       ir+nrootu*iv2))                            1d14s21
                         do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do iv1=0,nvirt(isbvc)-1                                1d12s21
                       nbot=iv1+1                                           1d13s21
                       nbot=nbot+jsw*(0-nbot)                               1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv2=nbot,nvirt(jsbv2)-1                          1d13s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,lj,jsb)*(iv2+nvirt(jsbv2)*( 1d14s21
     $                       ir+nrootu*iv1))                            1d14s21
                         do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(lj.eq.1.and.jsbv12.eq.1)then                     1d14s21
                      nv=nvirt(jsbv1)*nrootu                             1d14s21
                      do iv2=0,nvirt(isbvc)-1                            1d14s21
                       do ir=0,nrootm                                    1d14s21
                        iad=joff+iv2+nvirt(jsbv1)*ir                     1d14s21
                        jprod=iprod+nfdat(2,lj,jsb)*(iv2+nvirt(jsbv1)    1d14s21
     $                      *(ir+nrootu*iv2))                           1d14s21
                        do j=0,nfdat(2,lj,jsb)-1                         1d14s21
                         gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh       1d14s21
                        end do                                           1d14s21
                       end do                                            1d14s21
                      end do                                             1d14s21
                     end if                                              1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                1d14s21
                   call ilimts(nvirt(isbr),nvirt(isbl),mynprocg,            1d13s21
     $                   mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0)then                 3d19s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdat(2,li,isb)*nfdat(2,lj,jsb)                  1d15s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    if(ibcoff.lt.0)write(6,*)('enough6 '),intden,ibcoff,
     $                   mn,nhere
                    call enough('hcddjk. 27',bc,ibc)
                    fact=0d0                                                1d13s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then            1d14s21
                      i2eu=invk1(1,lsa,lsb,isbr,1)                      1d14s21
                      if(nokk(li,lj,lsa).gt.0)then                      1d14s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjk. 28',bc,ibc)
                       kint=kmats(i2eu)                                 1d15s21
                       do j=0,nokk(li,lj,lsa)-1                                 1d12s21
                        jj=ibc(ndenkf(li,lj,lsa)+j)                              1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=kint+i+nhere*jj                                 1d12s21
                         ji=itmpi+j+nokk(li,lj,lsa)*i                   1d14s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0, 1d14s21
     $                  bc(idenkf(li,lj,lsa)),mn,bc(itmpi),             1d14s21
     $                       nokk(li,lj,lsa),fact,bc(intden),mn,        1d14s21
     d' hcddjk. 16')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
c        intden
c       isbr  isbl   ibr ibl
c     1 isbv2 jsbv1   1   1
c     2 isbv1 jsbv1   2   1
c     3 isbv2 jsbv2   1   2
c     4 isbv1 jsbv2   2   2
                     nrow=nfdat(2,li,isb)*nvirt(isbr)                       1d13s21
                     nmul=nfdat(2,lj,jsb)*nhere2                            1d13s21
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjk. 29',bc,ibc)
                     do i=itmp1,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     i10=i1s                                                1d12s21
                     i1n=nvirt(isbr)                                        1d13s21
                     jntden=intden                                          1d12s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s                                            1d12s21
                      if(i2.eq.i2e)i1n=i1e                                  1d12s21
                      do i1=i10,i1n                                         1d12s21
                       i1m=i1-1                                             1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                             1d12s21
                        iad=itmp1+nfdat(2,li,isb)*(i1m+nvirt(isbr)*(i2n     1d13s21
     $                  +nhere2*j))                                     1d12s21
                        do i=0,nfdat(2,li,isb)-1                              1d12s21
                         bc(iad+i)=bc(jntden+i)                            1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdat(2,li,isb)                        1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjk. 30',bc,ibc)
                     nx=nhere2*nfdat(2,lj,jsb)*nrootu                        1d12s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibl.eq.2)then                                      1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+jsw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                         do iv1=0,ntop                                       1d12s21
                          irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+jsw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv1=i2-1                                              1d12s21
                       nbot=iv1+1                                            1d12s21
                       nbot=nbot+jsw*(0-nbot)                               1d13s21
                       do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                         do iv2=nbot,nvirt(jsbv2)-1                         1d13s21
                          irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+jsw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv2*nx)=vd(iuse+itri)                   1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(lj.eq.1.and.jsbv12.eq.1)then                     1d15s21
                      jjoff=jjoffp-nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)    1d15s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do j=0,nfdat(2,1,jsb)-1                           1d15s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=jjoff+iv1+nvirt(jsbv1)*(ir+nrootu*j)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,1,jsb)*(ir    1d15s21
     $                        +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjk. 31',bc,ibc)
                     tfu=tf*phs                                         1d15s21
                     call dgemm('n','n',nrow,ncol,nmul,tfu,             1d15s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 17')
c     prod
c       isbr  isbvc  ibr ibl
c     1 isbv2 isbv1   1   1
c     2 isbv1 isbv2   2   1
c     3 isbv2 isbv1   1   2
c     4 isbv1 isbv2   2   2
                     if(ibr.eq.2)then                                       1d13s21
                      do iv2=0,nvirt(isbvc)-1                                1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                   1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv1=0,ntop                                        1d12s21
                         irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+isw*(irec-itri)                         1d13s21
                         iad=iioffp+itri+mvv*ir                               1d13s21
                         jprod=iprod+nfdat(2,li,isb)*(iv1+nvirt(isbv1)* 1d14s21
     $                        (ir+nrootu*iv2))                          1d14s21
                         do i=0,nfdat(2,li,isb)-1                             1d13s21
                          gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do iv1=0,nvirt(isbvc)-1                                1d12s21
                       nbot=iv1+1                                           1d13s21
                       nbot=nbot+isw*(0-nbot)                               1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv2=nbot,nvirt(isbv2)-1                          1d13s21
                         irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+isw*(irec-itri)                         1d13s21
                         iad=iioffp+itri+mvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,li,isb)*(iv2+nvirt(isbv2)* 1d14s21
     $                        (ir+nrootu*iv1))                          1d14s21
                         do i=0,nfdat(2,li,isb)-1                             1d13s21
                          gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(li.eq.1.and.isbv12.eq.1)then                    1d14s21
                      nv=nvirt(isbv1)*nrootu                            1d14s21
                      do iv2=0,nvirt(isbvc)-1                           1d14s21
                       do ir=0,nrootm                                   1d14s21
                        iad=ioff+iv2+nvirt(isbv1)*ir                    1d14s21
                        jprod=iprod+nfdat(2,li,isb)*(iv2+nvirt(isbv1)   1d14s21
     $                       *(ir+nrootu*iv2))                          1d14s21
                        do i=0,nfdat(2,li,isb)-1                        1d14s21
                         gd(iad+i*nv)=gd(iad+i*nv)+bc(jprod+i)*srh      1d18s21
                        end do                                          1d14s21
                       end do                                           1d14s21
                      end do                                            1d14s21
                     end if                                             1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                   1d13s21
                   jjoffp=jjoffp+nnvv*nfdat(2,lj,jsb)                     1d14s21
                  end if                                                  1d14s21
                  if(jtf.ne.0)tfj=-1d0                                   1d15s21
                 end do                                                  1d14s21
                end if                                                   1d14s21
                iioffp=iioffp+mmvv*nfdat(2,li,isb)                       1d14s21
                if(itf.ne.0)tfi=-1d0                                     1d14s21
               end do                                                    1d14s21
               iioffp=ioffp                                               1d12s21
               jjoffp=joffp                                              1d13s21
               tf=1d0                                                     1d13s21
               do l=1,4                                                   1d12s21
                if(min(kase,nfdat(2,l,isb),nfdat(2,l,jsb)).gt.0)then     1d14s21
                 mn=nfdat(2,l,isb)*nfdat(2,l,jsb)                         1d12s21
                 call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                 nhere=ih+1-il                                            1d12s21
                 nhere2=i2e+1-i2s                                         1d12s21
                 if(min(nhere,nvirt(isbvc)).gt.0)then                    3d19s21
                  intden=ibcoff                                            1d12s21
                  ibcoff=intden+mn*nhere                                   1d12s21
                  call enough('hcddjk. 32',bc,ibc)
                  fact=0d0                                                1d13s21
                  do lsa=1,nsymb                                              1d8s21
                   lsb=multh(lsa,ijsb)                                        1d9s21
                   if(min(irefo(lsa),irefo(lsb)).gt.0.and.lsb.ge.lsa)
     $                  then                                            1d12s21
                    i2eu=inv(1,lsa,lsb,isbl)                              1d13s21
                    icase=inv(2,lsa,lsb,isbl)                             1d13s21
                    if(nokj(l,lsa).gt.0)then                              1d12s21
                     itmpi=ibcoff                                         1d12s21
                     ibcoff=itmpi+nokj(l,lsa)*nhere                       1d12s21
                     call enough('hcddjk. 33',bc,ibc)
                     jint=jmats(i2eu)                                     1d12s21
                     if(icase.eq.1)then                                   1d12s21
                      do j=0,nokj(l,lsa)-1                                 1d12s21
                       jj=ibc(ndenjf(l,lsa)+j)                              1d12s21
                       do i=0,nhere-1                                       1d12s21
                        ij=jint+i+nhere*jj                                 1d12s21
                        ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                        bc(ji)=bc(ij)                                      1d12s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     else                                                 1d12s21
                      do j=0,nokj(l,lsa)-1                                1d12s21
                       jj=ibc(ndenjf(l,lsa)+j)                            1d12s21
                       jb=jj/irefo(lsa)                                   1d12s21
                       ja=jj-irefo(lsa)*jb                                1d12s21
                       kk=jb+irefo(lsb)*ja                                1d12s21
                       do i=0,nhere-1                                       1d12s21
                        ij=jint+i+nhere*kk                                1d12s21
                        ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                        bc(ji)=bc(ij)                                      1d12s21
                       end do                                              1d12s21
                      end do                                              1d12s21
                     end if                                               1d12s21
                     call dgemm('n','n',mn,nhere,nokj(l,lsa),1d0,         1d12s21
     $                  bc(idenjf(l,lsa)),mn,bc(itmpi),nokj(l,lsa),fact,1d12s21
     $                  bc(intden),mn,                                  1d12s21
     d' hcddjk. 18')
                     fact=1d0                                             1d12s21
                     ibcoff=itmpi                                         1d13s21
                    end if                                                1d12s21
                   end if                                                 1d13s21
                  end do                                                  1d13s21
                  if(fact.gt.0.5d0)then                                    1d12s21
                   nrow=nfdat(2,l,jsb)*nvirt(isbl)                        1d13s21
                   nmul=nfdat(2,l,isb)*nhere2
                   itmp1=ibcoff                                            1d12s21
                   ibcoff=itmp1+nrow*nmul                                 1d12s21
                   call enough('hcddjk. 34',bc,ibc)
                   do i=itmp1,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
c     intden
c       isbl  isbr  ibl ibr
c     1 jsbv1 isbv2  1   1
c     2 jsbv1 isbv1  1   2
c     3 jsbv2 isbv2  2   1
c     4 jsbv2 isbv1  2   2
                   i10=i1s                                                1d12s21
                   i1n=nvirt(isbl)                                        1d13s21
                   jntden=intden                                          1d12s21
                   do i2=i2s,i2e                                          1d12s21
                    i2n=i2-i2s                                            1d12s21
                    if(i2.eq.i2e)i1n=i1e                                  1d12s21
                    do i1=i10,i1n                                         1d12s21
                     i1m=i1-1                                             1d12s21
                     do j=0,nfdat(2,l,jsb)-1                              1d12s21
                      do i=0,nfdat(2,l,isb)-1                             1d12s21
                       iad=itmp1+j+nfdat(2,l,jsb)*(i1m+nvirt(isbl)*(i2n   1d13s21
     $                  +nhere2*i))                                     1d12s21
                       bc(iad)=bc(jntden+i)                               1d13s21
                      end do                                              1d12s21
                      jntden=jntden+nfdat(2,l,isb)                        1d12s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                    i10=1                                                 1d12s21
                   end do                                                 1d12s21
                   ncol=nvirt(isbvc)*nrootu                               1d12s21
                   itmp2=ibcoff                                           1d12s21
                   ibcoff=itmp2+nmul*ncol                                 1d12s21
                   call enough('hcddjk. 35',bc,ibc)
                   nx=nhere2*nfdat(2,l,isb)*nrootu                        1d12s21
                   do i=itmp2,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
                   if(ibr.eq.1)then                                       1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv2=i2-1                                              1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                     do i=0,nfdat(2,l,isb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,l,isb)*ir)        1d13s21
                       do iv1=0,ntop                                       1d12s21
                        irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+isw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv1=i2-1                                              1d12s21
                     nbot=iv1+1                                            1d12s21
                     nbot=nbot+isw*(0-nbot)                               1d13s21
                     do i=0,nfdat(2,l,isb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,l,isb)*ir)        1d13s21
                       do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                        irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+isw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(l.eq.1.and.isbv12.eq.1)then                        1d15s21
                    iioff=iioffp-nvirt(isbv1)*nrootu*nfdat(2,1,isb)      1d15s21
                    do i2=i2s,i2e                                        1d15s21
                     iv1=i2-1                                            1d15s21
                     i2n=i2-i2s                                          1d15s21
                     do i=0,nfdat(2,1,isb)-1                              1d15s21
                      do ir=0,nrootm                                     1d15s21
                       iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)         1d15s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,1,isb)*(ir      1d15s21
     $                     +nrootu*iv1))                                1d15s21
                       bc(jtmp2)=vd(iuse)*srh                            1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                   iprod=ibcoff                                           1d12s21
                   ibcoff=iprod+nrow*ncol                                 1d12s21
                   call enough('hcddjk. 36',bc,ibc)
                   call dgemm('n','n',nrow,ncol,nmul,tf,                  1d13s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 19')
c        prod
c       isbl  isbvc ibl ibr
c     1 jsbv1 jsbv2  1   1
c     2 jsbv1 jsbv2  1   2
c     3 jsbv2 jsbv1  2   1
c     4 jsbv2 jsbv1  2   2
                   if(ibl.eq.1)then                                       1d13s21
                    do iv2=0,nvirt(isbvc)-1                                1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+jsw*(nvirt(jsbv1)-1-ntop)                   1d13s21
                     do ir=0,nrootm                                        1d12s21
                      do iv1=0,ntop                                        1d12s21
                       irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+jsw*(irec-itri)                          1d13s21
                       iad=jjoffp+itri+nvv*ir                            1d13s21
                       jprod=iprod+nfdat(2,l,jsb)*(iv1+nvirt(jsbv1)*(ir    1d13s21
     $                     +nrootu*iv2))                                  1d12s21
                       do j=0,nfdat(2,l,jsb)-1                             1d13s21
                        gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do iv1=0,nvirt(isbvc)-1                                1d12s21
                     nbot=iv1+1                                           1d13s21
                     nbot=nbot+jsw*(0-nbot)                               1d13s21
                     do ir=0,nrootm                                        1d12s21
                      do iv2=nbot,nvirt(jsbv2)-1                          1d13s21
                       irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+jsw*(irec-itri)                          1d13s21
                       iad=jjoffp+itri+nvv*ir                            1d13s21
                       jprod=iprod+nfdat(2,l,jsb)*(iv2+nvirt(jsbv2)*(ir   1d13s21
     $                     +nrootu*iv1))                                1d13s21
                       do j=0,nfdat(2,l,jsb)-1                             1d13s21
                        gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(l.eq.1.and.jsbv12.eq.1)then                        1d15s21
                    jjoff=jjoffp-nvirt(isbvc)*nrootu*nfdat(2,1,jsb)      1d15s21
                    nv=nrootu*nvirt(isbvc)                               1d15s21
                    do iv1=0,nvirt(isbvc)-1                              1d15s21
                     do ir=0,nrootm                                      1d15s21
                      iad=jjoff+iv1+nvirt(isbvc)*ir                      1d15s21
                      jprod=iprod+nfdat(2,l,jsb)*(iv1+nvirt(jsbv2)*(ir   1d15s21
     $                    +nrootu*iv1))                                 1d15s21
                      do j=0,nfdat(2,l,jsb)-1                            1d15s21
                       gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh         1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                  end if                                                   1d12s21
                  ibcoff=intden                                            1d12s21
                 end if                                                   1d13s21
                 if(isb.ne.jsb)then
                  call ilimts(nvirt(isbr),nvirt(isbl),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                  nhere=ih+1-il                                            1d12s21
                  nhere2=i2e+1-i2s                                         1d12s21
                  if(min(nhere,nvirt(isbvc)).gt.0)then                   3d19s21
                   intden=ibcoff                                            1d12s21
                   ibcoff=intden+mn*nhere                                   1d12s21
                   call enough('hcddjk. 37',bc,ibc)
                   fact=0d0                                                1d13s21
                   do lsa=1,nsymb                                              1d8s21
                    lsb=multh(lsa,ijsb)                                        1d9s21
                    if(min(irefo(lsa),irefo(lsb)).gt.0.and.              1d13s21
     $                  lsb.ge.lsa)then                                 1d13s21
                     i2eu=inv(1,lsa,lsb,isbr)                             1d13s21
                     icase=inv(2,lsa,lsb,isbr)                            1d13s21
                     if(nokj(l,lsa).gt.0)then                              1d12s21
                      itmpi=ibcoff                                         1d12s21
                      ibcoff=itmpi+nokj(l,lsa)*nhere                       1d12s21
                      call enough('hcddjk. 38',bc,ibc)
                      jint=jmats(i2eu)                                     1d12s21
                      if(icase.eq.1)then                                   1d12s21
                       do j=0,nokj(l,lsa)-1                                 1d12s21
                        jj=ibc(ndenjf(l,lsa)+j)                              1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=jint+i+nhere*jj                                 1d12s21
                         ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      else                                                 1d12s21
                       do j=0,nokj(l,lsa)-1                                1d12s21
                        jj=ibc(ndenjf(l,lsa)+j)                            1d12s21
                        jb=jj/irefo(lsa)                                   1d12s21
                        ja=jj-irefo(lsa)*jb                                1d12s21
                        kk=jb+irefo(lsb)*ja                                1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=jint+i+nhere*kk                                1d12s21
                         ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                              1d12s21
                      end if                                               1d12s21
                      call dgemm('n','n',mn,nhere,nokj(l,lsa),1d0,         1d12s21
     $                  bc(idenjf(l,lsa)),mn,bc(itmpi),nokj(l,lsa),fact,1d12s21
     $                  bc(intden),mn,                                  1d12s21
     d' hcddjk. 20')
                      fact=1d0                                             1d12s21
                      ibcoff=itmpi                                         1d13s21
                     end if                                                1d12s21
                    end if                                                 1d13s21
                   end do                                                  1d13s21
                   if(fact.gt.0.5d0)then                                    1d12s21
c     intden
c       isbr  isbl   ibr ibl
c     1 isbv2 jsbv1   1   1
c     2 isbv1 jsbv1   2   1
c     3 isbv2 jsbv2   1   2
c     4 isbv1 jsbv2   2   2
                    nrow=nfdat(2,l,isb)*nvirt(isbr)                       1d13s21
                    nmul=nfdat(2,l,jsb)*nhere2                            1d13s21
                    itmp1=ibcoff                                            1d12s21
                    ibcoff=itmp1+nrow*nmul                                 1d12s21
                    call enough('hcddjk. 39',bc,ibc)
                    do i=itmp1,ibcoff-1                                    1d12s21
                     bc(i)=0d0                                             1d12s21
                    end do                                                 1d12s21
                    i10=i1s                                                1d12s21
                    i1n=nvirt(isbr)                                        1d13s21
                    jntden=intden                                          1d12s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s                                            1d12s21
                     if(i2.eq.i2e)i1n=i1e                                  1d12s21
                     do i1=i10,i1n                                         1d12s21
                      i1m=i1-1                                             1d12s21
                      do j=0,nfdat(2,l,jsb)-1                             1d12s21
                       iad=itmp1+nfdat(2,l,isb)*(i1m+nvirt(isbr)*(i2n     1d13s21
     $                  +nhere2*j))                                     1d12s21
                       do i=0,nfdat(2,l,isb)-1                              1d12s21
                        bc(iad+i)=bc(jntden+i)                            1d13s21
                       end do                                              1d12s21
                       jntden=jntden+nfdat(2,l,isb)                        1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                     i10=1                                                 1d12s21
                    end do                                                 1d12s21
                    ncol=nvirt(isbvc)*nrootu                               1d12s21
                    itmp2=ibcoff                                           1d12s21
                    ibcoff=itmp2+nmul*ncol                                 1d12s21
                    call enough('hcddjk. 40',bc,ibc)
                    nx=nhere2*nfdat(2,l,jsb)*nrootu                        1d12s21
                    do i=itmp2,ibcoff-1                                    1d12s21
                     bc(i)=0d0                                             1d12s21
                    end do                                                 1d12s21
                    if(ibl.eq.2)then                                      1d13s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s
                      iv2=i2-1                                              1d12s21
                      ntop=iv2-1                                            1d12s21
                      ntop=ntop+jsw*(nvirt(isbvc)-1-ntop)                   1d12s21
                      do j=0,nfdat(2,l,jsb)-1                               1d12s21
                       do ir=0,nrootm                                       1d12s21
                        iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,l,jsb)*ir)        1d13s21
                        do iv1=0,ntop                                       1d12s21
                         irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    else                                                   1d13s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s
                      iv1=i2-1                                              1d12s21
                      nbot=iv1+1                                            1d12s21
                      nbot=nbot+jsw*(0-nbot)                               1d13s21
                      do j=0,nfdat(2,l,jsb)-1                               1d12s21
                       do ir=0,nrootm                                       1d12s21
                        iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,l,jsb)*ir)        1d13s21
                        do iv2=nbot,nvirt(jsbv2)-1                         1d13s21
                         irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         bc(jtmp2+iv2*nx)=vd(iuse+itri)                   1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    end if                                                 1d13s21
                    if(l.eq.1.and.jsbv12.eq.1)then                       1d15s21
                     jjoff=jjoffp-nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)     1d15s21
                     do i2=i2s,i2e                                       1d15s21
                      iv1=i2-1                                           1d15s21
                      i2n=i2-i2s                                         1d15s21
                      do j=0,nfdat(2,1,jsb)-1                             1d15s21
                       do ir=0,nrootm                                    1d15s21
                        iuse=jjoff+iv1+nvirt(jsbv1)*(ir+nrootu*j)        1d15s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,1,jsb)*(ir     1d15s21
     $                      +nrootu*iv1))                               1d15s21
                        bc(jtmp2)=vd(iuse)*srh                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end if                                               1d15s21
                    iprod=ibcoff                                           1d12s21
                    ibcoff=iprod+nrow*ncol                                 1d12s21
                    call enough('hcddjk. 41',bc,ibc)
                    call dgemm('n','n',nrow,ncol,nmul,tf,                 1d13s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 21')
c     prod
c       isbr  isbvc  ibr ibl
c     1 isbv2 isbv1   1   1
c     2 isbv1 isbv2   2   1
c     3 isbv2 isbv1   1   2
c     4 isbv1 isbv2   2   2
                    if(ibr.eq.2)then                                       1d13s21
                     do iv2=0,nvirt(isbvc)-1                                1d12s21
                      ntop=iv2-1                                            1d12s21
                      ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                   1d13s21
                      do ir=0,nrootm                                        1d12s21
                       do iv1=0,ntop                                        1d12s21
                        irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                        itri=itri+isw*(irec-itri)                         1d13s21
                        iad=iioffp+itri+mvv*ir                               1d13s21
                        jprod=iprod+nfdat(2,l,isb)*(iv1+nvirt(isbv1)*(ir    1d13s21
     $                     +nrootu*iv2))                                  1d12s21
                        do i=0,nfdat(2,l,isb)-1                             1d13s21
                         gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    else                                                   1d13s21
                     do iv1=0,nvirt(isbvc)-1                                1d12s21
                      nbot=iv1+1                                           1d13s21
                      nbot=nbot+isw*(0-nbot)                               1d13s21
                      do ir=0,nrootm                                        1d12s21
                       do iv2=nbot,nvirt(isbv2)-1                          1d13s21
                        irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                        itri=itri+isw*(irec-itri)                         1d13s21
                        iad=iioffp+itri+mvv*ir                            1d13s21
                        jprod=iprod+nfdat(2,l,isb)*(iv2+nvirt(isbv2)*(ir    1d13s21
     $                      +nrootu*iv1))                                  1d12s21
                        do i=0,nfdat(2,l,isb)-1                             1d13s21
                         gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    end if                                                 1d13s21
                    if(l.eq.1.and.isbv12.eq.1)then                        1d15s21
                     iioff=iioffp-nvirt(isbvc)*nrootu*nfdat(2,1,isb)      1d15s21
                     nv=nrootu*nvirt(isbvc)                               1d15s21
                     do iv1=0,nvirt(isbvc)-1                              1d15s21
                      do ir=0,nrootm                                      1d15s21
                       iad=iioff+iv1+nvirt(isbvc)*ir                      1d15s21
                       jprod=iprod+nfdat(2,l,isb)*(iv1+nvirt(isbv2)*(ir   1d15s21
     $                    +nrootu*iv1))                                 1d15s21
                       do i=0,nfdat(2,l,isb)-1                            1d15s21
                        gd(iad+i*nv)=gd(iad+i*nv)+bc(jprod+i)*srh        1d15s21
                       end do                                             1d15s21
                      end do                                              1d15s21
                     end do                                               1d15s21
                    end if                                                1d15s21
                   end if                                                   1d12s21
                   ibcoff=intden                                            1d12s21
                  end if                                                   1d13s21
                 end if                                                    1d13s21
                end if                                                      1d13s21
                if(mod(jtf+itf,2).ne.0)then                              1d15s21
                 tf=-1d0                                                 1d14s21
                end if                                                   1d14s21
                iioffp=iioffp+mmvv*nfdat(2,l,isb)                          1d13s21
                jjoffp=jjoffp+nnvv*nfdat(2,l,jsb)                            1d13s21
               end do                                                      1d13s21
              end do                                                     1d13s21
             end if                                                      1d13s21
             if(jsbv12.eq.1)joff=joff                                    1d12s21
     $           +nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)                    1d12s21
             do l=1,4                                                    1d12s21
              joff=joff+nnvv*nfdat(2,l,jsb)                              1d12s21
             end do                                                      1d12s21
            end if                                                       1d12s21
           end do                                                        1d12s21
           if(isbv12.eq.1)ioff=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)   1d12s21
           do l=1,4                                                      1d12s21
            ioff=ioff+mmvv*nfdat(2,l,isb)                                1d12s21
           end do                                                        1d12s21
          end if                                                         1d12s21
         end do                                                          1d12s21
         if(ijsb.eq.1)then                                               1d11s21
          ioff=ioffdnon                                                  1d11s21
          do isbv1=1,nsymb                                               1d11s21
           isbv2=multh(isbv1,isbv12)                                     1d11s21
           if(isbv2.ge.isbv1)then                                        1d11s21
            nvv=0                                                        4d28s21
            if(isbv1.eq.isbv2)then                                       1d11s21
             ioffs=ioff                                                  1d12s21
             nnn=nvirt(isbv1)*nrootu                                     1d11s21
             call ilimts(nvirt(isbv1),nrootu,mynprocg,mynowprog,il,ih,   1d11s21
     $           i1s,i1e,i2s,i2e)                                       1d11s21
             nhere=ih+1-il                                               1d11s21
             if(nhere.gt.0)then                                          1d11s21
              itmp=ibcoff                                                1d11s21
              ibcoff=itmp+nhere*nfdat(2,1,isb)                           1d11s21
              call enough('hcddjk. 42',bc,ibc)
              call dgemm('n','n',nhere,nfdat(2,1,isb),nfdat(2,1,isb),    1d11s21
     $            1d0,vd(ioff+il-1),nnn,bc(iden1ef(1)),nfdat(2,1,isb),  1d11s21
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjk. 22')
              do i=0,nfdat(2,1,isb)-1                                    1d11s21
               i10=i1s                                                   1d11s21
               i1n=nvirt(isbv1)                                          1d11s21
               jtmp=itmp+nhere*i                                         1d11s21
               do i2=i2s,i2e                                             1d11s21
                if(i2.eq.i2e)i1n=i1e                                     1d11s21
                iad=ioff-1+nvirt(isbv1)*(i2-1+nrootu*i)                  1d11s21
                do i1=i10,i1n                                            1d11s21
                 gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                 jtmp=jtmp+1                                             1d11s21
                end do                                                   1d11s21
                i10=1                                                    1d11s21
               end do                                                    1d11s21
              end do                                                     1d11s21
              ibcoff=itmp                                                1d11s21
             end if                                                      1d11s21
             if(nnn.gt.0)then                                            1d25s21
              ioff=ioff+nnn*nfdat(2,1,isb)                                1d11s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       1d11s21
              isw=0                                                       1d11s21
              itmp=ibcoff                                                 1d18s21
              ibcoff=itmp+nnn*nfdat(2,1,isb)                              1d18s21
              call enough('hcddjk. 43',bc,ibc)
              call dgemm('n','n',nnn,nfdat(2,1,isb),nfdat(2,1,isb),       1d18s21
     $              2d0,vd(ioffs),nnn,bc(idenhvvf(1)),nfdat(2,1,isb),   1d18s21
     $              0d0,bc(itmp),nnn,                                   1d18s21
     d' hcddjk. 23')
              do i=0,nfdat(2,1,isb)-1                                     1d18s21
               do ir=0,nrootm                                             1d18s21
                do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   1d18s21
                 iad=ioffs+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                 ih0=ih0av(isbv1)+(iv+irefo(isbv1))*(nh0av(isbv1)+1)      1d18s21
                 jtmp=itmp+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                 gd(iad)=gd(iad)+bc(ih0)*bc(jtmp)                         1d18s21
                end do                                                    1d18s21
               end do                                                     1d18s21
              end do                                                      1d18s21
              ibcoff=itmp                                                 1d18s21
             end if                                                      1d25s21
            else                                                         1d11s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               1d11s21
             isw=1                                                       1d11s21
            end if                                                       1d11s21
            nnn=nvv*nrootu                                               1d11s21
            tf=1d0                                                       1d11s21
            do l=1,4                                                     1d11s21
             if(nfdat(2,l,isb).gt.0)then                                 1d11s21
              call ilimts(nvv,nrootu,mynprocg,mynowprog,il,ih,            1d11s21
     $            i1s,i1e,i2s,i2e)                                       1d11s21
              nhere=ih+1-il                                               1d11s21
              if(nhere.gt.0)then                                          1d11s21
               itmp=ibcoff                                                1d11s21
               ibcoff=itmp+nhere*nfdat(2,l,isb)                           1d11s21
               if(isbv12.eq.1.and.l.eq.1)then                            1d12s21
                nrow=nvirt(isbv1)*nrootu                                 1d12s21
                ibcoff=max(ibcoff,itmp+nrow*nfdat(2,l,isb))              1d12s21
               end if                                                    1d12s21
               call enough('hcddjk. 44',bc,ibc)
               call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(2,l,isb),    1d11s21
     $            1d0,vd(ioff+il-1),nnn,bc(iden1ef(l)),nfdat(2,l,isb),  1d11s21
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjk. 24')
               do i=0,nfdat(2,l,isb)-1                                    1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                  1d11s21
                jtmp=itmp+nhere*i                                         1d11s21
                do i2=i2s,i2e                                             1d11s21
                 if(i2.eq.i2e)i1n=i1e                                     1d11s21
                 iad=ioff-1+nvv*(i2-1+nrootu*i)                          1d11s21
                 do i1=i10,i1n                                            1d11s21
                  gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                  jtmp=jtmp+1                                             1d11s21
                 end do                                                   1d11s21
                 i10=1                                                    1d11s21
                end do                                                    1d11s21
               end do                                                     1d11s21
c
c     Gvv'ri=dij Hvv" Vv'v"rj ... if sym v = sym v' = sym v" we get this
c     if v"=v' > Gvv'ri=dij Hvv' Vv'v'rj.
c     if v=v' > Gvvri dij Hvv" Vvv"rj
               if(isbv12.eq.1.and.l.eq.1.and.nrow.gt.0)then              4d28s21
                call dgemm('n','n',nrow,nfdat(2,1,isb),nfdat(2,1,isb),
     $              sr2,vd(ioffs),nrow,bc(idenhvvf(1)),nfdat(2,1,isb),  1d12s21
     $              0d0,bc(itmp),nrow,                                  1d12s21
     d' hcddjk. 25')
                do i=0,nfdat(2,l,isb)-1                                   1d11s21
                 i10=i1s                                                   1d11s21
                 i1n=nvv                                                   1d11s21
                 do i2=i2s,i2e                                             1d11s21
                  ir=i2-1                                                 1d11s21
                  if(i2.eq.i2e)i1n=i1e                                    1d11s21
                  do i1=i10,i1n                                           1d11s21
                   i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                   iv2=i1m/nvirt(isbv1)                                   1d11s21
                   iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                   try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                   iv2t=sqrt(try)+0.5d0                                   1d11s21
                   iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                   iv1=iv1t+isw*(iv1-iv1t)
                   iv2=iv2t+isw*(iv2-iv2t)
c     Gvv'ri=dij Hvv' Vv'v'rj
                   ih0=ih0av(isbv1)+iv1+irefo(isbv1)+nh0av(isbv1)*(      1d12s21
     $                irefo(isbv1)+iv2)                                 1d12s21
                   iad=ioff+i1m+nvv*(ir+nrootu*i)                        1d12s21
                   jtmp=itmp+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                   ktmp=itmp+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                   gd(iad)=gd(iad)+(bc(ktmp)+bc(jtmp))*bc(ih0)           1d12s21
                  end do                                                 1d12s21
                  i10=1                                                  4d21s21
                 end do                                                  1d12s21
                end do                                                   1d12s21
               end if                                                    1d12s21
c     Gvv'ri=dij Hvv" Vv"v'rj, Gvv'ri=dij Hv'v" Vvv"rj
c
               call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(2,l,isb),   1d11s21
     $             tf,vd(ioff+il-1),nnn,bc(idenhvvf(l)),nfdat(2,l,isb), 1d11s21
     $             0d0,bc(itmp),nhere,                                  1d11s21
     d' hcddjk. 26')
               if(l.eq.1.and.isbv12.eq.1)then                            1d18s21
                do i=0,nfdat(2,l,isb)-1                                   1d11s21
                 jtmp=itmp+nhere*i                                        1d11s21
                 i10=i1s                                                   1d11s21
                 i1n=nvv                                                   1d11s21
                 do i2=i2s,i2e                                             1d11s21
                  ir=i2-1                                                 1d11s21
                  if(i2.eq.i2e)i1n=i1e                                    1d11s21
                  ltmp=jtmp                                              1d18s21
                  do i1=i10,i1n                                           1d11s21
                   i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                   iv2=i1m/nvirt(isbv1)                                   1d11s21
                   iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                   try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                   iv2t=sqrt(try)+0.5d0                                   1d11s21
                   iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                   iv1=iv1t+isw*(iv1-iv1t)
                   iv2=iv2t+isw*(iv2-iv2t)
c     if v=v' > Gvvri dij Hvv" Vvv"rj
                   ih0=ih0av(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(      1d18s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                   iad=ioffs+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                   gd(iad)=gd(iad)+bc(ih0)*bc(jtmp)*sr2                  1d18s21
                   jtmp=jtmp+1                                           1d18s21
                  end do                                                 1d18s21
                  if(isbv1.eq.isbv2)then                                 1d18s21
                   jtmp=ltmp                                             1d18s21
                   do i1=i10,i1n                                           1d11s21
                    i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                    iv2=i1m/nvirt(isbv1)                                   1d11s21
                    iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                    try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                    iv2t=sqrt(try)+0.5d0                                   1d11s21
                    iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                    iv1=iv1t+isw*(iv1-iv1t)
                    iv2=iv2t+isw*(iv2-iv2t)
c     if v=v' > Gvvri dij Hvv" Vv"vrj
                    ih0=ih0av(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(      1d18s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                    iad=ioffs+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                    gd(iad)=gd(iad)+bc(ih0)*bc(jtmp)*sr2                  1d18s21
                    jtmp=jtmp+1                                           1d18s21
                   end do                                                 1d18s21
                  end if                                                 1d18s21
                  i10=1                                                  4d21s21
                 end do                                                  1d18s21
                end do                                                   1d18s21
               end if                                                    1d18s21
               do i=0,nfdat(2,l,isb)-1                                   1d11s21
                jtmp=itmp+nhere*i                                        1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                   1d11s21
                do i2=i2s,i2e                                             1d11s21
                 ir=i2-1                                                 1d11s21
                 if(i2.eq.i2e)i1n=i1e                                    1d11s21
                 do i1=i10,i1n                                           1d11s21
                  i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                  iv2=i1m/nvirt(isbv1)                                   1d11s21
                  iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                  try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                  iv2t=sqrt(try)+0.5d0                                   1d11s21
                  iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                  iv1=iv1t+isw*(iv1-iv1t)
                  iv2=iv2t+isw*(iv2-iv2t)
c     Gvv'ri=dij Hvv" Vv"v'rj
                  ntop=iv2-1                                             1d11s21
                  ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                    1d11s21
                  ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*(           1d12s21
     $                irefo(isbv1)+iv1)                                 1d12s21
                  iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                  do iv=0,ntop                                           1d11s21
                   irec=iv+nvirt(isbv1)*iv2                              1d11s21
                   itri=((iv2*(iv2-1))/2)+iv                             1d18s21
                   irow=itri+isw*(irec-itri)                             1d11s21
                   gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)         1d11s21
                  end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vvv"rj
                  nbot=iv1+1                                             1d11s21
                  nbot=nbot+isw*(0-nbot)                                 1d11s21
                  ih0=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*            1d12s21
     $                (irefo(isbv2)+iv2)                                1d12s21
                  do iv=nbot,nvirt(isbv2)-1                              1d11s21
                   irec=iv1+nvirt(isbv1)*iv                              1d11s21
                   itri=((iv*(iv-1))/2)+iv1                              1d18s21
                   irow=itri+isw*(irec-itri)                             1d11s21
                   gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)
                  end do                                                 1d11s21
                  if(isbv1.eq.isbv2)then
c     Gvv'ri=dij Hvv" Vv'v"rj
                   ntop=iv1-1                                             1d18s21
                   ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                    1d18s21
                   ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*(           1d12s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                   iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                   do iv=0,ntop                                           1d11s21
                    irec=iv+nvirt(isbv1)*iv1                              1d11s21
                    itri=((iv1*(iv1-1))/2)+iv                             1d18s21
                    irow=itri+isw*(irec-itri)                             1d11s21
                    gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)*tf      1d18s21
                   end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vv"vrj                                           1d18s21
                   nbot=iv2+1                                             1d18s21
                   nbot=nbot+isw*(0-nbot)                                 1d18s21
                   ih0=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*            1d18s21
     $                  (irefo(isbv2)+iv1)                                1d18s21
                   do iv=nbot,nvirt(isbv2)-1                              1d18s21
                    irec=iv2+nvirt(isbv1)*iv                              1d18s21
                    itri=((iv*(iv-1))/2)+iv2                              1d18s21
                    irow=itri+isw*(irec-itri)                             1d18s21
                    gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)*tf      1d18s21
                   end do                                                 1d11s21
                  end if                                                 1d18s21
                  jtmp=jtmp+1                                            1d11s21
                 end do                                                  1d11s21
                 i10=1                                                   1d11s21
                end do                                                   1d11s21
               end do                                                    1d11s21
               ibcoff=itmp                                               1d12s21
              end if                                                     1d11s21
              ioff=ioff+nnn*nfdat(2,l,isb)                               1d11s21
              ioffd=ioffd+nvv*nfdat(2,l,isb)                             1d11s21
             end if                                                      1d11s21
             tf=-1d0                                                     1d11s21
            end do                                                       1d11s21
           end if                                                        1d11s21
          end do                                                         1d11s21
         end if                                                          1d11s21
        end if                                                          5d10s21
        ibcoff=ibc0                                                     1d8s21
        do jsbv1=1,nsymb                                                1d8s21
         jsbv2=multh(jsbv1,jsbv12)                                      1d8s21
         if(jsbv2.ge.jsbv1)then                                         1d8s21
          if(jsbv12.eq.1)then                                           1d8s21
           joffdnon=joffdnon+nrootu*nvirt(jsbv1)*nfdat(2,1,jsb)         1d8s21
           nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                        1d8s21
          else                                                          1d8s21
           nvv=nvirt(jsbv1)*nvirt(jsbv2)                                1d8s21
          end if                                                        1d8s21
          nvv=nvv*nrootu                                                1d8s21
          do l=1,4                                                      1d8s21
           joffdnon=joffdnon+nvv*nfdat(2,l,jsb)                         1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
       end do                                                           1d8s21
       do isbv1=1,nsymb                                                 1d8s21
        isbv2=multh(isbv1,isbv12)                                       1d8s21
        if(isbv2.ge.isbv1)then                                          1d8s21
         if(isbv12.eq.1)then                                            1d8s21
          ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)          1d8s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         1d8s21
         else                                                           1d8s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 1d8s21
         end if                                                         1d8s21
         nvv=nvv*nrootu                                                 1d8s21
         do l=1,4                                                       1d8s21
          ioffdnon=ioffdnon+nvv*nfdat(2,l,isb)                          1d8s21
         end do                                                         1d8s21
        end if                                                          1d8s21
       end do                                                           1d8s21
      end do                                                            1d8s21
      return                                                            1d8s21
      end                                                               1d8s21
c mpec2.1 version zeta copyright u.s. government
      subroutine hcddjkd12(nff22,nfdat,gd,vd,nsymb,mdon,mdoo,nec,multh,    1d8s21
     $     isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,ismult,ixw1,ixw2,   1d8s21
     $     norb,nrootu,ih0av,nh0av,ioooo,jmats,kmats,ndoub,mdoub,shift, 2d18s21
     $     tdendd,tovr,sr2,srh,bc,ibc,iden,id4o,jmden,kmden,idcont,     9d28s23
     $     idorth,veco,nsdlk,isblk,nsdlkk,isblkk)                       10d3s23
      implicit real*8 (a-h,o-z)                                         1d8s21
c                                                                       1d8s21
c     0 and 2 virt contribution to gd=hdd*vd                            1d8s21
c                                                                       1d8s21
      logical lprt,lchoice,ltest                                        11d13s23
      include "common.store"                                            1d8s21
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d22s21
      integer*8 ipack8,itc,ito,jtc,jto,itestc,itesto,itesta,itestb,     1d17s21
     $     jtesta,jtestb,gandcc,gandco,gandcb                           11d1s22
      external second                                                   2d18s21
      equivalence (ipack8,ipack4)                                       1d8s21
      dimension nff22(mdoo+1,2,nsymb),nfdat(5,4,*),vd(*),gd(*),         1d8s21
     $     multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),irel(*),ism(*),       1d8s21
     $     irefo(*),ih0av(*),nh0av(*),ioooo(*),jmats(*),kmats(*),       1d8s21
     $     ipack4(2),nprei(4),mprei(4),nli(4),nprej(4),                 1d12s21
     $     mprej(4),nlj(4),nab4(2,3),idenj(4,8),ndenj(4,8),             1d9s21
     $     mdenj(4,8),idenk(4,8),ndenk(4,8),mdenk(4,8),loope(6),                 1d9s21
     $     idenhvv(4),mdenhvv(4),itest(32,3),iden1e(4),mden1e(4),       1d8s21
     $     ldenj(4,8),ldenk(4,8),iden1ef(4),idenhvvf(4),ivects(4),      11d3s22
     $     idenjf(4,8),idenkf(4,4,8),ndenjf(4,8),ndenkf(4,4,8),
     $     nokj(4,8),nokk(4,4,8),isorb(32),idorb(32),idenjn(4,8),       2d23s21
     $     idenkn(4,4,8),iden1en(4),idenhvvn(4),tdendd(*),ioxx(2),      9d28s23
     $     iden(*),id4o(*),jmden(*),kmden(*),idcont(*),idorth(4,*),     9d28s23
     $     iden1enh(4,8),veco(*),iden14o(4,512),isblk(4,*),idcj(4,8),   10d3s23
     $     isblkk(4,*),nokkk(8),iokk(8),iiden(8),idck(4,4,8),           10d11s23
     $     ipdc1(4),idchvv(4)                                           10d12s23
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/fnd2cm/inv(2,8,8,8)                                        9d2s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/4000000/
c     5573b
c     5574b
      igoul=1
      ngoul=2
      igoala=2
      igoalb=2141184
      igoalc=2135615
      lgoul=igoalc
      loop=0
      nnes=0
      ndelta2=0
      nsamex=0
      ntotcalc=0
      nnzznn=0
      mplanb=0
      ibcoffo=ibcoff
      norbx=norb+1                                                      1d8s21
      norbxx=norbx+1                                                    1d8s21
      norbxxx=norbxx+1                                                  1d8s21
      ioffdnon=1                                                        1d8s21
      nrootm=nrootu-1                                                   1d11s21
      ioffd=1                                                           1d11s21
      do isb=1,nsymb                                                    1d8s21
       isbv12=multh(isb,isymmrci)                                       1d8s21
       joffdnon=1                                                       1d8s21
       do jsb=1,isb                                                     1d8s21
        ijsb=multh(jsb,isb)                                             1d9s21
        jsbv12=multh(jsb,isymmrci)                                      1d8s21
        idcjb4=ibcoff                                                   10d2s23
        if(jsb.eq.isb)then                                              10d11s23
         do l=1,4                                                       10d11s23
          nln=nfdat(2,l,isb)*nfdat(2,l,isb)*nrootu                      10d12s23
          ipdc1(l)=ibcoff                                               10d11s23
          ibcoff=ipdc1(l)+nln
          idchvv(l)=ibcoff                                              10d12s23
          ibcoff=ibcoff+nln                                             10d12s23
         end do                                                         10d11s23
        end if                                                          10d11s23
        idcj0=ibcoff                                                    10d16s23
        do lsa=1,nsymb                                                  10d2s23
         lsb=multh(lsa,ijsb)                                            10d2s23
         if(min(irefo(lsa),irefo(lsb)).gt.0.and.lsb.ge.lsa)then         10d2s23
          if(lsa.eq.lsb)then                                            10d2s23
           mcol=(irefo(lsa)*(irefo(lsa)+1))/2                           10d2s23
          else                                                          10d2s23
           mcol=irefo(lsa)*irefo(lsb)                                   10d2s23
          end if                                                        10d2s23
          do l=1,4                                                      10d2s23
           idcj(l,lsa)=ibcoff                                           10d2s23
           ibcoff=ibcoff+nfdat(2,l,isb)*nfdat(2,l,jsb)*mcol*nrootu      10d2s23
          end do                                                        10d2s23
         end if                                                         10d2s23
        end do                                                          10d2s23
        do lsa=1,nsymb                                                  10d10s23
         lsb=multh(lsa,ijsb)                                            10d2s23
         if(min(irefo(lsa),irefo(lsb)).gt.0)then                        10d10s23
          mcol=irefo(lsa)*irefo(lsb)                                    10d10s23
          do lj=1,4                                                     10d10s23
           do li=1,4                                                    10d10s23
            idck(li,lj,lsa)=ibcoff                                      10d10s23
            ibcoff=ibcoff+nfdat(2,li,isb)*nfdat(2,lj,jsb)*mcol*nrootu   10d10s23
           end do                                                       10d10s23
          end do                                                        10d10s23
         end if                                                         10d10s23
        end do                                                          10d10s23
        ndcj=ibcoff-idcj0                                               10d16s23
        call enough('hcddjkd12.dcjb4',bc,ibc)                           10d2s23
        do iz=idcjb4,ibcoff-1                                           10d2s23
         bc(iz)=0d0                                                     10d2s23
        end do                                                          10d2s23
        ioff=ioffdnon                                                   1d12s21
        do isbv1=1,nsymb                                                10d2s23
         isbv2=multh(isbv1,isbv12)
         if(isbv2.ge.isbv1)then
          if(isbv12.eq.1)then                                           1d12s21
           isw=0                                                        1d12s21
           mvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        1d12s21
           if(isb.eq.jsb)then                                           10d11s23
            naddn=nvirt(isbv1)*nrootu                                   10d11s23
            nkn=nfdat(2,1,isb)**2                                       10d12s23
            ntri=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                      10d12s23
            do ir=0,nrootm                                               10d11s23
             do iv=0,nvirt(isbv1)-1                                      10d11s23
              ivp=iv+irefo(isbv1)                                       10d12s23
              itry=ioff+iv+nvirt(isbv1)*ir                              10d11s23
              jpdc1=ipdc1(1)+nkn*ir                                     10d12s23
              itri=(iv*(iv+1))/2                                        10d12s23
              jdchvv=idchvv(1)+nkn*ir                                   10d12s23
              ih0=ih0av(isbv1)+ivp+nh0av(isbv1)*ivp                     10d12s23
              xx=2d0*bc(ih0)
              do k=0,nfdat(2,1,isb)-1                                    10d11s23
               iad1=itry+naddn*k                                        10d11s23
               do kp=0,nfdat(2,1,isb)-1                                  10d11s23
                iad2=itry+naddn*kp                                      10d11s23
                prod=vd(iad1)*vd(iad2)*2d0                              10d16s23
                bc(jpdc1+kp)=bc(jpdc1+kp)+prod                          10d12s23
                bc(jdchvv+kp)=bc(jdchvv+kp)+prod*xx                     10d12s23
               end do                                                   10d11s23
               jpdc1=jpdc1+nfdat(2,1,isb)                               10d11s23
               jdchvv=jdchvv+nfdat(2,1,isb)                             10d12s23
              end do                                                    10d11s23
             end do                                                     10d11s23
            end do                                                      10d11s23
           end if                                                       10d11s23
           ioffp=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)
          else                                                          1d12s21
           isw=1                                                        1d12s21
           mvv=nvirt(isbv1)*nvirt(isbv2)                                1d12s21
           ioffp=ioff
          end if                                                        1d12s21
          if(isb.eq.jsb)then                                            10d11s23
           naddn=mvv*nrootu                                             10d11s23
           do ir=0,nrootm                                               10d11s23
            ioffpp=ioffp
            tf=1d0                                                       1d11s21
            do l=1,4                                                    10d11s23
             if(nfdat(2,l,isb).gt.0)then                                10d12s23
              do iv=0,mvv-1                                              10d11s23
               itry=ioffpp+iv+mvv*ir                                     10d11s23
               jpdc1=ipdc1(l)                                            10d11s23
               do k=0,nfdat(2,l,isb)-1                                    10d11s23
                iad1=itry+naddn*k                                        10d11s23
                do kp=0,nfdat(2,l,isb)-1                                  10d11s23
                 iad2=itry+naddn*kp                                      10d11s23
                 bc(jpdc1+kp)=bc(jpdc1+kp)+vd(iad1)*vd(iad2)*2d0        10d16s23
                end do                                                   10d11s23
                jpdc1=jpdc1+nfdat(2,l,isb)                               10d11s23
               end do                                                    10d11s23
              end do                                                     10d11s23
              if(isbv12.eq.1.and.l.eq.1)then                                          10d12s23
               nkn=nfdat(2,1,isb)**2                                    10d12s23
               nvn=nvirt(isbv1)*nrootu                                  10d12s23
               ntri=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                   10d12s23
               do iv2=0,nvirt(isbv2)-1                                  10d12s23
                iv2p=iv2+irefo(isbv2)                                   10d12s23
                do iv1=0,iv2-1                                          10d12s23
                 iv1p=iv1+irefo(isbv1)                                  10d12s23
                 ivv=((iv2*(iv2-1))/2)+iv1                              10d13s23
                 iadvnotv=ioffpp+ivv+mvv*ir                             10d12s23
                 iadvisv1=ioff+iv1+nvirt(isbv1)*ir                      10d12s23
                 iadvisv2=ioff+iv2+nvirt(isbv1)*ir                      10d12s23
                 ih0=ih0av(isbv1)+iv1p+nh0av(isbv1)*iv2p                10d12s23
                 xx=bc(ih0)*sr2*2d0                                     10d16s23
                 jdchvv=idchvv(l)+nkn*ir                                10d12s23
                 do k=0,nfdat(2,l,isb)-1                                10d12s23
                  do kp=0,nfdat(2,l,isb)-1                              10d12s23
                   il0=iadvnotv+naddn*k                                 10d12s23
                   ir0=iadvnotv+naddn*kp                                10d12s23
                   il1=iadvisv1+nvn*k                                   10d12s23
                   il2=iadvisv2+nvn*k                                   10d12s23
                   ir1=iadvisv1+nvn*kp                                  10d12s23
                   ir2=iadvisv2+nvn*kp                                  10d12s23
                   bc(jdchvv+kp)=bc(jdchvv+kp)+xx*(                     10d13s23
     $                  vd(il0)*(vd(ir1)+vd(ir2))                       10d13s23
     $                  +vd(ir0)*(vd(il1)+vd(il2)))                     10d13s23
                  end do                                                10d12s23
                  jdchvv=jdchvv+nfdat(2,l,isb)                          10d12s23
                 end do                                                 10d12s23
                end do                                                  10d12s23
               end do                                                   10d12s23
              end if                                                       10d12s23
              if(min(nfdat(2,l,isb),nvirt(isbv1),nvirt(isbv2)).gt.0)then7d11s24
               itmpva=ibcoff                                             10d12s23
               itmpvb=itmpva+nfdat(2,l,isb)*nvirt(isbv1)*nvirt(isbv2)    10d12s23
               itmppa=itmpvb+nfdat(2,l,isb)*nvirt(isbv1)*nvirt(isbv2)    10d12s23
               itmppb=itmppa+nfdat(2,l,isb)*nvirt(isbv1)*nvirt(isbv2)    10d12s23
               ibcoff=itmppb+nfdat(2,l,isb)*nvirt(isbv1)*nvirt(isbv2)    10d12s23
               call enough('hcddjkd12.va',bc,ibc)                        10d12s23
               do iz=itmpva,itmppa                                       10d12s23
                bc(iz)=0d0                                               10d12s23
               end do                                                    10d12s23
               do k=0,nfdat(2,l,isb)-1                                   10d12s23
                do iv2=0,nvirt(isbv2)-1                                   10d12s23
                 ntop=iv2-1+isw*(nvirt(isbv1)-iv2)                        10d12s23
                 do iv1=0,ntop                                            10d12s23
                  irec=iv1+nvirt(isbv1)*iv2                              10d12s23
                  itri=((iv2*(iv2-1))/2)+iv1                             10d12s23
                  itri=itri+isw*(irec-itri)                              10d12s23
                  iipp=ioffpp+itri+mvv*(ir+nrootu*k)                     10d12s23
                  iada=itmpva+iv1+nvirt(isbv1)*(iv2+nvirt(isbv2)*k)      10d12s23
                  iadb=itmpvb+iv2+nvirt(isbv2)*(iv1+nvirt(isbv1)*k)      10d12s23
                  bc(iada)=vd(iipp)                                      10d12s23
                  bc(iadb)=vd(iipp)                                      10d12s23
                 end do                                                  10d12s23
                end do                                                   10d12s23
               end do                                                    10d12s23
               ncol=nvirt(isbv2)*nfdat(2,l,isb)                          10d12s23
               ih0=ih0av(isbv1)+irefo(isbv1)*(nh0av(isbv1)+1)            10d12s23
               call dgemm('n','n',nvirt(isbv1),ncol,nvirt(isbv1),2d0*tf, 10d16s23
     $             bc(ih0),nh0av(isbv1),bc(itmpva),nvirt(isbv1),0d0,    10d12s23
     $             bc(itmppa),nvirt(isbv1),'hcddjkd12.pa')              10d12s23
               ncol=nvirt(isbv1)*nfdat(2,l,isb)                          10d12s23
               ih0=ih0av(isbv2)+irefo(isbv2)*(nh0av(isbv2)+1)            10d12s23
               call dgemm('n','n',nvirt(isbv2),ncol,nvirt(isbv2),2d0*tf, 10d16s23
     $             bc(ih0),nh0av(isbv2),bc(itmpvb),nvirt(isbv2),0d0,    10d12s23
     $             bc(itmppb),nvirt(isbv2),'hcddjkd12.pb')              10d12s23
               do k=0,nfdat(2,l,isb)-1                                   10d12s23
                do iv2=0,nvirt(isbv2)-1                                  10d12s23
                 ntop=iv2-1+isw*(nvirt(isbv1)-iv2)                       10d12s23
                 do iv1=0,ntop                                           10d12s23
                  irec=iv1+nvirt(isbv1)*iv2                              10d12s23
                  itri=((iv2*(iv2-1))/2)+iv1                             10d12s23
                  itri=itri+isw*(irec-itri)                              10d12s23
                  jad=ioffpp+itri+mvv*(ir+nrootu*k)                      10d12s23
                  do kp=0,nfdat(2,l,isb)-1                               10d12s23
                   lad=itmppa+iv1+nvirt(isbv1)*(iv2+nvirt(isbv2)*kp)      10d12s23
                   kad=idchvv(l)+k+nfdat(2,l,isb)*(kp
     $                  +nfdat(2,l,isb)*ir)                              10d12s23
                   bc(kad)=bc(kad)+bc(lad)*vd(jad)                       10d12s23
                   lad=itmppb+iv2+nvirt(isbv2)*(iv1+nvirt(isbv1)*kp)      10d12s23
                   bc(kad)=bc(kad)+bc(lad)*vd(jad)                       10d12s23
                  end do                                                 10d12s23
                  if(isbv1.eq.isbv2)then                                 10d12s23
                   do kp=0,nfdat(2,l,isb)-1                              10d12s23
                    lad=itmppb+iv1+nvirt(isbv1)*(iv2+nvirt(isbv2)*kp)     10d12s23
                    kad=idchvv(l)+k+nfdat(2,l,isb)*(kp
     $                  +nfdat(2,l,isb)*ir)                             10d12s23
                    bc(kad)=bc(kad)+bc(lad)*vd(jad)*tf                      10d12s23
                    lad=itmppa+iv2+nvirt(isbv2)*(iv1+nvirt(isbv1)*kp)     10d12s23
                    bc(kad)=bc(kad)+bc(lad)*vd(jad)*tf                      10d12s23
                   end do                                                10d12s23
                  end if                                                 10d12s23
                 end do                                                  10d12s23
                end do                                                   10d12s23
               end do                                                    10d12s23
               ibcoff=itmpva                                             10d12s23
              end if                                                    7d11s24
              ioffpp=ioffpp+mvv*nrootu*nfdat(2,l,isb)                    10d11s23
             end if                                                     10d12s23
             tf=-1d0                                                     1d11s21
            end do                                                      10d11s23
           end do                                                       10d11s23
          end if                                                        10d11s23
          mmvv=mvv*nrootu                                               1d12s21
          joff=joffdnon                                                 1d12s21
          do jsbv1=1,nsymb
           jsbv2=multh(jsbv1,jsbv12)
           if(jsbv2.ge.jsbv1)then
            if(jsbv12.eq.1)then                                         1d12s21
             jsw=0                                                      1d12s21
             nvv=(nvirt(jsbv1)*(nvirt(jsbv2)-1))/2                      1d12s21
             joffp=joff+nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)              1d13s21
            else                                                        1d12s21
             jsw=1                                                      1d12s21
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                              1d12s21
             joffp=joff                                                 1d13s21
            end if                                                      1d12s21
            nnvv=nvv*nrootu                                             1d12s21
            if(jsbv2.eq.isbv1.or.jsbv2.eq.isbv2.or.jsbv1.eq.isbv1.or.   1d13s21
     $         jsbv1.eq.isbv2)then                                      1d13s21
             do imatch=1,4                                              1d13s21
              kase=0                                                    1d14s21
              if(jsbv2.eq.isbv1.and.imatch.eq.1)then                    1d13s21
               isbvc=jsbv2                                               1d13s21
               isbl=jsbv1                                                1d13s21
               isbr=isbv2                                                1d13s21
               itf=0                                                     1d13s21
               jtf=0                                                     1d13s21
               ibl=1                                                     1d13s21
               ibr=1                                                     1d13s21
               kase=1
              end if                                                    1d13s21
              if(jsbv2.eq.isbv2.and.imatch.eq.2)then                    1d13s21
               isbvc=jsbv2                                               1d13s21
               isbl=jsbv1                                                1d13s21
               isbr=isbv1                                                1d13s21
               itf=1                                                     1d13s21
               jtf=0                                                    1d15s21
               ibl=1                                                     1d13s21
               ibr=2                                                     1d13s21
               kase=2
              end if                                                    1d13s21
              if(jsbv1.eq.isbv1.and.imatch.eq.3)then                    1d14s21
               isbvc=jsbv1                                               1d13s21
               isbl=jsbv2                                                1d13s21
               isbr=isbv2                                                1d13s21
               jtf=1                                                     1d13s21
               itf=0                                                    1d15s21
               ibl=2                                                     1d13s21
               ibr=1                                                     1d13s21
               kase=3
              end if                                                    1d13s21
              if(jsbv1.eq.isbv2.and.imatch.eq.4)then                    1d13s21
               isbvc=jsbv1                                               1d13s21
               isbl=jsbv2                                                1d13s21
               isbr=isbv1                                                1d13s21
               itf=1                                                     1d13s21
               jtf=1                                                     1d13s21
               ibl=2                                                     1d13s21
               ibr=2                                                     1d13s21
               kase=4                                                    1d13s21
              end if                                                     1d13s21
              if(kase.ne.0)then                                         10d10s23
              iioffp=ioffp                                               1d12s21
              tfi=1d0                                                     1d13s21
              do li=1,4                                                   1d12s21
               tfj=1d0
               jjoffp=joffp                                              1d13s21
               do lj=1,4                                                10d10s23
                if(min(nfdat(2,li,isb),nfdat(2,lj,jsb)).gt.0)then       10d10s23
                 tf=tfi*tfj                                             10d10s23
                 mn=nfdat(2,li,isb)*nfdat(2,lj,jsb)                     10d10s23
                 call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                 nhere=ih+1-il                                            1d12s21
                 nhere2=i2e+1-i2s                                         1d12s21
                 if(min(nhere,nvirt(isbvc)).gt.0)then                    3d19s21
                  if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                   phs=-1d0                                                   1d11s21
                  else                                                        1d11s21
                   phs=+1d0                                                   1d11s21
                  end if                                                      1d11s21
c     tmp3 is nvirt(isbl),nfdat(2,1,jsb),nvirt(isbvc),ir
                  nrow=nfdat(2,lj,jsb)*nvirt(isbl)                        1d13s21
                  nmul=nfdat(2,li,isb)*nhere2
                  ncol=nvirt(isbvc)*nrootu                               1d12s21
                  itmp2=ibcoff                                           1d12s21
                  ibcoff=itmp2+nmul*ncol                                 1d12s21
                  call enough('hcddjk. 35',bc,ibc)
                  nx=nhere2*nfdat(2,li,isb)*nrootu                        1d12s21
                  do i=itmp2,ibcoff-1                                    1d12s21
                   bc(i)=0d0                                             1d12s21
                  end do                                                 1d12s21
c
c     tmp2 is vr,ki,r,vc
c
                  if(ibr.eq.1)then                                       1d13s21
                   do i2=i2s,i2e                                          1d12s21
                    i2n=i2-i2s
                    iv2=i2-1                                              1d12s21
                    ntop=iv2-1                                            1d12s21
                    ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                    do i=0,nfdat(2,li,isb)-1                               1d12s21
                     do ir=0,nrootm                                       1d12s21
                      iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                      jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                      do iv1=0,ntop                                       1d12s21
                       irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                       itri=itri+isw*(irec-itri)                          1d13s21
                       bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                      end do                                              1d12s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                   end do                                                 1d12s21
                  else                                                   1d13s21
                   do i2=i2s,i2e                                          1d12s21
                    i2n=i2-i2s
                    iv1=i2-1                                              1d12s21
                    nbot=iv1+1                                            1d12s21
                    nbot=nbot+isw*(0-nbot)                               1d13s21
                    do i=0,nfdat(2,li,isb)-1                               1d12s21
                     do ir=0,nrootm                                       1d12s21
                      iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                      jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                      do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                       irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                       itri=itri+isw*(irec-itri)                          1d13s21
                       bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                      end do                                              1d12s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                   end do                                                 1d12s21
                  end if                                                 1d13s21
                  if(li.eq.1.and.isbv12.eq.1)then                         1d15s21
                   iioff=iioffp-nvirt(isbv1)*nrootu*nfdat(2,1,isb)       1d15s21
                   do i2=i2s,i2e                                         1d15s21
                    iv1=i2-1                                             1d15s21
                    i2n=i2-i2s                                           1d15s21
                    do i=0,nfdat(2,1,isb)-1                               1d15s21
                     do ir=0,nrootm                                      1d15s21
                      iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)          1d15s21
                      jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,1,isb)*(ir       1d15s21
     $                     +nrootu*iv1))                                 1d15s21
                      bc(jtmp2)=vd(iuse)*srh                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end do                                                1d15s21
                  end if                                                 1d15s21
                  itmp3=ibcoff                                           10d2s23
                  ibcoff=itmp3+nrow*ncol                                5d2s24
                  call enough('hcddjkd12.tmp3',bc,ibc)
                  do iz=itmp3,ibcoff-1                                   10d2s23
                   bc(iz)=0d0                                            10d2s23
                  end do                                                 10d2s23
c
c     tmp3 is vl,kj,vc,r
c
                  if(ibl.eq.1)then                                       1d13s21
                   do iv2=0,nvirt(isbvc)-1                                1d12s21
                    ntop=iv2-1                                            1d12s21
                    ntop=ntop+jsw*(nvirt(isbl)-1-ntop)                   1d13s21
                    do ir=0,nrootm                                        1d12s21
                     do iv1=0,ntop                                        1d12s21
                      irec=iv1+nvirt(isbl)*iv2                           10d2s23
                      itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                      itri=itri+jsw*(irec-itri)                          1d13s21
                      iad=jjoffp+itri+nvv*ir                             1d13s21
                      jtmp3=itmp3+iv1+nvirt(isbl)*nfdat(2,lj,jsb)*(iv2+
     $                      nvirt(isbvc)*ir)                             10d2s23
                      do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                       bc(jtmp3)=vd(iad+j*nnvv)                          10d2s23
                       jtmp3=jtmp3+nvirt(isbl)                           10d2s23
                      end do                                              1d12s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                   end do                                                 1d12s21
                  else                                                   1d13s21
                   do iv1=0,nvirt(isbvc)-1                                1d12s21
                    nbot=iv1+1                                           1d13s21
                    nbot=nbot+jsw*(0-nbot)                               1d13s21
                    do ir=0,nrootm                                         1d12s21
                     do iv2=nbot,nvirt(isbl)-1                           1d13s21
                      irec=iv1+nvirt(isbvc)*iv2                          10d2s23
                      itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                      itri=itri+jsw*(irec-itri)                          1d13s21
                      iad=jjoffp+itri+nvv*ir                             1d13s21
                      jtmp3=itmp3+iv2+nvirt(isbl)*nfdat(2,lj,jsb)*(iv1    10d2s23
     $                      +nvirt(isbvc)*ir)                           10d2s23
                      do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                       bc(jtmp3)=bc(jtmp3)+vd(iad+j*nnvv)                10d2s23
                       jtmp3=jtmp3+nvirt(isbl)                           10d2s23
                      end do                                              1d12s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                   end do                                                 1d12s21
                  end if                                                 1d13s21
                  if(lj.eq.1.and.jsbv12.eq.1)then                         1d15s21
                   jjoff=jjoffp-nvirt(isbvc)*nrootu*nfdat(2,1,jsb)       1d15s21
                   nv=nrootu*nvirt(isbvc)                                1d15s21
                   do iv1=0,nvirt(isbvc)-1                               1d15s21
                    do ir=0,nrootm                                       1d15s21
                     iad=jjoff+iv1+nvirt(isbvc)*ir                       1d15s21
                     jtmp3=itmp3+iv1+nvirt(isbl)*nfdat(2,1,jsb)*(iv1     10d2s23
     $                    +nvirt(isbvc)*ir)                              10d2s23
                     do j=0,nfdat(2,1,jsb)-1                             1d15s21
                      bc(jtmp3)=bc(jtmp3)+vd(iad+j*nv)*srh               10d2s23
                      jtmp3=jtmp3+nvirt(isbl)                            10d2s23
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end do                                                1d15s21
                  end if                                                 1d15s21
                  itmp4=ibcoff                                           10d2s23
                  itmp5=itmp4+nvirt(isbvc)*nhere2*nfdat(2,li,isb)         10d2s23
                  ibcoff=itmp5+nvirt(isbr)*nvirt(isbl)*nfdat(2,lj,jsb)    10d2s23
     $                *nfdat(2,li,isb)                                   10d2s23
                  call enough('hcddjkd12.tmp5',bc,ibc)                   10d2s23
                  do ir=0,nrootm                                         10d2s23
                   jtmp3=itmp3+nvirt(isbl)*nfdat(2,lj,jsb)*nvirt(isbvc)   10d2s23
     $                 *ir                                              10d2s23
                   do ivc=0,nvirt(isbvc)-1                               10d2s23
                    jtmp2=itmp2+nhere2*nfdat(2,li,isb)*(ir                10d2s23
     $                   +nrootu*ivc)                                   10d2s23
                    do k=0,nfdat(2,li,isb)-1                              10d2s23
                     do i2=i2s,i2e                                       10d2s23k
                      ivr=i2-i2s
                      jtmp4=itmp4+ivc+nvirt(isbvc)*(k                    10d2s23
     $                     +nfdat(2,li,isb)*ivr)                          10d2s23
                      bc(jtmp4)=bc(jtmp2+ivr)                            10d2s23
                     end do                                              10d2s23
                     jtmp2=jtmp2+nhere2                                  10d2s23
                    end do                                               10d2s23
                   end do                                                10d2s23
                   call dgemm('n','n',nrow,nmul,nvirt(isbvc),2d0,       10d16s23
     $                   bc(jtmp3),nrow,bc(itmp4),nvirt(isbvc),0d0,      10d2s23
     $                   bc(itmp5),nrow,'hcddjkd12.tmp5')                10d2s23
c     tmp5 is nvirt(isbl),nfdat(2,l,jsb),nfdat(2,l,isb),nhere2           10d2s23
                   itmp6=ibcoff                                          10d2s23
                   ibcoff=itmp6+nhere*mn                                 10d2s23
                   call enough('hcddjkd12.tmp6',bc,ibc)                  10d2s23
                   i10=i1s                                               10d2s23
                   i1n=nvirt(isbl)                                       10d2s23
                   jcol=0                                                10d2s23
                   do i2=i2s,i2e                                         10d2s23
                    i2n=i2-i2s                                           10d2s23
                    if(i2.eq.i2e)i1n=i1e                                 10d2s23
                    do i1=i10,i1n                                        10d2s23
                     i1m=i1-1                                            10d2s23
                     do kj=0,nfdat(2,lj,jsb)-1                            10d2s23
                      do ki=0,nfdat(2,li,isb)-1                           10d2s23
                       iad1=itmp6+ki+nfdat(2,li,isb)*(kj                  10d2s23
     $                     +nfdat(2,lj,jsb)*jcol)                        10d2s23
                       iad2=itmp5+i1m+nvirt(isbl)*(kj+nfdat(2,lj,jsb)      10d2s23
     $                       *(ki+nfdat(2,li,isb)*i2n))                    10d2s23
                       bc(iad1)=bc(iad2)                                  10d2s23
                      end do                                             10d2s23
                     end do                                              10d2s23
                     jcol=jcol+1                                         10d2s23
                    end do                                               10d2s23
                    i10=1                                                10d2s23
                   end do                                                10d2s23
c     tmp6 is now ki,kj,nhere.
c     mult by integrals to yield ki,kj,i,i'
                   tff=tf*phs                                           10d10s23
                   do lsa=1,nsymb                                       10d10s23
                    lsb=multh(lsa,ijsb)                                        1d9s21
                    if(min(irefo(lsa),irefo(lsb)).gt.0)then             1d14s21
                     mrow=irefo(lsa)*irefo(lsb)                         10d2s23
                     jans=idck(li,lj,lsa)+mn*mrow*ir                         10d2s23
c     need to swap lsa,lsb indices to match denk
                     i2eu=invk1(1,lsb,lsa,isbl,1)                       10d10s23
                     itmpi=ibcoff                                       10d11s23
                     ibcoff=itmpi+nhere*mrow                            10d11s23
                     call enough('hcddjkd12.tmpi',bc,ibc)               10d11s23
                     do ia=0,irefo(lsa)-1                               10d11s23
                      do ib=0,irefo(lsb)-1                              10d11s23
                       iba=kmats(i2eu)+nhere*(ib+irefo(lsb)*ia)         10d11s23
                       iab=itmpi+nhere*(ia+irefo(lsa)*ib)               10d11s23
                       do l=0,nhere-1                                   10d11s23
                        bc(iab+l)=bc(iba+l)                             10d11s23
                       end do                                           10d11s23
                      end do                                            10d11s23
                     end do                                             10d11s23
                     call dgemm('n','n',mn,mrow,nhere,tff,              10d2s23
     $                    bc(itmp6),mn,bc(itmpi),nhere,1d0,             10d11s23
     $                     bc(jans),mn,'hcddjkd12.jansb')               10d2s23
                     ibcoff=itmpi                                       10d11s23
                    end if                                              10d10s23
                   end do                                               10d10s23
                   if(li.eq.lj)then                                     10d10s23
                    do lsa=1,nsymb                                        10d2s23
                     lsb=multh(lsa,ijsb)                                  10d2s23
                     if(min(irefo(lsa),irefo(lsb)).gt.0.and.              10d2s23
     $                   lsb.ge.lsa)then                                 10d2s23
                      if(lsa.eq.lsb)then                                  10d2s23
                       mrow=(irefo(lsa)*(irefo(lsa)+1))/2                 10d2s23
                      else                                                10d2s23
                       mrow=irefo(lsa)*irefo(lsb)                         10d2s23
                      end if                                              10d2s23
                      jans=idcj(li,lsa)+mn*mrow*ir                         10d2s23
                      i2eu=inv(1,lsa,lsb,isbl)                            10d2s23
                      icase=inv(2,lsa,lsb,isbl)                           10d2s23
                      if(icase.eq.1)then                                  10d2s23
                       call dgemm('n','n',mn,mrow,nhere,tf,             10d10s23
     $                     bc(itmp6),mn,bc(jmats(i2eu)),nhere,1d0,      10d4s23
     $                     bc(jans),mn,'hcddjk12.jans')                  10d2s23
                      else                                                10d2s23
                       itmp7=ibcoff                                        10d2s23
                       ibcoff=itmp7+mn*mrow                                10d2s23
                       call enough('hcddjkd12.tmp7',bc,ibc)                10d2s23
                       call dgemm('n','n',mn,mrow,nhere,tf,             10d10s23
     $                    bc(itmp6),mn,bc(jmats(i2eu)),nhere,0d0,       10d2s23
     $                     bc(itmp7),mn,'hcddjk12.tmp7')                 10d2s23
                       do ja=0,irefo(lsa)-1                               10d2s23
                        do jb=0,irefo(lsb)-1                              10d2s23
                         jba=jans+mn*(ja+irefo(lsa)*jb)                 7d1s24
                         jab=itmp7+mn*(jb+irefo(lsb)*ja)                7d1s24
                         do i=0,mn-1                                        1d12s21
                          bc(jba+i)=bc(jba+i)+bc(jab+i)                   10d2s23
                         end do                                               1d12s21
                        end do                                               1d12s21
                       end do                                             10d2s23
                       ibcoff=itmp7                                       10d2s23
                      end if                                              10d2s23
                     end if                                               10d2s23
                    end do                                                10d2s23
                   end if                                               10d10s23
                  end do                                                 10d2s23
                  ibcoff=itmp2                                           10d2s23
                 end if                                                   1d13s21
                 if(isb.ne.jsb)then
                  call ilimts(nvirt(isbr),nvirt(isbl),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                  nhere=ih+1-il                                            1d12s21
                  nhere2=i2e+1-i2s                                         1d12s21
                  if(min(nhere,nvirt(isbvc)).gt.0)then                    3d19s21
                   if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                    phs=-1d0                                                   1d11s21
                   else                                                        1d11s21
                    phs=+1d0                                                   1d11s21
                   end if                                                      1d11s21
c       isbr  isbl   ibr ibl
c     1 isbv2 jsbv1   1   1
c     2 isbv1 jsbv1   2   1
c     3 isbv2 jsbv2   1   2
c     4 isbv1 jsbv2   2   2
                   nrow=nfdat(2,li,isb)*nvirt(isbr)                        1d13s21
                   nmul=nfdat(2,lj,jsb)*nhere2                             1d13s21
                   ncol=nvirt(isbvc)*nrootu                               1d12s21
                   itmp2=ibcoff                                           1d12s21
                   ibcoff=itmp2+nmul*ncol                                 1d12s21
                   call enough('hcddjk. 40',bc,ibc)
                   nx=nhere2*nfdat(2,lj,jsb)*nrootu                        1d12s21
                   do i=itmp2,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
c
c     tmp2 is vl,kj,r,vc
c
                   if(ibl.eq.2)then                                       1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv2=i2-1                                              1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+jsw*(nvirt(isbvc)-1-ntop)                   1d12s21
                     do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=jjoffp+nvv*(ir+nrootu*j)                      1d13s21
                       jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                       do iv1=0,ntop                                       1d12s21
                        irec=iv1+nvirt(isbvc)*iv2                         10d2s23
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+jsw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv1=i2-1                                              1d12s21
                     nbot=iv1+1                                            1d12s21
                     nbot=nbot+jsw*(0-nbot)                               1d13s21
                     do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=jjoffp+nvv*(ir+nrootu*j)                      1d13s21
                       jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                       do iv2=nbot,nvirt(jsbv2)-1                         1d13s21
                        irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+jsw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv2*nx)=vd(iuse+itri)                    1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(lj.eq.1.and.jsbv12.eq.1)then                         1d15s21
                    jjoff=jjoffp-nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)       1d15s21
                    do i2=i2s,i2e                                         1d15s21
                     iv1=i2-1                                             1d15s21
                     i2n=i2-i2s                                           1d15s21
                     do j=0,nfdat(2,1,jsb)-1                               1d15s21
                      do ir=0,nrootm                                      1d15s21
                       iuse=jjoff+iv1+nvirt(jsbv1)*(ir+nrootu*j)          1d15s21
                       jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,1,jsb)*(ir       1d15s21
     $                     +nrootu*iv1))                                 1d15s21
                       bc(jtmp2)=vd(iuse)*srh                             1d15s21
                      end do                                              1d15s21
                     end do                                               1d15s21
                    end do                                                1d15s21
                   end if                                                 1d15s21
                   itmp3=ibcoff                                           10d2s23
                   ibcoff=itmp3+nrow*ncol                               5d2s24
                   call enough('hcddjkd12.tmp3',bc,ibc)
                   do iz=itmp3,ibcoff-1                                   10d2s23
                    bc(iz)=0d0                                            10d2s23
                   end do                                                 10d2s23
c
c     tmp3 is vr,ki,vc,r
c
                   if(ibr.eq.2)then                                       1d13s21
                    do iv2=0,nvirt(isbvc)-1                                1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+isw*(nvirt(isbr)-1-ntop)                   10d2s23
                     do ir=0,nrootm                                        1d12s21
                      do iv1=0,ntop                                         1d12s21
                       irec=iv1+nvirt(isbr)*iv2                           10d2s23
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+isw*(irec-itri)                          1d13s21
                       iad=iioffp+itri+mvv*ir                               1d13s21
                       jtmp3=itmp3+iv1+nvirt(isbr)*nfdat(2,li,isb)*(iv2    10d2s23
     $                     +nvirt(isbvc)*ir)                             10d2s23
                       do i=0,nfdat(2,li,isb)-1                             1d13s21
                        bc(jtmp3)=vd(iad+i*mmvv)                          10d2s23
                        jtmp3=jtmp3+nvirt(isbr)                           10d2s23
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do iv1=0,nvirt(isbvc)-1                                1d12s21
                     nbot=iv1+1                                           1d13s21
                     nbot=nbot+isw*(0-nbot)                               1d13s21
                     do ir=0,nrootm                                        1d12s21
                      do iv2=nbot,nvirt(isbr)-1                           10d2s23
                       irec=iv1+nvirt(isbvc)*iv2                          10d2s23
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+isw*(irec-itri)                          1d13s21
                       iad=iioffp+itri+mvv*ir                             1d13s21
                       jtmp3=itmp3+iv2+nvirt(isbr)*nfdat(2,li,isb)*        10d2s23
     $                     (iv1+nvirt(isbvc)*ir)                        10d2s23
                       do i=0,nfdat(2,li,isb)-1                             1d13s21
                        bc(jtmp3)=vd(iad+i*mmvv)                          10d2s23
                        jtmp3=jtmp3+nvirt(isbr)                           10d2s23
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(li.eq.1.and.isbv12.eq.1)then                         1d15s21
                    iioff=iioffp-nvirt(isbvc)*nrootu*nfdat(2,1,isb)       1d15s21
                    nv=nrootu*nvirt(isbvc)                                1d15s21
                    do iv1=0,nvirt(isbvc)-1                               1d15s21
                     do ir=0,nrootm                                       1d15s21
                      iad=iioff+iv1+nvirt(isbvc)*ir                       1d15s21
                      jtmp3=itmp3+iv1+nvirt(isbv2)*nfdat(2,1,isb)         10d2s23
     $                     *(iv1+nvirt(isbvc)*ir)                         10d2s23
                      do i=0,nfdat(2,1,isb)-1                             1d15s21
                       bc(jtmp3)=vd(iad+i*nv)*srh                         10d2s23
                       jtmp3=jtmp3+nvirt(isbv2)                           10d2s23
                      end do                                              10d2s23
                     end do                                               1d15s21
                    end do                                                1d15s21
                   end if                                                 1d15s21
                   itmp4=ibcoff                                           10d2s23
                   itmp5=itmp4+nvirt(isbvc)*nhere2*nfdat(2,lj,jsb)         10d2s23
                   ibcoff=itmp5+nvirt(isbr)*nvirt(isbl)*nfdat(2,lj,jsb)    10d2s23
     $                  *nfdat(2,li,isb)                                 10d2s23
                   call enough('hcddjkd12.tmp5',bc,ibc)                   10d2s23
                   do ir=0,nrootm                                         10d2s23
                    jtmp3=itmp3+nvirt(isbr)*nfdat(2,li,isb)*nvirt(isbvc)   10d2s23
     $                   *ir                                            10d2s23
                    do ivc=0,nvirt(isbvc)-1                               10d2s23
                     jtmp2=itmp2+nhere2*nfdat(2,lj,jsb)*(ir                10d2s23
     $                   +nrootu*ivc)                                    10d2s23
                     do k=0,nfdat(2,lj,jsb)-1                              10d2s23
                      do i2=i2s,i2e                                       10d2s23k
                       ivl=i2-i2s
                       jtmp4=itmp4+ivc+nvirt(isbvc)*(k                    10d2s23
     $                      +nfdat(2,lj,jsb)*ivl)                          10d2s23
                       bc(jtmp4)=bc(jtmp2+ivl)                            10d2s23
                      end do                                              10d2s23
                      jtmp2=jtmp2+nhere2                                  10d2s23
                     end do                                               10d2s23
                    end do                                                10d2s23
                    call dgemm('n','n',nrow,nmul,nvirt(isbvc),2d0,      10d16s23
     $                 bc(jtmp3),nrow,bc(itmp4),nvirt(isbvc),0d0,       10d2s23
     $                   bc(itmp5),nrow,'hcddjkd12.tmp5')               10d2s23
c     tmp5 is nvirt(isbr),nfdat(2,l,isb),nfdat(2,l,jsb),nhere2          10d2s23
                    itmp6=ibcoff                                          10d2s23
                    ibcoff=itmp6+nhere*mn                                 10d2s23
                    call enough('hcddjkd12.tmp6',bc,ibc)                  10d2s23
                    i10=i1s                                               10d2s23
                    i1n=nvirt(isbr)                                       10d2s23
                    jcol=0                                                10d2s23
                    do i2=i2s,i2e                                         10d2s23
                     i2n=i2-i2s                                           10d2s23
                     if(i2.eq.i2e)i1n=i1e                                 10d2s23
                     do i1=i10,i1n                                        10d2s23
                      i1m=i1-1                                            10d2s23
                      do k=0,mn-1                                         10d2s23
                       iad1=itmp6+k+mn*jcol                               10d2s23
                       iad2=itmp5+i1m+nvirt(isbr)*(k+mn*i2n)              10d2s23
                       bc(iad1)=bc(iad2)                                  10d2s23
                      end do                                              10d2s23
                      jcol=jcol+1                                         10d2s23
                     end do                                               10d2s23
                     i10=1                                                10d2s23
                    end do                                                10d2s23
c     tmp6 is nfdat(2,l,isb),nfdat(2,l,jsb),nhere
c     mult by integrals ...
                    tff=tf*phs                                          10d10s23
                    do lsa=1,nsymb                                      10d10s23
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then             1d14s21
                      mrow=irefo(lsa)*irefo(lsb)                         10d2s23
                      jans=idck(li,lj,lsa)+mn*mrow*ir                         10d2s23
                      i2eu=invk1(1,lsa,lsb,isbr,1)                      10d10s23
                      call dgemm('n','n',mn,mrow,nhere,tff,              10d2s23
     $                     bc(itmp6),mn,bc(kmats(i2eu)),nhere,1d0,      10d2s23
     $                     bc(jans),mn,'hcddjkd12.jansb')               10d2s23
                     end if                                             10d10s23
                    end do                                              10d10s23
                    if(li.eq.lj)then                                    10d10s23
                     do lsa=1,nsymb                                        10d2s23
                      lsb=multh(lsa,ijsb)                                  10d2s23
                      if(min(irefo(lsa),irefo(lsb)).gt.0.and.              10d2s23
     $                    lsb.ge.lsa)then                                 10d2s23
                       if(lsa.eq.lsb)then                                  10d2s23
                        mrow=(irefo(lsa)*(irefo(lsa)+1))/2                 10d2s23
                       else                                                10d2s23
                        mrow=irefo(lsa)*irefo(lsb)                         10d2s23
                       end if                                              10d2s23
c     jans is nfdat(2,l,isb),nfdat(2,l,jsb),i,i',ir
                       jans=idcj(li,lsa)+mn*mrow*ir                         10d2s23
                       i2eu=inv(1,lsa,lsb,isbr)                            10d2s23
                       icase=inv(2,lsa,lsb,isbr)                           10d2s23
                       if(icase.eq.1)then                                  10d2s23
                        call dgemm('n','n',mn,mrow,nhere,tf,              10d2s23
     $                     bc(itmp6),mn,bc(jmats(i2eu)),nhere,1d0,      10d2s23
     $                     bc(jans),mn,'hcddjkd12.jansb')               10d2s23
                       else                                                10d2s23
                        itmp7=ibcoff                                        10d2s23
                        ibcoff=itmp7+mn*mrow                               10d2s23
                        call enough('hcddjkd12.tmp7',bc,ibc)                10d2s23
                        call dgemm('n','n',mn,mrow,nhere,tf,              10d2s23
     $                       bc(itmp6),mn,bc(jmats(i2eu)),nhere,0d0,       10d2s23
     $                      bc(itmp7),mn,'hcddjk12.tmp7b')                10d2s23
                        do jb=0,irefo(lsb)-1                               10d2s23
                         do ja=0,irefo(lsa)-1                              10d2s23
                          jab=itmp7+mn*(jb+irefo(lsb)*ja)               7d1s24
                          jba=jans+mn*(ja+irefo(lsa)*jb)                7d1s24
                          do i=0,mn-1                                      10d2s23
                           bc(jba+i)=bc(jba+i)+bc(jab+i)                       10d2s23
                          end do                                           10d2s23
                         end do                                            10d2s23
                        end do                                              1d12s21
                        ibcoff=itmp7                                       10d2s23
                       end if                                              10d2s23
                      end if                                               10d2s23
                     end do                                                10d2s23
                    end if                                              10d10s23
                   end do                                                 10d2s23
                   ibcoff=itmp2                                           10d2s23
                  end if                                                    1d13s21
                 end if                                                   10d2s23
                end if                                                   10d2s23
                if(jtf.ne.0)tfj=-1d0                                    10d10s23
                jjoffp=jjoffp+nnvv*nfdat(2,lj,jsb)                            1d13s21
c     end lj loop
               end do
               if(itf.ne.0)tfi=-1d0                                     1d14s21
               iioffp=iioffp+mmvv*nfdat(2,li,isb)                          1d13s21
c     end li loop
              end do                                                      1d13s21
 4114         format(3i5,5x,2i5)
              end if
             end do
            end if
            joff=joffp                                                  10d4s23
            do l=1,4                                                    10d4s23
             joff=joff+nvv*nfdat(2,l,jsb)*nrootu                        10d4s23
            end do                                                      10d4s23
           end if                                                       10d2s23
          end do
          ioff=ioffp                                                    10d4s23
          do l=1,4                                                      10d4s23
           ioff=ioff+mvv*nfdat(2,l,isb)*nrootu                          10d4s23
          end do                                                        10d4s23
         end if
        end do
        call dws_gsumf(bc(idcj0),ndcj)                                  10d16s23
        ibc0=ibcoff                                                     1d8s21
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    2d24s21
          do ll=1,4                                                     1d8s21
           if(nfdat(2,ll,jsb).gt.0)then                                 2d24s21
            nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                          2d24s21
            if(l.eq.ll)then                                             1d11s21
             iden1en(l)=ibcoff                                          2d24s21
             idenhvvn(l)=iden1en(l)+nll                                 2d24s21
             ibcoff=idenhvvn(l)+nll                                     2d24s21
             if(isb.eq.jsb)then                                         9d29s23
              do ksa=1,nsymb                                             9d28s23
               iden1enh(l,ksa)=ibcoff                                    9d28s23
               ntri=(irefo(ksa)*(irefo(ksa)+1))/2                        9d28s23
               ibcoff=ibcoff+nll*ntri                                    9d28s23
              end do                                                     9d28s23
              do is=1,nsdlk                                             9d29s23
               if(isblk(1,is).eq.isblk(2,is))then                       9d29s23
                nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2      9d29s23
                ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2      9d29s23
               else                                                     9d29s23
                nrow=irefo(isblk(1,is))*irefo(isblk(2,is))              9d29s23
                ncol=irefo(isblk(3,is))*irefo(isblk(4,is))              9d29s23
               end if                                                   9d29s23
               if(min(nrow,ncol).gt.0)then                              9d29s23
                iden14o(l,is)=ibcoff                                    9d29s23
                ibcoff=iden14o(l,is)+nrow*ncol*nll                      9d29s23
               end if                                                   9d29s23
              end do                                                    9d29s23
             end if                                                     9d29s23
            end if                                                      1d11s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(lsb.ge.lsa.and.l.eq.ll)then                             1d11s21
              if(lsb.eq.lsa)then                                        1d9s21
               nn=(irefo(lsa)*(irefo(lsa)+1))/2                         1d9s21
              else                                                      1d9s21
               nn=irefo(lsa)*irefo(lsb)                                 1d9s21
              end if                                                    1d9s21
              ndenjf(l,lsa)=ibcoff                                      2d25s21
              ibcoff=ndenjf(l,lsa)+nn                                    1d11s21
              idenjn(l,lsa)=ibcoff                                      1d11s21
              ibcoff=idenjn(l,lsa)+nn*nll                               2d24s21
             end if                                                     1d9s21
             nn=irefo(lsa)*irefo(lsb)                                   1d9s21
             ndenkf(l,ll,lsa)=ibcoff                                    2d25s21
             ibcoff=ndenkf(l,ll,lsa)+nn                                     1d11s21
             idenkn(l,ll,lsa)=ibcoff                                    1d9s21
             ibcoff=idenkn(l,ll,lsa)+nn*nll                             2d24s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call enough('hcddjk.  1',bc,ibc)
        nwds=ibcoff-ibc0                                                1d8s21
        do i=ibc0,ibcoff-1                                              1d8s21
         bc(i)=0d0                                                      1d8s21
        end do                                                          1d8s21
        ibclast=ibcoff                                                  1d17s21
        idoit=0                                                         2d23s21
        ldcont=idcont(isb)                                              10d4s23
        do ncloi=mdon,mdoo                                              1d8s21
         ncloip=ncloi+1                                                 1d8s21
         if(nff22(ncloip,1,isb).gt.0)then                               1d8s21
          iarg=ncloip-mdon                                              1d8s21
          nopeni=nec-2*ncloip                                           1d8s21
          nopenip=nopeni+2                                              1d8s21
          ibcst=ibcoff                                                  1d8s21
          ibcnd=ibcoff-1                                                1d8s21
          ivcv=nfdat(5,1,isb)+nff22(ncloip,2,isb)                       1d8s21
          do if=1,nff22(ncloip,1,isb)                                   1d8s21
           ipack8=ibc(ivcv)                                             1d8s21
           itc=ipack4(1)                                                1d8s21
           ito=ipack4(2)                                                1d8s21
           ito=ibset(ito,norbx)                                         1d8s21
           ito=ibset(ito,norbxx)                                        1d8s21
           nspacei=ibc(ivcv+1)                                          3d19s21
           ibcvects=ibcoff                                              11d3s22
           ndcont=0                                                     10d4s23
           do l=1,4                                                     3d19s21
            nli(l)=ibc(ivcv+1+l)                                        3d19s21
            ndcont=ndcont+nli(l)*ncsf2(l,iarg)                          10d4s23
            if(nli(l).gt.0)then                                         11d3s22
             iad1=ivcv+ibc(ivcv+5+l)                                    11d3s22
             iad2=iad1+nli(l)                                           11d3s22
             ivects(l)=ibcoff                                           11d3s22
             ibcoff=ivects(l)+nli(l)*ncsf2(l,iarg)                      11d3s22
             call enough('hcddjk.1.1',bc,ibc)
             do i=0,nli(l)-1                                            11d3s22
              do j=0,ncsf2(l,iarg)                                      11d3s22
               ji=iad2+j+ncsf2(l,iarg)*i                                11d3s22
               ij=ivects(l)+i+nli(l)*j                                  11d3s22
               bc(ij)=bc(ji)                                            11d3s22
              end do                                                    11d3s22
             end do                                                     11d3s22
            end if                                                      11d3s22
           end do                                                       3d19s21
           do i=1,norb                                                  1d8s21
            itest(i,1)=0                                                1d8s21
           end do                                                       1d8s21
           do i=1,norb                                                  1d8s21
            if(btest(itc,i))itest(i,1)=2                                1d8s21
            if(btest(ito,i))itest(i,1)=1                                1d8s21
           end do                                                       1d8s21
           itest(norbx,1)=1                                             1d8s21
           itest(norbxx,1)=1                                            1d8s21
           kdcont=idcont(jsb)                                           11d13s23
           do ncloj=mdon,max(ncloi-2,mdon)-1                            11d13s23
            nclojp=ncloj+1                                                1d8s21
            jarg=nclojp-mdon                                             1d8s21
            if(nff22(nclojp,1,jsb).gt.0)then                            1d14s21
             jvcv=nfdat(5,1,jsb)+nff22(nclojp,2,jsb)                      1d8s21
             do jf=1,nff22(nclojp,1,jsb)                                  1d8s21
              nspacej=ibc(jvcv+1)                                       3d19s21
              mdcont=0                                                  11d13s23
              do l=1,4                                                  11d13s23
               nlj(l)=ibc(jvcv+1+l)                                     11d13s23
               mdcont=mdcont+nlj(l)*ncsf2(l,jarg)                       11d13s23
              end do                                                    11d13s23
              jvcv=jvcv+nspacej                                         11d13s23
              kdcont=kdcont+mdcont                                      11d13s23
             end do                                                     11d13s23
            end if                                                      11d13s23
           end do                                                       11d13s23
           do ncloj=max(ncloi-2,mdon),min(ncloi+2,mdoo)                 1d8s21
            nclojp=ncloj+1                                                1d8s21
            if(nff22(nclojp,1,jsb).gt.0)then                            1d14s21
             if(mod(idoit,mynprocg).eq.mynowprog)then                   3d1s21
              jarg=nclojp-mdon                                             1d8s21
              nopenj=nec-2*nclojp                                          1d8s21
              nopenjp=nopenj+2                                           1d8s21
              jvcv=nfdat(5,1,jsb)+nff22(nclojp,2,jsb)                      1d8s21
              do jf=1,nff22(nclojp,1,jsb)                                  1d8s21
               ipack8=ibc(jvcv)                                            1d8s21
               jtc=ipack4(1)                                               1d8s21
               jto=ipack4(2)                                               1d8s21
               jto=ibset(jto,norbx)                                        1d8s21
               jto=ibset(jto,norbxx)                                       1d8s21
               nspacej=ibc(jvcv+1)                                       3d19s21
               mdcont=0                                                 11d13s23
               do l=1,4                                                  3d19s21
                nlj(l)=ibc(jvcv+1+l)                                     3d19s21
                mdcont=mdcont+nlj(l)*ncsf2(l,jarg)                      11d13s23
               end do                                                    3d19s21
               if(isb.eq.jsb)then                                        1d17s21
c                                                                       1d8s21
                gandcc=ieor(itc,jtc)                                     10d14s22
                gandco=ieor(ito,jto)                                     10d14s22
                gandcb=ior(gandcc,gandco)                                10d21s22
                ndifb=popcnt(gandcb)                                     10d21s22
                if(ndifb.le.4)then                                       10d21s22
                 ndifs=popcnt(gandco)                                     10d14s22
                 ndifd=popcnt(gandcc)                                     10d14s22
                 if(ndifd.eq.0.and.ndifs.eq.0)then                       10d21s22
                  ii=1                                                           11d21s20
                  do i=1,norb                                                    11d21s20
                   if(btest(ipack4(1),i))then                                    11d21s20
                    idorb(ii)=i                                                  11d21s20
                    ii=ii+1                                                      11d21s20
                   end if                                                        11d21s20
                  end do                                                         11d21s20
                  ii=1                                                           11d21s20
                  do i=1,norb                                                    11d21s20
                   if(btest(ipack4(2),i))then                                    11d21s20
                    isorb(ii)=i                                                  11d21s20
                    ii=ii+1                                                      11d21s20
                   end if                                                        11d21s20
                  end do                                                         11d21s20
                  sumx=0d0                                                 2d24s21
                  sumc=0d0                                                         9d23s19
                  do i=1,ncloi                                                   11d23s20
                   ig=irel(idorb(i))-1                                           11d23s20
                   is=ism(idorb(i))                                              11d23s20
                   iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
                   h0=bc(iad)
                   do j=1,ncloi                                                  11d23s20
                    jg=irel(idorb(j))-1                                          11d23s20
                    js=ism(idorb(j))                                             11d23s20
                    xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)      9d28s23
                    xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)      9d28s23
                    sumc=sumc+2d0*xj-xk                                            7d15s19
                   end do                                                          7d15s19
                   sumc=sumc+h0*2d0                                                7d15s19
                  end do                                                           7d15s19
                  sumo=0d0                                                         7d15s19
                  do i=1,nopeni                                                  11d23s20
                   ig=irel(isorb(i))-1                                           11d23s20
                   is=ism(isorb(i))                                              11d23s20
                   iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
                   h0=bc(iad)
                   sumo=sumo+h0                                                    7d15s19
                   do j=1,nopeni                                                 11d23s20
                    if(i.ne.j)then                                                 7d15s19
                     jg=irel(isorb(j))-1                                         11d23s20
                     js=ism(isorb(j))                                            11d23s20
                     xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)     9d28s23
                     sumo=sumo+xj*0.5d0                                            7d15s19
                    end if                                                         7d15s19
                   end do                                                          7d15s19
                  end do                                                           7d15s19
                  sum=sumc+sumo
                  do i=1,ncloi                                                   11d23s20
                   is=ism(idorb(i))                                              11d23s20
                   ig=irel(idorb(i))-1                                           11d23s20
                   do j=1,nopeni                                                 11d23s20
                    js=ism(isorb(j))                                             11d23s20
                    jg=irel(isorb(j))-1                                          11d23s20
                    xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)      9d28s23
                    xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)      9d28s23
                    sum=sum+xj*2d0-xk                                              7d15s19
                   end do                                                          7d15s19
                  end do                                                           7d15s19
                  sumx=sum                                               2d24s21
                  do i1=1,nopeni                                                 11d23s20
                   jsa=ism(isorb(i1))                                            11d22s20
                   jga=irel(isorb(i1))-1                                         11d22s20
                   do i2=i1+1,nopenip                                            11d19s20
                    if(i2.le.nopeni)then                                         11d23s20
                     jtesta=itc                                               11d19s20
                     jtestb=ito                                               11d19s20
                     nopenk=nopenip-2                                            11d19s20
                     karg=iarg+1                                                11d19s20
                     nqq=karg+mdon-1
                     xint=0d0
                     ksb=ism(isorb(i2))                                          11d22s20
                     kgb=irel(isorb(i2))-1                                       11d22s20
                     xint=getint(ioooo,jsa,ksb,ksb,jsa,jga,kgb,kgb,jga,  9d28s23
     $                   bc,ibc)                                        9d28s23
                     sumx=sumx-xint                                      2d24s21
                     if(nqq.ge.mdon.and.nqq.le.mdoo)then                         11d19s20
                      jtestb=ibclr(jtestb,isorb(i2))                      2d10s21
                      jtestb=ibclr(jtestb,isorb(i1))                      2d10s21
                      jtesta=ibset(jtesta,isorb(i1))                      2d10s21
                      call gandc(itc,ito,jtesta,jtestb,nopenip,nopenk,      11d19s20
     $             iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot,nab,iwpb,iwpk,       1d17s21
     $                   ncsfmid,bc,ibc)                                9d28s23
                      call gandc(jtesta,jtestb,itc,ito,nopenk,nopenip,   2d25s21
     $             karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,    2d25s21
     $                     iwpk1,ncsfmid1,bc,ibc)                       9d28s23
                      itmpid=ibcoff                                      10d11s23
                      jtmpid=itmpid+ncsf(karg)*ncsf(iarg)                10d11s23
                      ibcoff=jtmpid+ncsf(karg)*ncsf(iarg)                10d11s23
                      call enough('hcddjkd12.tmpid',bc,ibc)              10d11s23
                      call prodn(iwpb1,iwpk1,ncsf(karg),ncsf(iarg),      10d4s23
     $                  ncsfmid1,bc(jtmpid),bc,ibc,1d0,0d0)             10d11s23
                      do iii=0,ncsf(iarg)-1                                10d4s23
                       do j=0,ncsf(karg)-1                                 10d4s23
                        ji=jtmpid+j+ncsf(karg)*iii                         10d4s23
                        ij=itmpid+iii+ncsf(iarg)*j                         10d4s23
                        bc(ij)=bc(ji)                                      10d4s23
                       end do                                              10d4s23
                      end do                                               10d4s23
                      ibcoff=jtmpid                                        10d4s23
                      jtmpid=itmpid                                      10d11s23
                      jdcont=ldcont                                           10d4s23
                      iargo=1                                            2d24s21
                      do l=1,4                                           2d24s21
                       if(nli(l).gt.0)then                               2d24s21
                        iad1=ivcv+ibc(ivcv+5+l)                            3d19s21
                        iad2=iad1+nli(l)                                    3d19s21
                        itmp1=ibcoff                                     2d24s21
                        itmp2=itmp1+nli(l)*ncsf(karg)                    2d24s21
                        itmp3=itmp2+nli(l)*ncsf(karg)                    2d24s21
                        ibcoff=itmp3+nli(l)*nli(l)                       2d24s21
                        call enough('hcddjk.  2',bc,ibc)
                        call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1,    2d25s21
     $                       nli(l),iwpb1,iwpk1,bc(iad2),ncsf2(l,iarg), 2d25s21
     $                       bc(itmp1),ncsf(karg),1d0,0d0,iargo,        2d24s21
     $                       ncsf2(l,iarg),bc,ibc)                      9d28s23
                        do iii=0,nli(l)-1                                2d24s21
                         do k=0,ncsf(karg)-1                             2d24s21
                          ki=itmp1+k+ncsf(karg)*iii                      2d24s21
                          ik=itmp2+iii+nli(l)*k                          2d24s21
                          bc(ik)=bc(ki)                                  2d24s21
                         end do                                          2d24s21
                        end do                                           2d24s21
                        call dgemm('n','n',nli(l),nli(l),ncsf(karg),     2d24s21
     $                       1d0,bc(itmp2),nli(l),bc(itmp1),ncsf(karg), 9d29s23
     $                       0d0,bc(itmp3),nli(l),                      2d24s21
     d' hcddjk.  1')
                        jtmp3=itmp3                                      2d24s21
                        do j=0,nli(l)-1                                   2d23s21
                         jj=ibc(iad1+j)-1                                  1d8s21
                         jden=iden1en(l)+nfdat(2,l,isb)*jj-1             2d24s21
                         do k=0,nli(l)-1                                 2d24s21
                          kk=ibc(iad1+k)                                 2d24s21
                          bc(jden+kk)=bc(jden+kk)+bc(jtmp3+k)*xint       9d29s23
                         end do                                          2d24s21
                         jtmp3=jtmp3+nli(l)                              2d24s21
                        end do                                           2d24s21
                        itmpc0=ibcoff                                     10d4s23
                        ibcoff=itmpc0+ncsf2(l,iarg)*nli(l)               10d11s23
                        itmpc1=ibcoff                                     10d4s23
                        itmpc2=itmpc1+nli(l)*nli(l)                      10d11s23
                        ibcoff=itmpc2+ncsf2(l,iarg)*nli(l)               10d11s23
                        call enough('hcddjkd12.tmpmd',bc,ibc)             10d4s23
                        call dgemm('n','n',ncsf2(l,iarg),nli(l),         10d11s23
     $                     ncsf(karg),1d0,bc(jtmpid),ncsf(iarg),        10d4s23
     $                     bc(itmp1),ncsf(karg),0d0,bc(itmpc0),         10d4s23
     $                     ncsf2(l,iarg),'hcddjkd12.tmpmd')             10d11s23
                        do ir=0,nrootm                                    10d4s23
                         do iii=0,nli(l)-1                                 10d4s23
                          jj=ibc(iad1+iii)-1                               10d4s23
                          lad2=ipdc1(l)-1+nfdat(2,l,isb)*(jj             10d11s23
     $                        +nfdat(2,l,isb)*ir)                       10d11s23
                          lad1=itmpc1+iii                                 10d4s23
                          do j=0,nli(l)-1                                  10d4s23
                           kk=ibc(iad1+j)                                  10d4s23
                           bc(lad1)=bc(lad2+kk)                           10d4s23
                           lad1=lad1+nli(l)                              10d11s23
                          end do                                          10d4s23
                         end do                                           10d4s23
                         call dgemm('n','n',ncsf2(l,iarg),nli(l),        10d11s23
     $                        nli(l),xint,bc(itmpc0),ncsf2(l,iarg),     10d11s23
     $                        bc(itmpc1),nli(l),0d0,bc(itmpc2),         10d11s23
     $                        ncsf2(l,iarg),'hcddjkd12.tmpc2')          10d11s23
                         ltmpc2=itmpc2                                    10d4s23
                         do j=0,nli(l)-1                                 10d11s23
                          lad1=jdcont+ncsf2(l,iarg)*(j+nli(l)*ir)        10d11s23
                          do iii=0,ncsf2(l,iarg)-1                       10d11s23
                           bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)       10d4s23
                          end do                                          10d4s23
                          ltmpc2=ltmpc2+ncsf2(l,iarg)                    10d11s23
                         end do                                           10d4s23
                        end do                                            10d4s23
                        call den14oc(jsa,ksb,ksb,jsa,jga,kgb,kgb,jga,    9d29s23
     $                      irefo,itmp3,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,1d0,nli,iad1)                      10d17s23
                        ibcoff=itmp1                                     2d24s21
                       end if                                            2d24s21
                       iargo=iargo+ncsf2(l,iarg)                         2d24s21
                       jtmpid=jtmpid+ncsf2(l,iarg)                       10d11s23
                       jdcont=jdcont+nli(l)*ncsf2(l,iarg)*nrootu         10d11s23
                      end do                                             2d24s21
                     end if                                                      11d19s20
                    end if                                                       11d23s20
                   end do                                                1d18s21
                  end do                                                 1d18s21
                  iargo=1                                                2d24s21
                  jdcont=ldcont                                           10d4s23
                  do l=1,4                                               2d24s21
                   if(nli(l).gt.0)then                                   2d24s21
                    itmpt=ibcoff                                         2d24s21
                    itmpm=itmpt+ncsf2(l,iarg)*nli(l)                     2d24s21
                    ibcoff=itmpm+nli(l)*nli(l)                           2d24s21
                    call enough('hcddjk.  3',bc,ibc)
                    iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                    iad2=iad1+nli(l)                                     3d19s21
                    do iii=0,nli(l)-1                                    2d24s21
                     do j=0,ncsf2(l,iarg)-1                              2d24s21
                      ji=iad2+j+ncsf2(l,iarg)*iii                        2d24s21
                      ij=itmpt+iii+nli(l)*j                              2d24s21
                      bc(ij)=bc(ji)                                      2d24s21
                     end do                                              2d24s21
                    end do                                               2d24s21
                    call dgemm('n','n',nli(l),nli(l),ncsf2(l,iarg),1d0,  9d29s23
     $                   bc(itmpt),nli(l),bc(iad2),ncsf2(l,iarg),0d0,   2d24s21
     $                   bc(itmpm),nli(l),                              2d24s21
     d' hcddjk.  2')
                    do i=1,norb                                          9d29s23
                     fh0=0d0                                             9d29s23
                     if(btest(ipack4(1),i))then                          9d29s23
                      fh0=2d0                                            9d29s23
                     else if(btest(ipack4(2),i))then                     9d29s23
                      fh0=1d0                                            9d29s23
                     end if                                              9d29s23
                     if(fh0.gt.1d-5)then                                 9d29s23
                      jtmpm=itmpm                                          2d24s21
                      ksa=ism(i)                                         9d29s23
                      kga=irel(i)-1                                      9d29s23
                      ih0dcol=((kga*(kga+1))/2)+kga                      9d29s23
                      nh0dcol=(irefo(ksa)*(irefo(ksa)+1))/2              9d29s23
                      do iii=0,nlj(l)-1                                    2d24s21
                       jj=ibc(iad1+iii)-1                                  2d24s21
                       jjdend=iden1enh(l,ksa)-1+nfdat(2,l,isb)*(ih0dcol+   9d28s23
     $                   nh0dcol*jj)                                    9d28s23
                       do j=0,nli(l)-1                                     2d24s21
                        kk=ibc(iad1+j)                                     2d24s21
                        bc(jjdend+kk)=bc(jjdend+kk)+bc(jtmpm+j)*fh0      9d29s23
                       end do                                              2d24s21
                       jtmpm=jtmpm+nli(l)                                  2d24s21
                      end do                                               2d24s21
                     end if                                              9d29s23
                    end do                                               9d29s23
                    do i=1,ncloi                                                   11d23s20
                     ig=irel(idorb(i))-1                                           11d23s20
                     is=ism(idorb(i))                                              11d23s20
                     do j=1,ncloi                                                  11d23s20
                      jg=irel(idorb(j))-1                                          11d23s20
                      js=ism(idorb(j))                                             11d23s20
                      call den14oc(is,is,js,js,ig,ig,jg,jg,              9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,2d0,nli,iad1)                      10d17s23
                      call den14oc(is,js,js,is,ig,jg,jg,ig,              9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,-1d0,nli,iad1)                     10d17s23
                     end do                                              9d29s23
                    end do                                               9d29s23
                    do i=1,nopeni                                                  11d23s20
                     ig=irel(isorb(i))-1                                           11d23s20
                     is=ism(isorb(i))                                              11d23s20
                     do j=1,nopeni                                                 11d23s20
                      if(i.ne.j)then                                                 7d15s19
                       jg=irel(isorb(j))-1                                         11d23s20
                       js=ism(isorb(j))                                            11d23s20
                       call den14oc(is,is,js,js,ig,ig,jg,jg,              9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,0.5d0,nli,iad1)                    10d17s23
                      end if                                             9d29s23
                     end do                                              9d29s23
                    end do                                               9d29s23
                    do i=1,ncloi                                                   11d23s20
                     is=ism(idorb(i))                                              11d23s20
                     ig=irel(idorb(i))-1                                           11d23s20
                     do j=1,nopeni                                                 11d23s20
                      js=ism(isorb(j))                                             11d23s20
                      jg=irel(isorb(j))-1                                          11d23s20
                      call den14oc(is,is,js,js,ig,ig,jg,jg,              9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,2d0,nli,iad1)                      10d17s23
                      call den14oc(is,js,js,is,ig,jg,jg,ig,              9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,-1d0,nli,iad1)                     10d17s23
                     end do                                              9d29s23
                    end do                                               9d29s23
                    do i1=1,nopeni                                       9d29s23
                     jsa=ism(isorb(i1))                                  9d29s23
                     jga=irel(isorb(i1))-1                               9d29s23
                     do i2=i1+1,nopenip                                  9d29s23
                      if(i2.le.nopeni)then                               9d29s23
                       ksb=ism(isorb(i2))                                9d29s23
                       kgb=irel(isorb(i2))-1                             9d29s23
                       karg=iarg+1                                                11d19s20
                       nqq=karg+mdon-1
                       if(nqq.ge.mdon.and.nqq.le.mdoo)then               9d29s23
                        call den14oc(jsa,ksb,ksb,jsa,jga,kgb,kgb,jga,    9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,-1d0,nli,iad1)                     10d17s23
                       end if                                            9d29s23
                      end if
                     end do                                              9d29s23
                    end do                                               9d29s23
                    jtmpm=itmpm                                          2d24s21
                    do j=0,nli(l)-1                                      2d24s21
                     jj=ibc(iad1+j)-1                                    2d24s21
                     jden=iden1en(l)+nfdat(2,l,isb)*jj-1                 2d24s21
                     do k=0,nli(l)-1                                     2d24s21
                      kk=ibc(iad1+k)                                     2d24s21
                      bc(jden+kk)=bc(jden+kk)+bc(jtmpm+k)*sumx           9d29s23
                     end do                                              2d24s21
                     jtmpm=jtmpm+nli(l)                                  2d24s21
                    end do                                               2d24s21
                    itmpc1=ibcoff                                        10d11s23
                    itmpc2=itmpc1+nli(l)*nli(l)                          10d11s23
                    ibcoff=itmpc2+ncsf2(l,iarg)*nli(l)                   10d11s23
                    call enough('hcddjkd12.tmpc1',bc,ibc)                10d11s23
                    do ir=0,nrootm                                       10d11s23
                     do iii=0,nli(l)-1                                   10d11s23
                      jj=ibc(iad1+iii)-1                                 10d11s23
                      lad2=ipdc1(l)-1+nfdat(2,l,isb)*(jj                 10d11s23
     $                        +nfdat(2,l,isb)*ir)                       10d11s23
                      lad1=itmpc1+iii                                    10d11s23
                      do j=0,nli(l)-1                                    10d11s23
                       kk=ibc(iad1+j)                                    10d11s23
                       bc(lad1)=bc(lad2+kk)                              10d11s23
                       lad1=lad1+nli(l)                                  10d11s23
                      end do                                             10d11s23
                     end do                                              10d11s23
                     call dgemm('n','n',ncsf2(l,iarg),nli(l),            10d11s23
     $                        nli(l),sumx,bc(iad2),ncsf2(l,iarg),       10d11s23
     $                        bc(itmpc1),nli(l),0d0,bc(itmpc2),         10d11s23
     $                        ncsf2(l,iarg),'hcddjkd12.tmpc2')          10d11s23
                     ltmpc2=itmpc2                                       10d11s23
                     do j=0,nli(l)-1                                     10d11s23
                      lad1=jdcont+ncsf2(l,iarg)*(j+nli(l)*ir)            10d11s23
                      do iii=0,ncsf2(l,iarg)-1                           10d11s23
                       bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)          10d11s23
                      end do                                             10d11s23
                      ltmpc2=ltmpc2+ncsf2(l,iarg)                        10d11s23
                     end do                                              10d11s23
                    end do                                               10d11s23
                    ibcoff=itmpt                                         2d24s21
                   end if                                                2d24s21
                   iargo=iargo+ncsf2(l,iarg)                             2d24s21
                   jdcont=jdcont+nli(l)*ncsf2(l,iarg)*nrootu             10d11s23
                  end do                                                 2d24s21
                 else if(ndifs.eq.2.and.ndifb.eq.2)then                  10d21s22
                  do i=1,norbxx
                   if(btest(gandco,i))then                                          10d14s22
                    if((btest(ito,i).and..not.btest(jtc,i)).or.
     $                 (btest(itc,i).and.btest(jto,i)))then                           10d14s22
                     nab4(1,1)=i
                    else                                                           10d14s22
                     nab4(2,1)=i
                    end if
                   end if                                                          10d14s22
                  end do                                                           10d14s22
                  call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,  1d8s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $               ncsfmid1,bc,ibc)                                   9d28s23
                  call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,  2d24s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,  2d24s21
     $               ncsfmid1b,bc,ibc)                                  9d28s23
                  ksa=ism(nab1(1))                                       1d8s21
                  kga=irel(nab1(1))-1                                    1d8s21
                  ksb=ism(nab1(2))                                       1d8s21
                  kgb=irel(nab1(2))-1                                    1d8s21
                  ih0=ih0av(ksa)+kga+nh0av(ksa)*kgb                      1d8s21
                  ix=max(kga,kgb)                                        9d28s23
                  in=min(kga,kgb)                                        9d28s23
                  ih0dcol=((ix*(ix+1))/2)+in                             9d28s23
                  nh0dcol=(irefo(ksa)*(irefo(ksa)+1))/2                  9d28s23
                  sum=bc(ih0)                                            1d8s21
                  do i=1,norb                                            12d18s20
                   itest(i,2)=0                                          12d18s20
                  end do                                                 12d18s20
                  do i=1,norb                                            12d18s20
                   if(btest(jtc,i))then                                  12d18s20
                    itest(i,2)=2                                         12d18s20
                   end if                                                12d18s20
                   if(btest(jto,i))then                                  12d18s20
                    itest(i,2)=1                                         12d18s20
                   end if                                                12d18s20
                  end do                                                 12d18s20
                  nok=0                                                         11d13s20
                  do i=1,norb                                            1d8s21
                   ixn=min(itest(i,1),itest(i,2))
                   if(ixn.gt.0)then                                             11d13s20
                    nok=nok+1                                                   11d13s20
                    itest(nok,3)=ixn                                            11d13s20
                    itest(nok,2)=i                                              11d13s20
                   end if                                                       11d13s20
                  end do                                                        11d13s20
                  do i=1,nok                                             1d8s21
                   lsa=ism(itest(i,2))                                   1d8s21
                   lga=irel(itest(i,2))-1                                1d8s21
                   xj=getint(ioooo,lsa,lsa,ksa,ksb,lga,lga,kga,kgb,bc,   9d28s23
     $                 ibc)                                             9d28s23
                   sum=sum+xj                                            1d8s21
                   if(itest(i,3).eq.2)then                               1d8s21
                    xk=getint(ioooo,lsa,ksa,lsa,ksb,lga,kga,lga,kgb,bc,  9d28s23
     $                 ibc)                                             9d28s23
                    sum=sum+xj-xk                                        1d8s21
                   end if                                                1d8s21
                  end do                                                 1d8s21
                  if(iwpb1.gt.0)then                                     10d12s23
                   jtmpid=iwpb1                                          10d12s23
                  else                                                   10d12s23
                   itmpid=ibcoff                                          10d12s23
                   ibcoff=itmpid+ncsfmid1*ncsf(iarg)                     10d12s23
                   call enough('hcddjkd12.tmpid',bc,ibc)                 10d12s23
                   nusedi=ibc(-iwpb1)/2                                  10d12s23
                   if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1            10d12s23
                   icmp1=-iwpb1+1                                        10d12s23
                   icmp2=icmp1+ncsfmid1                                  10d12s23
                   icmp3=icmp2+nusedi                                    10d12s23
                   call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),        10d12s23
     $                 bc(itmpid),ncsf(iarg),ncsfmid1)                  10d12s23
                   jtmpid=itmpid                                         10d12s23
                  end if                                                 10d12s23
                  jdcont=ldcont                                           10d4s23
                  jargo=1                                                2d23s21
                  iargo=1                                                2d24s21
                  do l=1,4                                                1d8s21
                   if(min(nlj(l),nli(l)).gt.0)then                       2d25s21
                    jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                    jad2=jad1+nlj(l)                                     3d19s21
                    iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                    iad2=iad1+nli(l)                                     3d19s21
                    itmp1=ibcoff                                         2d23s21
                    itmp1b=itmp1+ncsfmid1*nlj(l)                         2d23s21
                    ibcoff=itmp1b+ncsfmid1*nli(l)                        2d23s21
                    call enough('hcddjk.  4',bc,ibc)
                    if(iwpk1.lt.0)then                                   2d23s21
                     nusedi=ibc(-iwpk1)/2                                2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1          2d23s21
                     icmp1=-iwpk1+1                                      2d23s21
                     icmp2=icmp1+ncsf(jarg)                              2d23s21
                     icmp3=icmp2+nusedi                                  2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   bc(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),   2d23s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                    else                                                              12d5s20
                     jwpk1=iwpk1+ncsfmid1*(jargo-1)                      2d23s21
                     call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                   1d0,bc(jwpk1),ncsfmid1,bc(jad2),ncsf2(l,jarg),
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjk.  3')
                    end if                                                            12d5s20
                    if(iwpk1b.lt.0)then                                   2d23s21
                     nusedi=ibc(-iwpk1b)/2                                2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                     icmp1=-iwpk1b+1                                      2d23s21
                     icmp2=icmp1+ncsf(iarg)                              2d23s21
                     icmp3=icmp2+nusedi                                  2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   bc(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),   2d23s21
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                    else                                                              12d5s20
                     jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                     call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                   1d0,bc(jwpk1b),ncsfmid1,bc(iad2),ncsf2(l,iarg),
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjk.  4')
                    end if                                                            12d5s20
                    itmpt=ibcoff                                         2d23s21
                    itmpm=itmpt+nli(l)*ncsfmid1                          2d24s21
                    ibcoff=itmpm+nli(l)*nlj(l)                           2d24s21
                    call enough('hcddjk.  5',bc,ibc)
                    do iii=0,nli(l)-1                                    2d23s21
                     do j=0,ncsfmid1-1                                   2d23s21
                      ji=itmp1b+j+ncsfmid1*iii                           2d23s21
                      ij=itmpt+iii+nli(l)*j                              2d23s21
                      bc(ij)=bc(ji)                                      2d23s21
                     end do                                              2d23s21
                    end do                                               2d23s21
                    call dgemm('n','n',nli(l),nlj(l),ncsfmid1,1d0,       9d28s23
     $                   bc(itmpt),nli(l),bc(itmp1),ncsfmid1,0d0,        2d23s21
     $                  bc(itmpm),nli(l),                               2d23s21
     d' hcddjk.  5')
                    jtmpm=itmpm                                          2d24s21
                    do iii=0,nlj(l)-1                                    2d24s21
                     jj=ibc(jad1+iii)-1                                  2d24s21
                     jjden=iden1en(l)+nfdat(2,l,isb)*jj-1                2d24s21
                     jjdend=iden1enh(l,ksa)-1+nfdat(2,l,isb)*(ih0dcol+   9d28s23
     $                   nh0dcol*jj)                                    9d28s23
                     do j=0,nli(l)-1                                     2d24s21
                      kk=ibc(iad1+j)                                     2d24s21
                      bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)*sum          9d28s23
                      bc(jjdend+kk)=bc(jjdend+kk)+bc(jtmpm+j)            9d28s23
                     end do                                              2d24s21
                     jtmpm=jtmpm+nli(l)                                  2d24s21
                    end do                                               2d24s21
                    itmpc0=ibcoff                                        10d12s23
                    ibcoff=itmpc0+ncsf2(l,iarg)*nlj(l)                   10d12s23
                    itmpc1=ibcoff                                        10d12s23
                    itmpc2=itmpc1+nlj(l)*nli(l)                          10d12s23
                    ibcoff=itmpc2+ncsf2(l,iarg)*nli(l)                   10d12s23
                    call enough('hcddjkd12.tmpmd',bc,ibc)                10d12s23
                    call dgemm('n','n',ncsf2(l,iarg),nlj(l),             10d12s23
     $                  ncsfmid1,1d0,bc(jtmpid),ncsf(iarg),             10d12s23
     $                     bc(itmp1),ncsfmid1,0d0,bc(itmpc0),           10d12s23
     $                     ncsf2(l,iarg),'hcddjkd12.tmpmd')             10d11s23
                    do ir=0,nrootm                                       10d12s23
                     do iii=0,nlj(l)-1                                   10d12s23
                      jj=ibc(jad1+iii)-1                                 10d12s23
                      lad2=ipdc1(l)-1+nfdat(2,l,isb)*(jj                 10d12s23
     $                     +nfdat(2,l,isb)*ir)                           10d12s23
                      lad1=itmpc1+iii                                    10d12s23
                      do j=0,nli(l)-1                                    10d12s23
                       kk=ibc(iad1+j)                                    10d12s23
                       bc(lad1)=bc(lad2+kk)                              10d12s23
                       lad1=lad1+nlj(l)                                  10d12s23
                      end do                                             10d12s23
                     end do                                              10d12s23
                     call dgemm('n','n',ncsf2(l,iarg),nli(l),            10d12s23
     $                        nlj(l),sum,bc(itmpc0),ncsf2(l,iarg),      10d12s23
     $                        bc(itmpc1),nlj(l),0d0,bc(itmpc2),         10d11s23
     $                        ncsf2(l,iarg),'hcddjkd12.tmpc2')          10d11s23
                     ltmpc2=itmpc2                                       10d12s23
                     do j=0,nli(l)-1                                     10d12s23
                      lad1=jdcont+ncsf2(l,iarg)*(j+nli(l)*ir)            10d12s23
                      do iii=0,ncsf2(l,iarg)-1                           10d12s23
                       bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)          10d12s23
                      end do                                             10d12s23
                      ltmpc2=ltmpc2+ncsf2(l,iarg)                        10d12s23
                     end do                                              10d12s23
                    end do                                               10d12s23
                    do i=1,nok                                           9d29s23
                     lsa=ism(itest(i,2))                                   1d8s21
                     lga=irel(itest(i,2))-1                                1d8s21
                     if(itest(i,3).eq.2)then                             9d29s23
                      call den14oc(lsa,lsa,ksa,ksb,lga,lga,kga,kgb,       9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,2d0,nlj,jad1)                      10d17s23
                      call den14oc(lsa,ksa,lsa,ksb,lga,kga,lga,kgb,       9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,-1d0,nlj,jad1)                     10d17s23
                     else                                                9d29s23
                      call den14oc(lsa,lsa,ksa,ksb,lga,lga,kga,kgb,       9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,1d0,nlj,jad1)                      10d17s23
                     end if
                    end do                                               9d29s23
                    ibcoff=itmp1                                         2d24s21
                   end if                                                 1d8s21
                   jtmpid=jtmpid+ncsf2(l,iarg)                           10d12s23
                   jdcont=jdcont+ncsf2(l,iarg)*nli(l)*nrootu             10d12s23
                   jargo=jargo+ncsf2(l,jarg)                             2d23s21
                   iargo=iargo+ncsf2(l,iarg)                             2d24s21
                  end do                                                  1d8s21
                  do i=1,nok                                             1d8s21
                   if(itest(i,3).eq.1)then                                  12d14s20
                    itestc=itc                                              12d14s20
                    itesto=ito                                              12d14s20
                    nopenk=nopenip                                                11d13s20
c
c     anihilate common
c
                    if(btest(itestc,itest(i,2)))then                             11d13s20
                     itestc=ibclr(itestc,itest(i,2))                             11d13s20
                     itesto=ibset(itesto,itest(i,2))                             11d13s20
                     karg=iarg-1                                                11d13s20
                     nopenk=nopenk+1                                             11d13s20
                    else                                                         11d13s20
                     itesto=ibclr(itesto,itest(i,2))                             11d13s20
                     karg=iarg                                                  11d13s20
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
                    nnot1=0                                              8d4s22
                    nnot2=0                                              8d4s22
                    if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                     call gandc(itestc,itesto,itc,ito,nopenk,nopenip,    2d24s21
     $                karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,    2d24s21
     $                iwpb1b,iwpk1b,ncsfmid1b,bc,ibc)                   9d28s23
                     nnot1=nnot1b                                        2d25s21
                     nab1(1)=nab1b(2)                                    2d25s21
                     nab1(2)=nab1b(1)                                    2d25s21
                     call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,        12d14s20
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            9d28s23
                    end if
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     lsa=ism(nab1(1))                                    1d8s21
                     lga=irel(nab1(1))-1                                 1d11s21
                     lsb=ism(nab1(2))                                    1d8s21
                     lgb=irel(nab1(2))-1                                 1d11s21
                     lsc=ism(nab2(1))                                    1d8s21
                     lgc=irel(nab2(1))-1                                 1d11s21
                     lsd=ism(nab2(2))                                    1d8s21
                     lgd=irel(nab2(2))-1                                 1d11s21
                     xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  9d28s23
     $                   bc,ibc)                                        9d28s23
                     itmpid=ibcoff                                       10d11s23
                     jtmpid=itmpid+ncsf(karg)*ncsf(iarg)                 10d11s23
                     ibcoff=jtmpid+ncsf(karg)*ncsf(iarg)                 10d11s23
                     call enough('hcddjkd12.tmpid',bc,ibc)               10d11s23
                     call prodn(iwpb1b,iwpk1b,ncsf(karg),ncsf(iarg),     10d11s23
     $                  ncsfmid1b,bc(jtmpid),bc,ibc,1d0,0d0)            10d11s23
                     do iii=0,ncsf(iarg)-1                                10d4s23
                      do j=0,ncsf(karg)-1                                 10d4s23
                       ji=jtmpid+j+ncsf(karg)*iii                         10d4s23
                       ij=itmpid+iii+ncsf(iarg)*j                         10d4s23
                       bc(ij)=bc(ji)                                      10d4s23
                      end do                                              10d4s23
                     end do                                               10d4s23
                     ibcoff=jtmpid                                        10d4s23
                     jtmpid=itmpid                                       10d11s23
                     jdcont=ldcont                                           10d4s23
                     jargo=1                                             2d23s21
                     iargo=1                                             2d24s21
                     do l=1,4                                             1d8s21
                      if(min(ncsf(karg),ncsf(jarg),ncsfmid2,             4d12s22
     $                     nli(l),nlj(l)).gt.0)then                      4d12s22
                       jad1=jvcv+ibc(jvcv+5+l)                           3d19s21
                       jad2=jad1+nlj(l)                                  3d19s21
                       iad1=ivcv+ibc(ivcv+5+l)                           3d19s21
                       iad2=iad1+nli(l)                                  3d19s21
                       itmpj=ibcoff                                      2d24s21
                       ibcoff=itmpj+ncsf(karg)*nlj(l)                    2d24s21
                       call enough('hcddjk.  6',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,bc(jad2),ncsf2(l,jarg),     2d24s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                       itmpi=ibcoff                                      2d24s21
                       itmpt=itmpi+ncsf(karg)*nli(l)                     2d24s21
                       ibcoff=itmpt+ncsf(karg)*nli(l)                    2d24s21
                       call enough('hcddjk.  7',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,bc(iad2),ncsf2(l,iarg),   2d24s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                       do iii=0,nli(l)-1                                 2d24s21
                        do k=0,ncsf(karg)-1                                2d23s21
                         ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                         ik=itmpi+iii+nli(l)*k                            2d23s21
                         bc(ik)=bc(ki)                                     2d23s21
                        end do                                             2d23s21
                       end do                                              2d23s21
                       itmpm=itmpt                                       2d24s21
                       ibcoff=itmpm+nli(l)*nlj(l)                        2d24s21
                       call enough('hcddjk.  8',bc,ibc)
                       call dgemm('n','n',nli(l),nlj(l),ncsf(karg),      2d25s21
     $                     1d0,                                         9d29s23
     $                      bc(itmpi),nli(l),bc(itmpj),ncsf(karg),0d0,  2d24s21
     $                      bc(itmpm),nli(l),                           2d24s21
     d' hcddjk.  6')
                       jtmpm=itmpm                                       2d24s21
                       do iii=0,nlj(l)-1                                   2d24s21
                        jj=ibc(jad1+iii)-1                                 2d24s21
                        jjden=iden1en(l)+nfdat(2,l,isb)*jj-1               2d24s21
                        do j=0,nli(l)-1                                    2d24s21
                         kk=ibc(iad1+j)                                    2d24s21
                         bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)*xint      9d29s23
                        end do                                             2d24s21
                        jtmpm=jtmpm+nli(l)                                 2d24s21
                       end do                                              2d24s21
                       itmpc0=ibcoff                                     10d4s23
                       ibcoff=itmpc0+ncsf2(l,iarg)*nlj(l)                10d11s23
                       itmpc1=ibcoff                                      10d4s23
                       itmpc2=itmpc1+nlj(l)*nli(l)                       10d11s23
                       ibcoff=itmpc2+ncsf2(l,iarg)*nli(l)                10d11s23
                       call enough('hcddjkd12.tmpmd',bc,ibc)             10d4s23
                       call dgemm('n','n',ncsf2(l,iarg),nlj(l),          10d11s23
     $                     ncsf(karg),1d0,bc(jtmpid),ncsf(iarg),        10d4s23
     $                     bc(itmpj),ncsf(karg),0d0,bc(itmpc0),         10d4s23
     $                     ncsf2(l,iarg),'hcddjkd12.tmpmd')             10d11s23
                       do ir=0,nrootm                                    10d4s23
                        do iii=0,nlj(l)-1                                 10d4s23
                         jj=ibc(jad1+iii)-1                               10d4s23
                         lad2=ipdc1(l)-1+nfdat(2,l,isb)*(jj              10d11s23
     $                        +nfdat(2,l,isb)*ir)                       10d11s23
                         lad1=itmpc1+iii                                 10d4s23
                         do j=0,nli(l)-1                                  10d4s23
                          kk=ibc(iad1+j)                                  10d4s23
                          bc(lad1)=bc(lad2+kk)                           10d4s23
                          lad1=lad1+nlj(l)                               10d11s23
                         end do                                          10d4s23
                        end do                                           10d4s23
                        call dgemm('n','n',ncsf2(l,iarg),nli(l),         10d11s23
     $                        nlj(l),xint,bc(itmpc0),ncsf2(l,iarg),     10d11s23
     $                        bc(itmpc1),nlj(l),0d0,bc(itmpc2),         10d11s23
     $                        ncsf2(l,iarg),'hcddjkd12.tmpc2')          10d11s23
                        ltmpc2=itmpc2                                    10d4s23
                        do j=0,nli(l)-1                                  10d11s23
                         lad1=jdcont+ncsf2(l,iarg)*(j+nli(l)*ir)         10d11s23
                         do iii=0,ncsf2(l,iarg)-1                        10d11s23
                          bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)       10d4s23
                         end do                                          10d4s23
                         ltmpc2=ltmpc2+ncsf2(l,iarg)                     10d11s23
                        end do                                           10d4s23
                       end do                                            10d4s23
                       call den14oc(lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,     9d29s23
     $                      irefo,itmpm,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,1d0,nlj,jad1)                      10d17s23
                       ibcoff=itmpj                                      2d24s21
                      end if                                              1d8s21
                      jargo=jargo+ncsf2(l,jarg)                          2d23s21
                      jtmpid=jtmpid+ncsf2(l,iarg)                        10d11s23
                      jdcont=jdcont+nli(l)*ncsf2(l,iarg)*nrootu          10d11s23
                      iargo=iargo+ncsf2(l,iarg)                          2d24s21
                     end do                                               1d8s21
                    end if                                                  12d14s20
                   end if                                                   12d14s20
                  end do                                                 1d8s21
                 else                                                    1d8s21
                  nnot=0                                                 10d21s22
                  if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s
                   nnot=4                                                          10d14s22
                   ioxx(1)=1                                                     10d17s22
                   ioxx(2)=1                                                     10d17s22
                   do i=1,norbxx                                                     10d17s22
                    if(btest(gandcb,i))then                                         10d14s22
                     if((btest(jtc,i).and.btest(ito,i)).or.                         10d17s22
     $                   (btest(jto,i).and..not.btest(itc,i)))then                   10d14s22
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
     $                   ((btest(itc,i).and..not.btest(jto,i)).or.                 10d14s22
     $                   (btest(jtc,i).and..not.btest(ito,i))))then                     10d14s22
                      if(btest(jtc,i))iswap=1                                        10d17s22
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
                    if(btest(itc,nab4(1,2)).and.
     $                   .not.btest(itc,nab4(1,1)))nbt=1                10d21s22
                   else                                                             10d17s22
                    nbt=0                                                           10d17s22
                    if(btest(jtc,nab4(2,2)).and.
     $                   .not.btest(jtc,nab4(2,1)))nbt=1                10d21s22
                   end if                                                           10d17s22
                   if(nbt.ne.0)then                                                 10d17s22
                    nab4(1,1)=nab4(1,2)                                             10d17s22
                    nab4(2,1)=nab4(2,2)                                             10d17s22
                   end if                                                           10d17s22
                  else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                   nnot=3
                   do i=1,norbxx
                    if(btest(gandcb,i))then                                         10d14s22
                     if(btest(itc,i))then
                      nab4(1,1)=i
                      nab4(1,2)=i
                     else                                                           10d14s22
                      nab4(2,1)=i
                      nab4(2,2)=i
                     end if
                    end if                                                          10d14s22
                   end do
                  end if                                                            10d14s22
                  ipssx=0                                                10d21s22
                  if(nnot.eq.3)then                                          12d8s20
                   ipssx=1                                                   12d8s20
                  else if(nnot.eq.4)then                                 11d1s22
                   ipssx=3                                                 12d18s20
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
                   itestc=itc                                              12d8s20
                   itesto=ito                                              12d8s20
                   if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                    itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                    itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                    nopenk=nopenip+1                                              11d13s20
                    karg=iarg-1                                             12d8s20
                   else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                    itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                    nopenk=nopenip-1                                              11d13s20
                    karg=iarg
                   else                                                          11d13s20
                    write(6,*)('bit not set for nab4(1,1) ='),           3d1s21
     $                   nab4(1,iu1)                                    3d1s21
                    stop 'nab4(1,1)'                                           11d27s20
                   end if                                                        11d13s20
                   if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                    itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                    itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                    nopenk=nopenk-1                                              11d13s20
                    karg=karg+1                                                  11d13s20
                   else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                    write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                 nab4(2,iu2)                                       12d18s20
                    stop 'nab4(2,1)'                                           11d27s20
                   else                                                          11d13s20
                    itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                    nopenk=nopenk+1                                              11d13s20
                   end if                                                        11d13s20
                   nqq=karg+mdon-1                                       4d20s21
                   nnot1=0                                               8d4s22
                   nnot2=0                                               8d4s22
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                   4d20s21
c
c     plan a: both contraction vectors untransposed, do
c     iwpb2*iwpk2*vecj, iwpb1b*iwpk1b*veci, then multiply the 2.
c     plan b: veci is transposed, form iwpb1bt*iwpb2 outside of l
c     loop, then multiply 5 matrices together via genmatn3
c     which is faster will depend on large ncsf(karg) is (a larger
c     value favoring plan b), and how small nli and nlj are
c     (small values favor plan a).
c
                    call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,       1d8s21
     $         karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d18s20
     $         ncsfmid2,bc,ibc)                                         9d28s23
                    rtb=dfloat(ncsf(karg))/dfloat(ncsfmid2)              11d3s22
                    rta=0d0                                              11d3s22
                    nls=0                                                11d3s22
                    do l=1,4                                             11d3s22
                     if(min(nli(l),nlj(l)).gt.0.and.
     $                   min(ncsf2(l,iarg),ncsf2(l,jarg)).gt.0)then          11d3s22
                      nls=nls+1                                          11d3s22
                      rta=max(rta,dfloat(ncsf2(l,jarg))/dfloat(nlj(l)))  11d3s22
                     end if                                              11d3s22
                    end do                                               11d3s22
                    nplanb=0
                    if(rtb.gt.-rta)then                                   11d3s22
                     call gandc(itc,ito,itestc,itesto,nopenip,nopenk,    11d3s22
     $         iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,     2d25s21
     $                  iwpk1b,ncsfmid1b,bc,ibc)                        9d28s23
                     nab1(1)=nab1b(1)                                    11d3s22
                     nab1(2)=nab1b(2)                                    11d3s22
                     lsa=ism(nab1(1))                                     1d8s21
                     lga=irel(nab1(1))-1                                  1d8s21
                     lsb=ism(nab1(2))                                     1d8s21
                     lgb=irel(nab1(2))-1                                  1d8s21
                     lsc=ism(nab2(1))                                     1d8s21
                     lgc=irel(nab2(1))-1                                  1d8s21
                     lsd=ism(nab2(2))                                     1d8s21
                     lgd=irel(nab2(2))-1                                  1d8s21
                     xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  9d28s23
     $                   bc,ibc)                                        9d28s23
                     fmul=1d0                                            9d29s23
                     if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2))     4d20s21
     $                   .and.                                          4d20s21
     $                max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))     1d12s21
     $                  fmul=0.5d0                                      9d29s23
                     xint=xint*fmul                                      9d29s23
                     if(abs(xint).gt.1d-14)then                          11d3s22
                      itmpid=ibcoff                                       10d11s23
                      ibcoff=itmpid+ncsf(karg)*ncsf(iarg)                10d12s23
                      call enough('hcddjkd12.tmpid',bc,ibc)               10d11s23
                      call prodn(iwpb1b,iwpk1b,ncsf(iarg),ncsf(karg),     10d11s23
     $                  ncsfmid1b,bc(itmpid),bc,ibc,1d0,0d0)            10d11s23
                      nplanb=1
                      nsz=ncsfmid1b*ncsfmid2                             11d3s22
                      iprod=ibcoff                                        11d3s22
                      ibcoff=iprod+nsz                                   11d3s22
                      call enough('hcddjk.8.1',bc,ibc)
                      call prodn(iwpk1b,iwpb2,ncsfmid1b,ncsfmid2,         11d3s22
     $                  ncsf(karg),bc(iprod),bc,ibc,1d0,0d0)            9d28s23
                      itrans=ibcoff
                      ibcoff=itrans+ncsfmid1b*ncsfmid2
                      call enough('hcddjk.8.3',bc,ibc)
                      do i=0,ncsfmid2-1
                       do j=0,ncsfmid1b-1
                        ji=iprod+j+ncsfmid1b*i
                        ij=itrans+i+ncsfmid2*j
                        bc(ij)=bc(ji)
                       end do
                      end do
                      iwprod=ibcoff                                      11d3s22
                      icmp1=iwprod+1                                     11d3s22
                      icmp2=icmp1+ncsfmid2                               11d3s22
                      icmp3=icmp2+nsz                                    11d3s22
                      ibcoff=icmp3+nsz                                   11d3s22
                      call enough('hcddjk.8.2',bc,ibc)
                      call cmpvcsf(ncsfmid2,ncsfmid1b,ibc(icmp1),        11d3s22
     $                    ibc(icmp2),bc(icmp3),bc(itrans),nkeep)         11d3s22
                      nkeeph=nkeep/2                                     11d3s22
                      if(2*nkeeph.ne.nkeep)nkeeph=nkeeph+1               11d3s22
                      if(ncsfmid1b+nkeeph+nkeep.lt.nsz)then              11d3s22
                       icmp3n=icmp2+nkeeph                               11d3s22
                       do i=0,nkeep-1                                    11d3s22
                        bc(icmp3n+i)=bc(icmp3+i)                         11d3s22
                       end do                                            11d3s22
                       ibc(iwprod)=nkeep                                 11d3s22
                       iwprod=-iwprod                                    11d3s22
                      else                                               11d3s22
                       iwprod=iprod                                      11d3s22
                      end if                                             11d3s22
                      jtmpid=itmpid                                       10d11s23
                      jdcont=ldcont                                           10d4s23
                      icsf0=1                                             11d3s22
                      jcsf0=1                                             11d3s22
                      do l=1,4                                            11d3s22
                       if(min(nlj(l),nli(l)).gt.0)then                    11d3s22
                        jad1=jvcv+ibc(jvcv+5+l)                           3d19s21
                        jad2=jad1+nlj(l)                                  3d19s21
                        iad1=ivcv+ibc(ivcv+5+l)                           3d19s21
                        iad2=iad1+nli(l)                                  3d19s21
                        itmp=ibcoff                                       11d3s22
                        ibcoff=itmp+nli(l)*nlj(l)                         11d3s22
                        call enough('hcddjk.8.2',bc,ibc)
c     nli,ncsf(iarg),ncsfmid1b,ncsf(karg),ncsfmid2,ncsf(jarg),nlj
c       ivects     iwpb1b    iprod     iwpk2    jad2
c      nli,ncsf(iarg),ncsfmid1b,ncsfmid2,ncsf(jarg),nlj
c       1    2           3         4         5       6
                        call genmatn3(nli(l),ncsfmid1b,ncsf(jarg),        11d3s22
     $                     ncsf(iarg),ivects(l),iwpb1b,ncsfmid2,        11d3s22
     $                     iwprod,iwpk2,bc(jad2),ncsf2(l,jarg),nlj(l),  11d3s22
     $                     bc(itmp),0d0,jcsf0,ncsf2(l,jarg),            11d3s22
     $                     icsf0,ncsf2(l,iarg),1d0,bc,ibc)              9d29s23
                        call den14oc(lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,       9d29s23
     $                      irefo,itmp,nli,iad1,nfdat(1,1,isb),iden14o, 9d29s23
     $                      l,bc,ibc,fmul,nlj,jad1)                     10d17s23
                        do j=1,nlj(l)
                         jm=j-1                                          11d4s22
                         jj=ibc(jad1+jm)-1                               11d4s22
                         jjden=iden1en(l)+nfdat(2,l,isb)*jj-1               2d24s21
                         do iii=1,nli(l)                                 11d4s22
                          iiim=iii-1                                     11d4s22
                          iad=itmp+iiim+nli(l)*jm                        11d4s22
                          kk=ibc(iad1+iiim)                              11d4s22
                          bc(jjden+kk)=bc(jjden+kk)+bc(iad)*xint         9d29s23
                         end do
                        end do
                        itmpj=ibcoff                                      2d24s21
                        ibcoff=itmpj+ncsf(karg)*nlj(l)                    2d24s21
                        call enough('hcddjk.  6',bc,ibc)
                        call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,bc(jad2),ncsf2(l,jarg),     2d24s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jcsf0,            10d12s23
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                        itmpc0=ibcoff                                     10d4s23
                        ibcoff=itmpc0+ncsf2(l,iarg)*nlj(l)                10d11s23
                        itmpc1=ibcoff                                      10d4s23
                        itmpc2=itmpc1+nlj(l)*nli(l)                       10d11s23
                        ibcoff=itmpc2+ncsf2(l,iarg)*nli(l)                10d11s23
                        call enough('hcddjkd12.tmpmd',bc,ibc)             10d4s23
                        call dgemm('n','n',ncsf2(l,iarg),nlj(l),          10d11s23
     $                      ncsf(karg),1d0,bc(jtmpid),ncsf(iarg),        10d4s23
     $                      bc(itmpj),ncsf(karg),0d0,bc(itmpc0),         10d4s23
     $                      ncsf2(l,iarg),'hcddjkd12.tmpmd')             10d11s23
                        do ir=0,nrootm                                    10d4s23
                         do iii=0,nlj(l)-1                                 10d4s23
                          jj=ibc(jad1+iii)-1                               10d4s23
                          lad2=ipdc1(l)-1+nfdat(2,l,isb)*(jj              10d11s23
     $                        +nfdat(2,l,isb)*ir)                       10d11s23
                          lad1=itmpc1+iii                                 10d4s23
                          do j=0,nli(l)-1                                  10d4s23
                           kk=ibc(iad1+j)                                  10d4s23
                           bc(lad1)=bc(lad2+kk)                           10d4s23
                           lad1=lad1+nlj(l)                               10d11s23
                          end do                                          10d4s23
                         end do                                           10d4s23
                         call dgemm('n','n',ncsf2(l,iarg),nli(l),         10d11s23
     $                        nlj(l),xint,bc(itmpc0),ncsf2(l,iarg),     10d11s23
     $                        bc(itmpc1),nlj(l),0d0,bc(itmpc2),         10d11s23
     $                        ncsf2(l,iarg),'hcddjkd12.tmpc2')          10d11s23
                         ltmpc2=itmpc2                                    10d4s23
                         do j=0,nli(l)-1                                  10d11s23
                          lad1=jdcont+ncsf2(l,iarg)*(j+nli(l)*ir)         10d11s23
                          do iii=0,ncsf2(l,iarg)-1                        10d11s23
                           bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)       10d4s23
                          end do                                          10d4s23
                          ltmpc2=ltmpc2+ncsf2(l,iarg)                     10d11s23
                         end do                                           10d4s23
                        end do                                            10d4s23
                        ibcoff=itmp                                       11d3s22
                       end if                                             11d3s22
                       jtmpid=jtmpid+ncsf2(l,iarg)                        10d11s23
                       jdcont=jdcont+nli(l)*ncsf2(l,iarg)*nrootu          10d11s23
                       icsf0=icsf0+ncsf2(l,iarg)                          11d3s22
                       jcsf0=jcsf0+ncsf2(l,jarg)                          11d3s22
                      end do                                              11d3s22
                      ibcoff=iprod                                       11d3s22
                      if(ipss.eq.2)go to 3                                   12d18s20
                     end if                                              11d3s22
                    else                                                 11d4s22
                     call gandc(itestc,itesto,itc,ito,nopenk,nopenip,      2d25s21
     $         karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,     2d25s21
     $                  iwpk1b,ncsfmid1b,bc,ibc)                        9d28s23
                     nnot1=nnot1b                                          2d25s21
                     nab1(1)=nab1b(2)                                      2d25s21
                     nab1(2)=nab1b(1)                                      2d25s21
                     if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                      lsa=ism(nab1(1))                                     1d8s21
                      lga=irel(nab1(1))-1                                  1d8s21
                      lsb=ism(nab1(2))                                     1d8s21
                      lgb=irel(nab1(2))-1                                  1d8s21
                      lsc=ism(nab2(1))                                     1d8s21
                      lgc=irel(nab2(1))-1                                  1d8s21
                      lsd=ism(nab2(2))                                     1d8s21
                      lgd=irel(nab2(2))-1                                  1d8s21
                      xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  9d28s23
     $                   bc,ibc)                                        9d28s23
                      fmul=1d0                                            9d29s23
                      if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2))     4d20s21
     $                   .and.                                          4d20s21
     $                max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))     1d12s21
     $                  fmul=0.5d0                                      9d29s23
                      xint=xint*fmul                                      9d29s23
                      jargo=1                                              2d23s21
                      iargo=1                                              2d24s21
                      do l=1,4                                               12d19s20
                       if(min(nli(l),nlj(l)).gt.0)then                    3d1s21
                        jad1=jvcv+ibc(jvcv+5+l)                           3d19s21
                        jad2=jad1+nlj(l)                                  3d19s21
                        iad1=ivcv+ibc(ivcv+5+l)                           3d19s21
                        iad2=iad1+nli(l)                                  3d19s21
                        itmpj=ibcoff                                      2d24s21
                        ibcoff=itmpj+ncsf(karg)*nlj(l)                    2d24s21
                        call enough('hcddjk.  9',bc,ibc)
                        call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,bc(jad2),ncsf2(l,jarg),     2d24s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                        itmpi=ibcoff                                      2d24s21
                        itmpt=itmpi+ncsf(karg)*nli(l)                     2d24s21
                        ibcoff=itmpt+ncsf(karg)*nli(l)                    2d24s21
                        call enough('hcddjk. 10',bc,ibc)
                        call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,bc(iad2),ncsf2(l,iarg),   2d24s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            9d29s23
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                        itmptt=ibcoff                                     9d29s23
                        itmpttt=itmptt+nli(l)*ncsf(karg)                  9d29s23
                        ibcoff=itmpttt+nli(l)*nlj(l)                      9d29s23
                        call enough('hcddjkd12.tmptt',bc,ibc)             9d29s23
                        do iii=0,nli(l)-1                                 9d29s23
                         do k=0,ncsf(karg)-1                              9d29s23
                          ki=itmpt+k+ncsf(karg)*iii                       9d29s23
                          ik=itmptt+iii+nli(l)*k                          9d29s23
                          bc(ik)=bc(ki)                                   9d29s23
                         end do                                           9d29s23
                        end do                                            9d29s23
                        call dgemm('n','n',nli(l),nlj(l),ncsf(karg),1d0,  9d29s23
     $                     bc(itmptt),nli(l),bc(itmpj),ncsf(karg),0d0,  9d29s23
     $                     bc(itmpttt),nli(l),'hcddjkd12.tmpttt')       9d29s23
                        call den14oc(lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,       9d29s23
     $                    irefo,itmpttt,nli,iad1,nfdat(1,1,isb),iden14o,9d29s23
     $                      l,bc,ibc,fmul,nlj,jad1)                     10d17s23
c     ?????
                        do iii=0,nlj(l)-1                                 11d2s22
                         iadj=itmpj+ncsf(karg)*iii                        11d2s22
                         jj=ibc(jad1+iii)-1                                 2d24s21
                         jjden=iden1en(l)+nfdat(2,l,isb)*jj-1               2d24s21
                         do j=0,nli(l)-1                                  11d2s22
                          iadi=itmpt+ncsf(karg)*j                         11d2s22
                          kk=ibc(iad1+j)                                    2d24s21
                          bc(jjden+kk)=bc(jjden+kk)+sum*xint              9d29s23
                         end do                                           11d2s22
                        end do                                            11d2s22
                        ibcoff=itmpj                                      2d24s21
                       end if                                             2d25s21
                       jargo=jargo+ncsf2(l,jarg)                           2d23s21
                       iargo=iargo+ncsf2(l,iarg)                           2d24s21
                      end do                                                 12d19s20
                      if(ipss.eq.2)go to 3                                   12d18s20
                     end if                                              11d4s22
                    end if                                                  12d18s20
                   end if                                                4d20s21
                  end do                                                   12d18s20
    3             continue                                                 12d18s20
                 end if                                                  1d8s21
                end if                                                   10d14s22
               end if                                                    1d8s21
               jto=ibclr(jto,norbx)                                      1d8s21
               jto=ibset(jto,norbxxx)                                    1d8s21
c                                                                       1d8s21
               gandcc=ieor(itc,jtc)                                      10d14s22
               gandco=ieor(ito,jto)                                      10d14s22
               gandcb=ior(gandcc,gandco)                                 10d21s22
               ndifb=popcnt(gandcb)                                      10d21s22
               if(ndifb.le.4)then                                        10d21s22
                ndifd=popcnt(gandcc)                                      10d14s22
                ndifs=popcnt(gandco)                                      10d14s22
                if(ndifs.eq.2.and.ndifb.eq.2)then                        10d21s22
                 do i=1,norbxxx                                           10d21s22
                  if(btest(gandco,i))then                                          10d14s22
                   if((btest(ito,i).and..not.btest(jtc,i)).or.
     $                 (btest(itc,i).and.btest(jto,i)))then                           10d14s22
                    nab4(1,1)=i
                   else                                                           10d14s22
                    nab4(2,1)=i
                   end if
                  end if                                                          10d14s22
                 end do                                                           10d14s22
                 if(isb.ne.jsb)then                                      1d9s21
                  write(6,*)('hey, what''s up with isb ne jsb ? ')       1d9s21
                  stop 'hcddjk:g'                                                   1d9s21
                 end if                                                  1d9s21
                 call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,   1d8s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $              ncsfmid1,bc,ibc)                                    9d28s23
                 call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,   2d23s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,      1d8s21
     $              ncsfmid1b,bc,ibc)                                   9d28s23
                 itmpn=ibcoff                                            2d23s21
                 nlij=0                                                  2d23s21
                 do l=1,4                                                2d23s21
                  nlij=nlij+nli(l)*nlj(l)                                2d23s21
                 end do                                                  2d23s21
                 itmpn=ibcoff                                            2d23s21
                 ibcoff=itmpn+nlij                                       2d23s21
                 call enough('hcddjk. 12',bc,ibc)
                 nok=0                                                    1d8s21
                 do i=1,norb                                             1d8s21
                  if(itest(i,1).gt.0)then                                1d8s21
                   nok=nok+1                                             1d8s21
                   itest(nok,3)=itest(i,1)                               1d8s21
                   itest(nok,2)=i                                        1d8s21
                  end if                                                 1d8s21
                 end do                                                  1d8s21
                 jtmpn=itmpn                                             2d23s21
                 jargo=1                                                  2d22s21
                 iargo=1                                                 2d23s21
                 jdcont=ldcont                                           10d4s23
                 do l=1,4                                                 1d8s21
                  if(min(nlj(l),nli(l)).gt.0)then                         2d23s21
                   iad1=ivcv+ibc(ivcv+5+l)                               3d19s21
                   iad2=iad1+nli(l)                                      3d19s21
                   jad1=jvcv+ibc(jvcv+5+l)                               3d19s21
                   jad2=jad1+nlj(l)                                      3d19s21
                   itmp1=ibcoff                                          2d23s21
                   itmp1b=itmp1+ncsfmid1*nlj(l)                          2d23s21
                   ibcoff=itmp1b+ncsfmid1*nli(l)                         2d23s21
                   call enough('hcddjk. 13',bc,ibc)
c     tmp1=k1*jvcv ncsfmid1*nlj(l)
                   if(iwpk1.lt.0)then                                    2d23s21
                    nusedi=ibc(-iwpk1)/2                                 2d23s21
                    if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1           2d23s21
                    icmp1=-iwpk1+1                                       2d23s21
                    icmp2=icmp1+ncsf(jarg)                               2d23s21
                    icmp3=icmp2+nusedi                                   2d23s21
                    call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),    2d23s21
     $                   bc(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),   2d23s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          9d28s23
                   else                                                              12d5s20
                    jwpk1=iwpk1+ncsfmid1*(jargo-1)                       2d23s21
                    call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                   1d0,bc(jwpk1),ncsfmid1,bc(jad2),ncsf2(l,jarg),
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjk.  8')
c     wpk1*jvcv
                   end if                                                            12d5s20
                   itmpdc1=ibcoff                                        10d4s23
                   itmpdc2=itmpdc1+ncsfmid1*nlj(l)                       10d4s23
                   itmpdc3=itmpdc2+ncsf2(l,iarg)*nlj(l)                  10d4s23
                   itmpdc4=itmpdc3+nli(l)*nlj(l)                         10d4s23
                   ibcoff=itmpdc4+nli(l)*ncsf2(l,iarg)                   10d4s23
                   call enough('hcddjkd12.tmpdc14',bc,ibc)               10d4s23
                   do iii=0,ncsfmid1-1                                   10d4s23
                    do j=0,nlj(l)-1                                      10d4s23
                     ij=itmp1+iii+ncsfmid1*j                             10d4s23
                     ji=itmpdc1+j+nlj(l)*iii                             10d4s23
                     bc(ji)=bc(ij)                                       10d4s23
                    end do                                               10d4s23
                   end do                                                10d4s23
c     tmp1b=k1b*ivcv ncsfmid1*nli
c     (tmp1b)T*tmp1=ivcvT*k1bT*tmp1
                   if(iwpk1b.lt.0)then                                   2d23s21
                    nusedi=ibc(-iwpk1b)/2                                2d23s21
                    if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                    icmp1=-iwpk1b+1                                      2d23s21
                    icmp2=icmp1+ncsf(iarg)                               2d23s21
                    icmp3=icmp2+nusedi                                   2d23s21
                    call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),    2d23s21
     $                   bc(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),   2d23s21
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          9d28s23
                    ixtmp=ibcoff                                         10d4s23
                    ibcoff=ixtmp+ncsfmid1*ncsf(iarg)                     10d4s23
                    call enough('hcddjkd12.xtmp',bc,ibc)                 10d4s23
                    call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),       10d4s23
     $                  bc(ixtmp),ncsfmid1,ncsf(iarg))                  10d4s23
                   else                                                              12d5s20
                    jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                    ixtmp=iwpk1b                                         10d4s23
c     jwpk1b*ivcv
                    call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                   1d0,bc(jwpk1b),ncsfmid1,bc(iad2),ncsf2(l,iarg),
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjk.  9')
                   end if                                                            12d5s20
                   jxtmp=ixtmp+ncsfmid1*(iargo-1)                        10d4s23
                   call dgemm('n','n',nlj(l),ncsf2(l,iarg),ncsfmid1,1d0, 10d4s23
     $                  bc(itmpdc1),nlj(l),bc(jxtmp),ncsfmid1,0d0,      10d4s23
     $                  bc(itmpdc2),nlj(l),'hcddjkd12.tmpdc21x')        10d4s23
c     tmpdc2 is nlj*ncsf2
                   do i=1,nok                                            10d4s23
                    lsa=ism(itest(i,2))                                    1d8s21
                    lga=irel(itest(i,2))-1                                 1d8s21
                    icolj=((lga*(lga+1))/2)+lga                            1d8s21
                    ncolj=(irefo(lsa)*(irefo(lsa)+1))/2                  10d4s23
                    ftj=1d0                                              10d4s23
                    if(itest(i,3).eq.2)then                              10d4s23
                     ftj=2d0                                             10d4s23
                     icolk=lga+irefo(lsa)*lga                            10d10s23
                     ncolk=irefo(lsa)*irefo(lsa)                         10d10s23
                     do ir=0,nrootm                                       10d4s23
                      jtmpdc3=itmpdc3                                     10d4s23
                      do iii=0,nlj(l)-1                                     10d4s23
                       jj=ibc(jad1+iii)-1                                   10d4s23
                       lad1=idck(l,l,lsa)-1+nfdat(2,l,isb)*(jj              10d4s23
     $                    +nfdat(2,l,jsb)*(icolk+ncolk*ir))             10d4s23
                       do j=0,nli(l)-1                                    10d4s23
                        kk=ibc(iad1+j)                                    10d4s23
                        bc(jtmpdc3+j)=bc(lad1+kk)                         10d4s23
                       end do                                             10d4s23
                       jtmpdc3=jtmpdc3+nli(l)                             10d4s23
                      end do                                              10d4s23
                      ftk=-1d0                                           10d10s23
                      call dgemm('n','n',nli(l),ncsf2(l,iarg),nlj(l),   10d16s23
     $                    ftk,bc(itmpdc3),nli(l),bc(itmpdc2),nlj(l),0d0,10d16s23
     $                   bc(itmpdc4),nli(l),'hcddjkd12.tmpdc4j')        10d4s23
                      lad1=jdcont+ncsf2(l,iarg)*nli(l)*ir                 10d4s23
                      do j=0,nli(l)-1                                     10d4s23
                       jtmpdc4=itmpdc4+j                                  10d4s23
                       do iii=0,ncsf2(l,iarg)-1                            10d4s23
                        bc(lad1+iii)=bc(lad1+iii)+bc(jtmpdc4)             10d4s23
                        jtmpdc4=jtmpdc4+nli(l)                            10d4s23
                       end do                                             10d4s23
                       lad1=lad1+ncsf2(l,iarg)                            10d4s23
                      end do                                              10d4s23
                     end do                                               10d4s23
                    end if                                               10d4s23
                    do ir=0,nrootm                                       10d4s23
                     jtmpdc3=itmpdc3                                     10d4s23
                     do iii=0,nlj(l)-1                                     10d4s23
                      jj=ibc(jad1+iii)-1                                   10d4s23
                      lad1=idcj(l,lsa)-1+nfdat(2,l,isb)*(jj              10d4s23
     $                    +nfdat(2,l,jsb)*(icolj+ncolj*ir))             10d4s23
                      do j=0,nli(l)-1                                    10d4s23
                       kk=ibc(iad1+j)                                    10d4s23
                       bc(jtmpdc3+j)=bc(lad1+kk)                         10d4s23
                      end do                                             10d4s23
                      jtmpdc3=jtmpdc3+nli(l)                             10d4s23
                     end do                                              10d4s23
                     call dgemm('n','n',nli(l),ncsf2(l,iarg),nlj(l),ftj, 10d4s23
     $                   bc(itmpdc3),nli(l),bc(itmpdc2),nlj(l),0d0,     10d4s23
     $                   bc(itmpdc4),nli(l),'hcddjkd12.tmpdc4j')        10d4s23
                     lad1=jdcont+ncsf2(l,iarg)*nli(l)*ir                 10d4s23
                     do j=0,nli(l)-1                                     10d4s23
                      jtmpdc4=itmpdc4+j                                  10d4s23
                      do iii=0,ncsf2(l,iarg)-1                            10d4s23
                       bc(lad1+iii)=bc(lad1+iii)+bc(jtmpdc4)             10d4s23
                       jtmpdc4=jtmpdc4+nli(l)                            10d4s23
                      end do                                             10d4s23
                      lad1=lad1+ncsf2(l,iarg)                            10d4s23
                     end do                                              10d4s23
                    end do                                               10d4s23
                   end do                                                10d4s23
                   do ir=0,nrootm                                        10d12s23
                    jtmpdc3=itmpdc3                                      10d12s23
                    do iii=0,nlj(l)-1                                     10d4s23
                     jj=ibc(jad1+iii)-1                                   10d4s23
                     lad1=idchvv(l)-1+nfdat(2,l,isb)*(jj                 10d12s23
     $                    +nfdat(2,l,jsb)*ir)                           10d12s23
                     do j=0,nli(l)-1                                     10d12s23
                      kk=ibc(iad1+j)                                     10d12s23
                      bc(jtmpdc3+j)=bc(lad1+kk)                          10d12s23
                     end do                                              10d12s23
                     jtmpdc3=jtmpdc3+nli(l)                              10d12s23
                    end do                                               10d12s23
                    call dgemm('n','n',nli(l),ncsf2(l,iarg),nlj(l),1d0,  10d12s23
     $                   bc(itmpdc3),nli(l),bc(itmpdc2),nlj(l),0d0,     10d4s23
     $                   bc(itmpdc4),nli(l),'hcddjkd12.tmpdc4j')        10d4s23
                    lad1=jdcont+ncsf2(l,iarg)*nli(l)*ir                  10d12s23
                    do j=0,nli(l)-1                                      10d12s23
                     jtmpdc4=itmpdc4+j                                   10d12s23
                     do iii=0,ncsf2(l,iarg)-1                            10d12s23
                      orig=bc(lad1+iii)
                      bc(lad1+iii)=bc(lad1+iii)+bc(jtmpdc4)              10d12s23
                      jtmpdc4=jtmpdc4+nli(l)                             10d12s23
                     end do                                              10d12s23
                     lad1=lad1+ncsf2(l,iarg)                             10d12s23
                    end do                                               10d12s23
                   end do                                                10d12s23
                   itmpt=ibcoff                                          2d23s21
                   ibcoff=itmpt+nli(l)*ncsfmid1                          2d23s21
                   call enough('hcddjk. 14',bc,ibc)
                   do iii=0,nli(l)-1                                     2d23s21
                    do j=0,ncsfmid1-1                                    2d23s21
                     ji=itmp1b+j+ncsfmid1*iii                            2d23s21
                     ij=itmpt+iii+nli(l)*j                               2d23s21
                     bc(ij)=bc(ji)                                       2d23s21
                    end do                                               2d23s21
                   end do                                                2d23s21
                   call dgemm('n','n',nli(l),nlj(l),ncsfmid1,1d0,        2d23s21
     $                  bc(itmpt),nli(l),bc(itmp1),ncsfmid1,0d0,         2d23s21
     $                  bc(jtmpn),nli(l),                                2d23s21
     d' hcddjk. 10')
                   ibcoff=itmp1                                          2d23s21
                   ktmpn=jtmpn                                           2d24s21
                   do iii=0,nlj(l)-1                                     2d24s21
                    jj=ibc(jad1+iii)-1                                   2d24s21
                    jjden=idenhvvn(l)+nfdat(2,l,isb)*jj-1                 2d24s21
                    do j=0,nli(l)-1                                      2d24s21
                     kk=ibc(iad1+j)                                      2d24s21
                     bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)               2d24s21
                    end do                                               2d24s21
                    ktmpn=ktmpn+nli(l)                                   2d24s21
                   end do                                                2d24s21
                   jtmpn=jtmpn+nli(l)*nlj(l)                             2d23s21
                  end if                                                 2d23s21
                  jargo=jargo+ncsf2(l,jarg)                              2d22s21
                  iargo=iargo+ncsf2(l,iarg)                              2d24s21
                  jdcont=jdcont+ncsf2(l,iarg)*nli(l)*nrootu              10d4s23
                 end do                                                  1d8s21
                 do i=1,nok                                              1d8s21
                  lsa=ism(itest(i,2))                                    1d8s21
                  lga=irel(itest(i,2))-1                                 1d8s21
                  icolj=((lga*(lga+1))/2)+lga                            1d8s21
                  icolk=lga+irefo(lsa)*lga                               1d8s21
                  if(itest(i,3).eq.2)then                                1d8s21
                   jtmpn=itmpn                                           2d23s21
                   do l=1,4                                              1d8s21
                    if(min(nli(l),nlj(l)).gt.0)then                      2d23s21
                     jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                     iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                     jden=idenjn(l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)    2d23s21
     $                    *icolj                                        2d23s21
                     kden=idenkn(l,l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)  2d24s21
     $                    *icolk                                        2d24s21
                     do iii=0,nlj(l)-1                                   2d23s21
                      jj=ibc(jad1+iii)-1                                 1d8s21
                      jjden=jden+nfdat(2,l,isb)*jj-1                     2d23s21
                      kkden=kden+nfdat(2,l,isb)*jj-1                     2d24s21
                      do j=0,nli(l)-1                                    2d23s21
                       kk=ibc(iad1+j)                                    2d23s21
                       bc(jjden+kk)=bc(jjden+kk)+bc(jtmpn+j)*2d0         2d23s21
                       bc(kkden+kk)=bc(kkden+kk)-bc(jtmpn+j)             2d24s21
                      end do                                             2d23s21
                      jtmpn=jtmpn+nli(l)                                 2d23s21
                     end do                                              2d23s21
                    end if                                               2d23s21
                   end do                                                1d8s21
                  else                                                   1d8s21
                   jtmpn=itmpn                                           2d23s21
                   do l=1,4                                              1d8s21
                    if(nlj(l).gt.0)then                                  1d8s21
                     jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                     iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                     jden=idenjn(l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)    2d23s21
     $                    *icolj                                        2d23s21
                     do iii=0,nlj(l)-1                                   2d23s21
                      jj=ibc(jad1+iii)-1                                 1d8s21
                      jjden=jden+nfdat(2,l,isb)*jj-1                     2d23s21
                      do j=0,nli(l)-1                                    2d23s21
                       kk=ibc(iad1+j)                                    2d23s21
                       bc(jjden+kk)=bc(jjden+kk)+bc(jtmpn+j)             2d23s21
                      end do                                             2d23s21
                      jtmpn=jtmpn+nli(l)                                 2d23s21
                     end do                                              2d23s21
                    end if                                               1d8s21
                   end do                                                1d8s21
                  end if                                                 1d8s21
                 end do                                                  1d8s21
                 ibcoff=itmpn                                            4d28s21
                 do i=1,nok                                              1d8s21
                  if(itest(i,3).eq.1)then                                  12d14s20
                   itestc=itc                                              12d14s20
                   itesto=ito                                              12d14s20
                   nopenk=nopenip                                                11d13s20
c
c     i is bra, j is ket.
c
c         I   J    K     thus  IK  KJ
c     b   n+1 n    n+1        g   g   (ck|bc)
c     k   m   m+1  m+1         ck  bc
c     c   p   p    p-1
c
c
c     anihilate common
c
                   if(btest(itestc,itest(i,2)))then                             11d13s20
                    itestc=ibclr(itestc,itest(i,2))                             11d13s20
                    itesto=ibset(itesto,itest(i,2))                             11d13s20
                    karg=iarg-1                                                11d13s20
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,2))                             11d13s20
                    karg=iarg                                                  11d13s20
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
                   nnot1=0                                               8d4s22
                   nnot2=0                                               8d4s22
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                    call gandc(itestc,itesto,itc,ito,nopenk,nopenip,     2d24s21
     $            karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b, 2d24s21
     $               iwpk1b,ncsfmid1b,bc,ibc)                           9d28s23
                    nnot1=nnot1b                                         2d25s21
                    nab1(1)=nab1b(2)                                     2d25s21
                    nab1(2)=nab1b(1)                                     2d25s21
                    call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,        12d14s20
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            9d28s23
                   end if
                   if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                    if(nab1(2).gt.nab2(1))then                           1d8s21
                     ia=nab1(1)
                     ib=nab2(2)
                    else
                     ia=nab2(2)
                     ib=nab1(1)
                    end if                                               1d8s21
                    lsa=ism(ia)                                          1d8s21
                    lga=irel(ia)-1                                       1d8s21
                    lsb=ism(ib)                                          1d8s21
                    lgb=irel(ib)-1                                       1d8s21
                    lsab=multh(lsa,lsb)                                  1d9s21
                    icol=lga+irefo(lsa)*lgb                              1d8s21
                    ncol=irefo(lsa)*irefo(lsb)                           10d11s23
                    if(lsab.ne.ijsb)then                                 1d9s21
                     write(6,*)('hey, lsab = '),lsab,(' ne ijsb = '),
     $                   ijsb
                     stop 'hcddjk:i'
                    end if                                               1d9s21
                    jargo=1                                               2d23s21
                    itmpj=ibcoff                                          2d23s21
                    jtmpj=itmpj                                           2d23s21
                    do lj=1,4                                             2d23s21
                     if(nlj(lj).gt.0)then                                 2d23s21
                      jad1=jvcv+ibc(jvcv+5+lj)                              3d19s21
                      jad2=jad1+nlj(lj)                                     3d19s21
                      ibcoff=jtmpj+ncsf(karg)*nlj(lj)                     2d23s21
                      call enough('hcddjk. 15',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                    nlj(lj),iwpb2,iwpk2,bc(jad2),ncsf2(lj,jarg),   2d23s21
     $                   bc(jtmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(lj,jarg),bc,ibc)                         9d28s23
                      jtmpj=jtmpj+ncsf(karg)*nlj(lj)                      2d23s21
                     end if                                               2d23s21
                     jargo=jargo+ncsf2(lj,jarg)                           2d23s21
                    end do                                                2d23s21
                    iargo=1                                               2d23s21
                    itmpid=ibcoff                                        10d4s23
                    ibcoff=itmpid+2*ncsf(karg)*ncsf(iarg)                10d4s23
                    call enough('hcddjkd12.tmpid',bc,ibc)                10d4s23
                    jtmpid=itmpid+ncsf(karg)*ncsf(iarg)                  10d4s23
                    call prodn(iwpb1b,iwpk1b,ncsf(karg),ncsf(iarg),      10d4s23
     $                  ncsfmid1b,bc(jtmpid),bc,ibc,1d0,0d0)            10d4s23
                    do iii=0,ncsf(iarg)-1                                10d4s23
                     do j=0,ncsf(karg)-1                                 10d4s23
                      ji=jtmpid+j+ncsf(karg)*iii                         10d4s23
                      ij=itmpid+iii+ncsf(iarg)*j                         10d4s23
                      bc(ij)=bc(ji)                                      10d4s23
                     end do                                              10d4s23
                    end do                                               10d4s23
                    ibcoff=jtmpid                                        10d4s23
                    itmpi=ibcoff                                          2d23s21
                    jtmpi=itmpi                                           2d23s21
                    do li=1,4                                             2d23s21
                     if(nli(li).gt.0)then                                 2d23s21
                      iad1=ivcv+ibc(ivcv+5+li)                              3d19s21
                      iad2=iad1+nli(li)                                     3d19s21
                      ibcoff=jtmpi+ncsf(karg)*nli(li)                     2d23s21
                      itmpt=ibcoff                                        2d23s21
                      ibcoff=itmpt+ncsf(karg)*nli(li)                     2d23s21
                      call enough('hcddjk. 16',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                    nli(li),iwpb1b,iwpk1b,bc(iad2),ncsf2(li,iarg), 2d23s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(li,iarg),bc,ibc)                         9d28s23
                      do iii=0,nli(li)-1                                  2d23s21
                       do k=0,ncsf(karg)-1                                2d23s21
                        ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                        ik=jtmpi+iii+nli(li)*k                            2d23s21
                        bc(ik)=bc(ki)                                     2d23s21
                       end do                                             2d23s21
                      end do                                              2d23s21
                      ibcoff=itmpt                                        2d23s21
                      jtmpi=jtmpi+ncsf(karg)*nli(li)                      2d23s21
                     end if                                               2d23s21
                     iargo=iargo+ncsf2(li,iarg)                            2d23s21
                    end do                                                2d23s21
                    jtmpj=itmpj                                          2d23s21
                    do lj=1,4                                             2d23s21
                     if(nlj(lj).gt.0)then                                2d24s21
                      jad1=jvcv+ibc(jvcv+5+lj)                              3d19s21
                      jtmpi=itmpi                                          2d23s21
                      jtmpid=itmpid                                      10d11s23
                      jdcont=ldcont                                      10d11s23
                      do li=1,4                                          2d24s21
                       if(nli(li).gt.0)then                              2d24s21
                        iad1=ivcv+ibc(ivcv+5+li)                         3d19s21
                        iad2=iad1+nli(li)
                        jden=idenkn(li,lj,lsa)+nfdat(2,li,isb)           2d24s21
     $                       *nfdat(2,lj,jsb)*icol                       2d24s21
                        itmpm=ibcoff                                       2d23s21
                        ibcoff=itmpm+nlj(lj)*nli(li)                     2d24s21
                        call enough('hcddjk. 17',bc,ibc)
                        call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),   2d24s21
     $                       1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),2d24s21
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' hcddjk. 11')
                        jtmpm=itmpm                                        2d23s21
                        do iii=0,nlj(lj)-1                               2d24s21
                         jj=ibc(jad1+iii)-1                                 1d8s21
                         jjden=jden+nfdat(2,li,isb)*jj-1                 2d24s21
                         do j=0,nli(li)-1                                2d24s21
                          kk=ibc(iad1+j)                                   2d23s21
                          bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)            2d23s21
                         end do                                            2d23s21
                         jtmpm=jtmpm+nli(li)                                2d23s21
                        end do                                             2d23s21
                        itmpc0=ibcoff                                     10d4s23
                        ibcoff=itmpc0+ncsf2(li,iarg)*nlj(lj)             10d10s23
                        itmpc1=ibcoff                                     10d4s23
                        itmpc2=itmpc1+nlj(lj)*nli(li)                    10d10s23
                        ibcoff=itmpc2+ncsf2(li,iarg)*nli(li)             10d10s23
                        call enough('hcddjkd12.tmpmd',bc,ibc)             10d4s23
                        call dgemm('n','n',ncsf2(li,iarg),nlj(lj),       10d10s23
     $                     ncsf(karg),1d0,bc(jtmpid),ncsf(iarg),        10d4s23
     $                     bc(jtmpj),ncsf(karg),0d0,bc(itmpc0),         10d4s23
     $                     ncsf2(li,iarg),'hcddjkd12.tmpmd')            10d10s23
                        do ir=0,nrootm                                    10d4s23
                         do iii=0,nlj(lj)-1                                 10d4s23
                          jj=ibc(jad1+iii)-1                               10d4s23
                          lad2=idck(li,lj,lsa)-1+nfdat(2,li,isb)*(jj     10d10s23
     $                      +nfdat(2,lj,jsb)*(icol+ncol*ir))            10d10s23
                          lad1=itmpc1+iii                                 10d4s23
                          do j=0,nli(li)-1                                  10d4s23
                           kk=ibc(iad1+j)                                  10d4s23
                           bc(lad1)=bc(lad2+kk)                           10d4s23
                           lad1=lad1+nlj(lj)                               10d4s23
                          end do                                          10d4s23
                         end do                                           10d4s23
                         call dgemm('n','n',ncsf2(li,iarg),nli(li),      10d10s23
     $                        nlj(lj),1d0,bc(itmpc0),ncsf2(li,iarg),    10d10s23
     $                        bc(itmpc1),nlj(lj),0d0,bc(itmpc2),        10d10s23
     $                        ncsf2(li,iarg),'hcddjkd12.tmpc2')         10d10s23
                         ltmpc2=itmpc2                                    10d4s23
                         do j=0,nli(li)-1                                  10d4s23
                          lad1=jdcont+ncsf2(li,iarg)*(j+nli(li)*ir)      10d10s23
                          do iii=0,ncsf2(li,iarg)-1                      10d10s23
                           bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)       10d4s23
                          end do                                          10d4s23
                          ltmpc2=ltmpc2+ncsf2(li,iarg)                   10d10s23
                         end do                                           10d4s23
                        end do                                            10d4s23
                        ibcoff=itmpm                                       2d23s21
                       end if                                              2d23s21
                       jtmpi=jtmpi+nli(li)*ncsf(karg)                    2d24s21
                       jtmpid=jtmpid+ncsf2(li,iarg)                      10d11s23
                       jdcont=jdcont+ncsf2(li,iarg)*nli(li)*nrootu       10d11s23
                      end do                                             2d24s21
                     end if                                              2d24s21
                     jtmpj=jtmpj+nlj(lj)*ncsf(karg)                      2d24s21
                    end do                                               2d23s21
                   end if                                                  12d14s20
                  end if                                                   12d14s20
                 end do                                                  12d18s20
                else                                                     1d8s21
                 ipssx=0                                                 10d21s22
                 nnot=0                                                  10d21s22
                 if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s2
                  nnot=4                                                          10d14s22
                  ioxx(1)=1                                                     10d17s22
                  ioxx(2)=1                                                     10d17s22
                  do i=1,norbxxx                                                     10d17s22
                   if(btest(gandcb,i))then                                         10d14s22
                    if((btest(jtc,i).and.btest(ito,i)).or.                         10d17s22
     $                  (btest(jto,i).and..not.btest(itc,i)))then                   10d14s22
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
                  do i=1,norbxxx
                   if(btest(gandcb,i))then                                         10d14s22
                    if(btest(gandcc,i).and.                                        10d14s22
     $        ((btest(itc,i).and..not.btest(jto,i)).or.                 10d14s22
     $         (btest(jtc,i).and..not.btest(ito,i))))then                     10d14s22
                     if(btest(jtc,i))iswap=1                                        10d17s22
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
                   if(btest(itc,nab4(1,2)).and.
     $                  .not.btest(itc,nab4(1,1)))nbt=1                 10d21s22
                  else                                                             10d17s22
                   nbt=0                                                           10d17s22
                   if(btest(jtc,nab4(2,2)).and.
     $                  .not.btest(jtc,nab4(2,1)))nbt=1                 10d21s22
                  end if                                                           10d17s22
                  if(nbt.ne.0)then                                                 10d17s22
                   nab4(1,1)=nab4(1,2)                                             10d17s22
                   nab4(2,1)=nab4(2,2)                                             10d17s22
                  end if                                                           10d17s22
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                  nnot=3
                  do i=1,norbxxx
                   if(btest(gandcb,i))then                                         10d14s22
                    if(btest(itc,i))then
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
                 else if(nnot.eq.4)then                                  11d1s22
                  ipssx=3                                                 12d18s20
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
                  itestc=itc                                              12d8s20
                  itesto=ito                                              12d8s20
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenip+1                                              11d13s20
                   karg=iarg-1                                             12d8s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenip-1                                              11d13s20
                   karg=iarg
                  else                                                          11d13s20
                   write(6,*)('bit not set for nab4(1,1) ='),nab4(1,iu1)       11d27s20
                   stop 'nab4(1,1)'                                           11d27s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                   karg=karg+1                                                  11d13s20
                  else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                   write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                  nab4(2,iu2)                                       12d18s20
                   stop 'nab4(2,1)'                                           11d27s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  nqq=karg+mdon-1                                        4d20s21
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                    4d20s21
                   call gandc(itestc,itesto,itc,ito,nopenk,nopenip,       2d23s21
     $         karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,    2d23s21
     $                iwpk1b,ncsfmid1b,bc,ibc)                          9d28s23
                   nnot1=nnot1b                                           2d25s21
                   nab1(1)=nab1b(2)                                       2d25s21
                   nab1(2)=nab1b(1)                                       2d25s21
                   call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,       1d8s21
     $         karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d18s20
     $         ncsfmid2,bc,ibc)                                         9d28s23
                   if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                    if(nab1(2).gt.norb)then                               1d8s21
                     if(nab2(1).lt.nab2(1))then                           1d8s21
                      ia=nab2(2)
                      ib=nab1(1)                                          1d8s21
                     else                                                 1d8s21
                      ia=nab1(1)                                          1d8s21
                      ib=nab2(2)                                          1d8s21
                     end if                                               1d8s21
                     lsa=ism(ia)                                           1d8s21
                     lsb=ism(ib)                                           1d8s21
                     lga=irel(ia)-1                                        1d8s21
                     lgb=irel(ib)-1                                        1d8s21
                     icol=lga+irefo(lsa)*lgb                              1d8s21
                     ncol=irefo(lsa)*irefo(lsb)                          10d10s23
                    else                                                  1d8s21
                     ia=min(nab1(1),nab1(2))                              1d8s21
                     ib=max(nab1(1),nab1(2))                              1d8s21
                     lsa=ism(ia)                                           1d8s21
                     lsb=ism(ib)                                           1d8s21
                     lga=irel(ia)-1                                        1d8s21
                     lgb=irel(ib)-1                                        1d8s21
                     if(lsa.eq.lsb)then                                   1d8s21
                      icol=((lgb*(lgb+1))/2)+lga                          1d10s21
                      ncol=(irefo(lsa)*(irefo(lsa)+1))/2                 10d4s23
                     else                                                 1d8s21
                      icol=lga+irefo(lsa)*lgb                             1d9s21
                      ncol=irefo(lsa)*irefo(lsb)                         10d4s23
                     end if                                               1d8s21
                    end if                                                1d8s21
                    lsab=multh(lsa,lsb)                                   1d9s21
                    jargo=1                                               2d23s21
                    itmpj=ibcoff                                          2d23s21
                    jtmpj=itmpj                                           2d23s21
                    do lj=1,4                                             2d23s21
                     if(nlj(lj).gt.0)then                                 2d23s21
                      jad1=jvcv+ibc(jvcv+5+lj)                           3d19s21
                      jad2=jad1+nlj(lj)                                     3d19s21
                      ibcoff=jtmpj+ncsf(karg)*nlj(lj)                     2d23s21
                      call enough('hcddjk. 18',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(lj),iwpb2,iwpk2,bc(jad2),ncsf2(lj,jarg),   2d23s21
     $                   bc(jtmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(lj,jarg),bc,ibc)                         9d28s23
                      jtmpj=jtmpj+ncsf(karg)*nlj(lj)                      2d23s21
                     end if                                               2d23s21
                     jargo=jargo+ncsf2(lj,jarg)                           2d23s21
                    end do                                                2d23s21
                    iargo=1                                               2d23s21
                    itmpjd=ibcoff                                       11d13s23
                    itmpid=itmpjd+ncsf(karg)*ncsf(jarg)                 11d13s23
                    ibcoff=itmpid+2*ncsf(karg)*ncsf(iarg)                10d4s23
                    call enough('hcddjkd12.tmpid',bc,ibc)                10d4s23
                    call prodn(iwpb2,iwpk2,ncsf(karg),ncsf(jarg),       11d13s23
     $                  ncsfmid2,bc(itmpjd),bc,ibc,1d0,0d0)             11d13s23
                    jtmpid=itmpid+ncsf(karg)*ncsf(iarg)                  10d4s23
                    call prodn(iwpb1b,iwpk1b,ncsf(karg),ncsf(iarg),      10d4s23
     $                  ncsfmid1b,bc(jtmpid),bc,ibc,1d0,0d0)            10d4s23
                    do iii=0,ncsf(iarg)-1                                10d4s23
                     do j=0,ncsf(karg)-1                                 10d4s23
                      ji=jtmpid+j+ncsf(karg)*iii                         10d4s23
                      ij=itmpid+iii+ncsf(iarg)*j                         10d4s23
                      bc(ij)=bc(ji)                                      10d4s23
                     end do                                              10d4s23
                    end do                                               10d4s23
                    ibcoff=jtmpid                                        10d4s23
                    itmpi=ibcoff                                          2d23s21
c     for dcont.
c     dcont is like iad2, i.e. ncsf2(li,iarg),nli(li)
c     here we multiply iwpb1b and iwpk1b to yield karg,nli(li).
c     i.e. in total we form
c     (b1b*k1b*ivl)T*(b2*k2*jvl')
c     =(ivl)T(b1b*k1b)T*(b2*k2*jvl')
                    jtmpi=itmpi                                          10d4s23
                    jtmpid=itmpid                                        10d4s23
                    do li=1,4                                             2d23s21
                     if(nli(li).gt.0)then                                 2d23s21
                      iad1=ivcv+ibc(ivcv+5+li)                              3d19s21
                      iad2=iad1+nli(li)                                     3d19s21
                      ibcoff=jtmpi+ncsf(karg)*nli(li)                     2d23s21
                      itmpt=ibcoff                                        2d23s21
                      ibcoff=itmpt+ncsf(karg)*nli(li)                     2d23s21
                      call enough('hcddjk. 19',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(li),iwpb1b,iwpk1b,bc(iad2),ncsf2(li,iarg), 2d23s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(li,iarg),bc,ibc)                         9d28s23
                      do iii=0,nli(li)-1                                  2d23s21
                       do k=0,ncsf(karg)-1                                2d23s21
                        ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                        ik=jtmpi+iii+nli(li)*k                            2d23s21
                        bc(ik)=bc(ki)                                     2d23s21
                       end do                                             2d23s21
                      end do                                              2d23s21
                      ibcoff=itmpt                                        2d23s21
                      jtmpi=jtmpi+ncsf(karg)*nli(li)                      2d23s21
                     end if                                               2d23s21
                     jtmpid=jtmpid+ncsf2(li,iarg)                        10d4s23
                     iargo=iargo+ncsf2(li,iarg)                            2d23s21
                    end do                                                2d23s21
                    if(nab1(2).gt.norb)then                               2d23s21
                     jtmpj=itmpj                                          2d23s21
                     jtmpjd=itmpjd                                      11d13s23
                     kkdcont=kdcont                                     11d13s23
                     do lj=1,4                                           2d24s21
                      if(nlj(lj).gt.0)then                               2d24s21
                       jad1=jvcv+ibc(jvcv+5+lj)                              3d19s21
                       jtmpi=itmpi                                       2d24s21
                       jtmpid=itmpid                                     10d10s23
                       jdcont=ldcont                                       10d4s23
                       do li=1,4                                         2d24s21
                        if(nli(li).gt.0)then                             2d24s21
                         iad1=ivcv+ibc(ivcv+5+li)                        3d19s21
                         iad2=iad1+nli(li)                               10d10s23
                         jden=idenkn(li,lj,lsa)+nfdat(2,li,isb)          2d24s21
     $                        *nfdat(2,lj,jsb)*icol                      2d24s21
                         itmpm=ibcoff                                       2d23s21
                         ibcoff=itmpm+nlj(lj)*nli(li)                    2d24s21
                         call enough('hcddjk. 20',bc,ibc)
                         call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),  2d24s21
     $                       1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),2d24s21
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' hcddjk. 12')
                         jtmpm=itmpm                                        2d23s21
                         do iii=0,nlj(lj)-1                                  2d23s21
                          jj=ibc(jad1+iii)-1                                 1d8s21
                          jjden=jden+nfdat(2,li,isb)*jj-1                2d24s21
                          do j=0,nli(li)-1                               2d24s21
                           kk=ibc(iad1+j)                                   2d23s21
                           bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)            2d23s21
                          end do                                            2d23s21
                          jtmpm=jtmpm+nli(li)                            2d24s21
                         end do                                             2d23s21
                         itmpd0=ibcoff                                    11d13s23
                         itmpd1=itmpd0+ncsf2(lj,jarg)*nli(li)           11d13s23
                         itmpd2=itmpd1+nlj(lj)*nli(li)                  11d13s23
                         ibcoff=itmpd2+ncsf2(lj,jarg)*nlj(lj)           11d13s23
                         itmpc0=ibcoff                                     10d4s23
                         ibcoff=itmpc0+ncsf2(li,iarg)*nlj(lj)            10d10s23
                         itmpc1=ibcoff                                     10d4s23
                         itmpc2=itmpc1+nlj(lj)*nli(li)                   10d10s23
                         ibcoff=itmpc2+ncsf2(li,iarg)*nli(li)            10d10s23
                         call enough('hcddjkd12.tmpmd',bc,ibc)             10d4s23
                         call dgemm('n','n',ncsf2(li,iarg),nlj(lj),      10d10s23
     $                     ncsf(karg),1d0,bc(jtmpid),ncsf(iarg),        10d4s23
     $                     bc(jtmpj),ncsf(karg),0d0,bc(itmpc0),         10d4s23
     $                     ncsf2(li,iarg),'hcddjkd12.tmpmd')            10d10s23
                         call dgemm('n','n',nli(li),ncsf2(lj,jarg),     11d13s23
     $                      ncsf(karg),1d0,bc(jtmpi),nli(li),bc(jtmpjd),11d13s23
     $                      ncsf(karg),0d0,bc(itmpd0),nli(li),          11d13s23
     $                      'hcddjkd12.tmpmdd')                         11d13s23
                         do ir=0,nrootm                                    10d4s23
                          do iii=0,nlj(lj)-1                                 10d4s23
                           jj=ibc(jad1+iii)-1                               10d4s23
                           lad2=idck(li,lj,lsa)-1+nfdat(2,li,isb)*(jj    10d10s23
     $                      +nfdat(2,lj,jsb)*(icol+ncol*ir))            10d10s23
                           lad1=itmpc1+iii                                 10d4s23
                           do j=0,nli(li)-1                                  10d4s23
                            kk=ibc(iad1+j)                                  10d4s23
                            bc(lad1)=bc(lad2+kk)                           10d4s23
                            lad1=lad1+nlj(lj)                               10d4s23
                           end do                                          10d4s23
                          end do                                           10d4s23
                          call dgemm('n','n',ncsf2(li,iarg),nli(li),     10d10s23
     $                        nlj(lj),0.5d0,bc(itmpc0),ncsf2(li,iarg),  11d13s23
     $                        bc(itmpc1),nlj(lj),0d0,bc(itmpc2),        10d10s23
     $                        ncsf2(li,iarg),'hcddjkd12.tmpc2')         10d10s23
                          ltmpc2=itmpc2                                    10d4s23
                          do j=0,nli(li)-1                                  10d4s23
                           lad1=jdcont+ncsf2(li,iarg)*(j+nli(li)*ir)     10d10s23
                           do iii=0,ncsf2(li,iarg)-1                     10d10s23
                            bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)       10d4s23
                           end do                                          10d4s23
                           ltmpc2=ltmpc2+ncsf2(li,iarg)                  10d10s23
                          end do                                           10d4s23
                          do iii=0,nli(li)-1                               11d13s23
                           jj=ibc(iad1+iii)-1                             11d13s23
                           lad2=idck(li,lj,lsa)+jj+nfdat(2,li,isb)*     11d13s23
     $                      nfdat(2,lj,jsb)*(icol+ncol*ir)              11d13s23
                           lad1=itmpd1+iii*nlj(lj)                      11d13s23
                           do j=0,nlj(lj)-1                             11d13s23
                            kk=ibc(jad1+j)-1                              11d13s23
                            bc(lad1+j)=bc(lad2+kk*nfdat(2,li,isb))      11d13s23
                           end do                                          10d4s23
                          end do                                           10d4s23
                          call dgemm('n','n',nlj(lj),ncsf2(lj,jarg),    11d13s23
     $                     nli(li),0.5d0,bc(itmpd1),nlj(lj),bc(itmpd0), 11d13s23
     $                      nli(li),0d0,bc(itmpd2),nlj(lj),             11d13s23
     $                      'hcddjkd12.tmpc2')                          10d4s23
                          do j=0,nlj(lj)-1                                 11d13s23
                           lad1=kkdcont+ncsf2(lj,jarg)*(j+nlj(lj)*ir)   11d13s23
                           ltmpd2=itmpd2+j                                11d13s23
                           do iii=0,ncsf2(lj,jarg)-1                    11d13s23
                            bc(lad1+iii)=bc(lad1+iii)+bc(ltmpd2)          11d13s23
                            ltmpd2=ltmpd2+nlj(lj)                       11d13s23
                           end do                                          10d4s23
                          end do                                           10d4s23
                         end do                                            10d4s23
                         ibcoff=itmpm                                       2d23s21
                        end if                                              2d23s21
                        jtmpi=jtmpi+nli(li)*ncsf(karg)                   10d10s23
                        jtmpid=jtmpid+ncsf2(li,iarg)                     10d10s23
                        jdcont=jdcont+ncsf2(li,iarg)*nli(li)*nrootu      10d10s23
                       end do                                            2d24s21
                      end if                                             2d24s21
                      jtmpj=jtmpj+nlj(lj)*ncsf(karg)                       2d23s21
                      jtmpjd=jtmpjd+ncsf(karg)*ncsf2(lj,jarg)           11d13s23
                      kkdcont=kkdcont+ncsf2(lj,jarg)*nlj(lj)*nrootu     11d13s23
                     end do                                               2d23s21
                    else                                                  2d23s21
                     jtmpj=itmpj                                          2d23s21
                     jtmpi=itmpi                                          2d23s21
                     jtmpid=itmpid                                           2d23s21
                     jdcont=ldcont                                       10d4s23
                     jtmpjd=itmpjd                                      11d13s23
                     kkdcont=kdcont                                     11d13s23
                     do l=1,4                                             2d23s21
                      if(min(nli(l),nlj(l)).gt.0)then                     2d23s21
                       jad1=jvcv+ibc(jvcv+5+l)                              3d19s21
                       iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                       iad2=iad1+nli(l)                                  10d4s23
                       jden=idenjn(l,lsa)+nfdat(2,l,isb)*nfdat(2,l,jsb)   2d23s21
     $                    *icol                                         2d23s21
                       itmpm=ibcoff                                       2d23s21
                       ibcoff=itmpm+nlj(l)*nli(l)                         2d23s21
                       call enough('hcddjk. 21',bc,ibc)
c     tmpi is nli*ncsf(karg)
c     tmpj is ncsf(karg)*nlj
                       call dgemm('n','n',nli(l),nlj(l),ncsf(karg),1d0,   2d23s21
     $                    bc(jtmpi),nli(l),bc(jtmpj),ncsf(karg),0d0,    2d23s21
     $                    bc(itmpm),nli(l),                             2d23s21
     d' hcddjk. 13')
                       jtmpm=itmpm                                        2d23s21
                       do iii=0,nlj(l)-1                                  2d23s21
                        jj=ibc(jad1+iii)-1                                 1d8s21
                        jjden=jden+nfdat(2,l,isb)*jj-1                    2d23s21
                        do j=0,nli(l)-1                                   2d23s21
                         kk=ibc(iad1+j)                                   2d23s21
                         bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)            2d23s21
                        end do                                            2d23s21
                        jtmpm=jtmpm+nli(l)                                2d23s21
                       end do                                             2d23s21
                       itmpd0=ibcoff                                    11d13s23
                       itmpd1=itmpd0+ncsf2(l,jarg)*nli(l)               11d13s23
                       itmpd2=itmpd1+nlj(l)*nli(l)                      11d13s23
                       ibcoff=itmpd2+ncsf2(l,jarg)*nlj(l)               11d13s23
                       itmpc0=ibcoff                                     10d4s23
                       ibcoff=itmpc0+ncsf2(l,iarg)*nlj(l)                 10d4s23
                       itmpc1=ibcoff                                     10d4s23
                       itmpc2=itmpc1+nlj(l)*nli(l)                       10d4s23
                       ibcoff=itmpc2+ncsf2(l,iarg)*nli(l)                10d4s23
                       call enough('hcddjkd12.tmpmd',bc,ibc)             10d4s23
                       call dgemm('n','n',ncsf2(l,iarg),nlj(l),          10d4s23
     $                     ncsf(karg),1d0,bc(jtmpid),ncsf(iarg),        10d4s23
     $                     bc(jtmpj),ncsf(karg),0d0,bc(itmpc0),         10d4s23
     $                     ncsf2(l,iarg),'hcddjkd12.tmpmd')             10d4s23
                       call dgemm('n','n',nli(l),ncsf2(l,jarg),         11d13s23
     $                      ncsf(karg),1d0,bc(jtmpi),nli(l),bc(jtmpjd), 11d13s23
     $                      ncsf(karg),0d0,bc(itmpd0),nli(l),           11d13s23
     $                      'hcddjkd12.tmpmdd')                         11d13s23
                       do ir=0,nrootm                                    10d4s23
                        do iii=0,nlj(l)-1                                 10d4s23
                         jj=ibc(jad1+iii)-1                               10d4s23
                         lad2=idcj(l,lsa)-1+nfdat(2,l,isb)*(jj           10d4s23
     $                      +nfdat(2,l,jsb)*(icol+ncol*ir))             10d4s23
                         lad1=itmpc1+iii                                 10d4s23
                         do j=0,nli(l)-1                                  10d4s23
                          kk=ibc(iad1+j)                                  10d4s23
                          bc(lad1)=bc(lad2+kk)                           10d4s23
                          lad1=lad1+nlj(l)                               10d4s23
                         end do                                          10d4s23
                        end do                                           10d4s23
                        call dgemm('n','n',ncsf2(l,iarg),nli(l),nlj(l),  10d4s23
     $                      0.5d0,bc(itmpc0),ncsf2(l,iarg),bc(itmpc1),  11d13s23
     $                      nlj(l),0d0,bc(itmpc2),ncsf2(l,iarg),        10d4s23
     $                      'hcddjkd12.tmpc2')                          10d4s23
                        ltmpc2=itmpc2                                    10d4s23
                        do j=0,nli(l)-1                                  10d4s23
                         lad1=jdcont+ncsf2(l,iarg)*(j+nli(l)*ir)         10d10s23
                         do iii=0,ncsf2(l,iarg)-1                        10d4s23
                          bc(lad1+iii)=bc(lad1+iii)+bc(ltmpc2+iii)       10d4s23
                         end do                                          10d4s23
                         ltmpc2=ltmpc2+ncsf2(l,iarg)                     10d4s23
                        end do                                           10d4s23
                        do iii=0,nli(l)-1                               11d13s23
                         jj=ibc(iad1+iii)-1                             11d13s23
                         lad2=idcj(l,lsa)+jj+nfdat(2,l,isb)*            11d13s23
     $                      nfdat(2,l,jsb)*(icol+ncol*ir)               11d13s23
                         lad1=itmpd1+iii*nlj(l)                         11d13s23
                         do j=0,nlj(l)-1                                  10d4s23
                          kk=ibc(jad1+j)-1                              11d13s23
                          bc(lad1+j)=bc(lad2+kk*nfdat(2,l,isb))         11d13s23
                         end do                                          10d4s23
                        end do                                           10d4s23
                        call dgemm('n','n',nlj(l),ncsf2(l,jarg),nli(l), 11d13s23
     $                      0.5d0,bc(itmpd1),nlj(l),bc(itmpd0),         11d13s23
     $                      nli(l),0d0,bc(itmpd2),nlj(l),               11d13s23
     $                      'hcddjkd12.tmpc2')                          10d4s23
                        do j=0,nlj(l)-1                                 11d13s23
                         lad1=kkdcont+ncsf2(l,jarg)*(j+nlj(l)*ir)       11d13s23
                         ltmpd2=itmpd2+j                                11d13s23
                         do iii=0,ncsf2(l,jarg)-1                       11d13s23
                          bc(lad1+iii)=bc(lad1+iii)+bc(ltmpd2)          11d13s23
                          ltmpd2=ltmpd2+nlj(l)                          11d13s23
                         end do                                          10d4s23
                        end do                                           10d4s23
                       end do                                            10d4s23
                       ibcoff=itmpm                                       2d23s21
                      end if                                              2d23s21
                      jtmpj=jtmpj+nlj(l)*ncsf(karg)                       2d23s21
                      jtmpi=jtmpi+nli(l)*ncsf(karg)                       2d23s21
                      jtmpid=jtmpid+ncsf2(l,iarg)                        10d4s23
                      jtmpjd=jtmpjd+ncsf(karg)*ncsf2(l,jarg)            11d13s23
                      jdcont=jdcont+ncsf2(l,iarg)*nli(l)*nrootu          10d4s23
                      kkdcont=kkdcont+ncsf2(l,jarg)*nlj(l)*nrootu       11d13s23
                     end do                                               2d23s21
                    end if                                                2d23s21
                    if(ipss.eq.2)go to 4                                   12d18s20
                   end if                                                  12d18s20
                  end if                                                 4d20s21
                 end do                                                   12d18s20
    4            continue                                                 12d18s20
                end if                                                   1d8s21
               end if                                                    10d14s22
               jvcv=jvcv+nspacej                                         3d19s21
               kdcont=kdcont+mdcont                                     11d13s23
              end do                                                       1d8s21
             else                                                       11d13s23
              if(nff22(nclojp,1,jsb).gt.0)then                            1d14s21
               jarg=nclojp-mdon                                         11d14s23
               jvcv=nfdat(5,1,jsb)+nff22(nclojp,2,jsb)                      1d8s21
               do jf=1,nff22(nclojp,1,jsb)                                  1d8s21
                nspacej=ibc(jvcv+1)                                       3d19s21
                mdcont=0                                                  11d13s23
                do l=1,4                                                  11d13s23
                 nlj(l)=ibc(jvcv+1+l)                                     11d13s23
                 mdcont=mdcont+nlj(l)*ncsf2(l,jarg)                       11d13s23
                end do                                                    11d13s23
                jvcv=jvcv+nspacej                                         11d13s23
                kdcont=kdcont+mdcont                                      11d13s23
               end do                                                     11d13s23
              end if                                                      11d13s23
             end if                                                     3d1s21
             idoit=idoit+1                                               3d1s21
            end if                                                      1d8s21
           end do                                                       1d8s21
           ibcoff=ibcvects                                              11d3s22
           ivcv=ivcv+nspacei                                            3d19s21
           ldcont=ldcont+ndcont                                         10d4s23
          end do                                                        1d8s21
          ibcoff=ibcst                                                  1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call dws_gsumf(bc(ibc0),nwds)                                   2d23s21
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    1d8s21
          do ll=1,4                                                     1d8s21
           if(nfdat(2,ll,jsb).gt.0)then                                 1d8s21
            nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                          1d8s21
            if(l.eq.ll.and.nll.gt.0)then                                2d24s21
             iden1ef(l)=iden1en(l)                                      2d25s21
             idenhvvf(l)=idenhvvn(l)                                    2d25s21
            end if                                                      1d11s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(lsb.ge.lsa.and.l.eq.ll.and.                             2d23s21
     $            min(irefo(lsa),irefo(lsb),nll).gt.0)then              2d23s21
              idenjf(l,lsa)=idenjn(l,lsa)                               2d25s21
             end if                                                     1d9s21
             nn=irefo(lsa)*irefo(lsb)                                   1d9s21
             if(min(nn,nll).gt.0)then                                   2d24s21
              idenkf(l,ll,lsa)=idenkn(l,ll,lsa)                         2d25s21
             end if                                                     2d24s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        nall=0                                                          5d10s21
        do ll=1,4                                                       1d8s21
         if(nfdat(2,ll,jsb).gt.0)then                                   1d8s21
          do l=1,4                                                      1d8s21
           if(nfdat(2,l,isb).gt.0)then                                  1d8s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             nokk(l,ll,lsa)=0                                           1d15s21
             if(min(irefo(lsa),irefo(lsb)).gt.0)then                    1d13s21
              nnn=nfdat(2,l,isb)*nfdat(2,ll,jsb)                        1d9s21
              if(lsb.ge.lsa.and.l.eq.ll)then                            1d13s21
               if(lsa.eq.lsb)then                                        1d8s21
                nn=(irefo(lsa)*(irefo(lsa)+1))/2                         1d8s21
               else                                                      1d8s21
                nn=irefo(lsa)*irefo(lsb)                                 1d8s21
               end if                                                    1d8s21
               jden=idenjf(l,lsa)                                        1d11s21
               jdent=idenjf(l,lsa)                                       1d11s21
               nokj(l,lsa)=0                                             1d12s21
               kden=idenjf(l,lsa)                                        1d12s21
               do icol=0,nn-1                                            1d9s21
                sz=0d0                                                    1d8s21
                do i=0,nnn-1                                             1d9s21
                 sz=sz+bc(jden+i)**2                                    1d13s21
                end do                                                    1d8s21
                sz=sqrt(sz/dfloat(nnn))                                  1d9s21
                if(sz.gt.1d-10)then                                       1d8s21
                 ibc(ndenjf(l,lsa)+nokj(l,lsa))=icol                     1d12s21
                 nokj(l,lsa)=nokj(l,lsa)+1                               1d12s21
                 do i=0,nnn-1                                            1d12s21
                  bc(kden+i)=bc(jden+i)                                  1d12s21
                 end do                                                  1d12s21
                 kden=kden+nnn                                           1d12s21
                 nall=nall+1                                            5d10s21
                end if                                                   1d13s21
                jden=jden+nnn                                            1d9s21
               end do                                                    1d9s21
              end if                                                    1d8s21
              nn=irefo(lsa)*irefo(lsb)                                  1d8s21
              jden=idenkf(l,ll,lsa)                                     1d9s21
              jdent=idenkf(ll,l,lsb)                                    1d10s21
              kden=idenkf(l,ll,lsa)                                     1d12s21
              do icol=0,nn-1                                            1d9s21
               ib=icol/irefo(lsa)
               ia=icol-irefo(lsa)*ib
               icolt=ib+irefo(lsb)*ia                                   1d10s21
               jdent=idenkf(ll,l,lsb)+nnn*icolt                         1d10s21
               sz=0d0                                                    1d8s21
               do i=0,nnn-1                                             1d9s21
                sz=sz+bc(jden+i)**2                                     1d9s21
               end do                                                    1d8s21
               sz=sqrt(sz/dfloat(nnn))                                  1d9s21
               if(sz.gt.1d-10)then                                       1d8s21
                ibc(ndenkf(l,ll,lsa)+nokk(l,ll,lsa))=icol               1d12s21
                nokk(l,ll,lsa)=nokk(l,ll,lsa)+1                         1d12s21
                do i=0,nnn-1                                            1d12s21
                 bc(kden+i)=bc(jden+i)                                  1d12s21
                end do                                                  1d12s21
                kden=kden+nnn                                           1d12s21
                nall=nall+1                                             5d10s21
               end if                                                    1d8s21
               jden=jden+nnn                                            1d9s21
              end do                                                    1d9s21
             end if                                                     1d8s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        if(nall.gt.0)then                                               5d10s21
         ioff=ioffdnon                                                   1d12s21
         do isbv1=1,nsymb                                                1d12s21
          isbv2=multh(isbv1,isbv12)                                      1d12s21
          if(isbv2.ge.isbv1)then                                         1d12s21
           if(isbv12.eq.1)then                                           1d12s21
            isw=0                                                        1d12s21
            mvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        1d12s21
            ioffp=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)
           else                                                          1d12s21
            isw=1                                                        1d12s21
            mvv=nvirt(isbv1)*nvirt(isbv2)                                1d12s21
            ioffp=ioff
           end if                                                        1d12s21
           mmvv=mvv*nrootu                                               1d12s21
           joff=joffdnon                                                 1d12s21
           do jsbv1=1,nsymb                                              1d12s21
            jsbv2=multh(jsbv1,jsbv12)                                    1d12s21
            if(jsbv2.ge.jsbv1)then                                       1d12s21
             if(jsbv12.eq.1)then                                         1d12s21
              jsw=0                                                      1d12s21
              nvv=(nvirt(jsbv1)*(nvirt(jsbv2)-1))/2                      1d12s21
              joffp=joff+nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)              1d13s21
             else                                                        1d12s21
              jsw=1                                                      1d12s21
              nvv=nvirt(jsbv1)*nvirt(jsbv2)                              1d12s21
              joffp=joff                                                 1d13s21
             end if                                                      1d12s21
             nnvv=nvv*nrootu                                             1d12s21
c     Gvvrj =djiab (ab|vv") Vvv"ri
c     Gvv'rj=djiab (ab|vv") Vv'v"ri
c     Gvv'rj=djiab (ab|vv") Vv"v'ri
c     Gv'vrj=djiab (ab|vv") Vv'v"ri
c     Gv'vrj=djiab (ab|vv") Vv"v'ri
             if(jsbv2.eq.isbv1.or.jsbv2.eq.isbv2.or.jsbv1.eq.isbv1.or.   1d13s21
     $         jsbv1.eq.isbv2)then                                      1d13s21
              do imatch=1,4                                              1d13s21
               kase=0                                                    1d14s21
               if(jsbv2.eq.isbv1.and.imatch.eq.1)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv2                                                1d13s21
                itf=0                                                     1d13s21
                jtf=0                                                     1d13s21
                ibl=1                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=1
               end if                                                    1d13s21
               if(jsbv2.eq.isbv2.and.imatch.eq.2)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=0                                                    1d15s21
                ibl=1                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=2
               end if                                                    1d13s21
               if(jsbv1.eq.isbv1.and.imatch.eq.3)then                    1d14s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv2                                                1d13s21
                jtf=1                                                     1d13s21
                itf=0                                                    1d15s21
                ibl=2                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=3
               end if                                                    1d13s21
               if(jsbv1.eq.isbv2.and.imatch.eq.4)then                    1d13s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=1                                                     1d13s21
                ibl=2                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=4                                                    1d13s21
               end if                                                     1d13s21
               iioffp=ioffp                                               1d12s21
               tfi=1d0                                                   1d14s21
               do li=1,4                                                 1d14s21
                if(min(kase,nfdat(2,li,isb)).gt.0)then                   1d14s21
                 jjoffp=joffp                                            1d14s21
                 tfj=1d0                                                 1d14s21
                 do lj=1,4                                               1d14s21
                  if(nfdat(2,lj,jsb).gt.0)then                           1d14s21
                   tf=tfi*tfj
                   call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0.and.isb.ne.jsb)then   3d19s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdat(2,li,isb)*nfdat(2,lj,jsb)                   1d14s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    call enough('hcddjk. 22',bc,ibc)
                    fact=0d0                                                1d13s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then             1d14s21
                      i2eu=invk1(1,lsb,lsa,isbl,1)                              1d13s21
                      if(nokk(li,lj,lsa).gt.0)then                              1d12s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjk. 23',bc,ibc)
                       kint=kmats(i2eu)                                     1d12s21
                       do j=0,nokk(li,lj,lsa)-1                                1d12s21
                        jj=ibc(ndenkf(li,lj,lsa)+j)                            1d12s21
                        jb=jj/irefo(lsa)                                   1d12s21
                        ja=jj-irefo(lsa)*jb                                1d12s21
                        kk=jb+irefo(lsb)*ja                                1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=kint+i+nhere*kk                                1d12s21
                         ji=itmpi+j+nokk(li,lj,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                              1d12s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0,   1d14s21
     $               bc(idenkf(li,lj,lsa)),mn,bc(itmpi),nokk(li,lj,lsa),1d14s21
     $                    fact,bc(intden),mn,                           1d14s21
     d' hcddjk. 14')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
                     nrow=nfdat(2,lj,jsb)*nvirt(isbl)                        1d13s21
                     nmul=nfdat(2,li,isb)*nhere2
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjk. 24',bc,ibc)
                     do i=itmp1,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     i10=i1s                                                1d12s21
                     i1n=nvirt(isbl)                                        1d13s21
                     jntden=intden                                          1d12s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s                                            1d12s21
                      if(i2.eq.i2e)i1n=i1e                                  1d12s21
                      do i1=i10,i1n                                         1d12s21
                       i1m=i1-1                                             1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                              1d12s21
                        do i=0,nfdat(2,li,isb)-1                             1d12s21
                         iad=itmp1+j+nfdat(2,lj,jsb)*(i1m+nvirt(isbl)*(  1d14s21
     $                       i2n+nhere2*i))                             1d14s21
                         bc(iad)=bc(jntden+i)                               1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdat(2,li,isb)                        1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjk. 25',bc,ibc)
                     nx=nhere2*nfdat(2,li,isb)*nrootu                        1d12s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibr.eq.1)then                                       1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do i=0,nfdat(2,li,isb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                         do iv1=0,ntop                                       1d12s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv1=i2-1                                              1d12s21
                       nbot=iv1+1                                            1d12s21
                       nbot=nbot+isw*(0-nbot)                               1d13s21
                       do i=0,nfdat(2,li,isb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                         do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(li.eq.1.and.isbv12.eq.1)then                     1d15s21
                      iioff=iioffp-nvirt(isbv1)*nrootu*nfdat(2,1,isb)    1d15s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do i=0,nfdat(2,1,isb)-1                           1d15s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,1,isb)*(ir    1d15s21
     $                       +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     nwant=0                                              10d3s23
                     do lsa=1,nsymb                                       10d3s23
                      lsb=multh(lsa,ijsb)                                 10d3s23
                      if(min(nokk(li,lj,lsa),irefo(lsa),irefo(lsb)).gt.0
     $                     )nwant=nwant+1
                      nokkk(lsa)=nokk(li,lj,lsa)                              10d3s23
                      iokk(lsa)=ndenkf(li,lj,lsa)
                      iiden(lsa)=idenkf(li,lj,lsa)
                     end do                                               10d3s23
                     tff=tf*phs                                          1d14s21
                     if(nwant.gt.0)then                                   10d3s23
                      call hcddjkd123(nrow,nmul,bc,ibc,ibl,nvirt,jsw,
     $                    isbl,isbvc,isbr,nrootm,jjoffp,nvv,
     $                    nfdat(2,lj,jsb),vd,nnvv,
     $                    lj.eq.1.and.jsbv12.eq.1,nfdat(2,li,isb),      10d3s23
     $                  nhere2,tff,i1s,i1e,i2s,i2e,1,nokkk,iokk,ijsb,    10d3s23
     $                   nsymb,.true.,iiden,invk1,kmden,kmats,nhere,    10d3s23
     $                  irefo,multh,itmp2,srh,sumkk,1)                  10d3s23
                     end if                                             10d3s23
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjk. 26',bc,ibc)
                     call dgemm('n','n',nrow,ncol,nmul,2d0*tff,         10d16s23
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 15')
                     if(ibl.eq.1)then                                       1d13s21
                      do iv2=0,nvirt(isbvc)-1                                1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+jsw*(nvirt(jsbv1)-1-ntop)                   1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv1=0,ntop                                        1d12s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,lj,jsb)*(iv1+nvirt(jsbv1)*( 1d14s21
     $                       ir+nrootu*iv2))                            1d14s21
                         do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do iv1=0,nvirt(isbvc)-1                                1d12s21
                       nbot=iv1+1                                           1d13s21
                       nbot=nbot+jsw*(0-nbot)                               1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv2=nbot,nvirt(jsbv2)-1                          1d13s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,lj,jsb)*(iv2+nvirt(jsbv2)*( 1d14s21
     $                       ir+nrootu*iv1))                            1d14s21
                         do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(lj.eq.1.and.jsbv12.eq.1)then                     1d14s21
                      nv=nvirt(jsbv1)*nrootu                             1d14s21
                      do iv2=0,nvirt(isbvc)-1                            1d14s21
                       do ir=0,nrootm                                    1d14s21
                        iad=joff+iv2+nvirt(jsbv1)*ir                     1d14s21
                        jprod=iprod+nfdat(2,lj,jsb)*(iv2+nvirt(jsbv1)    1d14s21
     $                      *(ir+nrootu*iv2))                           1d14s21
                        do j=0,nfdat(2,lj,jsb)-1                         1d14s21
                         gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh       1d14s21
                        end do                                           1d14s21
                       end do                                            1d14s21
                      end do                                             1d14s21
                     end if                                              1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                1d14s21
                   call ilimts(nvirt(isbr),nvirt(isbl),mynprocg,            1d13s21
     $                   mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0)then                 3d19s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdat(2,li,isb)*nfdat(2,lj,jsb)                  1d15s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    call enough('hcddjk. 27',bc,ibc)
                    fact=0d0                                                1d13s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then            1d14s21
                      i2eu=invk1(1,lsa,lsb,isbr,1)                      1d14s21
                      if(nokk(li,lj,lsa).gt.0)then                      1d14s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjk. 28',bc,ibc)
                       kint=kmats(i2eu)                                 1d15s21
                       do j=0,nokk(li,lj,lsa)-1                                 1d12s21
                        jj=ibc(ndenkf(li,lj,lsa)+j)                              1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=kint+i+nhere*jj                                 1d12s21
                         ji=itmpi+j+nokk(li,lj,lsa)*i                   1d14s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0, 1d14s21
     $                  bc(idenkf(li,lj,lsa)),mn,bc(itmpi),             1d14s21
     $                       nokk(li,lj,lsa),fact,bc(intden),mn,        1d14s21
     d' hcddjk. 16')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
c        intden
c       isbr  isbl   ibr ibl
c     1 isbv2 jsbv1   1   1
c     2 isbv1 jsbv1   2   1
c     3 isbv2 jsbv2   1   2
c     4 isbv1 jsbv2   2   2
                     nrow=nfdat(2,li,isb)*nvirt(isbr)                       1d13s21
                     nmul=nfdat(2,lj,jsb)*nhere2                            1d13s21
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjk. 29',bc,ibc)
                     do i=itmp1,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     i10=i1s                                                1d12s21
                     i1n=nvirt(isbr)                                        1d13s21
                     jntden=intden                                          1d12s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s                                            1d12s21
                      if(i2.eq.i2e)i1n=i1e                                  1d12s21
                      do i1=i10,i1n                                         1d12s21
                       i1m=i1-1                                             1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                             1d12s21
                        iad=itmp1+nfdat(2,li,isb)*(i1m+nvirt(isbr)*(i2n     1d13s21
     $                  +nhere2*j))                                     1d12s21
                        do i=0,nfdat(2,li,isb)-1                              1d12s21
                         bc(iad+i)=bc(jntden+i)                            1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdat(2,li,isb)                        1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjk. 30',bc,ibc)
                     nx=nhere2*nfdat(2,lj,jsb)*nrootu                        1d12s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibl.eq.2)then                                      1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+jsw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                         do iv1=0,ntop                                       1d12s21
                          irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+jsw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv1=i2-1                                              1d12s21
                       nbot=iv1+1                                            1d12s21
                       nbot=nbot+jsw*(0-nbot)                               1d13s21
                       do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                         do iv2=nbot,nvirt(jsbv2)-1                         1d13s21
                          irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+jsw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv2*nx)=vd(iuse+itri)                   1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(lj.eq.1.and.jsbv12.eq.1)then                     1d15s21
                      jjoff=jjoffp-nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)    1d15s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do j=0,nfdat(2,1,jsb)-1                           1d15s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=jjoff+iv1+nvirt(jsbv1)*(ir+nrootu*j)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,1,jsb)*(ir    1d15s21
     $                        +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     nwant=0                                              10d3s23
                     do lsa=1,nsymb                                       10d3s23
                      lsb=multh(lsa,ijsb)                                 10d3s23
                      if(min(nokk(li,lj,lsa),irefo(lsa),irefo(lsb)).gt.0
     $                     )nwant=nwant+1
                      nokkk(lsa)=nokk(li,lj,lsa)                              10d3s23
                      iokk(lsa)=ndenkf(li,lj,lsa)
                      iiden(lsa)=idenkf(li,lj,lsa)
                     end do                                               10d3s23
                     tfu=tf*phs                                         1d15s21
                     if(nwant.gt.0)then                                   10d3s23
                      call hcddjkd123(nrow,nmul,bc,ibc,3-ibr,nvirt,isw, 10d3s23
     $                    isbr,isbvc,isbl,nrootm,iioffp,mvv,            10d3s23
     $                    nfdat(2,li,isb),vd,mmvv,                      10d3s23
     $                    li.eq.1.and.isbv12.eq.1,nfdat(2,lj,jsb),      10d3s23
     $                  nhere2,tfu,i1s,i1e,i2s,i2e,0,nokkk,iokk,ijsb,    10d3s23
     $                   nsymb,.true.,iiden,invk1,kmden,kmats,nhere,    10d3s23
     $                  irefo,multh,itmp2,srh,sumkk,0)                  10d3s23
                     end if                                             10d3s23
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjk. 31',bc,ibc)
                     call dgemm('n','n',nrow,ncol,nmul,2d0*tfu,         10d16s23
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 17')
c     prod
c       isbr  isbvc  ibr ibl
c     1 isbv2 isbv1   1   1
c     2 isbv1 isbv2   2   1
c     3 isbv2 isbv1   1   2
c     4 isbv1 isbv2   2   2
                     if(ibr.eq.2)then                                       1d13s21
                      do iv2=0,nvirt(isbvc)-1                                1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                   1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv1=0,ntop                                        1d12s21
                         irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+isw*(irec-itri)                         1d13s21
                         iad=iioffp+itri+mvv*ir                               1d13s21
                         jprod=iprod+nfdat(2,li,isb)*(iv1+nvirt(isbv1)* 1d14s21
     $                        (ir+nrootu*iv2))                          1d14s21
                         do i=0,nfdat(2,li,isb)-1                             1d13s21
                          gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do iv1=0,nvirt(isbvc)-1                                1d12s21
                       nbot=iv1+1                                           1d13s21
                       nbot=nbot+isw*(0-nbot)                               1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv2=nbot,nvirt(isbv2)-1                          1d13s21
                         irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+isw*(irec-itri)                         1d13s21
                         iad=iioffp+itri+mvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,li,isb)*(iv2+nvirt(isbv2)* 1d14s21
     $                        (ir+nrootu*iv1))                          1d14s21
                         do i=0,nfdat(2,li,isb)-1                             1d13s21
                          gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(li.eq.1.and.isbv12.eq.1)then                    1d14s21
                      nv=nvirt(isbv1)*nrootu                            1d14s21
                      do iv2=0,nvirt(isbvc)-1                           1d14s21
                       do ir=0,nrootm                                   1d14s21
                        iad=ioff+iv2+nvirt(isbv1)*ir                    1d14s21
                        jprod=iprod+nfdat(2,li,isb)*(iv2+nvirt(isbv1)   1d14s21
     $                       *(ir+nrootu*iv2))                          1d14s21
                        do i=0,nfdat(2,li,isb)-1                        1d14s21
                         gd(iad+i*nv)=gd(iad+i*nv)+bc(jprod+i)*srh      1d18s21
                        end do                                          1d14s21
                       end do                                           1d14s21
                      end do                                            1d14s21
                     end if                                             1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                   1d13s21
                   jjoffp=jjoffp+nnvv*nfdat(2,lj,jsb)                     1d14s21
                  end if                                                  1d14s21
                  if(jtf.ne.0)tfj=-1d0                                   1d15s21
                 end do                                                  1d14s21
                end if                                                   1d14s21
                iioffp=iioffp+mmvv*nfdat(2,li,isb)                       1d14s21
                if(itf.ne.0)tfi=-1d0                                     1d14s21
               end do                                                    1d14s21
               iioffp=ioffp                                               1d12s21
               jjoffp=joffp                                              1d13s21
               tf=1d0                                                     1d13s21
               do l=1,4                                                   1d12s21
                if(min(kase,nfdat(2,l,isb),nfdat(2,l,jsb)).gt.0)then     1d14s21
                 mn=nfdat(2,l,isb)*nfdat(2,l,jsb)                         1d12s21
                 call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                 nhere=ih+1-il                                            1d12s21
                 nhere2=i2e+1-i2s                                         1d12s21
                 if(min(nhere,nvirt(isbvc)).gt.0)then                    3d19s21
                  intden=ibcoff                                            1d12s21
                  ibcoff=intden+mn*nhere                                   1d12s21
                  call enough('hcddjk. 32',bc,ibc)
                  fact=0d0                                                1d13s21
                  do lsa=1,nsymb                                              1d8s21
                   lsb=multh(lsa,ijsb)                                        1d9s21
                   if(min(irefo(lsa),irefo(lsb)).gt.0.and.lsb.ge.lsa)
     $                  then                                            1d12s21
                    i2eu=inv(1,lsa,lsb,isbl)                              1d13s21
                    icase=inv(2,lsa,lsb,isbl)                             1d13s21
                    if(nokj(l,lsa).gt.0)then                              1d12s21
                     itmpi=ibcoff                                         1d12s21
                     ibcoff=itmpi+nokj(l,lsa)*nhere                       1d12s21
                     call enough('hcddjk. 33',bc,ibc)
                     jint=jmats(i2eu)                                     1d12s21
                     if(icase.eq.1)then                                   1d12s21
                      do j=0,nokj(l,lsa)-1                                 1d12s21
                       jj=ibc(ndenjf(l,lsa)+j)                              1d12s21
                       do i=0,nhere-1                                       1d12s21
                        ij=jint+i+nhere*jj                                 1d12s21
                        ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                        bc(ji)=bc(ij)                                      1d12s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     else                                                 1d12s21
                      do j=0,nokj(l,lsa)-1                                1d12s21
                       jj=ibc(ndenjf(l,lsa)+j)                            1d12s21
                       jb=jj/irefo(lsa)                                   1d12s21
                       ja=jj-irefo(lsa)*jb                                1d12s21
                       kk=jb+irefo(lsb)*ja                                1d12s21
                       do i=0,nhere-1                                       1d12s21
                        ij=jint+i+nhere*kk                                1d12s21
                        ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                        bc(ji)=bc(ij)                                      1d12s21
                       end do                                              1d12s21
                      end do                                              1d12s21
                     end if                                               1d12s21
                     call dgemm('n','n',mn,nhere,nokj(l,lsa),1d0,         1d12s21
     $                  bc(idenjf(l,lsa)),mn,bc(itmpi),nokj(l,lsa),fact,1d12s21
     $                  bc(intden),mn,                                  1d12s21
     d' hcddjk. 18')
                     fact=1d0                                             1d12s21
                     ibcoff=itmpi                                         1d13s21
                    end if                                                1d12s21
                   end if                                                 1d13s21
                  end do                                                  1d13s21
                  if(fact.gt.0.5d0)then                                    1d12s21
                   nrow=nfdat(2,l,jsb)*nvirt(isbl)                        1d13s21
                   nmul=nfdat(2,l,isb)*nhere2
                   itmp1=ibcoff                                            1d12s21
                   ibcoff=itmp1+nrow*nmul                                 1d12s21
                   call enough('hcddjk. 34',bc,ibc)
                   do i=itmp1,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
c     intden
c       isbl  isbr  ibl ibr
c     1 jsbv1 isbv2  1   1
c     2 jsbv1 isbv1  1   2
c     3 jsbv2 isbv2  2   1
c     4 jsbv2 isbv1  2   2
                   i10=i1s                                                1d12s21
                   i1n=nvirt(isbl)                                        1d13s21
                   jntden=intden                                          1d12s21
                   do i2=i2s,i2e                                          1d12s21
                    i2n=i2-i2s                                            1d12s21
                    if(i2.eq.i2e)i1n=i1e                                  1d12s21
                    do i1=i10,i1n                                         1d12s21
                     i1m=i1-1                                             1d12s21
                     do j=0,nfdat(2,l,jsb)-1                              1d12s21
                      do i=0,nfdat(2,l,isb)-1                             1d12s21
                       iad=itmp1+j+nfdat(2,l,jsb)*(i1m+nvirt(isbl)*(i2n   1d13s21
     $                  +nhere2*i))                                     1d12s21
                       bc(iad)=bc(jntden+i)                               1d13s21
                      end do                                              1d12s21
                      jntden=jntden+nfdat(2,l,isb)                        1d12s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                    i10=1                                                 1d12s21
                   end do                                                 1d12s21
                   ncol=nvirt(isbvc)*nrootu                               1d12s21
                   itmp2=ibcoff                                           1d12s21
                   ibcoff=itmp2+nmul*ncol                                 1d12s21
                   call enough('hcddjk. 35',bc,ibc)
                   nx=nhere2*nfdat(2,l,isb)*nrootu                        1d12s21
                   do i=itmp2,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
c
c     tmp2 is vr,ki,r,vc
c
                   if(ibr.eq.1)then                                       1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv2=i2-1                                              1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                     do i=0,nfdat(2,l,isb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,l,isb)*ir)        1d13s21
                       do iv1=0,ntop                                       1d12s21
                        irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+isw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv1=i2-1                                              1d12s21
                     nbot=iv1+1                                            1d12s21
                     nbot=nbot+isw*(0-nbot)                               1d13s21
                     do i=0,nfdat(2,l,isb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,l,isb)*ir)        1d13s21
                       do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                        irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+isw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(l.eq.1.and.isbv12.eq.1)then                        1d15s21
                    iioff=iioffp-nvirt(isbv1)*nrootu*nfdat(2,1,isb)      1d15s21
                    do i2=i2s,i2e                                        1d15s21
                     iv1=i2-1                                            1d15s21
                     i2n=i2-i2s                                          1d15s21
                     do i=0,nfdat(2,1,isb)-1                              1d15s21
                      do ir=0,nrootm                                     1d15s21
                       iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)         1d15s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,1,isb)*(ir      1d15s21
     $                     +nrootu*iv1))                                1d15s21
                       bc(jtmp2)=vd(iuse)*srh                            1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                   nwant=0                                              10d3s23
                   do lsa=1,nsymb                                       10d3s23
                    lsb=multh(lsa,ijsb)                                 10d3s23
                    if(min(nokj(l,lsa),irefo(lsa),irefo(lsb)).gt.0.and.
     $                   lsb.ge.lsa)nwant=nwant+1
                    nokkk(lsa)=nokj(l,lsa)                              10d3s23
                    iokk(lsa)=ndenjf(l,lsa)
                    iiden(lsa)=idenjf(l,lsa)
                   end do                                               10d3s23
                   if(nwant.gt.0)then                                   10d3s23
                   call hcddjkd123(nrow,nmul,bc,ibc,ibl,nvirt,jsw,isbl, 10d3s23
     $                  isbvc,isbr,nrootm,jjoffp,nvv,nfdat(2,l,jsb),vd, 10d3s23
     $                  nnvv,l.eq.1.and.jsbv12.eq.1,nfdat(2,l,isb),     10d3s23
     $                  nhere2,tf,i1s,i1e,i2s,i2e,1,nokkk,iokk,ijsb,    10d3s23
     $                   nsymb,.false.,iiden,inv,jmden,jmats,nhere,     10d3s23
     $                  irefo,multh,itmp2,srh,sumjd,0)                    10d4s23
                   end if
                   iprod=ibcoff                                           1d12s21
                   ibcoff=iprod+nrow*ncol                                 1d12s21
                   call enough('hcddjk. 36',bc,ibc)
                   call dgemm('n','n',nrow,ncol,nmul,2d0*tf,            10d16s23
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 19')
c        prod
c       isbl  isbvc ibl ibr
c     1 jsbv1 jsbv2  1   1
c     2 jsbv1 jsbv2  1   2
c     3 jsbv2 jsbv1  2   1
c     4 jsbv2 jsbv1  2   2
                   if(ibl.eq.1)then                                       1d13s21
                    do iv2=0,nvirt(isbvc)-1                                1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+jsw*(nvirt(jsbv1)-1-ntop)                   1d13s21
                     do ir=0,nrootm                                        1d12s21
                      do iv1=0,ntop                                        1d12s21
                       irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+jsw*(irec-itri)                          1d13s21
                       iad=jjoffp+itri+nvv*ir                            1d13s21
                       jprod=iprod+nfdat(2,l,jsb)*(iv1+nvirt(jsbv1)*(ir    1d13s21
     $                     +nrootu*iv2))                                  1d12s21
                       do j=0,nfdat(2,l,jsb)-1                             1d13s21
                        gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do iv1=0,nvirt(isbvc)-1                                1d12s21
                     nbot=iv1+1                                           1d13s21
                     nbot=nbot+jsw*(0-nbot)                               1d13s21
                     do ir=0,nrootm                                        1d12s21
                      do iv2=nbot,nvirt(jsbv2)-1                          1d13s21
                       irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+jsw*(irec-itri)                          1d13s21
                       iad=jjoffp+itri+nvv*ir                            1d13s21
                       jprod=iprod+nfdat(2,l,jsb)*(iv2+nvirt(jsbv2)*(ir   1d13s21
     $                     +nrootu*iv1))                                1d13s21
                       do j=0,nfdat(2,l,jsb)-1                             1d13s21
                        gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(l.eq.1.and.jsbv12.eq.1)then                        1d15s21
                    jjoff=jjoffp-nvirt(isbvc)*nrootu*nfdat(2,1,jsb)      1d15s21
                    nv=nrootu*nvirt(isbvc)                               1d15s21
                    do iv1=0,nvirt(isbvc)-1                              1d15s21
                     do ir=0,nrootm                                      1d15s21
                      iad=jjoff+iv1+nvirt(isbvc)*ir                      1d15s21
                      jprod=iprod+nfdat(2,l,jsb)*(iv1+nvirt(jsbv2)*(ir   1d15s21
     $                    +nrootu*iv1))                                 1d15s21
                      do j=0,nfdat(2,l,jsb)-1                            1d15s21
                       gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh         1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                  end if                                                   1d12s21
                  ibcoff=intden                                            1d12s21
                 end if                                                   1d13s21
                 if(isb.ne.jsb)then
                  call ilimts(nvirt(isbr),nvirt(isbl),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                  nhere=ih+1-il                                            1d12s21
                  nhere2=i2e+1-i2s                                         1d12s21
                  if(min(nhere,nvirt(isbvc)).gt.0)then                   3d19s21
                   intden=ibcoff                                            1d12s21
                   ibcoff=intden+mn*nhere                                   1d12s21
                   call enough('hcddjk. 37',bc,ibc)
                   fact=0d0                                                1d13s21
                   do lsa=1,nsymb                                              1d8s21
                    lsb=multh(lsa,ijsb)                                        1d9s21
                    if(min(irefo(lsa),irefo(lsb)).gt.0.and.              1d13s21
     $                  lsb.ge.lsa)then                                 1d13s21
                     i2eu=inv(1,lsa,lsb,isbr)                             1d13s21
                     icase=inv(2,lsa,lsb,isbr)                            1d13s21
                     if(nokj(l,lsa).gt.0)then                              1d12s21
                      itmpi=ibcoff                                         1d12s21
                      ibcoff=itmpi+nokj(l,lsa)*nhere                       1d12s21
                      call enough('hcddjk. 38',bc,ibc)
                      jint=jmats(i2eu)                                     1d12s21
                      if(icase.eq.1)then                                   1d12s21
                       do j=0,nokj(l,lsa)-1                                 1d12s21
                        jj=ibc(ndenjf(l,lsa)+j)                              1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=jint+i+nhere*jj                                 1d12s21
                         ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      else                                                 1d12s21
                       do j=0,nokj(l,lsa)-1                                1d12s21
                        jj=ibc(ndenjf(l,lsa)+j)                            1d12s21
                        jb=jj/irefo(lsa)                                   1d12s21
                        ja=jj-irefo(lsa)*jb                                1d12s21
                        kk=jb+irefo(lsb)*ja                                1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=jint+i+nhere*kk                                1d12s21
                         ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                              1d12s21
                      end if                                               1d12s21
                      call dgemm('n','n',mn,nhere,nokj(l,lsa),1d0,         1d12s21
     $                  bc(idenjf(l,lsa)),mn,bc(itmpi),nokj(l,lsa),fact,1d12s21
     $                  bc(intden),mn,                                  1d12s21
     d' hcddjk. 20')
                      fact=1d0                                             1d12s21
                      ibcoff=itmpi                                         1d13s21
                     end if                                                1d12s21
                    end if                                                 1d13s21
                   end do                                                  1d13s21
                   if(fact.gt.0.5d0)then                                    1d12s21
c     intden
c       isbr  isbl   ibr ibl
c     1 isbv2 jsbv1   1   1
c     2 isbv1 jsbv1   2   1
c     3 isbv2 jsbv2   1   2
c     4 isbv1 jsbv2   2   2
                    nrow=nfdat(2,l,isb)*nvirt(isbr)                       1d13s21
                    nmul=nfdat(2,l,jsb)*nhere2                            1d13s21
                    itmp1=ibcoff                                            1d12s21
                    ibcoff=itmp1+nrow*nmul                                 1d12s21
                    call enough('hcddjk. 39',bc,ibc)
                    do i=itmp1,ibcoff-1                                    1d12s21
                     bc(i)=0d0                                             1d12s21
                    end do                                                 1d12s21
                    i10=i1s                                                1d12s21
                    i1n=nvirt(isbr)                                        1d13s21
                    jntden=intden                                          1d12s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s                                            1d12s21
                     if(i2.eq.i2e)i1n=i1e                                  1d12s21
                     do i1=i10,i1n                                         1d12s21
                      i1m=i1-1                                             1d12s21
                      do j=0,nfdat(2,l,jsb)-1                             1d12s21
                       iad=itmp1+nfdat(2,l,isb)*(i1m+nvirt(isbr)*(i2n     1d13s21
     $                  +nhere2*j))                                     1d12s21
                       do i=0,nfdat(2,l,isb)-1                              1d12s21
                        bc(iad+i)=bc(jntden+i)                            1d13s21
                       end do                                              1d12s21
                       jntden=jntden+nfdat(2,l,isb)                        1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                     i10=1                                                 1d12s21
                    end do                                                 1d12s21
                    ncol=nvirt(isbvc)*nrootu                               1d12s21
                    itmp2=ibcoff                                           1d12s21
                    ibcoff=itmp2+nmul*ncol                                 1d12s21
                    call enough('hcddjk. 40',bc,ibc)
                    nx=nhere2*nfdat(2,l,jsb)*nrootu                        1d12s21
                    do i=itmp2,ibcoff-1                                    1d12s21
                     bc(i)=0d0                                             1d12s21
                    end do                                                 1d12s21
c
c     tmp2 is vl,kj,r,vc
c
                    if(ibl.eq.2)then                                      1d13s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s
                      iv2=i2-1                                              1d12s21
                      ntop=iv2-1                                            1d12s21
                      ntop=ntop+jsw*(nvirt(isbvc)-1-ntop)                   1d12s21
                      do j=0,nfdat(2,l,jsb)-1                               1d12s21
                       do ir=0,nrootm                                       1d12s21
                        iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,l,jsb)*ir)        1d13s21
                        do iv1=0,ntop                                       1d12s21
                         irec=iv1+nvirt(isbvc)*iv2                      10d2s23
                         itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    else                                                   1d13s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s
                      iv1=i2-1                                              1d12s21
                      nbot=iv1+1                                            1d12s21
                      nbot=nbot+jsw*(0-nbot)                               1d13s21
                      do j=0,nfdat(2,l,jsb)-1                               1d12s21
                       do ir=0,nrootm                                       1d12s21
                        iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,l,jsb)*ir)        1d13s21
                        do iv2=nbot,nvirt(jsbv2)-1                         1d13s21
                         irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         bc(jtmp2+iv2*nx)=vd(iuse+itri)                   1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    end if                                                 1d13s21
                    if(l.eq.1.and.jsbv12.eq.1)then                       1d15s21
                     jjoff=jjoffp-nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)     1d15s21
                     do i2=i2s,i2e                                       1d15s21
                      iv1=i2-1                                           1d15s21
                      i2n=i2-i2s                                         1d15s21
                      do j=0,nfdat(2,1,jsb)-1                             1d15s21
                       do ir=0,nrootm                                    1d15s21
                        iuse=jjoff+iv1+nvirt(jsbv1)*(ir+nrootu*j)        1d15s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,1,jsb)*(ir     1d15s21
     $                      +nrootu*iv1))                               1d15s21
                        bc(jtmp2)=vd(iuse)*srh                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end if                                               1d15s21
                    nwant=0                                              10d3s23
                    do lsa=1,nsymb                                       10d3s23
                     lsb=multh(lsa,ijsb)                                 10d3s23
                     if(min(nokj(l,lsa),irefo(lsa),irefo(lsb)).gt.0.and.
     $                   lsb.ge.lsa)nwant=nwant+1
                     nokkk(lsa)=nokj(l,lsa)                              10d3s23
                     iokk(lsa)=ndenjf(l,lsa)
                     iiden(lsa)=idenjf(l,lsa)
                    end do                                               10d3s23
                    if(nwant.gt.0)then                                   10d3s23
                     call hcddjkd123(nrow,nmul,bc,ibc,3-ibr,nvirt,isw,
     $                   isbr,isbvc,isbl,nrootm,iioffp,mvv,
     $                   nfdat(2,l,isb),vd,mmvv,l.eq.1.and.isbv12.eq.1,
     $                   nfdat(2,l,jsb),nhere2,tf,i1s,i1e,i2s,i2e,0,
     $                   nokkk,iokk,ijsb,nsymb,.false.,iiden,inv,jmden,
     $                   jmats,nhere,irefo,multh,itmp2,srh,sumjd,0)     10d4s23
                    end if
                    iprod=ibcoff                                           1d12s21
                    ibcoff=iprod+nrow*ncol                                 1d12s21
                    call enough('hcddjk. 41',bc,ibc)
                    call dgemm('n','n',nrow,ncol,nmul,2d0*tf,           10d16s23
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 21')
c     prod
c       isbr  isbvc  ibr ibl
c     1 isbv2 isbv1   1   1
c     2 isbv1 isbv2   2   1
c     3 isbv2 isbv1   1   2
c     4 isbv1 isbv2   2   2
                    if(ibr.eq.2)then                                       1d13s21
                     do iv2=0,nvirt(isbvc)-1                                1d12s21
                      ntop=iv2-1                                            1d12s21
                      ntop=ntop+isw*(nvirt(isbr)-1-ntop)                10d2s23
                      do ir=0,nrootm                                        1d12s21
                       do iv1=0,ntop                                        1d12s21
                        irec=iv1+nvirt(isbr)*iv2                           1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                        itri=itri+isw*(irec-itri)                         1d13s21
                        iad=iioffp+itri+mvv*ir                               1d13s21
                        jprod=iprod+nfdat(2,l,isb)*(iv1+nvirt(isbr)*(ir 10d2s23
     $                     +nrootu*iv2))                                  1d12s21
                        do i=0,nfdat(2,l,isb)-1                             1d13s21
                         gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    else                                                   1d13s21
                     do iv1=0,nvirt(isbvc)-1                                1d12s21
                      nbot=iv1+1                                           1d13s21
                      nbot=nbot+isw*(0-nbot)                               1d13s21
                      do ir=0,nrootm                                        1d12s21
                       do iv2=nbot,nvirt(isbr)-1                        10d2s23
                        irec=iv1+nvirt(isbvc)*iv2                       10d2s23
                        itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                        itri=itri+isw*(irec-itri)                         1d13s21
                        iad=iioffp+itri+mvv*ir                            1d13s21
                        jprod=iprod+nfdat(2,l,isb)*(iv2+nvirt(isbr)*(ir 10d2s23
     $                      +nrootu*iv1))                                  1d12s21
                        do i=0,nfdat(2,l,isb)-1                             1d13s21
                         gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    end if                                                 1d13s21
                    if(l.eq.1.and.isbv12.eq.1)then                        1d15s21
                     iioff=iioffp-nvirt(isbvc)*nrootu*nfdat(2,1,isb)      1d15s21
                     nv=nrootu*nvirt(isbvc)                               1d15s21
                     do iv1=0,nvirt(isbvc)-1                              1d15s21
                      do ir=0,nrootm                                      1d15s21
                       iad=iioff+iv1+nvirt(isbvc)*ir                      1d15s21
                       jprod=iprod+nfdat(2,l,isb)*(iv1+nvirt(isbv2)*(ir   1d15s21
     $                    +nrootu*iv1))                                 1d15s21
                       do i=0,nfdat(2,l,isb)-1                            1d15s21
                        gd(iad+i*nv)=gd(iad+i*nv)+bc(jprod+i)*srh       10d16s23
                       end do                                             1d15s21
                      end do                                              1d15s21
                     end do                                               1d15s21
                    end if                                                1d15s21
                   end if                                                   1d12s21
                   ibcoff=intden                                            1d12s21
                  end if                                                   1d13s21
                 end if                                                    1d13s21
                end if                                                      1d13s21
                if(mod(jtf+itf,2).ne.0)then                              1d15s21
                 tf=-1d0                                                 1d14s21
                end if                                                   1d14s21
                iioffp=iioffp+mmvv*nfdat(2,l,isb)                          1d13s21
                jjoffp=jjoffp+nnvv*nfdat(2,l,jsb)                            1d13s21
               end do                                                      1d13s21
              end do                                                     1d13s21
             end if                                                      1d13s21
             if(jsbv12.eq.1)joff=joff                                    1d12s21
     $           +nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)                    1d12s21
             do l=1,4                                                    1d12s21
              joff=joff+nnvv*nfdat(2,l,jsb)                              1d12s21
             end do                                                      1d12s21
            end if                                                       1d12s21
           end do                                                        1d12s21
           if(isbv12.eq.1)ioff=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)   1d12s21
           do l=1,4                                                      1d12s21
            ioff=ioff+mmvv*nfdat(2,l,isb)                                1d12s21
           end do                                                        1d12s21
          end if                                                         1d12s21
         end do                                                          1d12s21
         if(ijsb.eq.1)then                                               1d11s21
          ioff=ioffdnon                                                  1d11s21
          do isbv1=1,nsymb                                               1d11s21
           isbv2=multh(isbv1,isbv12)                                     1d11s21
           if(isbv2.ge.isbv1)then                                        1d11s21
            nvv=0                                                        4d28s21
            if(isbv1.eq.isbv2)then                                       1d11s21
             ioffs=ioff                                                  1d12s21
             nnn=nvirt(isbv1)*nrootu                                     1d11s21
             call ilimts(nvirt(isbv1),nrootu,mynprocg,mynowprog,il,ih,   1d11s21
     $           i1s,i1e,i2s,i2e)                                       1d11s21
             nhere=ih+1-il                                               1d11s21
             if(nhere.gt.0)then                                          1d11s21
              itmp=ibcoff                                                1d11s21
              ibcoff=itmp+nhere*nfdat(2,1,isb)                           1d11s21
              call enough('hcddjk. 42',bc,ibc)
              call dgemm('n','n',nhere,nfdat(2,1,isb),nfdat(2,1,isb),    1d11s21
     $            2d0,vd(ioff+il-1),nnn,bc(iden1ef(1)),nfdat(2,1,isb),  10d16s23
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjkd12. 22')
              do i=0,nfdat(2,1,isb)-1                                    1d11s21
               i10=i1s                                                   1d11s21
               i1n=nvirt(isbv1)                                          1d11s21
               jtmp=itmp+nhere*i                                         1d11s21
               do i2=i2s,i2e                                             1d11s21
                if(i2.eq.i2e)i1n=i1e                                     1d11s21
                iad=ioff-1+nvirt(isbv1)*(i2-1+nrootu*i)                  1d11s21
                do i1=i10,i1n                                            1d11s21
                 gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                 jtmp=jtmp+1                                             1d11s21
                end do                                                   1d11s21
                i10=1                                                    1d11s21
               end do                                                    1d11s21
              end do                                                     1d11s21
              ibcoff=itmp                                                1d11s21
              do is=1,nsdlk                                             9d29s23
               if(isblk(1,is).eq.isblk(2,is))then                       9d29s23
                nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2      9d29s23
                ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2      9d29s23
               else                                                     9d29s23
                nrow=irefo(isblk(1,is))*irefo(isblk(2,is))              9d29s23
                ncol=irefo(isblk(3,is))*irefo(isblk(4,is))              9d29s23
               end if                                                   9d29s23
               if(min(nrow,ncol).gt.0)then                              9d29s23
                nadd=nrow*ncol                                          9d29s23
                nadd2=nfdat(2,1,isb)*nfdat(2,1,isb)                     9d29s23
                nncol=nfdat(2,1,isb)*nadd                               9d29s23
                itmp4o=ibcoff                                           9d29s23
                ibcoff=itmp4o+nhere*nncol                               9d29s23
                call enough('hcddjkd12.tmp4o',bc,ibc)                   9d29s23
                call dgemm('n','n',nhere,nncol,nfdat(2,1,isb),          9d29s23
     $               1d0,vd(ioff+il-1),nnn,bc(iden14o(1,is)),           9d29s23
     $               nfdat(2,1,isb),0d0,bc(itmp4o),nhere,               9d29s23
     $               'hcddjkd12.tmp4ovisv')                             9d29s23
                i10=i1s                                                 9d28s23
                i1n=nvirt(isbv1)                                        9d29s23
                jtmp4o=itmp4o                                           9d29s23
                ivd=ioff+il-1                                           9d28s23
                n4o=nhere*nfdat(2,1,isb)                                9d29s23
                do i2=i2s,i2e                                           9d28s23
                 ir=i2-1                                                9d28s23
                 if(i2.eq.i2e)i1n=i1e                                   9d28s23
                 iad2=id4o(is)+nadd*ir                                  9d29s23
                 do i1=i10,i1n                                          9d28s23
                  do k=0,nfdat(2,1,isb)-1                               9d28s23
                   iad1=jtmp4o+nhere*k
                   jvd=ivd+nnn*k                                        9d28s23
                   do i34=0,nadd-1                                       9d29s23
                    bc(iad2+i34)=bc(iad2+i34)+vd(jvd)*bc(iad1)           9d28s23
                    iad1=iad1+n4o                                       9d29s23
                   end do                                               9d28s23
                  end do                                                9d28s23
                  jtmp4o=jtmp4o+1                                       9d29s23
                  ivd=ivd+1                                             9d28s23
                 end do                                                 9d28s23
                 i10=1                                                  9d28s23
                end do                                                  9d28s23
                ibcoff=itmp4o                                           9d29s23
               end if                                                   9d29s23
              end do                                                    9d29s23
              do ksa=1,nsymb                                            9d28s23
               if(irefo(ksa).gt.0)then                                  9d28s23
                nh0dcol=(irefo(ksa)*(irefo(ksa)+1))/2                   9d28s23
                itmph0d=ibcoff                                            9d28s23
                ibcoff=itmph0d+nhere*nfdat(2,1,isb)*nh0dcol               9d28s23
                call enough('hcddjkd12.tmph0d',bc,ibc)                  9d28s23
                ncol=nfdat(2,1,isb)*nh0dcol                             9d28s23
                call dgemm('n','n',nhere,ncol,nfdat(2,1,isb),           9d28s23
     $            1d0,vd(ioff+il-1),nnn,bc(iden1enh(1,ksa)),
     $               nfdat(2,1,isb),0d0,bc(itmph0d),nhere,              9d28s23
     d' hcddjkd12. 22')
                i10=i1s                                                 9d28s23
                i1n=nvirt(isbv1)                                        9d29s23
                jtmph0d=itmph0d                                         9d28s23
                ivd=ioff+il-1                                           9d28s23
                do i2=i2s,i2e                                           9d28s23
                 ir=i2-1                                                9d28s23
                 if(i2.eq.i2e)i1n=i1e                                   9d28s23
                 iad2=iden(ksa)+nh0av(ksa)*nh0av(ksa)*ir                9d28s23
                 do i1=i10,i1n                                          9d28s23
                  do k=0,nfdat(2,1,isb)-1                               9d28s23
                   iad1=jtmph0d+nhere*nh0dcol*k                         9d28s23
                   jvd=ivd+nnn*k                                        9d28s23
                   do i4=0,irefo(ksa)-1                                  9d28s23
                    iad3=iad2+nh0av(ksa)*i4                              9d28s23
                    do i3=0,i4                                           9d28s23
                     bc(iad3+i3)=bc(iad3+i3)+vd(jvd)*bc(iad1)           9d28s23
                     iad1=iad1+nhere                                    9d28s23
                    end do                                               9d28s23
                   end do                                                9d28s23
                  end do                                                9d28s23
                  jtmph0d=jtmph0d+1                                     9d28s23
                  ivd=ivd+1                                             9d28s23
                 end do                                                 9d28s23
                 i10=1                                                  9d28s23
                end do                                                  9d28s23
                ibcoff=itmph0d                                          9d28s23
               end if                                                   9d28s23
              end do                                                    9d28s23
             end if                                                      1d11s21
             if(nnn.gt.0)then                                            1d25s21
              ioff=ioff+nnn*nfdat(2,1,isb)                                1d11s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       1d11s21
              isw=0                                                       1d11s21
              itmp=ibcoff                                                 1d18s21
              ibcoff=itmp+nnn*nfdat(2,1,isb)                              1d18s21
              call enough('hcddjk. 43',bc,ibc)
              call dgemm('n','n',nnn,nfdat(2,1,isb),nfdat(2,1,isb),       1d18s21
     $              2d0,vd(ioffs),nnn,bc(idenhvvf(1)),nfdat(2,1,isb),   1d18s21
     $              0d0,bc(itmp),nnn,                                   1d18s21
     d' hcddjk. 23')
              do i=0,nfdat(2,1,isb)-1                                     1d18s21
               do ir=0,nrootm                                             1d18s21
                do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   1d18s21
                 iad=ioffs+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                 ih0=ih0av(isbv1)+(iv+irefo(isbv1))*(nh0av(isbv1)+1)      1d18s21
                 jtmp=itmp+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                 term=bc(ih0)*bc(jtmp)*2d0                              10d16s23
                 gd(iad)=gd(iad)+term                                   9d28s23
                 ih0d=iden(isbv1)+(iv+irefo(isbv1))*(nh0av(isbv1)+1)    9d28s23
                 bc(ih0d)=bc(ih0d)+vd(iad)*bc(jtmp)                     9d28s23
                end do                                                    1d18s21
               end do                                                     1d18s21
              end do                                                      1d18s21
              ibcoff=itmp                                                 1d18s21
             end if                                                      1d25s21
            else                                                         1d11s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               1d11s21
             isw=1                                                       1d11s21
            end if                                                       1d11s21
            nnn=nvv*nrootu                                               1d11s21
            tf=1d0                                                       1d11s21
            do l=1,4                                                     1d11s21
             if(nfdat(2,l,isb).gt.0)then                                 1d11s21
              call ilimts(nvv,nrootu,mynprocg,mynowprog,il,ih,            1d11s21
     $            i1s,i1e,i2s,i2e)                                       1d11s21
              nhere=ih+1-il                                               1d11s21
              if(nhere.gt.0)then                                          1d11s21
               itmp=ibcoff                                                1d11s21
               ibcoff=itmp+nhere*nfdat(2,l,isb)                           1d11s21
               if(isbv12.eq.1.and.l.eq.1)then                            1d12s21
                nrow=nvirt(isbv1)*nrootu                                 1d12s21
                ibcoff=max(ibcoff,itmp+nrow*nfdat(2,l,isb))              1d12s21
               end if                                                    1d12s21
               call enough('hcddjk. 44',bc,ibc)
               call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(2,l,isb),    1d11s21
     $            2d0,vd(ioff+il-1),nnn,bc(iden1ef(l)),nfdat(2,l,isb),  10d16s23
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjk. 24')
               do i=0,nfdat(2,l,isb)-1                                    1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                  1d11s21
                jtmp=itmp+nhere*i                                         1d11s21
                do i2=i2s,i2e                                             1d11s21
                 if(i2.eq.i2e)i1n=i1e                                     1d11s21
                 iad=ioff-1+nvv*(i2-1+nrootu*i)                          1d11s21
                 do i1=i10,i1n                                            1d11s21
                  gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                  jtmp=jtmp+1                                             1d11s21
                 end do                                                   1d11s21
                 i10=1                                                    1d11s21
                end do                                                    1d11s21
               end do                                                     1d11s21
               do is=1,nsdlk                                             9d29s23
                if(isblk(1,is).eq.isblk(2,is))then                       9d29s23
                 nrw=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2      9d29s23
                 ncl=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2      9d29s23
                else                                                     9d29s23
                 nrw=irefo(isblk(1,is))*irefo(isblk(2,is))              9d29s23
                 ncl=irefo(isblk(3,is))*irefo(isblk(4,is))              9d29s23
                end if                                                   9d29s23
                if(min(nrw,ncl).gt.0)then                               9d29s23
                 nadd=nrw*ncl                                           10d13s23
                 nadd2=nfdat(2,l,isb)*nfdat(2,l,isb)                     9d29s23
                 nncol=nfdat(2,l,isb)*nadd                               9d29s23
                 itmp4o=ibcoff                                           9d29s23
                 ibcoff=itmp4o+nhere*nncol                               9d29s23
                 call enough('hcddjkd12.tmp4o',bc,ibc)                   9d29s23
                 call dgemm('n','n',nhere,nncol,nfdat(2,l,isb),          9d29s23
     $               1d0,vd(ioff+il-1),nnn,bc(iden14o(l,is)),           9d29s23
     $               nfdat(2,l,isb),0d0,bc(itmp4o),nhere,               9d29s23
     $               'hcddjkd12.tmp4ovisv')                             9d29s23
                 i10=i1s                                                 9d28s23
                 i1n=nvv                                                 9d29s23
                 jtmp4o=itmp4o                                           9d29s23
                 ivd=ioff+il-1                                           9d28s23
                 n4o=nhere*nfdat(2,l,isb)                                9d29s23
                 do i2=i2s,i2e                                           9d28s23
                  ir=i2-1                                                9d28s23
                  if(i2.eq.i2e)i1n=i1e                                   9d28s23
                  iad2=id4o(is)+nadd*ir                                  9d29s23
                  do i1=i10,i1n                                          9d28s23
                   do k=0,nfdat(2,l,isb)-1                               9d28s23
                    iad1=jtmp4o+nhere*k
                    jvd=ivd+nnn*k                                        9d28s23
                    do i34=0,nadd-1                                       9d29s23
                     bc(iad2+i34)=bc(iad2+i34)+vd(jvd)*bc(iad1)           9d28s23
                     iad1=iad1+n4o                                       9d29s23
                    end do                                               9d28s23
                   end do                                                9d28s23
                   jtmp4o=jtmp4o+1                                       9d29s23
                   ivd=ivd+1                                             9d28s23
                  end do                                                 9d28s23
                  i10=1                                                  9d28s23
                 end do                                                  9d28s23
                 ibcoff=itmp4o                                           9d29s23
                end if                                                   9d29s23
               end do                                                    9d29s23
               do ksa=1,nsymb                                            9d28s23
                if(irefo(ksa).gt.0)then                                  9d28s23
                 nh0dcol=(irefo(ksa)*(irefo(ksa)+1))/2                   9d28s23
                 itmph0d=ibcoff                                            9d28s23
                 ibcoff=itmph0d+nhere*nfdat(2,l,isb)*nh0dcol            9d29s23
                 call enough('hcddjkd12.tmph0d',bc,ibc)                  9d28s23
                 ncol=nfdat(2,l,isb)*nh0dcol                             9d28s23
                 call dgemm('n','n',nhere,ncol,nfdat(2,l,isb),           9d28s23
     $            1d0,vd(ioff+il-1),nnn,bc(iden1enh(l,ksa)),            9d29s23
     $               nfdat(2,l,isb),0d0,bc(itmph0d),nhere,              9d29s23
     d' hcddjkd12. 22')
                 i10=i1s                                                 9d28s23
                 i1n=nvv                                                9d29s23
                 jtmph0d=itmph0d                                         9d28s23
                 ivd=ioff+il-1                                           9d28s23
                 do i2=i2s,i2e                                           9d28s23
                  ir=i2-1                                                9d28s23
                  if(i2.eq.i2e)i1n=i1e                                   9d28s23
                  iad2=iden(ksa)+nh0av(ksa)*nh0av(ksa)*ir                9d28s23
                  do i1=i10,i1n                                          9d28s23
                   do k=0,nfdat(2,l,isb)-1                               9d28s23
                    iad1=jtmph0d+nhere*nh0dcol*k                         9d28s23
                    jvd=ivd+nnn*k                                        9d28s23
                    do i4=0,irefo(ksa)-1                                  9d28s23
                     iad3=iad2+nh0av(ksa)*i4                              9d28s23
                     do i3=0,i4                                           9d28s23
                      bc(iad3+i3)=bc(iad3+i3)+vd(jvd)*bc(iad1)           9d28s23
                      iad1=iad1+nhere                                    9d28s23
                     end do                                               9d28s23
                    end do                                                9d28s23
                   end do                                                9d28s23
                   jtmph0d=jtmph0d+1                                     9d28s23
                   ivd=ivd+1                                             9d28s23
                  end do                                                 9d28s23
                  i10=1                                                  9d28s23
                 end do                                                  9d28s23
                 ibcoff=itmph0d                                          9d28s23
                end if                                                   9d28s23
               end do                                                    9d28s23
c
c     Gvv'ri=dij Hvv" Vv'v"rj ... if sym v = sym v' = sym v" we get this
c     if v"=v' > Gvv'ri=dij Hvv' Vv'v'rj.
c     if v=v' > Gvvri dij Hvv" Vvv"rj
               if(isbv12.eq.1.and.l.eq.1.and.nrow.gt.0)then              4d28s21
                call dgemm('n','n',nrow,nfdat(2,1,isb),nfdat(2,1,isb),
     $              sr2,vd(ioffs),nrow,bc(idenhvvf(1)),nfdat(2,1,isb),  1d12s21
     $              0d0,bc(itmp),nrow,                                  1d12s21
     d' hcddjk. 25')
                do i=0,nfdat(2,l,isb)-1                                   1d11s21
                 i10=i1s                                                   1d11s21
                 i1n=nvv                                                   1d11s21
                 do i2=i2s,i2e                                             1d11s21
                  ir=i2-1                                                 1d11s21
                  if(i2.eq.i2e)i1n=i1e                                    1d11s21
                  do i1=i10,i1n                                           1d11s21
                   i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                   iv2=i1m/nvirt(isbv1)                                   1d11s21
                   iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                   try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                   iv2t=sqrt(try)+0.5d0                                   1d11s21
                   iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                   iv1=iv1t+isw*(iv1-iv1t)
                   iv2=iv2t+isw*(iv2-iv2t)
c     Gvv'ri=dij Hvv' Vv'v'rj
                   ih0=ih0av(isbv1)+iv1+irefo(isbv1)+nh0av(isbv1)*(      1d12s21
     $                irefo(isbv1)+iv2)                                 1d12s21
                   iad=ioff+i1m+nvv*(ir+nrootu*i)                        1d12s21
                   jtmp=itmp+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                   ktmp=itmp+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                   term=(bc(ktmp)+bc(jtmp))*bc(ih0)*2d0                 10d16s23
                   gd(iad)=gd(iad)+term                                 9d28s23
                   ih0d=iden(isbv1)+iv1+irefo(isbv1)+nh0av(isbv1)*(     9d28s23
     $                irefo(isbv1)+iv2)                                 1d12s21
                   bc(ih0d)=bc(ih0d)+vd(iad)*(bc(ktmp)+bc(jtmp))        9d28s23
                  end do                                                 1d12s21
                  i10=1                                                  4d21s21
                 end do                                                  1d12s21
                end do                                                   1d12s21
               end if                                                    1d12s21
c     Gvv'ri=dij Hvv" Vv"v'rj, Gvv'ri=dij Hv'v" Vvv"rj
c
               call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(2,l,isb),   1d11s21
     $             tf,vd(ioff+il-1),nnn,bc(idenhvvf(l)),nfdat(2,l,isb), 1d11s21
     $             0d0,bc(itmp),nhere,                                  1d11s21
     d              ' hcddjk. 26')
               if(l.eq.1.and.isbv12.eq.1)then                            1d18s21
                do i=0,nfdat(2,l,isb)-1                                   1d11s21
                 jtmp=itmp+nhere*i                                        1d11s21
                 i10=i1s                                                   1d11s21
                 i1n=nvv                                                   1d11s21
                 do i2=i2s,i2e                                             1d11s21
                  ir=i2-1                                                 1d11s21
                  if(i2.eq.i2e)i1n=i1e                                    1d11s21
                  ltmp=jtmp                                              1d18s21
                  do i1=i10,i1n                                           1d11s21
                   i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                   iv2=i1m/nvirt(isbv1)                                   1d11s21
                   iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                   try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                   iv2t=sqrt(try)+0.5d0                                   1d11s21
                   iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                   iv1=iv1t+isw*(iv1-iv1t)
                   iv2=iv2t+isw*(iv2-iv2t)
c     if v=v' > Gvvri dij Hvv" Vvv"rj
                   ih0=ih0av(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(      1d18s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                   iad=ioffs+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                   term=bc(ih0)*bc(jtmp)*sr2*2d0                        10d16s23
                   gd(iad)=gd(iad)+term                                 9d28s23
                   ih0d=iden(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(     9d28s23
     $                irefo(isbv1)+iv2)                                 1d18s21
                   bc(ih0d)=bc(ih0d)+vd(iad)*bc(jtmp)*sr2               9d28s23
                   jtmp=jtmp+1                                           1d18s21
                  end do                                                 1d18s21
                  if(isbv1.eq.isbv2)then                                 1d18s21
                   jtmp=ltmp                                             1d18s21
                   do i1=i10,i1n                                           1d11s21
                    i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                    iv2=i1m/nvirt(isbv1)                                   1d11s21
                    iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                    try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                    iv2t=sqrt(try)+0.5d0                                   1d11s21
                    iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                    iv1=iv1t+isw*(iv1-iv1t)
                    iv2=iv2t+isw*(iv2-iv2t)
c     if v=v' > Gvvri dij Hvv" Vv"vrj
                    ih0=ih0av(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(      1d18s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                    iad=ioffs+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                    term=bc(ih0)*bc(jtmp)*sr2*2d0                       10d16s23
                    gd(iad)=gd(iad)+term                                9d28s23
                    ih0d=iden(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(    9d28s23
     $                irefo(isbv1)+iv2)                                 9d28s23
                    bc(ih0d)=bc(ih0d)+vd(iad)*bc(jtmp)*sr2              9d28s23
                    jtmp=jtmp+1                                           1d18s21
                   end do                                                 1d18s21
                  end if                                                 1d18s21
                  i10=1                                                  4d21s21
                 end do                                                  1d18s21
                end do                                                   1d18s21
               end if                                                    1d18s21
               do i=0,nfdat(2,l,isb)-1                                   1d11s21
                jtmp=itmp+nhere*i                                        1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                   1d11s21
                do i2=i2s,i2e                                             1d11s21
                 ir=i2-1                                                 1d11s21
                 if(i2.eq.i2e)i1n=i1e                                    1d11s21
                 do i1=i10,i1n                                           1d11s21
                  i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                  iv2=i1m/nvirt(isbv1)                                   1d11s21
                  iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                  try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                  iv2t=sqrt(try)+0.5d0                                   1d11s21
                  iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                  iv1=iv1t+isw*(iv1-iv1t)
                  iv2=iv2t+isw*(iv2-iv2t)
c     Gvv'ri=dij Hvv" Vv"v'rj
                  ntop=iv2-1                                             1d11s21
                  ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                    1d11s21
                  ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*(           1d12s21
     $                irefo(isbv1)+iv1)                                 1d12s21
                  ih0d=iden(isbv1)+irefo(isbv1)+nh0av(isbv1)*(          9d28s23
     $                irefo(isbv1)+iv1)                                 1d12s21
                  iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                  do iv=0,ntop                                           1d11s21
                   irec=iv+nvirt(isbv1)*iv2                              1d11s21
                   itri=((iv2*(iv2-1))/2)+iv                             1d18s21
                   irow=itri+isw*(irec-itri)                             1d11s21
                   term=bc(ih0+iv)*bc(jtmp)*2d0                         10d16s23
                   gd(iad+irow)=gd(iad+irow)+term                       9d28s23
                   bc(ih0d+iv)=bc(ih0d+iv)+vd(iad+irow)*bc(jtmp)        9d28s23
                  end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vvv"rj
                  nbot=iv1+1                                             1d11s21
                  nbot=nbot+isw*(0-nbot)                                 1d11s21
                  ih0=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*            1d12s21
     $                (irefo(isbv2)+iv2)                                1d12s21
                  ih0d=iden(isbv2)+irefo(isbv2)+nh0av(isbv2)*           9d28s23
     $                (irefo(isbv2)+iv2)                                1d12s21
                  do iv=nbot,nvirt(isbv2)-1                              1d11s21
                   irec=iv1+nvirt(isbv1)*iv                              1d11s21
                   itri=((iv*(iv-1))/2)+iv1                              1d18s21
                   irow=itri+isw*(irec-itri)                             1d11s21
                   term=bc(ih0+iv)*bc(jtmp)*2d0                         10d16s23
                   gd(iad+irow)=gd(iad+irow)+term                       9d28s23
                   bc(ih0d+iv)=bc(ih0d+iv)+vd(iad+irow)*bc(jtmp)        9d28s23
                  end do                                                 1d11s21
                  if(isbv1.eq.isbv2)then
c     Gvv'ri=dij Hvv" Vv'v"rj
                   ntop=iv1-1                                             1d18s21
                   ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                    1d18s21
                   ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*(           1d12s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                   ih0d=iden(isbv1)+irefo(isbv1)+nh0av(isbv1)*(         9d28s23
     $                irefo(isbv1)+iv2)                                 1d18s21
                   iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                   do iv=0,ntop                                           1d11s21
                    irec=iv+nvirt(isbv1)*iv1                              1d11s21
                    itri=((iv1*(iv1-1))/2)+iv                             1d18s21
                    irow=itri+isw*(irec-itri)                             1d11s21
                    term=bc(ih0+iv)*bc(jtmp)*tf*2d0                     10d16s23
                    gd(iad+irow)=gd(iad+irow)+term                      9d28s23
                    bc(ih0d+iv)=bc(ih0d+iv)+vd(iad+irow)*bc(jtmp)*tf    9d28s23
                   end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vv"vrj                                           1d18s21
                   nbot=iv2+1                                             1d18s21
                   nbot=nbot+isw*(0-nbot)                                 1d18s21
                   ih0=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*            1d18s21
     $                  (irefo(isbv2)+iv1)                                1d18s21
                   ih0d=iden(isbv2)+irefo(isbv2)+nh0av(isbv2)*          9d28s23
     $                  (irefo(isbv2)+iv1)                                1d18s21
                   do iv=nbot,nvirt(isbv2)-1                              1d18s21
                    irec=iv2+nvirt(isbv1)*iv                              1d18s21
                    itri=((iv*(iv-1))/2)+iv2                              1d18s21
                    irow=itri+isw*(irec-itri)                             1d18s21
                    term=bc(ih0+iv)*bc(jtmp)*tf*2d0                     10d16s23
                    gd(iad+irow)=gd(iad+irow)+term                      9d28s23
                    bc(ih0d+iv)=bc(ih0d+iv)+vd(iad+irow)*bc(jtmp)*tf    9d28s23
                   end do                                                 1d11s21
                  end if                                                 1d18s21
                  jtmp=jtmp+1                                            1d11s21
                 end do                                                  1d11s21
                 i10=1                                                   1d11s21
                end do                                                   1d11s21
               end do                                                    1d11s21
               ibcoff=itmp                                               1d12s21
              end if                                                     1d11s21
              ioff=ioff+nnn*nfdat(2,l,isb)                               1d11s21
              ioffd=ioffd+nvv*nfdat(2,l,isb)                             1d11s21
             end if                                                      1d11s21
             tf=-1d0                                                     1d11s21
            end do                                                       1d11s21
           end if                                                        1d11s21
          end do                                                         1d11s21
         end if                                                          1d11s21
        end if                                                          5d10s21
        ibcoff=idcjb4                                                   10d2s23
        do jsbv1=1,nsymb                                                1d8s21
         jsbv2=multh(jsbv1,jsbv12)                                      1d8s21
         if(jsbv2.ge.jsbv1)then                                         1d8s21
          if(jsbv12.eq.1)then                                           1d8s21
           joffdnon=joffdnon+nrootu*nvirt(jsbv1)*nfdat(2,1,jsb)         1d8s21
           nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                        1d8s21
          else                                                          1d8s21
           nvv=nvirt(jsbv1)*nvirt(jsbv2)                                1d8s21
          end if                                                        1d8s21
          nvv=nvv*nrootu                                                1d8s21
          do l=1,4                                                      1d8s21
           joffdnon=joffdnon+nvv*nfdat(2,l,jsb)                         1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
       end do                                                           1d8s21
       do isbv1=1,nsymb                                                 1d8s21
        isbv2=multh(isbv1,isbv12)                                       1d8s21
        if(isbv2.ge.isbv1)then                                          1d8s21
         if(isbv12.eq.1)then                                            1d8s21
          ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)          1d8s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         1d8s21
         else                                                           1d8s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 1d8s21
         end if                                                         1d8s21
         nvv=nvv*nrootu                                                 1d8s21
         do l=1,4                                                       1d8s21
          ioffdnon=ioffdnon+nvv*nfdat(2,l,isb)                          1d8s21
         end do                                                         1d8s21
        end if                                                          1d8s21
       end do                                                           1d8s21
      end do                                                            1d8s21
      call dws_synca
      ioffn=0                                                           9d13s23
      ioffo=0                                                           9d13s23
      do isb=1,nsymb                                                    9d13s23
       isbv12=multh(isb,isymmrci)                                       9d13s23
       do isbv1=1,nsymb                                                 9d13s23
        isbv2=multh(isbv1,isbv12)                                       9d13s23
        if(isbv2.ge.isbv1)then                                          9d13s23
         if(isbv1.eq.isbv2)then                                         9d13s23
          do ifo=0,nfdat(3,1,isb)-1                                     9d13s23
           do ifn=0,nfdat(2,1,isb)-1                                    9d13s23
            do ir=0,nrootu-1                                            9d13s23
             sum=0d0                                                     9d13s23
             iadn=ioffn+nvirt(isbv1)*(ir+nrootu*ifn)                    9d13s23
             iado=ioffo+nvirt(isbv1)*(ir+nrootu*ifo)                    9d13s23
             do i=1,nvirt(isbv1)                                        9d13s23
              sum=sum+veco(iado+i)*gd(iadn+i)                           9d13s23
             end do                                                     9d13s23
             iad=idorth(1,isb)+ifn+nfdat(2,1,isb)*(ifo                  9d13s23
     $            +nfdat(3,1,isb)*ir)                                   9d13s23
             bc(iad)=bc(iad)+sum                                        9d13s23
            end do                                                      9d13s23
           end do                                                       9d13s23
          end do                                                        9d13s23
          ioffn=ioffn+nvirt(isbv1)*nrootu*nfdat(2,1,isb)                 9d13s23
          ioffo=ioffo+nvirt(isbv1)*nrootu*nfdat(3,1,isb)                 9d13s23
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         9d13s23
         else                                                           9d13s23
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 9d13s23
         end if                                                         9d13s23
         do l=1,4                                                       9d13s23
          do ifo=0,nfdat(3,l,isb)-1                                     9d13s23
           do ifn=0,nfdat(2,l,isb)-1                                    9d13s23
            do ir=0,nrootu-1                                             9d13s23
             sum=0d0                                                     9d13s23
             iadn=ioffn+nvv*(ir+nrootu*ifn)                             9d13s23
             iado=ioffo+nvv*(ir+nrootu*ifo)                             9d13s23
             do i=1,nvv                                                 9d13s23
              sum=sum+veco(iado+i)*gd(iadn+i)                           9d13s23
             end do                                                     9d13s23
             iad=idorth(l,isb)+ifn+nfdat(2,l,isb)*(ifo                  9d13s23
     $            +nfdat(3,l,isb)*ir)                                   9d13s23
             bc(iad)=bc(iad)+sum                                        9d13s23
            end do                                                      9d13s23
           end do                                                       9d13s23
          end do                                                        9d13s23
          ioffn=ioffn+nvv*nrootu*nfdat(2,l,isb)                         9d13s23
          ioffo=ioffo+nvv*nrootu*nfdat(3,l,isb)                         9d13s23
         end do                                                         9d13s23
        end if                                                          9d13s23
       end do                                                           9d13s23
      end do                                                            9d13s23
      return                                                            1d8s21
      end                                                               1d8s21
c     nfdat(2,l,jsb)
c     nok=nokj(l,lsa)
c     iok=ndenjf(l,lsa)
c     iden=idenjf(l,lsa)
c
      subroutine hcddjkd123(nrow,nmul,bc,ibc,ibl,nvirt,jsw,isbl,isbvc,  10d3s23
     $     isbr,nrootm,jjoffp,nvv,nfdat,vd,nnvv,lvisv,nfdati,nhere2,tf, 10d3s23
     $     i1s,i1e,i2s,i2e,itrans,nok,iok,ijsb,nsymb,lkmat,iden,inv,    10d3s23
     $     imden,imats,nhere,irefo,multh,itmp2,srh,sum,ilab)            10d3s23
      implicit real*8 (a-h,o-z)
      include "common.store"
      logical lvisv,lkmat                                                     10d3s23
      dimension nvirt(*),vd(*),iden(*),inv(2,8,8,8),nok(*),iok(*),
     $     imden(*),imats(*),irefo(*),multh(8,*)                        10d3s23
      igoal=1
      mn=nfdat*nfdati                                                   10d3s23
      nrootu=nrootm+1                                                   10d3s23
      itmp3=ibcoff                                                      10d3s23
      ibcoff=itmp3+nvirt(isbl)*nfdat*nvirt(isbvc)*nrootu                4d24s24
      call enough('hcddjkd123.tmp3',bc,ibc)                             10d3s23
      do iz=itmp3,ibcoff-1                                              10d3s23
       bc(iz)=0d0                                                       10d3s23
      end do                                                            10d3s23
      if(ibl.eq.1)then                                                  10d3s23
       do iv2=0,nvirt(isbvc)-1                                          10d3s23
        ntop=iv2-1                                                      10d3s23
        ntop=ntop+jsw*(nvirt(isbl)-1-ntop)                              10d3s23
        do ir=0,nrootm                                                  10d3s23
         do iv1=0,ntop                                                  10d3s23
          irec=iv1+nvirt(isbl)*iv2                                      10d3s23
          itri=((iv2*(iv2-1))/2)+iv1                                    10d3s23
          itri=itri+jsw*(irec-itri)                                     10d3s23
          iad=jjoffp+itri+nvv*ir                                        10d3s23
          jtmp3=itmp3+iv1+nvirt(isbl)*nfdat*(iv2+
     $                      nvirt(isbvc)*ir)                             10d2s23
          do j=0,nfdat-1                                                10d3s23
           bc(jtmp3)=vd(iad+j*nnvv)                                     10d3s23
           jtmp3=jtmp3+nvirt(isbl)                                      10d3s23
          end do                                                        10d3s23
         end do                                                         10d3s23
        end do                                                          10d3s23
       end do                                                           10d3s23
      else                                                              10d3s23
       do iv1=0,nvirt(isbvc)-1                                          10d3s23
        nbot=iv1+1                                                      10d3s23
        nbot=nbot+jsw*(0-nbot)                                          10d3s23
        do ir=0,nrootm                                                  10d3s23
         do iv2=nbot,nvirt(isbl)-1                                      10d3s23
          irec=iv1+nvirt(isbvc)*iv2                                     10d3s23
          itri=((iv2*(iv2-1))/2)+iv1                                    10d3s23
          itri=itri+jsw*(irec-itri)                                     10d3s23
          iad=jjoffp+itri+nvv*ir                                        10d3s23
          jtmp3=itmp3+iv2+nvirt(isbl)*nfdat*(iv1                        10d3s23
     $                      +nvirt(isbvc)*ir)                           10d2s23
          do j=0,nfdat-1                                                10d3s23
           bc(jtmp3)=bc(jtmp3)+vd(iad+j*nnvv)                           10d3s23
           jtmp3=jtmp3+nvirt(isbl)                                      10d3s23
          end do                                                        10d3s23
         end do                                                         10d3s23
        end do                                                          10d3s23
       end do                                                           10d3s23
      end if                                                            10d3s23
      if(lvisv)then                                                     10d3s23
       jjoff=jjoffp-nvirt(isbvc)*nrootu*nfdat                           10d3s23
       nv=nrootu*nvirt(isbvc)                                           10d3s23
       do iv1=0,nvirt(isbvc)-1                                          10d3s23
        do ir=0,nrootm                                                  10d3s23
         iad=jjoff+iv1+nvirt(isbvc)*ir                                  10d3s23
         jtmp3=itmp3+iv1+nvirt(isbl)*nfdat*(iv1+nvirt(isbvc)*ir)        10d3s23
         do j=0,nfdat-1                                                 10d3s23
          bc(jtmp3)=bc(jtmp3)+vd(iad+j*nv)*srh                          10d3s23
          jtmp3=jtmp3+nvirt(isbl)                                       10d3s23
         end do                                                         10d3s23
        end do                                                          10d3s23
       end do                                                           10d3s23
      end if                                                            10d3s23
      itmp4=ibcoff                                                      10d3s23
      itmp5=itmp4+nvirt(isbvc)*nhere2*nfdati                            10d3s23
      ibcoff=itmp5+nvirt(isbr)*nvirt(isbl)*nfdat*nfdati                 10d3s23
      call enough('hcddjkd123.tmp5',bc,ibc)                             10d3s23
      do ir=0,nrootm                                                    10d3s23
       jtmp3=itmp3+nvirt(isbl)*nfdat*nvirt(isbvc)*ir                    10d3s23
       do ivc=0,nvirt(isbvc)-1                                          10d3s23
        jtmp2=itmp2+nhere2*nfdati*(ir+nrootu*ivc)                       10d3s23
        do k=0,nfdati-1                                                 10d3s23
         do i2=i2s,i2e                                                  10d3s23
          ivr=i2-i2s                                                    10d3s23
          jtmp4=itmp4+ivc+nvirt(isbvc)*(k+nfdati*ivr)                   10d3s23
          bc(jtmp4)=bc(jtmp2+ivr)                                       10d3s23
         end do                                                         10d3s23
         jtmp2=jtmp2+nhere2                                             10d3s23
        end do                                                          10d3s23
       end do                                                           10d3s23
       call dgemm('n','n',nrow,nmul,nvirt(isbvc),tf,                    10d3s23
     $      bc(jtmp3),nrow,bc(itmp4),nvirt(isbvc),0d0,                  10d3s23
     $      bc(itmp5),nrow,'hcddjkd123.tmp5')                           10d3s23
       itmp6=ibcoff                                                     10d3s23
       ibcoff=itmp6+nhere*mn                                            10d3s23
       call enough('hcddjkd12.tmp6',bc,ibc)                             10d3s23
       i10=i1s                                                          10d3s23
       i1n=nvirt(isbl)                                                  10d3s23
       jtmp6=itmp6                                                      10d3s23
       do i2=i2s,i2e                                                    10d3s23
        i2n=i2-i2s                                                      10d3s23
        if(i2.eq.i2e)i1n=i1e                                            10d3s23
        do i1=i10,i1n                                                   10d3s23
         i1m=i1-1                                                       10d3s23
         do kj=0,nfdat-1                                                10d3s23
          do ki=0,nfdati-1                                              10d3s23
           kij=ki+nfdati*kj                                              10d3s23
           kji=kj+nfdat*ki                                              10d3s23
           kji=kij+itrans*(kji-kij)                                     10d3s23
           iad1=jtmp6+nhere*kij                                         10d3s23
           iad2=itmp5+i1m+nvirt(isbl)*(kji+mn*i2n)                      10d3s23
           bc(iad1)=bc(iad2)                                            10d3s23
          end do                                                        10d3s23
         end do                                                         10d3s23
         jtmp6=jtmp6+1                                                  10d3s23
        end do                                                          10d3s23
        i10=1                                                           10d3s23
       end do                                                           10d3s23
       do lsa=1,nsymb                                                   10d3s23
        lsb=multh(lsa,ijsb)                                             10d3s23
        if(min(nok(lsa),irefo(lsa),irefo(lsb)).gt.0.and.                     10d3s23
     $       (lsb.ge.lsa.or.lkmat))then                                 10d3s23
         if(lsa.eq.lsb.and..not.lkmat)then                               10d3s23
          mrow=(irefo(lsa)*(irefo(lsa)+1))/2                             10d3s23
         else                                                            10d3s23
          mrow=irefo(lsa)*irefo(lsb)                                     10d3s23
         end if                                                          10d3s23
         itmp7=ibcoff                                                    10d3s23
         ibcoff=itmp7+nhere*nok(lsa)                                          10d3s23
         call enough('hcddjkd123.tmp7',bc,ibc)                           10d3s23
         call dgemm('n','n',nhere,nok(lsa),mn,1d0,                            10d3s23
     $             bc(itmp6),nhere,bc(iden(lsa)),mn,0d0,                10d3s23
     $       bc(itmp7),nhere,'hcddjk123.tmp7')                          10d3s23
         if(ilab.eq.0)then                                              10d3s23
          i2eu=inv(1,lsa,lsb,isbl)                                        10d3s23
          icase=inv(2,lsa,lsb,isbl)                                       10d3s23
         else                                                           10d3s23
          i2eu=inv(1,lsb,lsa,isbl)                                      10d3s23
          icase=inv(2,lsb,lsa,isbl)                                     10d3s23
         end if                                                         10d3s23
         if(icase.eq.1.and.ilab.eq.0)then                               10d3s23
          do j=0,nok(lsa)-1                                             10d3s23
           jj=ibc(iok(lsa)+j)                                           10d3s23
           ij=imden(i2eu)+nhere*(jj+mrow*ir)                            10d3s23
           ll=imats(i2eu)+nhere*jj                                      10d3s23
           ji=itmp7+nhere*j                                             10d3s23
           do i=0,nhere-1                                               10d3s23
                   orig=bc(igoal)
            bc(ij+i)=bc(ij+i)+bc(ji+i)                                  10d3s23
           end do                                                       10d3s23
          end do                                                        10d3s23
         else                                                           10d3s23
          do j=0,nok(lsa)-1                                             10d3s23
           jj=ibc(iok(lsa)+j)                                           10d3s23
           jb=jj/irefo(lsa)                                             10d3s23
           ja=jj-irefo(lsa)*jb                                          10d3s23
           kk=jb+irefo(lsb)*ja                                          10d3s23
           ij=imden(i2eu)+nhere*(kk+mrow*ir)                            10d3s23
           ll=imats(i2eu)+nhere*kk                                      10d3s23
           ji=itmp7+nhere*j                                             10d3s23
           do i=0,nhere-1                                               10d3s23
                   orig=bc(igoal)
            bc(ij+i)=bc(ij+i)+bc(ji+i)                                  10d3s23
           end do                                                       10d3s23
          end do                                                        10d3s23
         end if                                                         10d3s23
         ibcoff=itmp7                                                   10d3s23
        end if                                                          10d3s23
       end do                                                           10d3s23
      end do                                                            10d3s23
      ibcoff=itmp3                                                      10d3s23
      return
      end                                                               10d3s23
