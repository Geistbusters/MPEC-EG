c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prehcddjk(nff22,nfdat,nsymb,mdon,mdoo,nec,multh,ncsf,  11d7s22
     $     ncsf2,irel,ism,irefo,ismult,ixw1,ixw2,norb,ih0av,nh0av,ioooo,11d7s22
     $     shift,iden1e,idenhvv,idenj,idenk,bc,ibc,idenput,idenorg)     5d17s23
      implicit real*8 (a-h,o-z)                                         1d8s21
c                                                                       1d8s21
c     precompute densities for 0 and 2 virt contribution to gd=hdd*vd   11d7s22
c                                                                       1d8s21
      logical lprt,lchoice,lq                                              3d19s21
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d22s21
      integer*8 ipack8,itc,ito,jtc,jto,itestc,itesto,itesta,itestb,     1d17s21
     $     jtesta,jtestb,gandcc,gandco,gandcb,idenput(nsymb,nsymb,2)    5d17s23
      external second                                                   2d18s21
      equivalence (ipack8,ipack4)                                       1d8s21
      dimension nff22(mdoo+1,2,nsymb),nfdat(5,4,*),multh(8,8),ncsf(*),
     $     ncsf2(4,*),irel(*),ism(*),irefo(*),ipack4(2),nprei(4),       11d7s22
     $     mprei(4),nli(4),nprej(4),ih0av(*),nh0av(*),ioooo(*),         11d7s22
     $     mprej(4),nlj(4),nab4(2,3),idenj(4,8,8,8),ndenj(4,8),ioxx(2), 11d7s22
     $     mdenj(4,8),idenk(4,4,8,8,8),ndenk(4,8),mdenk(4,8),loope(6),    11d7s22
     $     idenhvv(4,8,*),mdenhvv(4),itest(32,3),iden1e(4,*),mden1e(4), 11d7s22
     $     ldenj(4,8),ldenk(4,8),iden1ef(4),idenhvvf(4),ivects(4),      11d3s22
     $     idenjf(4,8),idenkf(4,4,8),ndenjf(4,8),ndenkf(4,4,8),
     $     nokj(4,8),nokk(4,4,8),isorb(32),idorb(32),idenjn(4,8),       2d23s21
     $     idenkn(4,4,8),iden1en(4),idenhvvn(4)                         11d7s22
      include "common.store"                                            1d8s21
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/fnd2cm/inv(2,8,8,8)                                        9d2s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      ibcoffo=ibcoff
      igoal=ibcoffo+3628+3*3
      norbx=norb+1                                                      1d8s21
      norbxx=norbx+1                                                    1d8s21
      norbxxx=norbxx+1                                                  1d8s21
      ioffd=1                                                           1d11s21
c
c     ibc0
c     iden1e l,isb
c     idenhvv l,isb,jsb
c     idenj l,lsa,isb,jsb
c     idenk l,ll,lsa,isb,jsb
c     ibctopr
      do isb=1,nsymb                                                    1d8s21
       do jsb=1,isb                                                     12d9s22
        ijsb=multh(jsb,isb)                                             1d9s21
        ibc0=ibcoff                                                     1d8s21
        njc0=0
        njc1=0
        nkc0=0
        nkc1=0
        do l=1,4                                                        11d7s22
         if(nfdat(2,l,isb).gt.0)then                                    11d7s22
          if(isb.eq.jsb)then                                            11d7s22
           iden1e(l,isb)=ibcoff                                         11d7s22
           ibcoff=ibcoff+nfdat(2,l,isb)*nfdat(2,l,isb)                  11d7s22
          end if                                                        11d7s22
          idenhvv(l,isb,jsb)=ibcoff                                     11d7s22
          ibcoff=idenhvv(l,isb,jsb)+nfdat(2,l,isb)*nfdat(2,l,jsb)       12d10s22
         end if                                                         11d7s22
        end do                                                          11d7s22
        ibctopr=ibcoff                                                  11d7s22
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    2d24s21
          do ll=1,4                                                     1d8s21
           if(nfdat(2,ll,jsb).gt.0)then                                 2d24s21
            nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                          2d24s21
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
        call enough('prehcddjk.  1',bc,ibc)
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
             call enough('prehcddjk.1.1',bc,ibc)
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
              jftop=nff22(nclojp,1,jsb)                                 11d7s22
             do jf=1,jftop
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
                   xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)      11d15s22
                   xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)      11d15s22
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
                    xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)     11d15s22
                    sumo=sumo+xj*0.5d0                                            7d15s19
                   end if                                                         7d15s19
                  end do                                                          7d15s19
                 end do                                                           7d15s19
                 sum=sumc+sumo+shift
                 do i=1,ncloi                                                   11d23s20
                  is=ism(idorb(i))                                              11d23s20
                  ig=irel(idorb(i))-1                                           11d23s20
                  do j=1,nopeni                                                 11d23s20
                   js=ism(isorb(j))                                             11d23s20
                   jg=irel(isorb(j))-1                                          11d23s20
                   xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)      11d15s22
                   xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)      11d15s22
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
                    xint=getint(ioooo,jsa,ksb,ksb,jsa,jga,kgb,kgb,jga,  11d15s22
     $                   bc,ibc)                                        11d15s22
                    sumx=sumx-xint                                      2d24s21
                    if(nqq.ge.mdon.and.nqq.le.mdoo)then                         11d19s20
                     jtestb=ibclr(jtestb,isorb(i2))                      2d10s21
                     jtestb=ibclr(jtestb,isorb(i1))                      2d10s21
                     jtesta=ibset(jtesta,isorb(i1))                      2d10s21
                     call gandc(itc,ito,jtesta,jtestb,nopenip,nopenk,      11d19s20
     $             iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot,nab,iwpb,iwpk,       1d17s21
     $                   ncsfmid,bc,ibc)                                11d14s22
                     call gandc(jtesta,jtestb,itc,ito,nopenk,nopenip,   2d25s21
     $             karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,    2d25s21
     $                     iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                     iargo=1                                            2d24s21
                     do l=1,4                                           2d24s21
                      if(nli(l).gt.0)then                               2d24s21
                       iad1=ivcv+ibc(ivcv+5+l)                            3d19s21
                       iad2=iad1+nli(l)                                    3d19s21
                       itmp1=ibcoff                                     2d24s21
                       itmp2=itmp1+nli(l)*ncsf(karg)                    2d24s21
                       itmp3=itmp2+nli(l)*ncsf(karg)                    2d24s21
                       ibcoff=itmp3+nli(l)*nli(l)                       2d24s21
                       call enough('prehcddjk.  2',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1,    2d25s21
     $                       nli(l),iwpb1,iwpk1,bc(iad2),ncsf2(l,iarg), 2d25s21
     $                       bc(itmp1),ncsf(karg),1d0,0d0,iargo,        2d24s21
     $                       ncsf2(l,iarg),bc,ibc)                      11d10s22
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
     d' prehcddjk.  1')
                       jtmp3=itmp3                                      2d24s21
                       do j=0,nli(l)-1                                   2d23s21
                        jj=ibc(iad1+j)-1                                  1d8s21
                        jden=iden1e(l,isb)+nfdat(2,l,isb)*jj-1          11d7s22
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
                   call enough('prehcddjk.  3',bc,ibc)
                   iad1=ivcv+ibc(ivcv+5+l)                              3d19s21
                   iad2=iad1+nli(l)                                     3d19s21
                   do iii=0,nli(l)-1                                    2d24s21
                    do j=0,ncsf2(l,iarg)-1                              2d24s21
                     ji=iad2+j+ncsf2(l,iarg)*iii                        2d24s21
                     ij=itmpt+iii+nli(l)*j                              2d24s21
                     bc(ij)=bc(ji)                                      2d24s21
                    end do                                              2d24s21
                   end do                                               2d24s21
                   call dgemm('n','n',nli(l),nli(l),ncsf2(l,iarg),sumx, 11d7s22
     $                   bc(itmpt),nli(l),bc(iad2),ncsf2(l,iarg),0d0,   2d24s21
     $                   bc(itmpm),nli(l),                              2d24s21
     d' prehcddjk.  2')
                   jtmpm=itmpm                                          2d24s21
                   do j=0,nli(l)-1                                      2d24s21
                    jj=ibc(iad1+j)-1                                    2d24s21
                    jden=iden1e(l,isb)+nfdat(2,l,isb)*jj-1              11d7s22
                    do k=0,nli(l)-1                                     2d24s21
                     kk=ibc(iad1+k)                                     2d24s21
                     bc(jden+kk)=bc(jden+kk)+bc(jtmpm+k)                2d24s21
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
     $                 (btest(itc,i).and.btest(jto,i)))then                           10d14s22
                    nab4(1,1)=i
                   else                                                           10d14s22
                    nab4(2,1)=i
                   end if
                  end if                                                          10d14s22
                 end do                                                           10d14s22
                 call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,  1d8s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $               ncsfmid1,bc,ibc)                                   11d14s22
                 call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,  2d24s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,  2d24s21
     $               ncsfmid1b,bc,ibc)                                  11d14s22
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
                  xj=getint(ioooo,lsa,lsa,ksa,ksb,lga,lga,kga,kgb,bc,   11d15s22
     $                 ibc)                                             11d15s22
                  sum=sum+xj                                            1d8s21
                  if(itest(i,3).eq.2)then                               1d8s21
                   xk=getint(ioooo,lsa,ksa,lsa,ksb,lga,kga,lga,kgb,bc,  11d15s22
     $                 ibc)                                             11d15s22
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
                   call enough('prehcddjk.  4',bc,ibc)
                   if(iwpk1.lt.0)then                                   2d23s21
                    nusedi=ibc(-iwpk1)/2                                2d23s21
                    if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1          2d23s21
                    icmp1=-iwpk1+1                                      2d23s21
                    icmp2=icmp1+ncsf(jarg)                              2d23s21
                    icmp3=icmp2+nusedi                                  2d23s21
                    call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   bc(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),   2d23s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d14s22
                   else                                                              12d5s20
                    jwpk1=iwpk1+ncsfmid1*(jargo-1)                      2d23s21
                    call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                   1d0,bc(jwpk1),ncsfmid1,bc(jad2),ncsf2(l,jarg),
     $                   0d0,bc(itmp1),ncsfmid1,
     d' prehcddjk.  3')
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
     $                   ncsf2(l,iarg),bc,ibc)                          11d14s22
                   else                                                              12d5s20
                    jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                    call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                   1d0,bc(jwpk1b),ncsfmid1,bc(iad2),ncsf2(l,iarg),
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' prehcddjk.  4')
                   end if                                                            12d5s20
                   itmpt=ibcoff                                         2d23s21
                   itmpm=itmpt+nli(l)*ncsfmid1                          2d24s21
                   ibcoff=itmpm+nli(l)*nlj(l)                           2d24s21
                   call enough('prehcddjk.  5',bc,ibc)
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
     d' prehcddjk.  5')
                   jtmpm=itmpm                                          2d24s21
                   do iii=0,nlj(l)-1                                    2d24s21
                    jj=ibc(jad1+iii)-1                                  2d24s21
                    jjden=iden1e(l,isb)+nfdat(2,l,isb)*jj-1             11d7s22
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
     $                iwpb1b,iwpk1b,ncsfmid1b,bc,ibc)                   11d14s22
                    nnot1=nnot1b                                        2d25s21
                    nab1(1)=nab1b(2)                                    2d25s21
                    nab1(2)=nab1b(1)                                    2d25s21
                    call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,        12d14s20
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            11d14s22
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
                    xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  11d15s22
     $                   bc,ibc)                                        11d15s22
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
                      call enough('prehcddjk.  6',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,bc(jad2),ncsf2(l,jarg),     2d24s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d10s22
                      itmpi=ibcoff                                      2d24s21
                      itmpt=itmpi+ncsf(karg)*nli(l)                     2d24s21
                      ibcoff=itmpt+ncsf(karg)*nli(l)                    2d24s21
                      call enough('prehcddjk.  7',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,bc(iad2),ncsf2(l,iarg),   2d24s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          11d10s22
                      do iii=0,nli(l)-1                                 2d24s21
                       do k=0,ncsf(karg)-1                                2d23s21
                        ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                        ik=itmpi+iii+nli(l)*k                            2d23s21
                        bc(ik)=bc(ki)                                     2d23s21
                       end do                                             2d23s21
                      end do                                              2d23s21
                      itmpm=itmpt                                       2d24s21
                      ibcoff=itmpm+nli(l)*nlj(l)                        2d24s21
                      call enough('prehcddjk.  8',bc,ibc)
                      call dgemm('n','n',nli(l),nlj(l),ncsf(karg),      2d25s21
     $                     xint,                                        2d25s21
     $                      bc(itmpi),nli(l),bc(itmpj),ncsf(karg),0d0,  2d24s21
     $                      bc(itmpm),nli(l),                           2d24s21
     d' prehcddjk.  6')
                      jtmpm=itmpm                                       2d24s21
                      do iii=0,nlj(l)-1                                   2d24s21
                       jj=ibc(jad1+iii)-1                                 2d24s21
                       jjden=iden1e(l,isb)+nfdat(2,l,isb)*jj-1          11d7s22
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
                 if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s2
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
     $         ncsfmid2,bc,ibc)                                         11d14s22
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
     $                  iwpk1b,ncsfmid1b,bc,ibc)                        11d14s22
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
                    xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  11d15s22
     $                   bc,ibc)                                        11d15s22
                    if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2))     4d20s21
     $                   .and.                                          4d20s21
     $                max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))     1d12s21
     $                  xint=xint*0.5d0                                 1d12s21
                    if(abs(xint).gt.1d-14)then                          11d3s22
                     nplanb=1
                     nsz=ncsfmid1b*ncsfmid2                             11d3s22
                     iprod=ibcoff                                        11d3s22
                     ibcoff=iprod+nsz                                   11d3s22
                     call enough('prehcddjk.8.1',bc,ibc)
                     call prodn(iwpk1b,iwpb2,ncsfmid1b,ncsfmid2,         11d3s22
     $                  ncsf(karg),bc(iprod),bc,ibc,1d0,0d0)            2d13s23
                     itrans=ibcoff
                     ibcoff=itrans+ncsfmid1b*ncsfmid2
                     call enough('prehcddjk.8.3',bc,ibc)
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
                     call enough('prehcddjk.8.2',bc,ibc)
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
                       call enough('prehcddjk.8.2',bc,ibc)
c     nli,ncsf(iarg),ncsfmid1b,ncsf(karg),ncsfmid2,ncsf(jarg),nlj
c       ivects     iwpb1b    iprod     iwpk2    jad2
c      nli,ncsf(iarg),ncsfmid1b,ncsfmid2,ncsf(jarg),nlj
c       1    2           3         4         5       6
                       call genmatn3(nli(l),ncsfmid1b,ncsf(jarg),        11d3s22
     $                     ncsf(iarg),ivects(l),iwpb1b,ncsfmid2,        11d3s22
     $                     iwprod,iwpk2,bc(jad2),ncsf2(l,jarg),nlj(l),  11d3s22
     $                     bc(itmp),0d0,jcsf0,ncsf2(l,jarg),            11d3s22
     $                     icsf0,ncsf2(l,iarg),xint,bc,ibc)             11d14s22
                       do j=1,nlj(l)
                        jm=j-1                                          11d4s22
                        jj=ibc(jad1+jm)-1                               11d4s22
                        jjden=iden1e(l,isb)+nfdat(2,l,isb)*jj-1         11d7s22
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
     $                  iwpk1b,ncsfmid1b,bc,ibc)                        11d14s22
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
                    xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  11d15s22
     $                   bc,ibc)                                        11d15s22
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
                      call enough('prehcddjk.  9',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,bc(jad2),ncsf2(l,jarg),     2d24s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d10s22
                      itmpi=ibcoff                                      2d24s21
                      itmpt=itmpi+ncsf(karg)*nli(l)                     2d24s21
                      ibcoff=itmpt+ncsf(karg)*nli(l)                    2d24s21
                      call enough('prehcddjk. 10',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,bc(iad2),ncsf2(l,iarg),   2d24s21
     $                   bc(itmpt),ncsf(karg),xint,0d0,iargo,           11d2s22
     $                   ncsf2(l,iarg),bc,ibc)                          11d10s22
                      do iii=0,nlj(l)-1                                 11d2s22
                       iadj=itmpj+ncsf(karg)*iii                        11d2s22
                       jj=ibc(jad1+iii)-1                                 2d24s21
                       jjden=iden1e(l,isb)+nfdat(2,l,isb)*jj-1          11d7s22
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
    3            continue                                                 12d18s20
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
                 stop 'prehcddjk:g'                                                   1d9s21
                end if                                                  1d9s21
                call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,   1d8s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $              ncsfmid1,bc,ibc)                                    11d14s22
                call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,   2d23s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,      1d8s21
     $              ncsfmid1b,bc,ibc)                                   11d14s22
                itmpn=ibcoff                                            2d23s21
                nlij=0                                                  2d23s21
                do l=1,4                                                2d23s21
                 nlij=nlij+nli(l)*nlj(l)                                2d23s21
                end do                                                  2d23s21
                itmpn=ibcoff                                            2d23s21
                ibcoff=itmpn+nlij                                       2d23s21
                call enough('prehcddjk. 12',bc,ibc)
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
                  call enough('prehcddjk. 13',bc,ibc)
                  if(iwpk1.lt.0)then                                    2d23s21
                   nusedi=ibc(-iwpk1)/2                                 2d23s21
                   if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1           2d23s21
                   icmp1=-iwpk1+1                                       2d23s21
                   icmp2=icmp1+ncsf(jarg)                               2d23s21
                   icmp3=icmp2+nusedi                                   2d23s21
                   call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),    2d23s21
     $                   bc(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),   2d23s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d14s22
                  else                                                              12d5s20
                   jwpk1=iwpk1+ncsfmid1*(jargo-1)                       2d23s21
                   call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                   1d0,bc(jwpk1),ncsfmid1,bc(jad2),ncsf2(l,jarg),
     $                   0d0,bc(itmp1),ncsfmid1,
     d' prehcddjk.  8')
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
     $                   ncsf2(l,iarg),bc,ibc)                          11d14s22
                  else                                                              12d5s20
                   jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                   call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                   1d0,bc(jwpk1b),ncsfmid1,bc(iad2),ncsf2(l,iarg),
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' prehcddjk.  9')
                  end if                                                            12d5s20
                  itmpt=ibcoff                                          2d23s21
                  ibcoff=itmpt+nli(l)*ncsfmid1                          2d23s21
                  call enough('prehcddjk. 14',bc,ibc)
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
     d' prehcddjk. 10')
                  ibcoff=itmp1                                          2d23s21
                  ktmpn=jtmpn                                           2d24s21
                  do iii=0,nlj(l)-1                                     2d24s21
                   jj=ibc(jad1+iii)-1                                   2d24s21
                   jjden=idenhvv(l,isb,jsb)+nfdat(2,l,isb)*jj-1         11d7s22
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
                      bc(kkden+kk)=bc(kkden+kk)-bc(jtmpn+j)             12d9s22
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
     $               iwpk1b,ncsfmid1b,bc,ibc)                           11d14s22
                   nnot1=nnot1b                                         2d25s21
                   nab1(1)=nab1b(2)                                     2d25s21
                   nab1(2)=nab1b(1)                                     2d25s21
                   call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,        12d14s20
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            11d14s22
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
                    stop 'prehcddjk:i'
                   end if                                               1d9s21
                   jargo=1                                               2d23s21
                   itmpj=ibcoff                                          2d23s21
                   jtmpj=itmpj                                           2d23s21
                   do lj=1,4                                             2d23s21
                    if(nlj(lj).gt.0)then                                 2d23s21
                     jad1=jvcv+ibc(jvcv+5+lj)                              3d19s21
                     jad2=jad1+nlj(lj)                                     3d19s21
                     ibcoff=jtmpj+ncsf(karg)*nlj(lj)                     2d23s21
                     call enough('prehcddjk. 15',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                    nlj(lj),iwpb2,iwpk2,bc(jad2),ncsf2(lj,jarg),   2d23s21
     $                   bc(jtmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(lj,jarg),bc,ibc)                         11d10s22
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
                     call enough('prehcddjk. 16',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                    nli(li),iwpb1b,iwpk1b,bc(iad2),ncsf2(li,iarg), 2d23s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(li,iarg),bc,ibc)                         11d10s22
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
                       call enough('prehcddjk. 17',bc,ibc)
                       call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),   2d24s21
     $                     1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),  11d7s22
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' prehcddjk. 11')
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
     $                iwpk1b,ncsfmid1b,bc,ibc)                          11d14s22
                  nnot1=nnot1b                                           2d25s21
                  nab1(1)=nab1b(2)                                       2d25s21
                  nab1(2)=nab1b(1)                                       2d25s21
                  call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,       1d8s21
     $         karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d18s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
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
                     call enough('prehcddjk. 18',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(lj),iwpb2,iwpk2,bc(jad2),ncsf2(lj,jarg),   2d23s21
     $                   bc(jtmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(lj,jarg),bc,ibc)                         11d10s22
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
                     call enough('prehcddjk. 19',bc,ibc)
                     call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(li),iwpb1b,iwpk1b,bc(iad2),ncsf2(li,iarg), 2d23s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(li,iarg),bc,ibc)                         11d10s22
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
                        call enough('prehcddjk. 20',bc,ibc)
                        call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),  2d24s21
     $                       1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),2d24s21
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' prehcddjk. 12')
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
                      call enough('prehcddjk. 21',bc,ibc)
                      call dgemm('n','n',nli(l),nlj(l),ncsf(karg),1d0,   2d23s21
     $                    bc(jtmpi),nli(l),bc(jtmpj),ncsf(karg),0d0,    2d23s21
     $                    bc(itmpm),nli(l),                             2d23s21
     d' prehcddjk. 13')
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
        do l=1,4
         if(min(nfdat(2,l,isb),nfdat(2,l,jsb)).gt.0)then
          if(jsb.eq.isb)then
          end if
         end if
         do ll=1,4
          if(min(nfdat(2,l,isb),nfdat(2,ll,jsb)).gt.0)then
           nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                           2d24s21
           do lsa=1,nsymb
            lsb=multh(lsa,ijsb)
            if(lsb.ge.lsa.and.l.eq.ll)then                              11d7s22
             if(lsb.eq.lsa)then                                           11d7s22
              nn=(irefo(lsa)*(irefo(lsa)+1))/2                            11d7s22
             else                                                         11d7s22
              nn=irefo(lsa)*irefo(lsb)                                    11d7s22
             end if                                                       11d7s22
             if(nn.gt.0)then
              itmp=ibcoff
              ibcoff=itmp+nll*nn
              call enough('prehcddjk. 22',bc,ibc)
              do i=0,nn-1
               do j=0,nll-1
                ji=idenjn(l,lsa)+j+nll*i
                ij=itmp+i+nn*j
                bc(ij)=bc(ji)
               end do
              end do
              npos=nn*nll
              icmp1=ibcoff
              icmp2=icmp1+nn
              icmp3=icmp2+npos
              ibcoff=icmp3+npos
              call enough('prehcddjk. 23',bc,ibc)
              call cmpvcsf(nn,nll,ibc(icmp1),ibc(icmp2),bc(icmp3),
     $             bc(itmp),nkeep)
              nkeeph=nkeep/2                                            11d7s22
              if(mod(nkeep,2).ne.0)nkeeph=nkeeph+1                      11d7s22
              ntot=nn+nkeeph+nkeep+1                                    11d7s22
              if(ntot.lt.npos)then
               njc1=njc1+1
               idenj(l,lsa,isb,jsb)=-ibctopr                            11d7s22
               ibc(ibctopr)=nkeep                                       11d7s22
               ibctopr=ibctopr+1                                        11d7s22
               do i=0,nn-1                                              11d7s22
                ibc(ibctopr+i)=ibc(icmp1+i)                             11d7s22
                ipack8=ibc(icmp1+i)
               end do                                                   11d7s22
               ibctopr=ibctopr+nn                                       11d7s22
               do i=0,nkeeph-1                                          11d7s22
                ibc(ibctopr+i)=ibc(icmp2+i)                             11d7s22
                ipack8=ibc(icmp2+i)
               end do                                                   11d7s22
               ibctopr=ibctopr+nkeeph                                   11d7s22
               do i=0,nkeep-1                                           11d7s22
                bc(ibctopr+i)=bc(icmp3+i)                               11d7s22
               end do                                                   11d7s22
               ibctopr=ibctopr+nkeep                                    11d7s22
              else
               njc0=njc0+1
               idenj(l,lsa,isb,jsb)=ibctopr                             11d7s22
               do i=0,npos-1
                bc(ibctopr+i)=bc(idenjn(l,lsa)+i)                       11d7s22
               end do                                                   11d7s22
               ibctopr=ibctopr+npos                                     11d7s22
              end if
              ibcoff=itmp                                               12d10s22
             end if
            end if
            nn=irefo(lsa)*irefo(lsb)                                    11d7s22
            if(nn.gt.0)then
             if(isb.eq.jsb.and.lsb.eq.lsa)then                          11d7s22
              if(l.le.ll)then
              end if
             end if
             itmp=ibcoff
             ibcoff=itmp+nll*nn
             call enough('prehcddjk. 24',bc,ibc)
             do i=0,nn-1
              do j=0,nll-1
               ji=idenkn(l,ll,lsa)+j+nll*i
               ij=itmp+i+nn*j
               bc(ij)=bc(ji)
              end do
             end do
             npos=nn*nll
             icmp1=ibcoff
             icmp2=icmp1+nn
             icmp3=icmp2+npos
             ibcoff=icmp3+npos
             call enough('prehcddjk. 25',bc,ibc)
             call cmpvcsf(nn,nll,ibc(icmp1),ibc(icmp2),bc(icmp3),
     $             bc(itmp),nkeep)
             nkeeph=nkeep/2                                             11d7s22
             if(mod(nkeep,2).ne.0)nkeeph=nkeeph+1                       11d7s22
             ntot=nn+nkeeph+nkeep+1                                     11d7s22
             if(ntot.lt.npos)then
              nkc1=nkc1+1
              idenk(l,ll,lsa,isb,jsb)=-ibctopr                          11d7s22
              ibc(ibctopr)=nkeep                                        11d7s22
              ibctopr=ibctopr+1                                         11d7s22
              do i=0,nn-1                                               11d7s22
               ibc(ibctopr+i)=ibc(icmp1+i)                              11d7s22
              end do                                                    11d7s22
              ibctopr=ibctopr+nn                                        11d7s22
              do i=0,nkeeph-1                                           11d7s22
               ibc(ibctopr+i)=ibc(icmp2+i)                              11d7s22
              end do                                                    11d7s22
              ibctopr=ibctopr+nkeeph                                    11d7s22
              do i=0,nkeep-1                                            11d7s22
               bc(ibctopr+i)=bc(icmp3+i)                                11d7s22
              end do                                                    11d7s22
              ibctopr=ibctopr+nkeep                                     11d7s22
             else
              nkc0=nkc0+1
              idenk(l,ll,lsa,isb,jsb)=ibctopr                           11d7s22
              do i=0,npos-1
               bc(ibctopr+i)=bc(idenkn(l,ll,lsa)+i)                     11d7s22
              end do                                                    11d7s22
              ibctopr=ibctopr+npos                                      11d7s22
             end if
            end if
           end do
          end if
         end do
        end do
        ibcoff=ibctopr                                                  11d7s22
        nwdstot=ibcoff-ibc0                                             5d17s23
        idenput(isb,jsb,2)=nwdstot                                      5d17s23
        call ilimts(1,nwdstot,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e) 5d17s23
        itc=1                                                           5d17s23
        call ddi_create(bc,ibc,itc,idenput(isb,jsb,2),                  5d17s23
     $       idenput(isb,jsb,1))                                        5d17s23
        jtc=il                                                          5d17s23
        jto=ih                                                          5d17s23
        ist=ibc0+il-1                                                   5d17s23
        call ddi_put(bc,ibc,idenput(isb,jsb,1),itc,itc,jtc,jto,bc(ist)) 5d17s23
        ibcoff=ibc0                                                     5d17s23
       end do                                                           1d8s21
      end do                                                            1d8s21
      idenorg=ibcoff                                                    5d17s23
      return                                                            1d8s21
      end                                                               1d8s21
c mpec2.1 version zeta copyright u.s. government
