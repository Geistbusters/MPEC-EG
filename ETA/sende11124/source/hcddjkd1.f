c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcddjkd1(nff22,nfdat,vd,nsymb,mdon,mdoo,nec,multh,     3d9s21
     $     isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,ixw1,ixw2,          3d9s21
     $     norb,nrootu,iden,nh0av,sr2,srh,iff22,ff22,bc,ibc,iroff)      1      11d10s22
      implicit real*8 (a-h,o-z)                                         1d8s21
c                                                                       1d8s21
c     1 part density
c                                                                       1d8s21
      logical lprt                                                      1d12s21
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d22s21
      integer*8 ipack8,itc,ito,jtc,jto,itestc,itesto,itesta,itestb,     1d17s21
     $     jtesta,jtestb,iff22,gandcc,gandco,gandcb                     2d6s23
      external second                                                   2d18s21
      equivalence (ipack8,ipack4)                                       1d8s21
      dimension nff22(mdoo+1,2,nsymb),nfdat(5,4,*),vd(*),               3d9s21
     $     multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),irel(*),ism(*),       1d8s21
     $     irefo(*),iden(*),nh0av(*),                                   3d9s21
     $     ipack4(2),nprei(4),mprei(4),nli(4),nprej(4),                 1d12s21
     $     mprej(4),nlj(4),nab4(2,3),idenj(4,8),ndenj(4,8),             1d9s21
     $     mdenj(4,8),idenk(4,8),ndenk(4,8),mdenk(4,8),                 1d9s21
     $     idenhvv(4),mdenhvv(4),itest(32,3),iden1e(4),mden1e(4),       1d8s21
     $     ldenj(4,8),ldenk(4,8),iden1ef(4),idenhvvf(4),                1d11s21
     $     idenjf(4,8),idenkf(4,4,8),ndenjf(4,8),ndenkf(4,4,8),
     $     nokj(4,8),nokk(4,4,8),isorb(32),idorb(32),idenjn(4,8),       2d23s21
     $     idenkn(4,4,8),iden1en(4,8),idenhvvn(4),                      3d9s21
     $     idenhvvd(4),iff22(*),ff22(*)                                 3d21s22
      include "common.store"                                            1d8s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      todate=0d0
      ion1h=1
      ionj=1
      ionk=1
      loop=0                                                            1d8s21
      loopx=12000000
      ibcoffo=ibcoff
      igoal=iden(1)+2*(1+nh0av(1))
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
        ibc0=ibcoff                                                     1d8s21
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    2d24s21
          nll=nfdat(2,l,isb)*nfdat(2,l,jsb)                             3d9s21
          do isc=1,nsymb                                                3d9s21
           nvv=(irefo(isc)*(irefo(isc)+1))/2                            3d9s21
           iden1en(l,isc)=ibcoff                                        3d9s21
           ibcoff=iden1en(l,isc)+nll*nvv                                3d9s21
          end do                                                        1d8s21
          idenhvvn(l)=ibcoff                                            3d9s21
          ibcoff=idenhvvn(l)+nll                                        3d9s21
          idenhvvd(l)=ibcoff                                            3d9s21
          ibcoff=idenhvvd(l)+nll                                        3d9s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call enough('hcddjkd1.  1',bc,ibc)
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
           ipack8=iff22(ivcv)                                           3d21s22
           itc=ipack4(1)                                                1d8s21
           ito=ipack4(2)                                                1d8s21
           ito=ibset(ito,norbx)                                         1d8s21
           ito=ibset(ito,norbxx)                                        1d8s21
           nspacei=iff22(ivcv+1)                                        3d21s22
           do l=1,4                                                       11d20s20
            nli(l)=iff22(ivcv+1+l)                                      3d21s22
           end do
           do i=1,norb                                                  1d8s21
            itest(i,1)=0                                                1d8s21
           end do                                                       1d8s21
           do i=1,norb                                                  1d8s21
            if(btest(itc,i))itest(i,1)=2                                1d8s21
            if(btest(ito,i))itest(i,1)=1                                1d8s21
           end do                                                       1d8s21
           itest(norbx,1)=1                                             1d8s21
           itest(norbxx,1)=1                                            1d8s21
           do ncloj=max(ncloi-1,mdon),min(ncloi+1,mdoo)                 3d9s21
            nclojp=ncloj+1                                                1d8s21
            if(nff22(nclojp,1,jsb).gt.0)then                            1d14s21
             if(mod(idoit,mynprocg).eq.mynowprog)then                   3d1s21
             jarg=nclojp-mdon                                             1d8s21
             nopenj=nec-2*nclojp                                          1d8s21
             nopenjp=nopenj+2                                           1d8s21
             jvcv=nfdat(5,1,jsb)+nff22(nclojp,2,jsb)                      1d8s21
             do jf=1,nff22(nclojp,1,jsb)                                  1d8s21
              ipack8=iff22(jvcv)                                        3d21s22
              jtc=ipack4(1)                                               1d8s21
              jto=ipack4(2)                                               1d8s21
              jto=ibset(jto,norbx)                                        1d8s21
              jto=ibset(jto,norbxx)                                       1d8s21
              nspacej=iff22(jvcv+1)                                     3d21s22
              do l=1,4                                                       11d20s20
               nlj(l)=iff22(jvcv+1+l)                                   3d21s22
              end do
              if(isb.eq.jsb)then                                        1d17s21
c                                                                       1d8s21
               gandcc=ieor(itc,jtc)                                     2d6s23
               gandco=ieor(ito,jto)                                     2d6s23
               gandcb=ior(gandcc,gandco)                                10d21s22
               ndifb=popcnt(gandcb)                                     10d21s22
               if(ndifb.le.2)then                                       2d6s23
                ndifs=popcnt(gandco)                                     10d14s22
                ndifd=popcnt(gandcc)                                     10d14s22
                 if(ndifs.eq.0.and.ndifd.eq.0)then                      2d6s23
                  iargo=1                                               2d24s21
                  do l=1,4                                              2d24s21
                   if(nli(l).gt.0)then                                  2d24s21
                    itmpt=ibcoff                                        2d24s21
                    itmpm=itmpt+ncsf2(l,iarg)*nli(l)                    2d24s21
                    ibcoff=itmpm+nli(l)*nli(l)                          2d24s21
                    call enough('hcddjkd1.  2',bc,ibc)
                    iad1=ivcv+iff22(ivcv+5+l)                           3d21s22
                    iad2=iad1+nli(l)                                     3d19s21
                    do iii=0,nli(l)-1                                   2d24s21
                     do j=0,ncsf2(l,iarg)-1                             2d24s21
                      ji=iad2+j+ncsf2(l,iarg)*iii                       2d24s21
                      ij=itmpt+iii+nli(l)*j                             2d24s21
                      bc(ij)=ff22(ji)                                   3d21s22
                     end do                                             2d24s21
                    end do                                              2d24s21
                    call dgemm('n','n',nli(l),nli(l),ncsf2(l,iarg),1d0, 3d9s21
     $                   bc(itmpt),nli(l),ff22(iad2),ncsf2(l,iarg),0d0, 3d21s22
     $                   bc(itmpm),nli(l),                              2d24s21
     d' hcddjkd1.  1')
                    nll=nfdat(2,l,isb)*nfdat(2,l,jsb)                   3d9s21
                    do i=1,norb                                         3d9s21
                     fact=0d0                                           3d9s21
                     if(btest(ipack4(1),i))then                         3d9s21
                      fact=2d0                                          3d9s21
                     else if(btest(ipack4(2),i))then                    3d9s21
                      fact=1d0                                          3d9s21
                     end if                                             3d9s21
                     if(fact.ne.0d0)then                                3d9s21
                      is=ism(i)                                         3d9s21
                      jg=irel(i)-1                                       3d9s21
                      jg=((jg*(jg+1))/2)+jg                             3d9s21
                      jtmpm=itmpm                                         2d24s21
                      do j=0,nli(l)-1                                     2d24s21
                       jj=iff22(iad1+j)-1                               3d21s22
                       jden=iden1en(l,is)+nfdat(2,l,isb)*jj-1+nll*jg    3d9s21
                       do k=0,nli(l)-1                                    2d24s21
                        kk=iff22(iad1+k)                                3d21s22
                        bc(jden+kk)=bc(jden+kk)+bc(jtmpm+k)*fact        3d9s21
                       end do                                             2d24s21
                       jtmpm=jtmpm+nli(l)                                 2d24s21
                      end do                                              2d24s21
                     end if                                             3d9s21
                    end do                                              3d9s21
                    jtmpm=itmpm                                         2d24s21
                    do j=0,nli(l)-1                                     2d24s21
                     jj=iff22(iad1+j)-1                                 3d21s22
                     jden=idenhvvd(l)+nfdat(2,l,isb)*jj-1               3d9s21
                     do k=0,nli(l)-1                                    2d24s21
                      kk=iff22(iad1+k)                                  3d21s22
                      bc(jden+kk)=bc(jden+kk)+bc(jtmpm+k)               3d9s21
                     end do                                             2d24s21
                     jtmpm=jtmpm+nli(l)                                 2d24s21
                    end do                                              2d24s21
                    ibcoff=itmpt                                        2d24s21
                   end if                                               2d24s21
                   iargo=iargo+ncsf2(l,iarg)                            2d24s21
                  end do                                                2d24s21
                 else if(ndifs.eq.2.and.ndifb.eq.2)then                 2d6s23
                  do i=1,norbxx                                         2d6s23
                   if(btest(gandco,i))then                              2d6s23
                    if((btest(ito,i).and..not.btest(jtc,i)).or.         2d6s23
     $                 (btest(itc,i).and.btest(jto,i)))then             2d6s23
                     nab4(1,1)=i                                        2d6s23
                    else                                                2d6s23
                     nab4(2,1)=i                                        2d6s23
                    end if                                              2d6s23
                   end if                                               2d6s23
                  end do                                                2d6s23
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
                  if(ksa.ne.ksb)then                                     1d9s21
                   write(6,*)('need ksa=ksb, but have '),ksa,ksb
                   stop 'hcddjk:d'
                  end if                                                 1d9s21
                  kgx=max(kga,kgb)                                      3d9s21
                  kgn=min(kga,kgb)                                      3d9s21
                  kgg=((kgx*(kgx+1))/2)+kgn                             3d9s21
                  jargo=1                                                2d23s21
                  iargo=1                                                2d24s21
                  do l=1,4                                                1d8s21
                   if(min(nlj(l),nli(l)).gt.0)then                       2d25s21
                    jad1=jvcv+iff22(jvcv+5+l)                           3d21s22
                    jad2=jad1+nlj(l)                                     3d19s21
                    iad1=ivcv+iff22(ivcv+5+l)                           3d21s22
                    iad2=iad1+nli(l)                                     3d19s21
                    itmp1=ibcoff                                         2d23s21
                    itmp1b=itmp1+ncsfmid1*nlj(l)                         2d23s21
                    ibcoff=itmp1b+ncsfmid1*nli(l)                        2d23s21
                    call enough('hcddjkd1.  3',bc,ibc)
                    if(iwpk1.lt.0)then                                   2d23s21
                     nusedi=ibc(-iwpk1)/2                                2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1          2d23s21
                     icmp1=-iwpk1+1                                      2d23s21
                     icmp2=icmp1+ncsf(jarg)                              2d23s21
                     icmp3=icmp2+nusedi                                  2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   ff22(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1), 3d21s22
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d14s22
                    else                                                              12d5s20
                     jwpk1=iwpk1+ncsfmid1*(jargo-1)                      2d23s21
                     call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                  1d0,bc(jwpk1),ncsfmid1,ff22(jad2),ncsf2(l,jarg),3d21s22
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjkd1.  2')
                    end if                                                            12d5s20
                    if(iwpk1b.lt.0)then                                   2d23s21
                     nusedi=ibc(-iwpk1b)/2                                2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                     icmp1=-iwpk1b+1                                      2d23s21
                     icmp2=icmp1+ncsf(iarg)                              2d23s21
                     icmp3=icmp2+nusedi                                  2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   ff22(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),3d21s22
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          11d14s22
                    else                                                              12d5s20
                     jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                     call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                 1d0,bc(jwpk1b),ncsfmid1,ff22(iad2),ncsf2(l,iarg),3d21s22
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjkd1.  3')
                    end if                                                            12d5s20
                    itmpt=ibcoff                                         2d23s21
                    itmpm=itmpt+nli(l)*ncsfmid1                         2d24s21
                    ibcoff=itmpm+nli(l)*nlj(l)                          2d24s21
                    call enough('hcddjkd1.  4',bc,ibc)
                    do iii=0,nli(l)-1                                    2d23s21
                     do j=0,ncsfmid1-1                                   2d23s21
                      ji=itmp1b+j+ncsfmid1*iii                           2d23s21
                      ij=itmpt+iii+nli(l)*j                              2d23s21
                      bc(ij)=bc(ji)                                      2d23s21
                     end do                                              2d23s21
                    end do                                               2d23s21
c                                                                       3d10s21
c     half here because we only store triangle of density and           3d10s21
c     we will finally have square density by reflecting lower triangle  3d10s21
c     to upper triangle.                                                3d10s21
c                                                                       3d10s21
                    if(kgx.eq.kgn)then                                  3d10s21
                     ffr=1d2                                            3d10s21
                    else                                                3d10s21
                     ffr=0.5d0                                          3d10s21
                    end if                                              3d10s21
                    call dgemm('n','n',nli(l),nlj(l),ncsfmid1,ffr,      3d10s21
     $                   bc(itmpt),nli(l),bc(itmp1),ncsfmid1,0d0,        2d23s21
     $                   bc(itmpm),nli(l),                               2d23s21
     d' hcddjkd1.  4')
                    nll=nfdat(2,l,isb)*nfdat(2,l,jsb)                   3d9s21
                    jtmpm=itmpm                                         2d24s21
                    do iii=0,nlj(l)-1                                   2d24s21
                     jj=iff22(jad1+iii)-1                               3d21s22
                     jjden=iden1en(l,ksa)+nfdat(2,l,isb)*jj-1+nll*kgg   3d9s21
                     do j=0,nli(l)-1                                    2d24s21
                      kk=iff22(iad1+j)                                  3d21s22
                      bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)             2d24s21
                     end do                                             2d24s21
                     jtmpm=jtmpm+nli(l)                                 2d24s21
                    end do                                              2d24s21
                    ibcoff=itmp1                                        2d24s21
                   end if                                                 1d8s21
                   jargo=jargo+ncsf2(l,jarg)                             2d23s21
                   iargo=iargo+ncsf2(l,iarg)                             2d24s21
                  end do                                                  1d8s21
                 end if                                                  1d8s21
                end if                                                   1d8s21
              end if                                                    1d8s21
              jto=ibclr(jto,norbx)                                      1d8s21
              jto=ibset(jto,norbxxx)                                    1d8s21
              gandcc=ieor(itc,jtc)                                      2d6s23
              gandco=ieor(ito,jto)                                      2d6s23
              gandcb=ior(gandcc,gandco)                                 2d6s23
              ndifb=popcnt(gandcb)                                      2d6s23
c                                                                       1d8s21
              if(ndifb.le.2)then                                        2d6s23
               ndifs=popcnt(gandco)                                     2d6s23
               if(ndifs.eq.2.and.ndifb.eq.2)then                        2d6s23
                do i=1,norbxxx                                          2d6s23
                 if(btest(gandco,i))then                                2d6s23
                  if((btest(ito,i).and..not.btest(jtc,i)).or.           2d6s23
     $                (btest(itc,i).and.btest(jto,i)))then               2d6s23
                   nab4(1,1)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,1)=i                                          2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
                 if(isb.ne.jsb)then                                      1d9s21
                  write(6,*)('hey, what''s up with isb ne jsb ? ')       1d9s21
                  stop 'hcddjk:g'                                                   1d9s21
                 end if                                                  1d9s21
                 call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,   1d8s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $              ncsfmid1,bc,ibc)                                    11d14s22
                 call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,   2d23s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,      1d8s21
     $              ncsfmid1b,bc,ibc)                                   11d14s22
                 nlij=0                                                 2d23s21
                 do l=1,4                                               2d23s21
                  nlij=nlij+nli(l)*nlj(l)                               2d23s21
                 end do                                                 2d23s21
                 itmpn=ibcoff                                           2d23s21
                 ibcoff=itmpn+nlij                                      2d23s21
                 call enough('hcddjkd1.  5',bc,ibc)
                 jtmpn=itmpn                                            2d23s21
                 jargo=1                                                 2d22s21
                 iargo=1                                                2d23s21
                 do l=1,4                                                1d8s21
                  if(min(nlj(l),nli(l)).gt.0)then                        2d23s21
                   iad1=ivcv+iff22(ivcv+5+l)                            3d21s22
                   iad2=iad1+nli(l)                                     3d19s21
                   jad1=jvcv+iff22(jvcv+5+l)                            3d21s22
                   jad2=jad1+nlj(l)                                     3d19s21
                   itmp1=ibcoff                                         2d23s21
                   itmp1b=itmp1+ncsfmid1*nlj(l)                         2d23s21
                   ibcoff=itmp1b+ncsfmid1*nli(l)                        2d23s21
                   call enough('hcddjkd1.  6',bc,ibc)
                   if(iwpk1.lt.0)then                                   2d23s21
                    nusedi=ibc(-iwpk1)/2                                2d23s21
                    if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1          2d23s21
                    icmp1=-iwpk1+1                                      2d23s21
                    icmp2=icmp1+ncsf(jarg)                              2d23s21
                    icmp3=icmp2+nusedi                                  2d23s21
                    call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   ff22(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1), 3d21s22
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d14s22
                   else                                                              12d5s20
                    jwpk1=iwpk1+ncsfmid1*(jargo-1)                      2d23s21
                    call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                  1d0,bc(jwpk1),ncsfmid1,ff22(jad2),ncsf2(l,jarg),3d21s22
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjkd1.  5')
                   end if                                                            12d5s20
                   if(iwpk1b.lt.0)then                                   2d23s21
                    nusedi=ibc(-iwpk1b)/2                                2d23s21
                    if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                    icmp1=-iwpk1b+1                                      2d23s21
                    icmp2=icmp1+ncsf(iarg)                              2d23s21
                    icmp3=icmp2+nusedi                                  2d23s21
                    call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   ff22(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),3d21s22
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          11d14s22
                   else                                                              12d5s20
                    jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                    call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                 1d0,bc(jwpk1b),ncsfmid1,ff22(iad2),ncsf2(l,iarg),3d21s22
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjkd1.  6')
                   end if                                                            12d5s20
                   itmpt=ibcoff                                         2d23s21
                   ibcoff=itmpt+nli(l)*ncsfmid1                         2d23s21
                   call enough('hcddjkd1.  7',bc,ibc)
                   do iii=0,nli(l)-1                                    2d23s21
                    do j=0,ncsfmid1-1                                   2d23s21
                     ji=itmp1b+j+ncsfmid1*iii                           2d23s21
                     ij=itmpt+iii+nli(l)*j                              2d23s21
                     bc(ij)=bc(ji)                                      2d23s21
                    end do                                              2d23s21
                   end do                                               2d23s21
                   call dgemm('n','n',nli(l),nlj(l),ncsfmid1,1d0,       2d23s21
     $                  bc(itmpt),nli(l),bc(itmp1),ncsfmid1,0d0,        2d23s21
     $                  bc(jtmpn),nli(l),                               2d23s21
     d' hcddjkd1.  7')
                   ibcoff=itmp1                                         2d23s21
                   ktmpn=jtmpn                                          2d24s21
                   do iii=0,nlj(l)-1                                    2d24s21
                    jj=iff22(jad1+iii)-1                                3d21s22
                    jjden=idenhvvn(l)+nfdat(2,l,isb)*jj-1                2d24s21
                    do j=0,nli(l)-1                                     2d24s21
                     kk=iff22(iad1+j)                                   3d21s22
                     bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)              2d24s21
                    end do                                              2d24s21
                    ktmpn=ktmpn+nli(l)                                  2d24s21
                   end do                                               2d24s21
                   jtmpn=jtmpn+nli(l)*nlj(l)                            2d23s21
                  end if                                                2d23s21
                  jargo=jargo+ncsf2(l,jarg)                              2d22s21
                  iargo=iargo+ncsf2(l,iarg)                             2d24s21
                 end do                                                  1d8s21
               end if                                                   1d8s21
              end if                                                    1d8s21
              jvcv=jvcv+nspacej                                         3d19s21
             end do                                                       1d8s21
             end if                                                     3d1s21
             idoit=idoit+1                                               3d1s21
            end if                                                      1d8s21
           end do                                                       1d8s21
           ivcv=ivcv+nspacei                                            3d19s21
          end do                                                        1d8s21
          ibcoff=ibcst                                                  1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call dws_gsumf(bc(ibc0),nwds)                                   2d23s21
        if(ijsb.eq.1.and.ion1h.ne.0)then                                               1d11s21
         ioff=ioffdnon                                                  1d11s21
         do isbv1=1,nsymb                                               1d11s21
          isbv2=multh(isbv1,isbv12)                                     1d11s21
          if(isbv2.ge.isbv1)then                                        1d11s21
           if(isbv1.eq.isbv2)then                                       1d11s21
            ioffs=ioff                                                  1d12s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       1d30s23
            isw=0                                                       1d30s23
            nnn=nvirt(isbv1)*nrootu                                     1d11s21
            call ilimts(nvirt(isbv1),nrootu,mynprocg,mynowprog,il,ih,   1d11s21
     $           i1s,i1e,i2s,i2e)                                       1d11s21
            nhere=ih+1-il                                               1d11s21
            if(min(nfdat(2,1,isb),nhere).gt.0)then                      11d17s21
             do isc=1,nsymb                                             3d9s21
              nrr=(irefo(isc)*(irefo(isc)+1))/2                         3d9s21
              if(nrr.gt.0)then                                          3d9s21
               mcol=nfdat(2,1,isb)*nrr                                  3d9s21
               itmp=ibcoff                                                1d11s21
               ibcoff=itmp+nhere*mcol                                   3d9s21
               call enough('hcddjkd1.  8',bc,ibc)
               call dgemm('n','n',nhere,mcol,nfdat(2,1,isb),1d0,        3d9s21
     $          vd(ioff+il-1),nnn,bc(iden1en(1,isc)),nfdat(2,1,isb),0d0,3d9s21
     $              bc(itmp),nhere,                                     3d9s21
     d' hcddjkd1.  8')
               do jka=0,irefo(isc)-1                                    3d9s21
                do jkb=0,jka                                            3d9s21
                 jkab=((jka*(jka+1))/2)+jkb                             3d9s21
                 do i=0,nfdat(2,1,isb)-1                                    1d11s21
                  i10=i1s                                                   1d11s21
                  i1n=nvirt(isbv1)                                          1d11s21
                  jtmp=itmp+nhere*(i+nfdat(2,1,isb)*jkab)               3d9s21
                  do i2=i2s,i2e                                             1d11s21
                   if(i2.eq.i2e)i1n=i1e                                     1d11s21
                   i2m=i2-1                                             4d5s21
                   iad=ioff-1+nvirt(isbv1)*(i2m+nrootu*i)               4d5s21
                   sum=0d0                                              3d22s22
                   do i1=i10,i1n                                            1d11s21
                    orig=sum
                    sum=sum+bc(jtmp)*vd(iad+i1)                         3d9s21
                    jtmp=jtmp+1                                             1d11s21
                   end do                                                   1d11s21
                   i10=1                                                    1d11s21
                   iad=iden(isc)+jka+nh0av(isc)*(jkb+nh0av(isc)*i2m)    3d22s22
                   bc(iad)=bc(iad)+sum                                  3d22s22
                  end do                                                    1d11s21
                 end do                                                     1d11s21
                end do                                                  3d9s21
               end do                                                   3d9s21
               ibcoff=itmp                                                1d11s21
              end if                                                    3d9s21
             end do                                                     3d9s21
            end if                                                      1d11s21
            if(min(nnn,nfdat(2,1,isb)).gt.0)then                        11d17s21
             ioff=ioff+nnn*nfdat(2,1,isb)                                1d11s21
             itmp=ibcoff                                                 1d18s21
             ibcoff=itmp+nnn*nfdat(2,1,isb)                              1d18s21
             call enough('hcddjkd1.  9',bc,ibc)
             call dgemm('n','n',nnn,nfdat(2,1,isb),nfdat(2,1,isb),       1d18s21
     $              2d0,vd(ioffs),nnn,bc(idenhvvn(1)),nfdat(2,1,isb),   1d18s21
     $              0d0,bc(itmp),nnn,                                   1d18s21
     d' hcddjkd1.  9')
             do i=0,nfdat(2,1,isb)-1                                     1d18s21
              do ir=0,nrootm                                             1d18s21
               iru=ir*iroff                                             1d27s23
               do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   1d18s21
                iad=ioffs+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                ih0=iden(isbv1)+(iv+irefo(isbv1))*(nh0av(isbv1)+1)      3d9s21
     $               +nh0av(isbv1)*nh0av(isbv1)*iru                     1d27s23
                jtmp=itmp+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                bc(ih0)=bc(ih0)+vd(iad)*bc(jtmp)                        3d9s21
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
              do isc=1,nsymb                                            3d9s21
               nrr=(irefo(isc)*(irefo(isc)+1))/2                        3d9s21
               if(nrr.gt.0)then                                         3d9s21
                mcol=nfdat(2,l,isb)*nrr                                 3d9s21
                itmp=ibcoff                                                1d11s21
                ibcoff=itmp+nhere*mcol                                  3d9s21
                call enough('hcddjkd1. 10',bc,ibc)
                call dgemm('n','n',nhere,mcol,nfdat(2,l,isb),1d0,       3d9s21
     $          vd(ioff+il-1),nnn,bc(iden1en(l,isc)),nfdat(2,l,isb),0d0,3d10s21
     $               bc(itmp),nhere,                                    3d9s21
     d' hcddjkd1. 10')
                do kga=0,irefo(isc)-1                                   3d9s21
                 do kgb=0,kga                                           3d9s21
                  kgab=((kga*(kga+1))/2)+kgb                            3d9s21
                  sum=0d0                                               3d9s21
                  do i=0,nfdat(2,l,isb)-1                                    1d11s21
                   i10=i1s                                                   1d11s21
                   i1n=nvv                                                  1d11s21
                   jtmp=itmp+nhere*(i+nfdat(2,l,isb)*kgab)              3d9s21
                   do i2=i2s,i2e                                             1d11s21
                    i2m=i2-1                                            4d5s21
                    i2mu=i2m*iroff                                      1d27s23
                    sum=0d0                                             3d22s22
                    if(i2.eq.i2e)i1n=i1e                                     1d11s21
                    iad=ioff-1+nvv*(i2m+nrootu*i)                       4d5s21
                    do i1=i10,i1n                                            1d11s21
                     sum=sum+vd(iad+i1)*bc(jtmp)                        3d9s21
                     jtmp=jtmp+1                                             1d11s21
                    end do                                                   1d11s21
                    i10=1                                                    1d11s21
                    iad=iden(isc)+kga+nh0av(isc)*(kgb+nh0av(isc)*i2mu)  1d27s23
                    bc(iad)=bc(iad)+sum                                 3d22s22
                   end do                                                    1d11s21
                  end do                                                     1d11s21
                 end do                                                 3d9s21
                end do                                                  3d9s21
                ibcoff=itmp                                             3d9s21
               end if                                                   3d9s21
              end do                                                    3d9s21
c
c     Gvv'ri=dij Hvv" Vv'v"rj ... if sym v = sym v' = sym v" we get this
c     if v"=v' > Gvv'ri=dij Hvv' Vv'v'rj.
c     if v=v' > Gvvri dij Hvv" Vvv"rj
              if(isbv12.eq.1.and.l.eq.1.and.nvirt(isbv1).gt.0)then      4d28s21
               itmp=ibcoff                                              3d9s21
               nrow=nvirt(isbv1)*nrootu                                 1d12s21
               ibcoff=itmp+nrow*nfdat(2,1,isb)                          3d9s21
               if(nfdat(2,1,isb).gt.0)then                              1d30s23
                call dgemm('n','n',nrow,nfdat(2,1,isb),nfdat(2,1,isb),
     $              sr2,vd(ioffs),nrow,bc(idenhvvn(1)),nfdat(2,1,isb),  3d9s21
     $              0d0,bc(itmp),nrow,                                  1d12s21
     d' hcddjkd1. 11')
               end if                                                   1d30s23
               do i=0,nfdat(2,l,isb)-1                                   1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                   1d11s21
                do i2=i2s,i2e                                             1d11s21
                 ir=i2-1                                                 1d11s21
                 iru=ir*iroff                                           1d27s23
                 if(i2.eq.i2e)i1n=i1e                                    1d11s21
                 do i1=i10,i1n                                           1d11s21
                  i1m=i1-1                                               1d11s21
                  iv2=i1m/nvirt(isbv1)                                   1d11s21
                  iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                  try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                  iv2t=sqrt(try)+0.5d0                                   1d11s21
                  iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                  iv1=iv1t+isw*(iv1-iv1t)
                  iv2=iv2t+isw*(iv2-iv2t)
c     Gvv'ri=dij Hvv' Vv'v'rj
                  iad=ioff+i1m+nvv*(ir+nrootu*i)                        1d12s21
                  jtmp=itmp+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                  ktmp=itmp+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                  ih0=iden(isbv1)+iv1+irefo(isbv1)+nh0av(isbv1)*(       3d9s21
     $                irefo(isbv1)+iv2+nh0av(isbv1)*iru)                1d27s23
                  avg=0.5d0*vd(iad)*(bc(ktmp)+bc(jtmp))
                  bc(ih0)=bc(ih0)+avg                                   3d10s21
                  ih0=iden(isbv1)+iv2+irefo(isbv1)+nh0av(isbv1)*(       3d10s21
     $                irefo(isbv1)+iv1+nh0av(isbv1)*iru)                1d27s23
                  bc(ih0)=bc(ih0)+avg                                   3d10s21
                 end do                                                 1d12s21
                 i10=1                                                  4d21s21
                end do                                                  1d12s21
               end do                                                   1d12s21
               ibcoff=itmp                                              3d9s21
              end if                                                    1d12s21
c     Gvv'ri=dij Hvv" Vv"v'rj, Gvv'ri=dij Hv'v" Vvv"rj
c
              itmp=ibcoff                                               3d9s21
              ibcoff=itmp+nhere*nfdat(2,l,isb)                          3d9s21
              call enough('hcddjkd1. 11',bc,ibc)
              call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(2,l,isb),   1d11s21
     $             tf,vd(ioff+il-1),nnn,bc(idenhvvn(l)),nfdat(2,l,isb), 1d11s21
     $             0d0,bc(itmp),nhere,                                  1d11s21
     d' hcddjkd1. 12')
              if(l.eq.1.and.isbv12.eq.1)then                            1d18s21
               do i=0,nfdat(2,l,isb)-1                                   1d11s21
                jtmp=itmp+nhere*i                                        1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                   1d11s21
                do i2=i2s,i2e                                             1d11s21
                 ir=i2-1                                                 1d11s21
                 iru=ir*iroff                                           1d27s23
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
                  ih0=iden(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(       3d9s21
     $                irefo(isbv1)+iv2+nh0av(isbv1)*iru)                1d27s23
                  iad=ioffs+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                  val=0.5d0*vd(iad)*bc(jtmp)*sr2                        3d10s21
                  bc(ih0)=bc(ih0)+val                                   3d10s21
                  ih0=iden(isbv1)+irefo(isbv1)+iv2+nh0av(isbv1)*(       3d10s21
     $                irefo(isbv1)+iv1+nh0av(isbv1)*iru)                1d27s23
                  bc(ih0)=bc(ih0)+val                                   3d10s21
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
                   ih0=iden(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(      3d9s21
     $                irefo(isbv1)+iv2+nh0av(isbv1)*iru)                1d27s23
                   iad=ioffs+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                   val=0.5d0*vd(iad)*bc(jtmp)*sr2                       3d10s21
                   bc(ih0)=bc(ih0)+val                                  3d10s21
                   ih0=iden(isbv1)+irefo(isbv1)+iv2+nh0av(isbv1)*(      3d10s21
     $                irefo(isbv1)+iv1+nh0av(isbv1)*iru)                1d27s23
                   bc(ih0)=bc(ih0)+val                                  3d10s21
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
                iru=ir*iroff                                            1d27s23
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
                 ih0=iden(isbv1)+irefo(isbv1)+nh0av(isbv1)*(            3d9s21
     $                irefo(isbv1)+iv1+nh0av(isbv1)*iru)                1d27s23
                 iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                 do iv=0,ntop                                           1d11s21
                  irec=iv+nvirt(isbv1)*iv2                              1d11s21
                  itri=((iv2*(iv2-1))/2)+iv                             1d18s21
                  irow=itri+isw*(irec-itri)                             1d11s21
                  bc(ih0+iv)=bc(ih0+iv)+vd(iad+irow)*bc(jtmp)           3d9s21
                 end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vvv"rj
                 nbot=iv1+1                                             1d11s21
                 nbot=nbot+isw*(0-nbot)                                 1d11s21
                 ih0=iden(isbv2)+irefo(isbv2)+nh0av(isbv2)*             3d9s21
     $                (irefo(isbv2)+iv2+nh0av(isbv2)*iru)               1d27s23
                 do iv=nbot,nvirt(isbv2)-1                              1d11s21
                  irec=iv1+nvirt(isbv1)*iv                              1d11s21
                  itri=((iv*(iv-1))/2)+iv1                              1d18s21
                  irow=itri+isw*(irec-itri)                             1d11s21
                  bc(ih0+iv)=bc(ih0+iv)+vd(iad+irow)*bc(jtmp)           3d9s21
                 end do                                                 1d11s21
                 if(isbv1.eq.isbv2)then
c     Gvv'ri=dij Hvv" Vv'v"rj
                  ntop=iv1-1                                             1d18s21
                  ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                    1d18s21
                  ih0=iden(isbv1)+irefo(isbv1)+nh0av(isbv1)*(           3d9s21
     $                irefo(isbv1)+iv2+nh0av(isbv1)*iru)                1d27s23
                  iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                  do iv=0,ntop                                           1d11s21
                   irec=iv+nvirt(isbv1)*iv1                              1d11s21
                   itri=((iv1*(iv1-1))/2)+iv                             1d18s21
                   irow=itri+isw*(irec-itri)                             1d11s21
                   bc(ih0+iv)=bc(ih0+iv)+vd(iad+irow)*bc(jtmp)*tf       3d9s21
                  end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vv"vrj                                           1d18s21
                  nbot=iv2+1                                             1d18s21
                  nbot=nbot+isw*(0-nbot)                                 1d18s21
                  ih0=iden(isbv2)+irefo(isbv2)+nh0av(isbv2)*            3d9s21
     $                (irefo(isbv2)+iv1+nh0av(isbv2)*iru)               1d27s23
                  do iv=nbot,nvirt(isbv2)-1                              1d18s21
                   irec=iv2+nvirt(isbv1)*iv                              1d18s21
                   itri=((iv*(iv-1))/2)+iv2                              1d18s21
                   irow=itri+isw*(irec-itri)                             1d18s21
                   bc(ih0+iv)=bc(ih0+iv)+vd(iad+irow)*bc(jtmp)*tf       3d9s21
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
            end if                                                      1d11s21
            tf=-1d0                                                     1d11s21
           end do                                                       1d11s21
          end if                                                        1d11s21
         end do                                                         1d11s21
        end if                                                          1d11s21
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
