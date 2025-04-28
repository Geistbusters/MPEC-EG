c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcddjkbk(nff22b,nfdatb,nff22k,nfdatk,gd,vd,nsymb,mdon, 8d26s21
     $     mdoop,nec,multh,isymbra,isymket,nvirt,ncsf,ncsf2,irel,ism,   8d26s21
     $     irefo,ixw1,ixw2,norb,nrootu,ixmtf,lxmt,phase1,n2e,ixmt,      8d26s21
     $     isymop,i2eop,phase2,shift,sr2,srh,idoubo,nbasdws,izero,      8d26s21
     $     iff22b,ff22b,iff22k,ff22k,bc,ibc)                            11d10s22
      implicit real*8 (a-h,o-z)                                         1d8s21
c                                                                       1d8s21
c     0 and 2 virt contribution to gd=hdd*vd                            1d8s21
c                                                                       1d8s21
      logical lprt,lchoice                                              3d19s21
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d22s21
      integer*8 ipack8,itc,ito,jtc,jto,itestc,itesto,itesta,itestb,     1d17s21
     $     jtesta,jtestb,iff22b(*),iff22k(*),gandcc,gandco,gandcb       2d6s23
      external second                                                   2d18s21
      equivalence (ipack8,ipack4)                                       1d8s21
      dimension nff22b(mdoop,2,nsymb),nfdatb(5,4,*),vd(*),gd(*),          8d9s21
     $     multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),irel(*),ism(*),       1d8s21
     $     irefo(*),ixmtf(*),ixmt(8,*),isymop(*),i2eop(2,3),idoubo(*),  8d9s21
     $     ipack4(2),nprei(4),mprei(4),nli(4),nprej(4),nbasdws(*),      8d9s21
     $     mprej(4),nlj(4),nab4(2,3),idenj(4,8),ndenj(4,8),             1d9s21
     $     mdenj(4,8),idenk(4,8),ndenk(4,8),mdenk(4,8),loope(6),                 1d9s21
     $     idenhvv(4),mdenhvv(4),itest(32,3),iden1e(4),mden1e(4),       1d8s21
     $     ldenj(4,8),ldenk(4,8),iden1ef(4),idenhvvf(4),                1d11s21
     $     idenjf(4,8),idenkf(4,4,8),ndenjf(4,8),ndenkf(4,4,8),
     $     nokj(4,8),nokk(4,4,8),isorb(32),idorb(32),idenjn(4,8),       2d23s21
     $     idenkn(4,4,8),iden1en(4),idenhvvn(4),ff22b(*),itmpl(4),      8d26s21
     $     nff22k(mdoop,2,nsymb),nfdatk(5,4,*),ff22k(*),ioxx(2)         2d6s23
      data loopx/10000000/
      include "common.store"                                            1d8s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data icall/0/
      save icall
      loop=0
      if(izero.ne.0)then                                                8d10s21
       ioffg=0                                                          8d10s21
       do isb=1,nsymb                                                   8d10s21
        isbv12=multh(isb,isymbra)                                       8d10s21
        do isbv1=1,nsymb                                                8d10s21
         isbv2=multh(isbv12,isbv1)                                      8d10s21
         if(isbv2.ge.isbv1)then                                         8d10s21
          if(isbv1.eq.isbv2)then                                        8d10s21
           nrow=nvirt(isbv1)*nrootu                                     8d10s21
           nn=nrow*nfdatb(2,1,isb)                                      8d26s21
           do iz=1,nn                                                   8d10s21
            gd(ioffg+iz)=0d0                                            8d10s21
           end do                                                       8d10s21
           ioffg=ioffg+nn                                               8d10s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        8d10s21
          else                                                          8d10s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                8d10s21
          end if                                                        8d10s21
          nrow=nvv*nrootu                                               8d10s21
          do l=1,4                                                      8d10s21
           nn=nrow*nfdatb(2,l,isb)                                      8d26s21
           do iz=1,nn                                                   8d10s21
            gd(ioffg+iz)=0d0                                            8d10s21
           end do                                                       8d10s21
           ioffg=ioffg+nn                                               8d10s21
          end do                                                        8d10s21
         end if                                                         8d10s21
        end do                                                          8d10s21
       end do                                                           8d10s21
      end if                                                            8d10s21
      mdoo=mdoop-1                                                      8d9s21
      ibcoffo=ibcoff
      norbx=norb+1                                                      1d8s21
      norbxx=norbx+1                                                    1d8s21
      norbxxx=norbxx+1                                                  1d8s21
      ioffdnon=1                                                        1d8s21
      nrootm=nrootu-1                                                   1d11s21
      do isb=1,nsymb                                                    1d8s21
       isbv12=multh(isb,isymket)                                        8d9s21
       joffdnon=1                                                       1d8s21
       do jsb=1,nsymb                                                   8d10s21
        ijsb=multh(jsb,isb)                                             1d9s21
        jsbv12=multh(jsb,isymbra)                                       8d9s21
        ibc0=ibcoff                                                     1d8s21
        do l=1,4                                                        1d8s21
         if(nfdatk(2,l,isb).gt.0)then                                   8d26s21
          do ll=1,4                                                     1d8s21
           if(nfdatb(2,ll,jsb).gt.0)then                                8d26s21
            nll=nfdatk(2,l,isb)*nfdatb(2,ll,jsb)                        8d26s21
            if(l.eq.ll)then                                             1d11s21
             iden1en(l)=ibcoff                                          2d24s21
             idenhvvn(l)=iden1en(l)+nll                                 2d24s21
             ibcoff=idenhvvn(l)+nll                                     2d24s21
            end if                                                      1d11s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(l.eq.ll)then                                            8d10s21
               nn=irefo(lsa)*irefo(lsb)                                 1d9s21
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
        call enough('hcddjkbk.  1',bc,ibc)
        nwds=ibcoff-ibc0                                                1d8s21
        do i=ibc0,ibcoff-1                                              1d8s21
         bc(i)=0d0                                                      1d8s21
        end do                                                          1d8s21
        ibclast=ibcoff                                                  1d17s21
        idoit=0                                                         2d23s21
        nlprt=0
        do ncloi=mdon,mdoo                                              1d8s21
         ncloip=ncloi+1                                                 1d8s21
         if(nff22k(ncloip,1,isb).gt.0)then                              8d26s21
          iarg=ncloip-mdon                                              1d8s21
          nopeni=nec-2*ncloip                                           1d8s21
          nopenip=nopeni+2                                              1d8s21
          ibcst=ibcoff                                                  1d8s21
          ibcnd=ibcoff-1                                                1d8s21
          ivcv=nfdatk(5,1,isb)+nff22k(ncloip,2,isb)                     8d26s21
          do if=1,nff22k(ncloip,1,isb)                                  8d26s21
           ipack8=iff22k(ivcv)                                          8d26s21
           itc=ipack4(1)                                                1d8s21
           ito=ipack4(2)                                                1d8s21
           ito=ibset(ito,norbx)                                         1d8s21
           ito=ibset(ito,norbxx)                                        1d8s21
           nspacei=iff22k(ivcv+1)                                       8d26s21
           do l=1,4                                                     3d19s21
            nli(l)=iff22k(ivcv+1+l)                                     8d26s21
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
            if(nff22b(nclojp,1,jsb).gt.0)then                           8d26s21
             if(mod(idoit,mynprocg).eq.mynowprog)then                   3d1s21
              jarg=nclojp-mdon                                             1d8s21
              nopenj=nec-2*nclojp                                          1d8s21
              nopenjp=nopenj+2                                           1d8s21
              jvcv=nfdatb(5,1,jsb)+nff22b(nclojp,2,jsb)                  8d26s21
              do jf=1,nff22b(nclojp,1,jsb)                               8d26s21
               ipack8=iff22b(jvcv)                                       8d26s21
               jtc=ipack4(1)                                               1d8s21
               jto=ipack4(2)                                               1d8s21
               jto=ibset(jto,norbx)                                        1d8s21
               jto=ibset(jto,norbxx)                                       1d8s21
               nspacej=iff22b(jvcv+1)                                    8d26s21
               do l=1,4                                                  3d19s21
                nlj(l)=iff22b(jvcv+1+l)                                  8d26s21
               end do                                                    3d19s21
c     j is bra, i is ket
               if(isbv12.eq.jsbv12)then                                  8d9s21
c                                                                       1d8s21
                gandcc=ieor(itc,jtc)                                      10d13s22
                gandco=ieor(ito,jto)                                      10d13s22
                gandcb=ior(gandcc,gandco)                                 10d20s22
                ndifb=popcnt(gandcb)                                      10d20s22
                if(ndifb.le.4)then                                        10d20s22
                 ndifs=popcnt(gandco)                                      10d13s22
                 ndifd=popcnt(gandcc)                                      10d13s22
                 if(ndifs.eq.0.and.ndifd.eq.0)then
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
                   is=ism(idorb(i))                                              11d23s20
                   ig=irel(idorb(i))-1+idoubo(is)                        8d9s21
                   if(isymop(lxmt).eq.1)then                             8d9s21
                    iad=ixmtf(is)+ig*(nbasdws(is)+1)                     8d9s21
                    h0=bc(iad)*phase1                                    8d9s21
                   else                                                  8d9s21
                    h0=0d0                                               8d9s21
                   end if                                                8d9s21
                   do j=1,ncloi                                                   11d23s20
                    js=ism(idorb(j))                                              11d23s20
                    jg=irel(idorb(j))-1+idoubo(js)                       8d9s21
                    xj=0d0                                               8d9s21
                    xk=0d0                                               8d9s21
                    do ii2e=1,n2e                                        8d9s21
                     if(isymop(i2eop(1,ii2e)).eq.1)then                  8d9s21
                      iss=ixmt(is,i2eop(1,ii2e))+ig*(nbasdws(is)+1)      8d9s21
                      jss=ixmt(js,i2eop(2,ii2e))+jg*(nbasdws(js)+1)      8d9s21
                      xj=xj+bc(iss)*bc(jss)*phase2                       8d9s21
                     end if                                              8d9s21
                     if(multh(is,js).eq.isymop(i2eop(1,ii2e)))then       8d9s21
                      isjs=ixmt(is,i2eop(1,ii2e))+ig+nbasdws(is)*jg      8d9s21
                      jsis=ixmt(js,i2eop(2,ii2e))+jg+nbasdws(js)*ig      8d9s21
                      xk=xk+bc(isjs)*bc(jsis)*phase2                     8d9s21
                     end if                                              8d9s21
                    end do                                               8d9s21
                    orig=sumc
                    sumc=sumc+2d0*xj-xk                                            7d15s19
                   end do                                                          7d15s19
                   orig=sumc
                   sumc=sumc+h0*2d0                                                7d15s19
                  end do                                                           7d15s19
                  sumo=0d0                                                         7d15s19
                  do i=1,nopeni                                                  11d23s20
                   is=ism(isorb(i))                                              11d23s20
                   ig=irel(isorb(i))-1+idoubo(is)                        8d9s21
                   if(isymop(lxmt).eq.1)then                             8d9s21
                    iad=ixmtf(is)+ig*(nbasdws(is)+1)                     8d9s21
                    h0=bc(iad)*phase1                                    8d9s21
                   else                                                  8d9s21
                    h0=0d0                                               8d9s21
                   end if                                                8d9s21
                   orig=sumo
                   sumo=sumo+h0                                                    7d15s19
                   do j=1,nopeni                                                 11d23s20
                    if(i.ne.j)then                                                 7d15s19
                     js=ism(isorb(j))                                            11d23s20
                     jg=irel(isorb(j))-1+idoubo(js)                      8d9s21
                     xj=0d0                                              8d9s21
                     do ii2e=1,n2e                                       8d9s21
                      if(isymop(i2eop(1,ii2e)).eq.1)then                 8d9s21
                       iss=ixmt(is,i2eop(1,ii2e))+ig*(nbasdws(is)+1)     8d9s21
                       jss=ixmt(js,i2eop(2,ii2e))+jg*(nbasdws(js)+1)     8d9s21
                       xj=xj+phase2*bc(iss)*bc(jss)                      8d9s21
                      end if                                             8d9s21
                     end do                                              8d9s21
                     orig=sumo
                     sumo=sumo+xj*0.5d0                                            7d15s19
                    end if                                                         7d15s19
                   end do                                                          7d15s19
                  end do                                                           7d15s19
                  sum=sumc+sumo
                  do i=1,ncloi                                                   11d23s20
                   is=ism(idorb(i))                                              11d23s20
                   ig=irel(idorb(i))-1+idoubo(is)                        8d9s21
                   do j=1,nopeni                                                  11d23s20
                    js=ism(isorb(j))                                              11d23s20
                    jg=irel(isorb(j))-1+idoubo(js)                       8d9s21
                    xj=0d0                                               8d9s21
                    xk=0d0                                               8d9s21
                    do ii2e=1,n2e                                        8d9s21
                     if(isymop(i2eop(1,ii2e)).eq.1)then                  8d9s21
                      iss=ixmt(is,i2eop(1,ii2e))+ig*(nbasdws(is)+1)      8d9s21
                      jss=ixmt(js,i2eop(2,ii2e))+jg*(nbasdws(js)+1)      8d9s21
                      xj=xj+phase2*bc(iss)*bc(jss)                       8d9s21
                     end if                                              8d9s21
                     if(multh(is,js).eq.isymop(i2eop(1,ii2e)))then       8d9s21
                      isjs=ixmt(is,i2eop(1,ii2e))+ig+nbasdws(is)*jg      8d9s21
                      jsis=ixmt(js,i2eop(2,ii2e))+jg+nbasdws(js)*ig      8d9s21
                      orig=xk
                      xk=xk+phase2*bc(isjs)*bc(jsis)                     8d9s21
                     end if                                              8d9s21
                    end do                                               8d9s21
                    orig=sum
                    sum=sum+xj*2d0-xk                                              7d15s19
                   end do                                                          7d15s19
                  end do                                                           7d15s19
                  sumx=sum                                               2d24s21
                  if(n2e.gt.0)then                                       2d6s23
                   do i1=1,nopeni                                                 11d23s20
                    jsa=ism(isorb(i1))                                            11d22s20
                    jga=irel(isorb(i1))-1+idoubo(jsa)                    8d9s21
                    do i2=i1+1,nopenip                                            11d19s20
                     if(i2.le.nopeni)then                                         11d23s20
                      jtesta=itc                                               11d19s20
                      jtestb=ito                                               11d19s20
                      nopenk=nopenip-2                                            11d19s20
                      karg=iarg+1                                                11d19s20
                      nqq=karg+mdon-1
                      xint=0d0
                      ksb=ism(isorb(i2))                                          11d22s20
                      kgb=irel(isorb(i2))-1+idoubo(ksb)                  8d9s21
                      xint=0d0                                           8d9s21
                      do ii2e=1,n2e                                      8d9s21
                       if(multh(jsa,ksb).eq.isymop(i2eop(1,ii2e)))then   8d9s21
                        ijk=ixmt(jsa,i2eop(1,ii2e))+jga+nbasdws(jsa)*kgb 8d9s21
                        ikj=ixmt(ksb,i2eop(2,ii2e))+kgb+nbasdws(ksb)*jga 8d9s21
                        xint=xint+phase2*bc(ijk)*bc(ikj)                 8d9s21
                       end if                                            8d9s21
                      end do                                             8d9s21
                      orig=sumx
                      sumx=sumx-xint                                     2d24s21
                      if(nqq.ge.mdon.and.nqq.le.mdoo)then                         11d19s20
                       jtestb=ibclr(jtestb,isorb(i2))                      2d10s21
                       jtestb=ibclr(jtestb,isorb(i1))                      2d10s21
                       jtesta=ibset(jtesta,isorb(i1))                      2d10s21
                       call gandc(itc,ito,jtesta,jtestb,nopenip,nopenk,      11d19s20
     $             iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot,nab,iwpb,iwpk,       1d17s21
     $                      ncsfmid,bc,ibc)                                11d14s22
                       call gandc(jtesta,jtestb,itc,ito,nopenk,nopenip,  2d25s21
     $             karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,    2d25s21
     $                     iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                       iargo=1                                           2d24s21
                       do l=1,4                                          2d24s21
                        if(nli(l).gt.0)then                              2d24s21
                         iad1=ivcv+iff22k(ivcv+5+l)                      8d26s21
                         iad2=iad1+nli(l)                                    3d19s21
                         itmp1=ibcoff                                    2d24s21
                         itmp2=itmp1+nli(l)*ncsf(karg)                   2d24s21
                         itmp3=itmp2+nli(l)*ncsf(karg)                   2d24s21
                         ibcoff=itmp3+nli(l)*nli(l)                      2d24s21
                         call enough('hcddjkbk.  2',bc,ibc)
                         call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1,   2d25s21
     $                     nli(l),iwpb1,iwpk1,ff22k(iad2),ncsf2(l,iarg),8d26s21
     $                       bc(itmp1),ncsf(karg),1d0,0d0,iargo,        2d24s21
     $                       ncsf2(l,iarg),bc,ibc)                      11d10s22
                         do iii=0,nli(l)-1                               2d24s21
                          do k=0,ncsf(karg)-1                            2d24s21
                           ki=itmp1+k+ncsf(karg)*iii                     2d24s21
                           ik=itmp2+iii+nli(l)*k                         2d24s21
                           bc(ik)=bc(ki)                                 2d24s21
                          end do                                         2d24s21
                         end do                                          2d24s21
                         call dgemm('n','n',nli(l),nli(l),ncsf(karg),    2d24s21
     $                       xint,bc(itmp2),nli(l),bc(itmp1),ncsf(karg),2d24s21
     $                       0d0,bc(itmp3),nli(l),                      2d24s21
     d' hcddjkbk.  1')
                         jtmp3=itmp3                                     2d24s21
                         do j=0,nli(l)-1                                  2d23s21
                          jj=iff22k(iad1+j)-1                            8d26s21
                          jden=iden1en(l)+nfdatk(2,l,isb)*jj-1           8d26s21
                          do k=0,nli(l)-1                                2d24s21
                           kk=iff22k(iad1+k)                             8d26s21
                           bc(jden+kk)=bc(jden+kk)+bc(jtmp3+k)           2d24s21
                          end do                                         2d24s21
                          jtmp3=jtmp3+nli(l)                             2d24s21
                         end do                                          2d24s21
                         ibcoff=itmp1                                    2d24s21
                        end if                                           2d24s21
                        iargo=iargo+ncsf2(l,iarg)                        2d24s21
                       end do                                            2d24s21
                      end if                                                      11d19s20
                     end if                                                       11d23s20
                    end do                                                1d18s21
                   end do                                                 1d18s21
                  end if                                                 2d6s23
                  iargo=1                                                2d24s21
                  do l=1,4                                               2d24s21
                   if(nli(l).gt.0)then                                   2d24s21
                    itmpt=ibcoff                                         2d24s21
                    itmpm=itmpt+ncsf2(l,iarg)*nli(l)                     2d24s21
                    ibcoff=itmpm+nli(l)*nli(l)                           2d24s21
                    call enough('hcddjkbk.  3',bc,ibc)
                    iad1=ivcv+iff22k(ivcv+5+l)                           8d26s21
                    iad2=iad1+nli(l)                                     3d19s21
                    do iii=0,nli(l)-1                                    2d24s21
                     do j=0,ncsf2(l,iarg)-1                              2d24s21
                      ji=iad2+j+ncsf2(l,iarg)*iii                        2d24s21
                      ij=itmpt+iii+nli(l)*j                              2d24s21
                      bc(ij)=ff22k(ji)                                   8d26s21
                     end do                                              2d24s21
                    end do                                               2d24s21
                    call dgemm('n','n',nli(l),nli(l),ncsf2(l,iarg),sumx, 2d24s21
     $                  bc(itmpt),nli(l),ff22k(iad2),ncsf2(l,iarg),0d0, 8d10s21
     $                   bc(itmpm),nli(l),                              2d24s21
     d' hcddjkbk.  2')
                    jtmpm=itmpm                                          2d24s21
                    do j=0,nli(l)-1                                      2d24s21
                     jj=iff22k(iad1+j)-1                                 8d26s21
                     jden=iden1en(l)+nfdatk(2,l,isb)*jj-1                8d26s21
                     do k=0,nli(l)-1                                     2d24s21
                      kk=iff22k(iad1+k)                                  8d26s21
                      bc(jden+kk)=bc(jden+kk)+bc(jtmpm+k)                2d24s21
                     end do                                              2d24s21
                     jtmpm=jtmpm+nli(l)                                  2d24s21
                    end do                                               2d24s21
                    ibcoff=itmpt                                         2d24s21
                   end if                                                2d24s21
                   iargo=iargo+ncsf2(l,iarg)                             2d24s21
                  end do                                                 2d24s21
                 else if(ndifs.eq.2.and.ndifb.eq.2)then                  2d6s23
                  do i=1,norbxx
                   if(btest(gandco,i))then
                    if((btest(itc,i).and.btest(jto,i)).or.                 10d21s22
     $                   (btest(ito,i).and..not.btest(jtc,i)))then            10d21s22
                     nab4(1,1)=i
                    else
                     nab4(2,1)=i
                    end if
                   end if
                  end do
                  call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,  1d8s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $               ncsfmid1,bc,ibc)                                   11d14s22
                  call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,  2d24s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,  2d24s21
     $               ncsfmid1b,bc,ibc)                                  11d14s22
                  ksa=ism(nab1(1))                                       1d8s21
                  kga=irel(nab1(1))-1+idoubo(ksa)                        8d9s21
                  ksb=ism(nab1(2))                                        1d8s21
                  kgb=irel(nab1(2))-1+idoubo(ksb)                        8d9s21
                  if(multh(ksa,ksb).eq.isymop(lxmt))then                 8d9s21
                   ih0=ixmtf(ksb)+kgb+nbasdws(ksb)*kga                   8d9s21
                   sum=bc(ih0)*phase1                                    8d11s21
                  else                                                   8d9s21
                   sum=0d0                                               8d9s21
                  end if                                                 8d9s21
                  do i=1,norb                                            12d18s20
                   itest(i,2)=0                                          12d18s20
                  end do                                                 12d18s20
                  do i=1,norb                                            12d18s20
                   if(btest(jtc,i))then                                  12d18s20
                    itest(i,2)=2                                         12d18s20
                   end if                                                12d18s20
                   if(btest(jto,i))then                                  12d18s20
                    itest(i,2)=1                                          12d18s20
                   end if                                                 12d18s20
                  end do                                                  12d18s20
                  nok=0                                                         11d13s20
                  if(n2e.gt.0)then                                        2d6s23
                   do i=1,norb                                            1d8s21
                    ixn=min(itest(i,1),itest(i,2))
                    if(ixn.gt.0)then                                             11d13s20
                     nok=nok+1                                                   11d13s20
                     itest(nok,3)=ixn                                            11d13s20
                     itest(nok,2)=i                                              11d13s20
                    end if                                                       11d13s20
                   end do                                                        11d13s20
                  end if                                                  2d6s23
                  do i=1,nok                                              1d8s21
                   lsa=ism(itest(i,2))                                    1d8s21
                   lga=irel(itest(i,2))-1+idoubo(lsa)                     8d9s21
                   xj=0d0                                                 8d9s21
                   do ii2e=1,n2e                                          8d9s21
                    if(isymop(i2eop(1,ii2e)).eq.1)then                    8d9s21
                     ill=ixmt(lsa,i2eop(1,ii2e))+lga+nbasdws(lsa)*lga     8d9s21
                     ikk=ixmt(ksb,i2eop(1,ii2e))+kgb+nbasdws(ksb)*kga     8d9s21
                     xj=xj+phase2*bc(ill)*bc(ikk)                         8d9s21
                    end if                                                8d9s21
                   end do                                                 8d9s21
                   sum=sum+xj                                             1d8s21
                   if(itest(i,3).eq.2)then                                1d8s21
                    xk=0d0                                                8d9s21
                    do ii2e=1,n2e                                         8d9s21
                     if(multh(lsa,ksa).eq.isymop(i2eop(1,ii2e)))then      8d9s21
                      ilka=ixmt(lsa,i2eop(1,ii2e))+lga+nbasdws(lsa)*kga   8d9s21
                      iklb=ixmt(ksb,i2eop(2,ii2e))+kgb+nbasdws(ksb)*lga   8d9s21
                      xk=xk+phase2*bc(ilka)*bc(iklb)                      8d9s21
                     end if                                               8d9s21
                    end do                                                8d9s21
                    sum=sum+xj-xk                                         1d8s21
                   end if                                                 1d8s21
                  end do                                                  1d8s21
                  jargo=1                                                 2d23s21
                  iargo=1                                                 2d24s21
                  do l=1,4                                                1d8s21
                   if(min(nlj(l),nli(l)).gt.0)then                        2d25s21
                    jad1=jvcv+iff22b(jvcv+5+l)                            8d26s21
                    jad2=jad1+nlj(l)                                      3d19s21
                    iad1=ivcv+iff22k(ivcv+5+l)                            8d26s21
                    iad2=iad1+nli(l)                                      3d19s21
                    itmp1=ibcoff                                          2d23s21
                    itmp1b=itmp1+ncsfmid1*nlj(l)                          2d23s21
                    ibcoff=itmp1b+ncsfmid1*nli(l)                         2d23s21
                    call enough('hcddjkbk.  4',bc,ibc)
                    if(iwpk1.lt.0)then                                    2d23s21
                     nusedi=ibc(-iwpk1)/2                                 2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1           2d23s21
                     icmp1=-iwpk1+1                                       2d23s21
                     icmp2=icmp1+ncsf(jarg)                               2d23s21
                     icmp3=icmp2+nusedi                                   2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),    2d23s21
     $                   ff22b(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),8d26s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d14s22
                    else                                                              12d5s20
                     jwpk1=iwpk1+ncsfmid1*(jargo-1)                       2d23s21
                     call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                 1d0,bc(jwpk1),ncsfmid1,ff22b(jad2),ncsf2(l,jarg),8d26s21
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjkbk.  3')
                    end if                                                            12d5s20
                    if(iwpk1b.lt.0)then                                   2d23s21
                     nusedi=ibc(-iwpk1b)/2                                2d23s21
                     if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                     icmp1=-iwpk1b+1                                      2d23s21
                     icmp2=icmp1+ncsf(iarg)                               2d23s21
                     icmp3=icmp2+nusedi                                   2d23s21
                     call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),    2d23s21
     $                  ff22k(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),8d26s21
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          11d14s22
                    else                                                              12d5s20
                     jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                     call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                1d0,bc(jwpk1b),ncsfmid1,ff22k(iad2),ncsf2(l,iarg),8d26s21
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjkbk.  4')
                    end if                                                            12d5s20
                    itmpt=ibcoff                                          2d23s21
                    itmpm=itmpt+nli(l)*ncsfmid1                           2d24s21
                    ibcoff=itmpm+nli(l)*nlj(l)                            2d24s21
                    call enough('hcddjkbk.  5',bc,ibc)
                    do iii=0,nli(l)-1                                     2d23s21
                     do j=0,ncsfmid1-1                                    2d23s21
                      ji=itmp1b+j+ncsfmid1*iii                            2d23s21
                      ij=itmpt+iii+nli(l)*j                               2d23s21
                      bc(ij)=bc(ji)                                       2d23s21
                     end do                                               2d23s21
                    end do                                                2d23s21
                    call dgemm('n','n',nli(l),nlj(l),ncsfmid1,sum,        2d24s21
     $                   bc(itmpt),nli(l),bc(itmp1),ncsfmid1,0d0,        2d23s21
     $                   bc(itmpm),nli(l),                               2d23s21
     d' hcddjkbk.  5')
                    jtmpm=itmpm                                           2d24s21
                    do iii=0,nlj(l)-1                                     2d24s21
                     jj=iff22b(jad1+iii)-1                                8d26s21
                     jjden=iden1en(l)+nfdatk(2,l,isb)*jj-1                8d26s21
                     do j=0,nli(l)-1                                      2d24s21
                      kk=iff22k(iad1+j)                                   8d26s21
                      bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)               2d24s21
                     end do                                               2d24s21
                     jtmpm=jtmpm+nli(l)                                   2d24s21
                    end do                                                2d24s21
                    ibcoff=itmp1                                          2d24s21
                   end if                                                 1d8s21
                   jargo=jargo+ncsf2(l,jarg)                              2d23s21
                   iargo=iargo+ncsf2(l,iarg)                              2d24s21
                  end do                                                  1d8s21
                  do i=1,nok                                              1d8s21
                   if(itest(i,3).eq.1)then                                2d6s23
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
                    nnot1=0                                               8d4s22
                    nnot2=0                                               8d4s22
                    if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                     call gandc(itestc,itesto,itc,ito,nopenk,nopenip,     2d24s21
     $                karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,    2d24s21
     $                iwpb1b,iwpk1b,ncsfmid1b,bc,ibc)                   11d14s22
                     nnot1=nnot1b                                         2d25s21
                     nab1(1)=nab1b(2)                                     2d25s21
                     nab1(2)=nab1b(1)                                     2d25s21
                     call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,        12d14s20
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            11d14s22
                    end if
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     lsa=ism(nab1(1))                                     1d8s21
                     lga=irel(nab1(1))-1+idoubo(lsa)                      8d9s21
                     lsb=ism(nab1(2))                                     1d8s21
                     lgb=irel(nab1(2))-1+idoubo(lsb)                      8d9s21
                     lsc=ism(nab2(1))                                     1d8s21
                     lgc=irel(nab2(1))-1+idoubo(lsc)                      8d9s21
                     lsd=ism(nab2(2))                                     1d8s21
                     lgd=irel(nab2(2))-1+idoubo(lsd)                      8d9s21
                     xint=0d0                                             8d9s21
                     do ii2e=1,n2e                                        8d9s21
                      if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)).and.
     $                   multh(lsc,lsd).eq.isymop(i2eop(2,ii2e)))then   8d9s21
                       iab=ixmt(lsa,i2eop(1,ii2e))+lga+nbasdws(lsa)*lgb   8d9s21
                       icd=ixmt(lsc,i2eop(2,ii2e))+lgc+nbasdws(lsc)*lgd   8d9s21
                       xint=xint+phase2*bc(iab)*bc(icd)                   8d9s21
                      end if                                              8d9s21
                     end do                                               8d9s21
                     jargo=1                                              2d23s21
                     iargo=1                                              2d24s21
                     do l=1,4                                             1d8s21
                      if(min(nli(l),nlj(l)).gt.0)then                     3d1s21
                       jad1=jvcv+iff22b(jvcv+5+l)                         8d26s21
                       jad2=jad1+nlj(l)                                   3d19s21
                       iad1=ivcv+iff22k(ivcv+5+l)                         8d26s21
                       iad2=iad1+nli(l)                                   3d19s21
                       itmpj=ibcoff                                       2d24s21
                       ibcoff=itmpj+ncsf(karg)*nlj(l)                     2d24s21
                       call enough('hcddjkbk.  6',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                    nlj(l),iwpb2,iwpk2,ff22b(jad2),ncsf2(l,jarg),  8d26s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d10s22
                       itmpi=ibcoff                                       2d24s21
                       itmpt=itmpi+ncsf(karg)*nli(l)                      2d24s21
                       ibcoff=itmpt+ncsf(karg)*nli(l)                     2d24s21
                       call enough('hcddjkbk.  7',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,ff22k(iad2),ncsf2(l,iarg),8d26s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          11d10s22
                       do iii=0,nli(l)-1                                  2d24s21
                        do k=0,ncsf(karg)-1                                2d23s21
                         ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                         ik=itmpi+iii+nli(l)*k                            2d23s21
                         bc(ik)=bc(ki)                                     2d23s21
                        end do                                             2d23s21
                       end do                                              2d23s21
                       itmpm=itmpt                                        2d24s21
                       ibcoff=itmpm+nli(l)*nlj(l)                         2d24s21
                       call enough('hcddjkbk.  8',bc,ibc)
                       call dgemm('n','n',nli(l),nlj(l),ncsf(karg),xint,  2d6s23
     $                      bc(itmpi),nli(l),bc(itmpj),ncsf(karg),0d0,  2d24s21
     $                      bc(itmpm),nli(l),                           2d24s21
     d' hcddjkbk.  6')
                       jtmpm=itmpm                                        2d24s21
                       do iii=0,nlj(l)-1                                   2d24s21
                        jj=iff22b(jad1+iii)-1                             8d26s21
                        jjden=iden1en(l)+nfdatk(2,l,isb)*jj-1             8d26s21
                        do j=0,nli(l)-1                                    2d24s21
                         kk=iff22k(iad1+j)                                8d26s21
                         bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)             2d24s21
                        end do                                             2d24s21
                        jtmpm=jtmpm+nli(l)                                 2d24s21
                       end do                                              2d24s21
                       ibcoff=itmpj                                       2d24s21
                      end if                                              1d8s21
                      jargo=jargo+ncsf2(l,jarg)                           2d23s21
                      iargo=iargo+ncsf2(l,iarg)                           2d24s21
                     end do                                               1d8s21
                    end if                                                  12d14s20
                   end if                                                   12d14s20
                  end do                                                  1d8s21
                 else if(n2e.gt.0)then                                    2d6s23
                  nnot=0                                                  2d6s23
                  ipssx=0                                                 2d6s23
                  if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s
                   nnot=4                                                          10d14s22
                   ioxx(1)=1                                                     10d17s22
                   ioxx(2)=1                                                     10d17s22
                   do i=1,norbxx                                                     10d17s22
                    if(btest(gandcb,i))then                                         10d14s22
                     if((btest(jtc,i).and.btest(ito,i)).or.                         10d17s22
     $                (btest(jto,i).and..not.btest(itc,i)))then                   10d14s22
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
     $                 ((btest(itc,i).and..not.btest(jto,i)).or.                 10d14s22
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
                    if(btest(itc,nab4(1,2)).and.                        2d6s23
     $                   .not.btest(itc,nab4(1,1)))nbt=1                2d6s23
                   else                                                             10d17s22
                    nbt=0                                                           10d17s22
                    if(btest(jtc,nab4(2,2)).and.                        2d7s23
     $                   .not.btest(jtc,nab4(2,1)))nbt=1                2d7s23
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
                  if(nnot.eq.3)then                                          12d8s20
                   ipssx=1                                                   12d8s20
                  else if(nnot.eq.4)then                                  2d6s23
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
                    write(6,*)('bit not set for nab4(1,1) ='),            3d1s21
     $                   nab4(1,iu1)                                    3d1s21
                    stop 'nab4(1,1)a'                                           11d27s20
                   end if                                                        11d13s20
                   if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                    itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                    itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                    nopenk=nopenk-1                                              11d13s20
                    karg=karg+1                                                  11d13s20
                   else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                    write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                   nab4(2,iu2)                                       12d18s20
                    stop 'nab4(2,1)'                                           11d27s20
                   else                                                          11d13s20
                    itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                    nopenk=nopenk+1                                              11d13s20
                   end if                                                        11d13s20
                   nqq=karg+mdon-1                                        4d20s21
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                    4d20s21
                    call gandc(itestc,itesto,itc,ito,nopenk,nopenip,      2d25s21
     $         karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,     2d25s21
     $                  iwpk1b,ncsfmid1b,bc,ibc)                        11d14s22
                    nnot1=nnot1b                                          2d25s21
                    nab1(1)=nab1b(2)                                      2d25s21
                    nab1(2)=nab1b(1)                                      2d25s21
                    call gandc(itestc,itesto,jtc,jto,nopenk,nopenjp,       1d8s21
     $         karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d18s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                    if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                     lsa=ism(nab1(1))                                     1d8s21
                     lga=irel(nab1(1))-1+idoubo(lsa)                      8d9s21
                     lsb=ism(nab1(2))                                     1d8s21
                     lgb=irel(nab1(2))-1+idoubo(lsb)                      8d9s21
                     lsc=ism(nab2(1))                                     1d8s21
                     lgc=irel(nab2(1))-1+idoubo(lsc)                      8d9s21
                     lsd=ism(nab2(2))                                     1d8s21
                     lgd=irel(nab2(2))-1+idoubo(lsd)                      8d9s21
                     xint=0d0                                             8d9s21
                     do ii2e=1,n2e                                        8d9s21
                      if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)).and.     8d9s21
     $                     multh(lsc,lsd).eq.isymop(i2eop(2,ii2e)))then   8d9s21
                       iab=ixmt(lsa,i2eop(1,ii2e))+lga+nbasdws(lsa)*lgb   8d9s21
                       icd=ixmt(lsc,i2eop(2,ii2e))+lgc+nbasdws(lsc)*lgd   8d9s21
                       xint=xint+phase2*bc(iab)*bc(icd)                   8d9s21
                      end if                                              8d9s21
                     end do                                               8d9s21
                     if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2))      4d20s21
     $                    .and.                                         4d20s21
     $                max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))     1d12s21
     $                  xint=xint*0.5d0                                 1d12s21
                     jargo=1                                              2d23s21
                     iargo=1                                              2d24s21
                     do l=1,4                                               12d19s20
                      if(min(nli(l),nlj(l)).gt.0)then                     3d1s21
                       jad1=jvcv+iff22b(jvcv+5+l)                         8d26s21
                       jad2=jad1+nlj(l)                                   3d19s21
                       iad1=ivcv+iff22k(ivcv+5+l)                         8d26s21
                       iad2=iad1+nli(l)                                   3d19s21
                       itmpj=ibcoff                                       2d24s21
                       ibcoff=itmpj+ncsf(karg)*nlj(l)                     2d24s21
                       call enough('hcddjkbk.  9',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(l),iwpb2,iwpk2,ff22b(jad2),ncsf2(l,jarg),  8d26s21
     $                   bc(itmpj),ncsf(karg),1d0,0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d10s22
                       itmpi=ibcoff                                       2d24s21
                       itmpt=itmpi+ncsf(karg)*nli(l)                      2d24s21
                       ibcoff=itmpt+ncsf(karg)*nli(l)                     2d24s21
                       call enough('hcddjkbk. 10',bc,ibc)
                       call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                   nli(l),iwpb1b,iwpk1b,ff22k(iad2),ncsf2(l,iarg),8d26s21
     $                   bc(itmpt),ncsf(karg),1d0,0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          11d10s22
                       do iii=0,nli(l)-1                                  2d24s21
                        do k=0,ncsf(karg)-1                                2d23s21
                         ki=itmpt+k+ncsf(karg)*iii                         2d23s21
                         ik=itmpi+iii+nli(l)*k                            2d23s21
                         bc(ik)=bc(ki)                                     2d23s21
                        end do                                             2d23s21
                       end do                                              2d23s21
                       itmpm=itmpt                                        2d24s21
                       ibcoff=itmpm+nli(l)*nlj(l)                         2d24s21
                       call enough('hcddjkbk. 11',bc,ibc)
                       call dgemm('n','n',nli(l),nlj(l),ncsf(karg),xint,  2d24s21
     $                      bc(itmpi),nli(l),bc(itmpj),ncsf(karg),0d0,  2d24s21
     $                      bc(itmpm),nli(l),                           2d24s21
     d' hcddjkbk.  7')
                       jtmpm=itmpm                                        2d24s21
                       do iii=0,nlj(l)-1                                   2d24s21
                        jj=iff22b(jad1+iii)-1                             8d26s21
                        jjden=iden1en(l)+nfdatk(2,l,isb)*jj-1             8d26s21
                        do j=0,nli(l)-1                                    2d24s21
                         kk=iff22k(iad1+j)                                8d26s21
                         bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)             2d24s21
                        end do                                             2d24s21
                        jtmpm=jtmpm+nli(l)                                 2d24s21
                       end do                                              2d24s21
                       ibcoff=itmpj                                       2d24s21
                      end if                                              2d25s21
                      jargo=jargo+ncsf2(l,jarg)                           2d23s21
                      iargo=iargo+ncsf2(l,iarg)                           2d24s21
                     end do                                                 12d19s20
                     if(ipss.eq.2)go to 3                                   12d18s20
                    end if                                                  12d18s20
                   end if                                                 4d20s21
                  end do                                                   12d18s20
    3             continue                                                 12d18s20
                 end if                                                   1d8s21
                end if                                                    1d8s21
               end if                                                    1d8s21
               jto=ibclr(jto,norbx)                                      1d8s21
               jto=ibset(jto,norbxxx)                                    1d8s21
c                                                                       1d8s21
               gandcc=ieor(itc,jtc)                                      10d13s22
               gandco=ieor(ito,jto)                                      10d13s22
               gandcb=ior(gandcc,gandco)                                 10d20s22
               ndifb=popcnt(gandcb)                                      10d20s22
               if(ndifb.le.4)then                                        10d20s22
                ndifs=popcnt(gandco)                                      10d13s22
                ndifd=popcnt(gandcc)                                      10d13s22
                if(ndifs.eq.2.and.ndifb.eq.2)then                       2d7s23
                 do i=1,norbxxx                                         2d7s23
                  if(btest(gandco,i))then                               2d7s23
                   if((btest(itc,i).and.btest(jto,i)).or.               2d7s23
     $             (btest(ito,i).and..not.btest(jtc,i)))then             2d7s23
                    nab4(1,1)=i                                         2d7s23
                   else                                                 2d7s23
                    nab4(2,1)=i                                         2d7s23
                   end if                                               2d7s23
                  end if                                                2d7s23
                 end do                                                 2d7s23
                 if(isb.ne.jsb)then                                      1d9s21
                  write(6,*)('hey, what''s up with isb ne jsb ? ')       1d9s21
                  stop 'hcddjkbk:g'                                                   1d9s21
                 end if                                                  1d9s21
                 call gandc(itc,ito,jtc,jto,nopenip,nopenjp,iarg,jarg,   1d8s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $              ncsfmid1,bc,ibc)                                    11d14s22
                 call gandc(jtc,jto,itc,ito,nopenjp,nopenip,jarg,iarg,   2d23s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1b,nab1b,iwpb1b,iwpk1b,      1d8s21
     $              ncsfmid1b,bc,ibc)                                   11d14s22
                 itmpn=ibcoff                                            2d23s21
                 nlij=0                                                 2d23s21
                 do l=1,4                                               2d23s21
                  nlij=nlij+nli(l)*nlj(l)                               2d23s21
                 end do                                                 2d23s21
                 itmpn=ibcoff                                           2d23s21
                 ibcoff=itmpn+nlij                                      2d23s21
                 call enough('hcddjkbk. 12',bc,ibc)
                 jtmpn=itmpn                                            2d23s21
                 jargo=1                                                 2d22s21
                 iargo=1                                                2d23s21
                 do l=1,4                                                1d8s21
                  if(min(nlj(l),nli(l)).gt.0)then                        2d23s21
                   iad1=ivcv+iff22k(ivcv+5+l)                           8d26s21
                   iad2=iad1+nli(l)                                     3d19s21
                   jad1=jvcv+iff22b(jvcv+5+l)                           8d26s21
                   jad2=jad1+nlj(l)                                     3d19s21
                   itmp1=ibcoff                                         2d23s21
                   itmp1b=itmp1+ncsfmid1*nlj(l)                         2d23s21
                   ibcoff=itmp1b+ncsfmid1*nli(l)                        2d23s21
                   call enough('hcddjkbk. 13',bc,ibc)
                   if(iwpk1.lt.0)then                                   2d23s21
                    nusedi=ibc(-iwpk1)/2                                2d23s21
                    if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1          2d23s21
                    icmp1=-iwpk1+1                                      2d23s21
                    icmp2=icmp1+ncsf(jarg)                              2d23s21
                    icmp3=icmp2+nusedi                                  2d23s21
                    call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                   ff22b(jad2),ncsf2(l,jarg),ncsf(jarg),bc(itmp1),8d26s21
     $                   ncsfmid1,ncsfmid1,nlj(l),0d0,jargo,            2d23s21
     $                   ncsf2(l,jarg),bc,ibc)                          11d14s22
                   else                                                              12d5s20
                    jwpk1=iwpk1+ncsfmid1*(jargo-1)                      2d23s21
                    call dgemm('n','n',ncsfmid1,nlj(l),ncsf2(l,jarg),
     $                 1d0,bc(jwpk1),ncsfmid1,ff22b(jad2),ncsf2(l,jarg),8d26s21
     $                   0d0,bc(itmp1),ncsfmid1,
     d' hcddjkbk.  8')
                   end if                                                            12d5s20
                   if(iwpk1b.lt.0)then                                   2d23s21
                    nusedi=ibc(-iwpk1b)/2                                2d23s21
                    if(2*nusedi.ne.ibc(-iwpk1b))nusedi=nusedi+1          2d23s21
                    icmp1=-iwpk1b+1                                      2d23s21
                    icmp2=icmp1+ncsf(iarg)                              2d23s21
                    icmp3=icmp2+nusedi                                  2d23s21
                    call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),   2d23s21
     $                  ff22k(iad2),ncsf2(l,iarg),ncsf(iarg),bc(itmp1b),8d26s21
     $                   ncsfmid1,ncsfmid1,nli(l),0d0,iargo,            2d23s21
     $                   ncsf2(l,iarg),bc,ibc)                          11d14s22
                   else                                                              12d5s20
                    jwpk1b=iwpk1b+ncsfmid1*(iargo-1)                      2d23s21
                    call dgemm('n','n',ncsfmid1,nli(l),ncsf2(l,iarg),
     $                1d0,bc(jwpk1b),ncsfmid1,ff22k(iad2),ncsf2(l,iarg),8d26s21
     $                   0d0,bc(itmp1b),ncsfmid1,
     d' hcddjkbk.  9')
                   end if                                                            12d5s20
                   itmpt=ibcoff                                         2d23s21
                   ibcoff=itmpt+nli(l)*ncsfmid1                         2d23s21
                   call enough('hcddjkbk. 14',bc,ibc)
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
     d' hcddjkbk. 10')
                   ibcoff=itmp1                                         2d23s21
                   ktmpn=jtmpn                                          2d24s21
                   do iii=0,nlj(l)-1                                    2d24s21
                    jj=iff22b(jad1+iii)-1                               8d26s21
                    jjden=idenhvvn(l)+nfdatk(2,l,isb)*jj-1              8d26s21
                    do j=0,nli(l)-1                                     2d24s21
                     kk=iff22k(iad1+j)                                  8d26s21
                     bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)              2d24s21
                    end do                                              2d24s21
                    ktmpn=ktmpn+nli(l)                                  2d24s21
                   end do                                               2d24s21
                   jtmpn=jtmpn+nli(l)*nlj(l)                            2d23s21
                  end if                                                2d23s21
                  jargo=jargo+ncsf2(l,jarg)                              2d22s21
                  iargo=iargo+ncsf2(l,iarg)                             2d24s21
                 end do                                                  1d8s21
                 nok=0                                                    1d8s21
                 if(n2e.gt.0)then                                       8d21s21
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
                   icolk=lga+irefo(lsa)*lga                               1d8s21
                   icolj=icolk                                           8d10s21
                   if(itest(i,3).eq.2)then                                1d8s21
                    jtmpn=itmpn                                           2d23s21
                    do l=1,4                                              1d8s21
                     if(min(nli(l),nlj(l)).gt.0)then                      2d23s21
                      jad1=jvcv+iff22b(jvcv+5+l)                        8d26s21
                      iad1=ivcv+iff22k(ivcv+5+l)                        8d26s21
                      jden=idenjn(l,lsa)                                2d7s23
     $                      +nfdatk(2,l,isb)*nfdatb(2,l,jsb)*icolj      2d7s23
                      kden=idenkn(l,l,lsa)+nfdatk(2,l,isb)              8d26s21
     $                     *nfdatb(2,l,jsb)*icolk                       8d26s21
                      do iii=0,nlj(l)-1                                  2d23s21
                       jj=iff22b(jad1+iii)-1                            8d26s21
                       jjden=jden+nfdatk(2,l,isb)*jj-1                  8d26s21
                       kkden=kden+nfdatk(2,l,isb)*jj-1                  8d26s21
                       do j=0,nli(l)-1                                   2d23s21
                        kk=iff22k(iad1+j)                               8d26s21
                        bc(jjden+kk)=bc(jjden+kk)+bc(jtmpn+j)*2d0        2d23s21
                        bc(kkden+kk)=bc(kkden+kk)-bc(jtmpn+j)            2d24s21
                       end do                                            2d23s21
                       jtmpn=jtmpn+nli(l)                                2d23s21
                      end do                                             2d23s21
                     end if                                               2d23s21
                    end do                                                1d8s21
                   else                                                   1d8s21
                    jtmpn=itmpn                                           2d23s21
                    do l=1,4                                              1d8s21
                     if(nlj(l).gt.0)then                                  1d8s21
                      jad1=jvcv+iff22b(jvcv+5+l)                        8d26s21
                      iad1=ivcv+iff22k(ivcv+5+l)                        8d26s21
                      jden=idenjn(l,lsa)+nfdatk(2,l,isb)                2d7s23
     $                     *nfdatb(2,l,jsb)*icolj                       2d7s23
                      do iii=0,nlj(l)-1                                  2d23s21
                       jj=iff22b(jad1+iii)-1                            8d26s21
                       jjden=jden+nfdatk(2,l,isb)*jj-1                  8d26s21
                       do j=0,nli(l)-1                                   2d23s21
                        kk=iff22k(iad1+j)                               8d26s21
                        bc(jjden+kk)=bc(jjden+kk)+bc(jtmpn+j)            2d23s21
                       end do                                            2d23s21
                       jtmpn=jtmpn+nli(l)                                2d23s21
                      end do                                             2d23s21
                     end if                                               1d8s21
                    end do                                                1d8s21
                   end if                                                 1d8s21
                  end do                                                  1d8s21
                 end if                                                 8d21s21
                 ibcoff=itmpn                                           4d28s21
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
                   nnot1=0                                              8d4s22
                   nnot2=0                                              8d4s22
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
                    ia=nab2(2)
                    ib=nab1(1)
                    lsa=ism(ia)                                          1d8s21
                    lga=irel(ia)-1                                       1d8s21
                    lsb=ism(ib)                                          1d8s21
                    lgb=irel(ib)-1                                       1d8s21
                    lsab=multh(lsa,lsb)                                  1d9s21
                    icol=lga+irefo(lsa)*lgb                              1d8s21
                    if(lsab.ne.ijsb)then                                 1d9s21
                     write(6,*)('hey, lsab = '),lsab,(' ne ijsb = '),
     $                   ijsb
                     stop 'hcddjkbk:i'
                    end if                                               1d9s21
                    jargo=1                                               2d23s21
                    itmpj=ibcoff                                          2d23s21
                    jtmpj=itmpj                                           2d23s21
                    do lj=1,4                                             2d23s21
                     if(nlj(lj).gt.0)then                                 2d23s21
                      jad1=jvcv+iff22b(jvcv+5+lj)                       8d26s21
                      jad2=jad1+nlj(lj)                                     3d19s21
                      ibcoff=jtmpj+ncsf(karg)*nlj(lj)                     2d23s21
                      call enough('hcddjkbk. 15',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(lj),iwpb2,iwpk2,ff22b(jad2),ncsf2(lj,jarg),8d26s21
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
                      iad1=ivcv+iff22k(ivcv+5+li)                       8d26s21
                      iad2=iad1+nli(li)                                     3d19s21
                      ibcoff=jtmpi+ncsf(karg)*nli(li)                     2d23s21
                      itmpt=ibcoff                                        2d23s21
                      ibcoff=itmpt+ncsf(karg)*nli(li)                     2d23s21
                      call enough('hcddjkbk. 16',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                 nli(li),iwpb1b,iwpk1b,ff22k(iad2),ncsf2(li,iarg),8d26s21
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
                     if(nlj(lj).gt.0)then                               2d24s21
                      jad1=jvcv+iff22b(jvcv+5+lj)                       8d26s21
                      jtmpi=itmpi                                          2d23s21
                      do li=1,4                                         2d24s21
                       if(nli(li).gt.0)then                             2d24s21
                        iad1=ivcv+iff22k(ivcv+5+li)                     8d26s21
                        jden=idenkn(li,lj,lsa)+nfdatk(2,li,isb)         8d26s21
     $                       *nfdatb(2,lj,jsb)*icol                     8d26s21
                        itmpm=ibcoff                                       2d23s21
                        ibcoff=itmpm+nlj(lj)*nli(li)                    2d24s21
                        call enough('hcddjkbk. 17',bc,ibc)
                        call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),  2d24s21
     $                       1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),2d24s21
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' hcddjkbk. 11')
                        jtmpm=itmpm                                        2d23s21
                        do iii=0,nlj(lj)-1                              2d24s21
                         jj=iff22b(jad1+iii)-1                          8d26s21
                         jjden=jden+nfdatk(2,li,isb)*jj-1               8d26s21
                         do j=0,nli(li)-1                               2d24s21
                          kk=iff22k(iad1+j)                             8d26s21
                          bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)            2d23s21
                         end do                                            2d23s21
                         jtmpm=jtmpm+nli(li)                                2d23s21
                        end do                                             2d23s21
                        ibcoff=itmpm                                       2d23s21
                       end if                                              2d23s21
                       jtmpi=jtmpi+nli(li)*ncsf(karg)                   2d24s21
                      end do                                            2d24s21
                     end if                                             2d24s21
                     jtmpj=jtmpj+nlj(lj)*ncsf(karg)                     2d24s21
                    end do                                               2d23s21
                   end if                                                  12d14s20
                  end if                                                   12d14s20
                 end do                                                  12d18s20
                else if(n2e.gt.0)then                                   2d7s23
                 nnot=0                                                 2d7s23
                 if(ndifs.eq.4.and.ndifb.eq.4)then                      2d7s23
                  nnot=4                                                2d7s23
                  ioxx(1)=1                                             2d7s23
                  ioxx(2)=1                                             2d7s23
                  do i=1,norbxxx                                        2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if((btest(jtc,i).and.btest(ito,i)).or.              2d7s23
     $                (btest(jto,i).and..not.btest(itc,i)))then         2d7s23
                     nab4(2,ioxx(2))=i                                  2d7s23
                     ioxx(2)=ioxx(2)+1                                  2d7s23
                    else                                                2d7s23
                     nab4(1,ioxx(1))=i                                  2d7s23
                     ioxx(1)=ioxx(1)+1                                  2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                 else if(ndifb.eq.3)then                                2d7s23
                  nnot=3                                                2d7s23
                  ioxx(1)=1                                             2d7s23
                  ioxx(2)=1                                             2d7s23
                  iswap=0                                               2d7s23
                  do i=1,norbxxx                                        2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if(btest(gandcc,i).and.                             2d7s23
     $        ((btest(itc,i).and..not.btest(jto,i)).or.                 2d7s23
     $         (btest(jtc,i).and..not.btest(ito,i))))then               2d7s23
                     if(btest(i1c,i))iswap=1                            2d7s23
                     nab4(1,1)=i                                        2d7s23
                     nab4(1,2)=i                                        2d7s23
                    else                                                2d7s23
                     nab4(2,ioxx(2))=i                                  2d7s23
                     ioxx(2)=ioxx(2)+1                                  2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                  if(iswap.ne.0)then                                    2d7s23
                   icpy=nab4(1,1)                                       2d7s23
                   nab4(1,1)=nab4(2,1)                                  2d7s23
                   nab4(2,1)=icpy                                       2d7s23
                   icpy=nab4(1,2)                                       2d7s23
                   nab4(1,2)=nab4(2,2)                                  2d7s23
                   nab4(2,2)=icpy                                       2d7s23
                   nbt=0                                                2d7s23
                   if(btest(itc,nab4(1,2)).and.                         2d7s23
     $                  .not.btest(itc,nab4(1,1)))nbt=1                 2d7s23
                  else                                                  2d7s23
                   nbt=0                                                2d7s23
                   if(btest(i1c,nab4(2,2)).and.                         2d7s23
     $                  .not.btest(i1c,nab4(2,1)))nbt=1                 2d7s23
                  end if                                                2d7s23
                  if(nbt.ne.0)then                                      2d7s23
                   nab4(1,1)=nab4(1,2)                                  2d7s23
                   nab4(2,1)=nab4(2,2)                                  2d7s23
                  end if                                                2d7s23
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                 2d7s23
                  nnot=3                                                2d7s23
                  do i=1,norbxxx                                        2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if(btest(itc,i))then                                2d7s23
                     nab4(1,1)=i                                        2d7s23
                     nab4(1,2)=i                                        2d7s23
                    else                                                2d7s23
                     nab4(2,1)=i                                        2d7s23
                     nab4(2,2)=i                                        2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                 end if                                                 2d7s23
                 ipssx=0                                                2d7s23
                 if(nnot.eq.3)then                                          12d8s20
                  ipssx=1                                                   12d8s20
                 else if(nnot.eq.4)then                                 2d7s23
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
                   write(6,*)
     $                   ('bit not set for nab4(1,1) ='),nab4(1,iu1)       11d27s20
                   stop 'nab4(1,1)b'                                           11d27s20
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
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                   4d20s21
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
                     ia=nab1(1)                                          1d8s21
                     ib=nab2(2)                                          1d8s21
                     lsa=ism(ia)                                           1d8s21
                     lsb=ism(ib)                                           1d8s21
                     lga=irel(ia)-1                                        1d8s21
                     lgb=irel(ib)-1                                        1d8s21
                     icol=lga+irefo(lsa)*lgb                              1d8s21
                    else                                                  1d8s21
                     ia=nab1(2)                                         8d10s21
                     ib=nab1(1)                                         8d10s21
                     lsa=ism(ia)                                           1d8s21
                     lsb=ism(ib)                                           1d8s21
                     lga=irel(ia)-1                                        1d8s21
                     lgb=irel(ib)-1                                        1d8s21
                     icol=lga+irefo(lsa)*lgb                             1d9s21
                    end if                                                1d8s21
                    lsab=multh(lsa,lsb)                                   1d9s21
                    jargo=1                                               2d23s21
                    itmpj=ibcoff                                          2d23s21
                    jtmpj=itmpj                                           2d23s21
                    do lj=1,4                                             2d23s21
                     if(nlj(lj).gt.0)then                                 2d23s21
                      jad1=jvcv+iff22b(jvcv+5+lj)                       8d26s21
                      jad2=jad1+nlj(lj)                                     3d19s21
                      ibcoff=jtmpj+ncsf(karg)*nlj(lj)                     2d23s21
                      call enough('hcddjkbk. 18',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(jarg),ncsfmid2,       2d23s21
     $                   nlj(lj),iwpb2,iwpk2,ff22b(jad2),ncsf2(lj,jarg),8d26s21
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
                      iad1=ivcv+iff22k(ivcv+5+li)                       8d26s21
                      iad2=iad1+nli(li)                                     3d19s21
                      ibcoff=jtmpi+ncsf(karg)*nli(li)                     2d23s21
                      itmpt=ibcoff                                        2d23s21
                      ibcoff=itmpt+ncsf(karg)*nli(li)                     2d23s21
                      call enough('hcddjkbk. 19',bc,ibc)
                      call xtimesn2(ncsf(karg),ncsf(iarg),ncsfmid1b,      2d23s21
     $                 nli(li),iwpb1b,iwpk1b,ff22k(iad2),ncsf2(li,iarg),8d26s21
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
                       jad1=jvcv+iff22b(jvcv+5+lj)                      8d26s21
                       jtmpi=itmpi                                       2d24s21
                       do li=1,4                                         2d24s21
                        if(nli(li).gt.0)then                             2d24s21
                         iad1=ivcv+iff22k(ivcv+5+li)                    8d26s21
                         jden=idenkn(li,lj,lsa)+nfdatk(2,li,isb)        8d26s21
     $                        *nfdatb(2,lj,jsb)*icol                     8d26s21
                         itmpm=ibcoff                                       2d23s21
                         ibcoff=itmpm+nlj(lj)*nli(li)                    2d24s21
                         call enough('hcddjkbk. 20',bc,ibc)
                         call dgemm('n','n',nli(li),nlj(lj),ncsf(karg),  2d24s21
     $                       1d0,bc(jtmpi),nli(li),bc(jtmpj),ncsf(karg),2d24s21
     $                       0d0,bc(itmpm),nli(li),                     2d24s21
     d' hcddjkbk. 12')
                         jtmpm=itmpm                                        2d23s21
                         do iii=0,nlj(lj)-1                                  2d23s21
                          jj=iff22b(jad1+iii)-1                         8d26s21
                          jjden=jden+nfdatk(2,li,isb)*jj-1              8d26s21
                          do j=0,nli(li)-1                               2d24s21
                           kk=iff22k(iad1+j)                            8d26s21
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
                       jad1=jvcv+iff22b(jvcv+5+l)                       8d26s21
                       iad1=ivcv+iff22k(ivcv+5+l)                       8d26s21
                       jden=idenjn(l,lsa)                               8d26s21
     $                      +nfdatk(2,l,isb)*nfdatb(2,l,jsb)*icol       8d26s21
                       itmpm=ibcoff                                       2d23s21
                       ibcoff=itmpm+nlj(l)*nli(l)                         2d23s21
                       call enough('hcddjkbk. 21',bc,ibc)
                       call dgemm('n','n',nli(l),nlj(l),ncsf(karg),1d0,   2d23s21
     $                    bc(jtmpi),nli(l),bc(jtmpj),ncsf(karg),0d0,    2d23s21
     $                    bc(itmpm),nli(l),                             2d23s21
     d' hcddjkbk. 13')
                       jtmpm=itmpm                                        2d23s21
                       do iii=0,nlj(l)-1                                  2d23s21
                        jj=iff22b(jad1+iii)-1                           8d26s21
                        jjden=jden+nfdatk(2,l,isb)*jj-1                 8d26s21
                        do j=0,nli(l)-1                                   2d23s21
                         kk=iff22k(iad1+j)                              8d26s21
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
                  end if                                                4d20s21
                 end do                                                   12d18s20
    4            continue                                                 12d18s20
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
        do l=1,4                                                        1d8s21
         if(nfdatk(2,l,isb).gt.0)then                                   8d26s21
          do ll=1,4                                                     1d8s21
           if(nfdatb(2,ll,jsb).gt.0)then                                 1d8s21
            nll=nfdatk(2,l,isb)*nfdatb(2,ll,jsb)                        8d26s21
            if(l.eq.ll.and.nll.gt.0)then                                2d24s21
             iden1ef(l)=iden1en(l)                                      2d25s21
             idenhvvf(l)=idenhvvn(l)                                    2d25s21
            end if                                                      1d11s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(l.eq.ll.and.                                            8d10s21
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
         if(nfdatb(2,ll,jsb).gt.0)then                                  8d26s21
          do l=1,4                                                      1d8s21
           if(nfdatk(2,l,isb).gt.0)then                                 8d26s21
            if(l.eq.ll)then                                             8d21s21
             sz=0d0                                                     8d21s21
             do i=0,nfdatk(2,l,isb)*nfdatb(2,ll,jsb)-1                  8d26s21
              sz=sz+bc(iden1ef(l)+i)**2                                 8d21s21
              sz=sz+bc(idenhvvf(l)+i)**2                                8d21s21
             end do                                                     8d21s21
             sz=sqrt(sz/dfloat(nfdatk(2,l,isb)*nfdatb(2,ll,jsb)))       8d26s21
             if(sz.gt.1d-10)nall=nall+1                                 8d21s21
            end if                                                      8d21s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             nokk(l,ll,lsa)=0                                           1d15s21
             if(min(irefo(lsa),irefo(lsb)).gt.0)then                    1d13s21
              nnn=nfdatk(2,l,isb)*nfdatb(2,ll,jsb)                      8d26s21
              if(l.eq.ll)then                                           8d10s21
                nn=irefo(lsa)*irefo(lsb)                                 1d8s21
               jden=idenjf(l,lsa)                                        1d11s21
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
              kden=idenkf(l,ll,lsa)                                     1d12s21
              do icol=0,nn-1                                            1d9s21
               ib=icol/irefo(lsa)
               ia=icol-irefo(lsa)*ib
               icolt=ib+irefo(lsb)*ia                                   1d10s21
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
            ioffp=ioff+nvirt(isbv1)*nrootu*nfdatk(2,1,isb)              8d26s21
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
              joffp=joff+nvirt(jsbv1)*nrootu*nfdatb(2,1,jsb)            8d26s21
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
                if(min(kase,nfdatk(2,li,isb),n2e).gt.0)then             8d26s21
                 jjoffp=joffp                                            1d14s21
                 tfj=1d0                                                 1d14s21
                 do lj=1,4                                               1d14s21
                  if(nfdatb(2,lj,jsb).gt.0)then                         8d26s21
                   tf=tfi*tfj
                   call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0)then                 8d11s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdatk(2,li,isb)*nfdatb(2,lj,jsb)                8d26s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    call enough('hcddjkbk. 22',bc,ibc)
                    fact=0d0                                                1d13s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then             1d14s21
                      if(nokk(li,lj,lsa).gt.0)then                              1d12s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjkbk. 23',bc,ibc)
                       do iz=itmpi,ibcoff-1                             8d10s21
                        bc(iz)=0d0                                      8d10s21
                       end do                                           8d10s21
                       do ii2e=1,n2e                                    8d10s21
                        if(multh(lsb,isbr).eq.isymop(i2eop(1,ii2e)).and.8d10s21
     $                     multh(lsa,isbl).eq.isymop(i2eop(2,ii2e)))then8d10s21
                         do j=0,nokk(li,lj,lsa)-1                                1d12s21
                          jj=ibc(ndenkf(li,lj,lsa)+j)                            1d12s21
                          jb=jj/irefo(lsa)                                   1d12s21
                          ja=jj-irefo(lsa)*jb                                1d12s21
                          jb=jb+idoubo(lsb)                             8d10s21
                          ja=ja+idoubo(lsa)                             8d10s21
                          i10=i1s                                       8d10s21
                          i1n=nvirt(isbl)                               8d10s21
                          jtmpi=itmpi+j                                 8d10s21
                          irb=ixmt(isbr,i2eop(1,ii2e))+nbasdws(isbr)*jb 8d10s21
                          ial=ixmt(lsa,i2eop(2,ii2e))+ja                8d10s21
                          do i2=i2s,i2e                                 8d10s21
                           ivr=i2-1+idoubo(isbr)+irefo(isbr)            8d10s21
                           if(i2.eq.i2e)i1n=i1e                         8d10s21
                           do i1=i10,i1n                                8d10s21
                            ivl=i1-1+idoubo(isbl)+irefo(isbl)           8d10s21
                            bc(jtmpi)=bc(jtmpi)+phase2*bc(irb+ivr)      8d10s21
     $                           *bc(ial+nbasdws(lsa)*ivl)              8d10s21
                            jtmpi=jtmpi+nokk(li,lj,lsa)                 8d10s21
                           end do                                       8d10s21
                           i10=1                                        8d10s21
                          end do                                        8d10s21
                         end do                                         8d10s21
                        end if                                          8d10s21
                       end do                                           8d10s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0,   1d14s21
     $               bc(idenkf(li,lj,lsa)),mn,bc(itmpi),nokk(li,lj,lsa),1d14s21
     $                    fact,bc(intden),mn,                           1d14s21
     d' hcddjkbk. 14')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
                     nrow=nfdatb(2,lj,jsb)*nvirt(isbl)                  8d26s21
                     nmul=nfdatk(2,li,isb)*nhere2                       8d26s21
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjkbk. 24',bc,ibc)
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
                       do j=0,nfdatb(2,lj,jsb)-1                        8d26s21
                        do i=0,nfdatk(2,li,isb)-1                       8d26s21
                         iad=itmp1+j+nfdatb(2,lj,jsb)*(i1m+nvirt(isbl)*(8d26s21
     $                       i2n+nhere2*i))                             1d14s21
                         bc(iad)=bc(jntden+i)                               1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdatk(2,li,isb)                  8d26s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjkbk. 25',bc,ibc)
                     nx=nhere2*nfdatk(2,li,isb)*nrootu                  8d26s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibr.eq.1)then                                       1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do i=0,nfdatk(2,li,isb)-1                        8d26s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,li,isb)*ir) 8d26s21
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
                       do i=0,nfdatk(2,li,isb)-1                        8d26s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,li,isb)*ir) 8d26s21
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
                      iioff=iioffp-nvirt(isbv1)*nrootu*nfdatk(2,1,isb)  8d26s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do i=0,nfdatk(2,1,isb)-1                         8d26s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,1,isb)*(ir  8d26s21
     $                       +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjkbk. 26',bc,ibc)
                     tff=tf*phs                                          1d14s21
                     call dgemm('n','n',nrow,ncol,nmul,tff,              1d14s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjkbk. 15')
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
                         jprod=iprod+nfdatb(2,lj,jsb)*(iv1              8d26s21
     $                        +nvirt(jsbv1)*(ir+nrootu*iv2))            8d26s21
                         do j=0,nfdatb(2,lj,jsb)-1                      8d26s21
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
                         jprod=iprod+nfdatb(2,lj,jsb)*(iv2              8d26s21
     $                        +nvirt(jsbv2)*(ir+nrootu*iv1))            8d26s21
                         do j=0,nfdatb(2,lj,jsb)-1                      8d26s21
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
                        jprod=iprod+nfdatb(2,lj,jsb)*(iv2+nvirt(jsbv1)  8d26s21
     $                      *(ir+nrootu*iv2))                           1d14s21
                        do j=0,nfdatb(2,lj,jsb)-1                       8d26s21
                         gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh       1d14s21
                        end do                                           1d14s21
                       end do                                            1d14s21
                      end do                                             1d14s21
                     end if                                              1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                1d14s21
                   jjoffp=jjoffp+nnvv*nfdatb(2,lj,jsb)                  8d26s21
                  end if                                                  1d14s21
                  if(jtf.ne.0)tfj=-1d0                                   1d15s21
                 end do                                                  1d14s21
                end if                                                   1d14s21
                iioffp=iioffp+mmvv*nfdatk(2,li,isb)                     8d26s21
                if(kase.ne.0)then                                       1d13s23
                 if(itf.ne.0)tfi=-1d0                                     1d14s21
                end if                                                  1d13s23
               end do                                                    1d14s21
               iioffp=ioffp                                               1d12s21
               jjoffp=joffp                                              1d13s21
               tf=1d0                                                     1d13s21
               do l=1,4                                                   1d12s21
                if(min(kase,nfdatk(2,l,isb),nfdatb(2,l,jsb)).gt.0)then  8d26s21
                 mn=nfdatk(2,l,isb)*nfdatb(2,l,jsb)                     8d26s21
                 call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                 nhere=ih+1-il                                            1d12s21
                 nhere2=i2e+1-i2s                                         1d12s21
                 if(min(nhere,nvirt(isbvc),n2e).gt.0)then               8d10s21
                  intden=ibcoff                                            1d12s21
                  ibcoff=intden+mn*nhere                                   1d12s21
                  call enough('hcddjkbk. 27',bc,ibc)
                  fact=0d0                                                1d13s21
                  do lsa=1,nsymb                                              1d8s21
                   lsb=multh(lsa,ijsb)                                        1d9s21
                   if(min(irefo(lsa),irefo(lsb)).gt.0)then              8d10s21
                    if(nokj(l,lsa).gt.0)then                              1d12s21
                     itmpi=ibcoff                                         1d12s21
                     ibcoff=itmpi+nokj(l,lsa)*nhere                       1d12s21
                     call enough('hcddjkbk. 28',bc,ibc)
                     do iz=itmpi,ibcoff-1                               8d10s21
                      bc(iz)=0d0                                        8d10s21
                     end do                                             8d10s21
                     do ii2e=1,n2e                                      8d10s21
                      if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)).and.   8d10s21
     $                   multh(isbl,isbr).eq.isymop(i2eop(2,ii2e)))then 8d10s21
                       do j=0,nokj(l,lsa)-1                                1d12s21
                        jj=ibc(ndenjf(l,lsa)+j)                            1d12s21
                        jb=jj/irefo(lsa)                                   1d12s21
                        ja=jj-irefo(lsa)*jb                                1d12s21
                        jb=jb+idoubo(lsb)                               8d10s21
                        ja=ja+idoubo(lsa)                               8d10s21
                        iab=ixmt(lsa,i2eop(1,ii2e))+ja+nbasdws(lsa)*jb  8d10s21
                        jtmpi=itmpi+j                                   8d10s21
                        i10=i1s                                         8d10s21
                        i1n=nvirt(isbl)                                 8d10s21
                        fab=phase2*bc(iab)                              8d10s21
                        do i2=i2s,i2e                                   8d10s21
                         ivr=i2-1+idoubo(isbr)+irefo(isbr)              8d10s21
                         if(i2.eq.i2e)i1n=i1e                           8d10s21
                         do i1=i10,i1n                                  8d10s21
                          ivl=i1-1+idoubo(isbl)+irefo(isbl)             8d10s21
                          ilr=ixmt(isbl,i2eop(2,ii2e))+ivl              8d10s21
     $                         +nbasdws(isbl)*ivr                       8d10s21
                          bc(jtmpi)=bc(jtmpi)+fab*bc(ilr)               8d10s21
                          jtmpi=jtmpi+nokj(l,lsa)                       8d10s21
                         end do                                         8d10s21
                         i10=1                                          8d10s21
                        end do                                          8d10s21
                       end do                                           8d10s21
                      end if                                            8d10s21
                     end do                                             8d10s21
                     call dgemm('n','n',mn,nhere,nokj(l,lsa),1d0,         1d12s21
     $                  bc(idenjf(l,lsa)),mn,bc(itmpi),nokj(l,lsa),fact,1d12s21
     $                  bc(intden),mn,                                  1d12s21
     d' hcddjkbk. 16')
                     fact=1d0                                             1d12s21
                     ibcoff=itmpi                                         1d13s21
                    end if                                                1d12s21
                   end if                                                 1d13s21
                  end do                                                  1d13s21
                  if(fact.gt.0.5d0)then                                    1d12s21
                   nrow=nfdatb(2,l,jsb)*nvirt(isbl)                     8d26s21
                   nmul=nfdatk(2,l,isb)*nhere2                          8d26s21
                   itmp1=ibcoff                                            1d12s21
                   ibcoff=itmp1+nrow*nmul                                 1d12s21
                   call enough('hcddjkbk. 29',bc,ibc)
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
                     do j=0,nfdatb(2,l,jsb)-1                           8d26s21
                      do i=0,nfdatk(2,l,isb)-1                          8d26s21
                       iad=itmp1+j+nfdatb(2,l,jsb)*(i1m+nvirt(isbl)*(i2n8d26s21
     $                  +nhere2*i))                                     1d12s21
                       bc(iad)=bc(jntden+i)                               1d13s21
                      end do                                              1d12s21
                      jntden=jntden+nfdatk(2,l,isb)                     8d26s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                    i10=1                                                 1d12s21
                   end do                                                 1d12s21
                   ncol=nvirt(isbvc)*nrootu                               1d12s21
                   itmp2=ibcoff                                           1d12s21
                   ibcoff=itmp2+nmul*ncol                                 1d12s21
                   call enough('hcddjkbk. 30',bc,ibc)
                   nx=nhere2*nfdatk(2,l,isb)*nrootu                     8d26s21
                   do i=itmp2,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
                   if(ibr.eq.1)then                                       1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv2=i2-1                                              1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                     do i=0,nfdatk(2,l,isb)-1                           8d26s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,l,isb)*ir)    8d26s21
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
                     do i=0,nfdatk(2,l,isb)-1                           8d26s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,l,isb)*ir)    8d26s21
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
                    iioff=iioffp-nvirt(isbv1)*nrootu*nfdatk(2,1,isb)    8d26s21
                    do i2=i2s,i2e                                        1d15s21
                     iv1=i2-1                                            1d15s21
                     i2n=i2-i2s                                          1d15s21
                     do i=0,nfdatk(2,1,isb)-1                           8d26s21
                      do ir=0,nrootm                                     1d15s21
                       iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)         1d15s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,1,isb)*(ir    8d26s21
     $                     +nrootu*iv1))                                1d15s21
                       bc(jtmp2)=vd(iuse)*srh                            1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                   iprod=ibcoff                                           1d12s21
                   ibcoff=iprod+nrow*ncol                                 1d12s21
                   call enough('hcddjkbk. 31',bc,ibc)
                   call dgemm('n','n',nrow,ncol,nmul,tf,                  1d13s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjkbk. 17')
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
                       jprod=iprod+nfdatb(2,l,jsb)*(iv1+nvirt(jsbv1)*(ir8d26s21
     $                     +nrootu*iv2))                                  1d12s21
                       do j=0,nfdatb(2,l,jsb)-1                         8d26s21
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
                       jprod=iprod+nfdatb(2,l,jsb)*(iv2+nvirt(jsbv2)*(ir8d26s21
     $                     +nrootu*iv1))                                1d13s21
                       do j=0,nfdatb(2,l,jsb)-1                         8d26s21
                        gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(l.eq.1.and.jsbv12.eq.1)then                        1d15s21
                    jjoff=jjoffp-nvirt(isbvc)*nrootu*nfdatb(2,1,jsb)    8d26s21
                    nv=nrootu*nvirt(isbvc)                               1d15s21
                    do iv1=0,nvirt(isbvc)-1                              1d15s21
                     do ir=0,nrootm                                      1d15s21
                      iad=jjoff+iv1+nvirt(isbvc)*ir                      1d15s21
                      jprod=iprod+nfdatb(2,l,jsb)*(iv1+nvirt(jsbv2)*(ir 8d26s21
     $                    +nrootu*iv1))                                 1d15s21
                      do j=0,nfdatb(2,l,jsb)-1                          8d26s21
                       gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh         1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                  end if                                                   1d12s21
                  ibcoff=intden                                            1d12s21
                 end if                                                   1d13s21
                end if                                                      1d13s21
                if(kase.ne.0)then                                       1d13s23
                 if(mod(jtf+itf,2).ne.0)then                              1d15s21
                  tf=-1d0                                                 1d14s21
                 end if                                                   1d14s21
                end if                                                  1d13s23
                iioffp=iioffp+mmvv*nfdatk(2,l,isb)                      8d26s21
                jjoffp=jjoffp+nnvv*nfdatb(2,l,jsb)                      8d26s21
               end do                                                      1d13s21
              end do                                                     1d13s21
             end if                                                      1d13s21
             if(jsbv12.eq.1)joff=joff                                    1d12s21
     $           +nvirt(jsbv1)*nrootu*nfdatb(2,1,jsb)                   8d26s21
             do l=1,4                                                    1d12s21
              joff=joff+nnvv*nfdatb(2,l,jsb)                            8d26s21
             end do                                                      1d12s21
            end if                                                       1d12s21
           end do                                                        1d12s21
           if(isbv12.eq.1)ioff=ioff+nvirt(isbv1)*nrootu*nfdatk(2,1,isb) 8d26s21
           do l=1,4                                                      1d12s21
            ioff=ioff+mmvv*nfdatk(2,l,isb)                              8d26s21
           end do                                                        1d12s21
          end if                                                         1d12s21
         end do                                                          1d12s21
c     we are here if isbv1=jsbv1 and isbv2=jsbv2. This occurs when
c     isbv12=jsbv12=isb*isymket=jsb*isymbra.
c     ijsb=isb*jsb
         lprt=nlprt.ne.0
         if(isbv12.eq.jsbv12)then                                       8d11s21
          ioff=ioffdnon                                                  1d11s21
          joff=joffdnon                                                 8d11s21
          do isbv1=1,nsymb                                               1d11s21
           isbv2=multh(isbv1,isbv12)                                     1d11s21
           if(isbv2.ge.isbv1)then                                        1d11s21
            nvv=0                                                        4d28s21
            if(isbv1.eq.isbv2)then                                       1d11s21
             ioffs=ioff                                                  1d12s21
             joffs=joff                                                 8d11s21
             nnn=nvirt(isbv1)*nrootu                                     1d11s21
             call ilimts(nvirt(isbv1),nrootu,mynprocg,mynowprog,il,ih,   1d11s21
     $           i1s,i1e,i2s,i2e)                                       1d11s21
             nhere=ih+1-il                                               1d11s21
             if(min(nfdatb(2,1,jsb),nfdatk(2,1,isb),nhere).gt.0)then    8d26s21
              itmp=ibcoff                                                1d11s21
              ibcoff=itmp+nhere*nfdatb(2,1,jsb)                         8d26s21
              call enough('hcddjkbk. 32',bc,ibc)
              call dgemm('n','n',nhere,nfdatb(2,1,jsb),nfdatk(2,1,isb), 8d26s21
     $            1d0,vd(ioff+il-1),nnn,bc(iden1ef(1)),nfdatk(2,1,isb), 8d26s21
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjkbk. 18')
              do i=0,nfdatb(2,1,jsb)-1                                   8d11s21
               i10=i1s                                                   1d11s21
               i1n=nvirt(isbv1)                                          1d11s21
               jtmp=itmp+nhere*i                                         1d11s21
               do i2=i2s,i2e                                             1d11s21
                if(i2.eq.i2e)i1n=i1e                                     1d11s21
                iad=joff-1+nvirt(isbv1)*(i2-1+nrootu*i)                  1d11s21
                do i1=i10,i1n                                            1d11s21
                 gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                 jtmp=jtmp+1                                             1d11s21
                end do                                                   1d11s21
                i10=1                                                    1d11s21
               end do                                                    1d11s21
              end do                                                     1d11s21
              ibcoff=itmp                                                1d11s21
             end if                                                      1d11s21
             ioff=ioff+nnn*nfdatk(2,1,isb)                              8d26s21
             joff=joff+nnn*nfdatb(2,1,jsb)                              8d26s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       1d11s21
            else                                                         1d11s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               1d11s21
            end if                                                       1d11s21
            nnn=nvv*nrootu                                               1d11s21
            tf=1d0                                                       1d11s21
            do l=1,4                                                     1d11s21
             if(min(nfdatb(2,l,jsb),nfdatk(2,l,isb)).gt.0)then          8d26s21
              call ilimts(nvv,nrootu,mynprocg,mynowprog,il,ih,            1d11s21
     $            i1s,i1e,i2s,i2e)                                       1d11s21
              nhere=ih+1-il                                               1d11s21
              if(nhere.gt.0)then                                          1d11s21
               itmp=ibcoff                                                1d11s21
               ibcoff=itmp+nhere*nfdatb(2,l,jsb)                        8d26s21
               call enough('hcddjkbk. 33',bc,ibc)
               call dgemm('n','n',nhere,nfdatb(2,l,jsb),nfdatk(2,l,isb),8d26s21
     $            1d0,vd(ioff+il-1),nnn,bc(iden1ef(l)),nfdatk(2,l,isb), 8d26s21
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjkbk. 19')
               do i=0,nfdatb(2,l,jsb)-1                                 8d26s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                  1d11s21
                jtmp=itmp+nhere*i                                         1d11s21
                do i2=i2s,i2e                                             1d11s21
                 if(i2.eq.i2e)i1n=i1e                                     1d11s21
                 iad=joff-1+nvv*(i2-1+nrootu*i)                         8d11s21
                 do i1=i10,i1n                                            1d11s21
                  gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                  jtmp=jtmp+1                                             1d11s21
                 end do                                                   1d11s21
                 i10=1                                                    1d11s21
                end do                                                    1d11s21
               end do                                                     1d11s21
               ibcoff=itmp                                               1d12s21
              end if                                                     1d11s21
             end if                                                      1d11s21
             ioff=ioff+nnn*nfdatk(2,l,isb)                              8d26s21
             joff=joff+nnn*nfdatb(2,l,jsb)                              8d26s21
             tf=-1d0                                                     1d11s21
            end do                                                       1d11s21
           end if                                                        1d11s21
          end do                                                         1d11s21
         end if                                                          1d11s21
c
c     for xmtf, we need isbvn*jsbvm=isymop(lxmt) and
c     isbvn'=jsbvm', where n,m are 1 or 2 and n'=3-n and m'=3-m.        8d11s21
c     since isbvn'=jsbvm', then isbvn'*jsbvm'=1 and                     8d11s21
c     isbvn*jsbvm=isbv12*jsbv12=isymop(lxmt).                           8d11s21
c
c     densities were computed with order bra v v" ket v v'
c
         if(multh(isbv12,jsbv12).eq.isymop(lxmt))then                   8d11s21
          ioff=ioffdnon                                                  1d11s21
          do isbv1=1,nsymb                                               1d11s21
           isbv2=multh(isbv1,isbv12)                                     1d11s21
           if(isbv2.ge.isbv1)then                                        1d11s21
            ioffs=ioff                                                  1d12s21
            ibctop=ibcoff                                               8d11s21
            if(isbv12.eq.1)then                                         8d11s21
             ioff=ioff+nrootu*nfdatk(2,1,isb)*nvirt(isbv1)              8d26s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      8d11s21
             isw=0                                                      8d11s21
             itmpvisv=ibcoff                                            8d11s21
             if(min(nvirt(isbv1),nfdatb(2,1,jsb),nfdatk(2,1,isb))       11d16s22
     $            .gt.0)then                                            11d16s22
              nnn=nvirt(isbv1)*nrootu                                    8d11s21
              ibcoff=itmpvisv+nnn*nfdatb(2,1,jsb)                        8d26s21
              call enough('hcddjkbk. 34',bc,ibc)
              fff=2d0*phase1                                             8d11s21
              call dgemm('n','n',nnn,nfdatb(2,1,jsb),nfdatk(2,1,isb),    8d26s21
     $              fff,vd(ioffs),nnn,bc(idenhvvf(1)),nfdatk(2,1,isb),  8d26s21
     $              0d0,bc(itmpvisv),nnn,                               8d11s21
     d' hcddjkbk. 20')
             end if                                                     8d27s21
            else                                                        8d11s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                              8d11s21
             isw=1                                                      8d11s21
            end if                                                      8d11s21
            call ilimts(nvv,nrootu,mynprocg,mynowprog,il,ih,i1s,i1e,    8d11s21
     $              i2s,i2e)                                            8d11s21
            nhere=ih+1-il                                               8d11s21
            if(nhere.gt.0)then                                          8d11s21
             nnn=nvv*nrootu                                             8d11s21
             do l=1,4                                                   8d11s21
              itmpl(l)=ibcoff                                           8d11s21
              ibcoff=itmpl(l)+nhere*nfdatb(2,l,jsb)                      8d11s21
              call enough('hcddjkbk. 35',bc,ibc)
              if(min(nfdatb(2,l,jsb),nfdatk(2,l,isb)).gt.0)then         8d26s21
               call dgemm('n','n',nhere,nfdatb(2,l,jsb),                8d26s21
     $                nfdatk(2,l,isb),phase1,vd(ioff+il-1),nnn,         8d26s21
     $                bc(idenhvvf(l)),nfdatk(2,l,isb),0d0,bc(itmpl(l))  8d26s21
     $                ,nhere,                                           8d11s21
     d' hcddjkbk. 21')
              end if                                                    8d11s21
              ioff=ioff+nfdatk(2,l,isb)*nvv*nrootu                      8d26s21
             end do                                                     8d11s21
            else                                                        12d12s22
             do l=1,4                                                   12d12s22
              ioff=ioff+nfdatk(2,l,isb)*nvv*nrootu                      8d26s21
             end do                                                     12d12s22
            end if                                                      8d11s21
            joff=joffdnon                                               8d11s21
            do jsbv1=1,nsymb                                            8d11s21
             jsbv2=multh(jsbv1,jsbv12)                                  8d11s21
             if(jsbv2.ge.jsbv1)then                                     8d11s21
c
c     densities were computed with order bra v v" ket v v'
c
c     cases                                        tf
c     Gvv'rj= dji Xv'v" Vvv"ri    jsbv1=isbv1      no
c     Gvv'rj= dji Xv'v" Vv"vri    jsbv1=isbv2      yes
c     Gvv'rj= dji Xvv"  Vv"v'ri   jsbv2=isbv2      no
c     Gvv'rj= dji Xvv"  Vv'v"ri   jsbv2=isbv1      yes
c
c     ket visv:
c     Gvv'rj= dji Xv'v Vvvri
c     Gvv'rj= dji Xv'v Vvvri        same as above
c     Gvv'rj= dji Xvv' Vv'v'ri
c     Gvv'rj= dji Xvv' Vv'v'ri      same as above
c
c     bra visv:
c     Gvvrj= dji Xvv" Vvv"ri
c     Gvvrj= dji Xvv" Vv"vri
c     Gvvrj= dji Xvv" Vv"vri        same as above
c     Gvvrj= dji Xvv" Vvv"ri        same as top
c
              joffs=joff                                                8d11s21
              if(jsbv12.eq.1)then                                       8d11s21
               joff=joff+nfdatb(2,1,jsb)*nrootu*nvirt(jsbv1)             8d11s21
               mvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                    8d11s21
               jsw=0                                                    8d11s21
              else                                                      8d11s21
               mvv=nvirt(jsbv1)*nvirt(jsbv2)                            8d11s21
               jsw=1                                                    8d11s21
              end if                                                    8d11s21
              if(isymop(lxmt).eq.1.and.jsbv12.eq.1.and.jsbv1.eq.isbv1)  8d11s21
     $             then                                                 8d11s21
               do i=0,nfdatb(2,1,jsb)-1                                     1d18s21
                do ir=0,nrootm                                             1d18s21
                 do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   1d18s21
                  ivp=iv+idoubo(isbv1)+irefo(isbv1)                      8d10s21
                  iad=joffs+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                  ih0=ixmtf(isbv1)+ivp*(nbasdws(isbv1)+1)                8d10s21
                  jtmp=itmpvisv+iv+nvirt(isbv1)*(ir+nrootu*i)           8d11s21
                  gd(iad)=gd(iad)+bc(ih0)*bc(jtmp)                      8d11s21
                 end do                                                    1d18s21
                end do                                                     1d18s21
               end do                                                      1d18s21
              end if                                                    8d11s21
              nnn=nvv*nrootu                                               1d11s21
              tf=1d0                                                       1d11s21
              do l=1,4                                                     1d11s21
               if(min(nfdatb(2,l,jsb),nfdatk(2,l,isb)).gt.0)then          8d11s21
                if(l.eq.1)then                                          8d11s21
                 if(isbv12.eq.1.and.(jsbv1.eq.isbv1.or.jsbv2.eq.isbv1)) 8d11s21
     $               then                                               8d11s21
c     ket visv:
c     Gvv'rj= dji Xv'v Vvvri
c     Gvv'rj= dji Xvv' Vv'v'ri
                  do i=0,nfdatb(2,1,jsb)-1                              8d26s21
                   do ir=0,nrootm                                             1d18s21
                    iad=joff+mvv*(ir+nrootu*i)                          8d11s21
                    do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   1d18s21
                     ivp=iv+1                                           8d11s21
                     iadv=itmpvisv+iv+nvirt(isbv1)*(ir+nrootu*i)        8d11s21
                     ivpp=iv+idoubo(isbv1)+irefo(isbv1)                 8d11s21
                     ix1=ixmtf(jsbv1)+idoubo(jsbv1)+irefo(jsbv1)         8d11s21
     $                   +nbasdws(jsbv1)*ivpp                           8d11s21
                     ix2=ixmtf(jsbv2)+idoubo(jsbv2)+irefo(jsbv2)         8d11s21
     $                   +nbasdws(jsbv2)*ivpp                           8d11s21
                     if(jsbv1.eq.isbv1)then                             8d11s21
c     Gvv'rj= dji Xv'v Vvvri
                      ibot=ivp*(1-jsw)                                   8d11s21
                      do iv2=ibot,nvirt(jsbv2)-1                         8d11s21
                       irec=iv+nvirt(jsbv1)*iv2                         8d11s21
                       itri=((iv2*(iv2-1))/2)+iv                        8d11s21
                       irow=itri+jsw*(irec-itri)                        8d11s21
                       gd(iad+irow)=gd(iad+irow)+bc(ix2+iv2)*bc(iadv)   8d12s21
     $                      *srh                                        8d12s21
                      end do                                             8d11s21
                     end if                                             8d11s21
                     if(jsbv2.eq.isbv1)then                             8d11s21
c     Gvv'rj= dji Xvv' Vv'v'ri
                      itop=(iv+jsw*(nvirt(jsbv1)-iv))-1                  8d11s21
                      do iv1=0,itop                                      8d11s21
                       irec=iv1+nvirt(jsbv1)*iv                         8d11s21
                       itri=((iv*(iv-1))/2)+iv1                         8d11s21
                       irow=itri+jsw*(irec-itri)                        8d11s21
                       gd(iad+irow)=gd(iad+irow)+bc(ix1+iv1)*bc(iadv)   8d12s21
     $                      *srh                                        8d12s21
                      end do                                             8d11s21
                     end if                                             8d11s21
                    end do                                              8d11s21
                   end do                                               8d11s21
                  end do                                                8d11s21
                 end if                                                 8d11s21
                 if(nhere.gt.0.and.jsbv12.eq.1.and.                     8d11s21
     $                (jsbv1.eq.isbv1.or.jsbv1.eq.isbv2))then           8d11s21
c     bra visv:
c     Gvvrj= dji Xvv" Vvv"ri
c     Gvvrj= dji Xvv" Vv"vri
                  do i=0,nfdatb(2,1,jsb)-1                              8d26s21
                   i10=i1s                                              8d11s21
                   i1n=nvv                                              8d11s21
                   jtmp=itmpl(l)+nhere*i                                8d11s21
                   do i2=i2s,i2e                                        8d11s21
                    if(i2.eq.i2e)i1n=i1e                                8d11s21
                    ir=i2-1
                    iad=joffs+nvirt(jsbv1)*(ir+nrootu*i)                8d11s21
                    if(jsbv1.eq.isbv1)then                              8d11s21
                     ltmp=jtmp                                          8d12s21
c     Gvvrj= dji Xvv" Vvv"ri
                     do i1=i10,i1n                                       8d11s21
                      i1m=i1-1                                           8d11s21
                      iv2=i1m/nvirt(isbv1)                              8d21s21
                      iv1=i1m-nvirt(isbv1)*iv2                          8d21s21
                      try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                      iv2t=sqrt(try)+0.5d0                                   1d11s21
                      iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                      iv1=iv1t+isw*(iv1-iv1t)                           8d21s21
                      iv2=iv2t+isw*(iv2-iv2t)                           8d21s21
                      ix=ixmtf(isbv1)+iv1+idoubo(isbv1)+irefo(isbv1)    8d11s21
     $                  +nbasdws(isbv1)*(iv2+idoubo(isbv2)+irefo(isbv2))8d11s21
                      gd(iad+iv1)=gd(iad+iv1)+bc(ix)*bc(ltmp)*sr2       8d11s21
                      ltmp=ltmp+1                                        8d11s21
                     end do                                              8d11s21
                    end if                                              8d11s21
                    if(jsbv1.eq.isbv2)then                              8d11s21
                     ltmp=jtmp                                          8d12s21
c     Gvvrj= dji Xvv" Vv"vri
                     do i1=i10,i1n                                       8d11s21
                      i1m=i1-1                                           8d11s21
                      iv2=i1m/nvirt(isbv1)                              8d21s21
                      iv1=i1m-nvirt(isbv1)*iv2                          8d21s21
                      try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                      iv2t=sqrt(try)+0.5d0                                   1d11s21
                      iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                      iv1=iv1t+isw*(iv1-iv1t)                           8d21s21
                      iv2=iv2t+isw*(iv2-iv2t)                           8d21s21
                      ix=ixmtf(isbv2)+iv2+idoubo(isbv2)+irefo(isbv2)    8d11s21
     $                  +nbasdws(isbv2)*(iv1+idoubo(isbv1)+irefo(isbv1))8d11s21
                      gd(iad+iv2)=gd(iad+iv2)+bc(ix)*bc(ltmp)*sr2       8d11s21
                      ltmp=ltmp+1                                        8d11s21
                     end do                                              8d11s21
                    end if                                              8d11s21
                    jtmp=jtmp+i1n+1-i10                                 2d11s22
                    i10=1                                               8d11s21
                   end do                                               8d11s21
                  end do                                                8d11s21
                 end if                                                 8d11s21
                end if                                                  8d11s21
                if(nhere.gt.0.and.(jsbv1.eq.isbv1.or.jsbv1.eq.isbv2.or. 8d11s21
     $               jsbv2.eq.isbv1.or.jsbv2.eq.isbv2))then             8d11s21
                 do i=0,nfdatb(2,l,jsb)-1                               8d26s21
                  jtmp=itmpl(l)+nhere*i                                 8d11s21
                  i10=i1s                                               8d11s21
                  i1n=nvv                                               8d11s21
                  do i2=i2s,i2e                                         8d11s21
                   if(i2.eq.i2e)i1n=i1e                                 8d11s21
                   ir=i2-1                                              8d11s21
                   iad=joff+mvv*(ir+nrootu*i)                           8d11s21
                   do i1=i10,i1n                                        8d11s21
                    i1m=i1-1                                               1d11s21
                    iv2=i1m/nvirt(isbv1)                                8d11s21
                    iv1=i1m-nvirt(isbv1)*iv2                            8d11s21
                    try=2d0*dfloat(i1m)+0.25d0                          8d11s21
                    iv2t=sqrt(try)+0.5d0                                8d11s21
                    iv1t=i1m-((iv2t*(iv2t-1))/2)                        8d11s21
                    iv1=iv1t+isw*(iv1-iv1t)                             8d11s21
                    iv2=iv2t+isw*(iv2-iv2t)                             8d11s21
                    iv1p=iv1+1                                          8d11s21
                    iv2p=iv2+1                                          8d11s21
                    iv1pp=iv1+idoubo(isbv1)+irefo(isbv1)                8d11s21
                    iv2pp=iv2+idoubo(isbv2)+irefo(isbv2)                8d11s21
                    fff0=bc(jtmp)                                       8d11s21
                    fff1=bc(jtmp)*tf                                    8d11s21
                    if(jsbv1.eq.isbv1)then                              8d11s21
c     Gvv'rj= dji Xv'v" Vvv"ri  no
                     ibot=iv1p*(1-jsw)                                  8d11s21
                     ix=ixmtf(jsbv2)+idoubo(jsbv2)+irefo(jsbv2)         8d11s21
     $                    +nbasdws(jsbv2)*iv2pp                         8d11s21
                     do iv=ibot,nvirt(jsbv2)-1                          8d11s21
                      irec=iv1+nvirt(jsbv1)*iv                          8d11s21
                      itri=((iv*(iv-1))/2)+iv1                          8d11s21
                      irow=itri+jsw*(irec-itri)                         8d11s21
                      gd(iad+irow)=gd(iad+irow)+bc(ix+iv)*fff1          8d12s21
                     end do                                             8d11s21
                    end if                                              8d11s21
                    if(jsbv1.eq.isbv2)then                              8d11s21
c     Gvv'rj= dji Xv'v" Vv"vri  yes
                     ibot=iv2p*(1-jsw)                                  8d11s21
                     ix=ixmtf(jsbv2)+idoubo(jsbv2)+irefo(jsbv2)         8d11s21
     $                    +nbasdws(jsbv2)*iv1pp                         8d11s21
                     do iv=ibot,nvirt(jsbv2)-1                          8d11s21
                      irec=iv2+nvirt(jsbv1)*iv                          8d11s21
                      itri=((iv*(iv-1))/2)+iv2                          8d11s21
                      irow=itri+jsw*(irec-itri)                         8d11s21
                      gd(iad+irow)=gd(iad+irow)+bc(ix+iv)*fff0          8d11s21
                     end do                                             8d11s21
                    end if                                              8d11s21
                    if(jsbv2.eq.isbv2)then                              8d11s21
c     Gvv'rj= dji Xvv"  Vv"v'ri no
                     itop=(iv2+jsw*(nvirt(jsbv1)-iv2))-1                8d11s21
                     ix=ixmtf(jsbv1)+idoubo(jsbv1)+irefo(jsbv1)         8d11s21
     $                    +nbasdws(jsbv1)*iv1pp                         8d11s21
                     do iv=0,itop                                       8d11s21
                      irec=iv+nvirt(jsbv1)*iv2                          8d11s21
                      itri=((iv2*(iv2-1))/2)+iv                         8d11s21
                      irow=itri+jsw*(irec-itri)                         8d11s21
                      gd(iad+irow)=gd(iad+irow)+bc(ix+iv)*fff1          8d11s21
                     end do                                             8d11s21
                    end if                                              8d11s21
                    if(jsbv2.eq.isbv1)then                              8d11s21
c     Gvv'rj= dji Xvv"  Vv'v"ri yes
                     itop=(iv1+jsw*(nvirt(jsbv1)-iv1))-1                8d11s21
                     ix=ixmtf(jsbv1)+idoubo(jsbv1)+irefo(jsbv1)         8d11s21
     $                    +nbasdws(jsbv1)*iv2pp                         8d11s21
                     do iv=0,itop                                       8d11s21
                      irec=iv+nvirt(jsbv1)*iv1                          8d11s21
                      itri=((iv1*(iv1-1))/2)+iv                         8d11s21
                      irow=itri+jsw*(irec-itri)                         8d11s21
                      gd(iad+irow)=gd(iad+irow)+bc(ix+iv)*fff0          8d11s21
                     end do                                             8d11s21
                    end if                                              8d11s21
                    jtmp=jtmp+1                                         8d11s21
                   end do                                               8d11s21
                   i10=1                                                8d11s21
                  end do                                                8d11s21
                 end do                                                 8d11s21
                end if                                                  8d11s21
               end if                                                   8d11s21
               joff=joff+mvv*nrootu*nfdatb(2,l,jsb)                     8d26s21
               tf=-1d0                                                  8d11s21
              end do                                                    8d11s21
             end if                                                     8d11s21
            end do                                                      8d11s21
            ibcoff=ibctop                                               8d11s21
           end if                                                       8d11s21
          end do                                                        8d11s21
         end if                                                          1d11s21
        end if                                                          5d10s21
        ibcoff=ibc0                                                     1d8s21
        do jsbv1=1,nsymb                                                1d8s21
         jsbv2=multh(jsbv1,jsbv12)                                      1d8s21
         if(jsbv2.ge.jsbv1)then                                         1d8s21
          if(jsbv12.eq.1)then                                           1d8s21
           joffdnon=joffdnon+nrootu*nvirt(jsbv1)*nfdatb(2,1,jsb)        8d26s21
           nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                        1d8s21
          else                                                          1d8s21
           nvv=nvirt(jsbv1)*nvirt(jsbv2)                                1d8s21
          end if                                                        1d8s21
          nvv=nvv*nrootu                                                1d8s21
          do l=1,4                                                      1d8s21
           joffdnon=joffdnon+nvv*nfdatb(2,l,jsb)                        8d26s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
       end do                                                           1d8s21
       do isbv1=1,nsymb                                                 1d8s21
        isbv2=multh(isbv1,isbv12)                                       1d8s21
        if(isbv2.ge.isbv1)then                                          1d8s21
         if(isbv12.eq.1)then                                            1d8s21
          ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdatk(2,1,isb)         8d26s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         1d8s21
         else                                                           1d8s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 1d8s21
         end if                                                         1d8s21
         nvv=nvv*nrootu                                                 1d8s21
         do l=1,4                                                       1d8s21
          ioffdnon=ioffdnon+nvv*nfdatk(2,l,isb)                         8d26s21
         end do                                                         1d8s21
        end if                                                          1d8s21
       end do                                                           1d8s21
      end do                                                            1d8s21
      return                                                            1d8s21
      end                                                               1d8s21
