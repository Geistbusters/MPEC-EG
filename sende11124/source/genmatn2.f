c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genmatn2(ncsfb,ncsfm,ncsfk,ncsf1,iwpb1,iwpk1,ncsf2,    2d23s21
     $     iwpb2,iwpk2,vec,ndim,mcol,xout,ff,icsf0,mcsf,fmul,bc,ibc)    11d10s22
      implicit real*8 (a-h,o-z)                                         11d13s20
c
c     like genmatn when icsf0=1 and mcsf=ncsfk.                         2d23s21
c     multiply iwpb1*iwpk1*iwpb2*iwpk2*vec and add to xout              11d26s20
c                                                                       11d13s20
      dimension vec(ndim,*),xout(ncsfb,*)
      include "common.store"                                            11d13s20
      nn=min(ncsfb,ncsfm,ncsfk,ncsf1,ncsf2,mcol)                        9d30s22
      nx=max(ncsfb,ncsfm,ncsfk,ncsf1,ncsf2,mcol)                        9d30s22
      ibcoffo=ibcoff                                                    10s3s22
      if(nx.le.5)then                                                   10d14s22
c
c     iwbp1*(iwpk1*(iwpb2*(iwpk2*vec))                                        11d13s20
c
      itmp1=ibcoff
      itmp2=itmp1+mcol*ncsf2
      itmp3=itmp2+mcol*ncsfm
      ibcoff=itmp3+mcol*ncsf1
      call enough('genmatn2.  1',bc,ibc)
      if(iwpk2.lt.0)then                                                12d5s20
       nusedi=ibc(-iwpk2)/2                                             12d5s20
       if(2*nusedi.ne.ibc(-iwpk2))nusedi=nusedi+1                       12d5s20
       icmp1=-iwpk2+1                                                   12d5s20
       icmp2=icmp1+ncsfk                                                12d5s20
       icmp3=icmp2+nusedi                                               12d5s20
       call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),vec,            2d23s21
     $       ndim,ncsfk,bc(itmp1),ncsf2,ncsf2,mcol,0d0,icsf0,mcsf,bc,   11d14s22
     $      ibc)                                                        11d14s22
      else                                                              12d5s20
       jwpk2=iwpk2+ncsf2*(icsf0-1)                                      2d23s21
       call dgemm('n','n',ncsf2,mcol,mcsf,1d0,bc(jwpk2),ncsf2,          2d23s21
     $     vec,ndim,0d0,bc(itmp1),ncsf2,
     d' genmatn2.  1')
      end if                                                            12d5s20
      if(iwpb2.lt.0)then                                                12d5s20
       nusedi=ibc(-iwpb2)/2                                             12d5s20
       if(2*nusedi.ne.ibc(-iwpb2))nusedi=nusedi+1                       12d5s20
       icmp1=-iwpb2+1                                                   12d5s20
       icmp2=icmp1+ncsf2                                                12d5s20
       icmp3=icmp2+nusedi                                               12d5s20
       call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmp1),       12d5s20
     $       ncsf2,ncsf2,bc(itmp2),ncsfm,ncsfm,mcol,0d0,1d0)            3d25s22
      else                                                              12d5s20
       call dgemm('n','n',ncsfm,mcol,ncsf2,1d0,bc(iwpb2),ncsfm,          11d26s20
     $     bc(itmp1),ncsf2,0d0,bc(itmp2),ncsfm,                         12d5s20
     d' genmatn2.  2')
      end if                                                            12d5s20
      if(iwpk1.lt.0)then                                                12d5s20
       nusedi=ibc(-iwpk1)/2                                             12d5s20
       if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
       icmp1=-iwpk1+1                                                   12d5s20
       icmp2=icmp1+ncsfm                                                12d5s20
       icmp3=icmp2+nusedi                                               12d5s20
       call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmp2),       12d5s20
     $       ncsfm,ncsfm,bc(itmp3),ncsf1,ncsf1,mcol,0d0,1d0)            3d25s22
      else                                                              12d5s20
       call dgemm('n','n',ncsf1,mcol,ncsfm,1d0,bc(iwpk1),ncsf1,
     $     bc(itmp2),ncsfm,0d0,bc(itmp3),ncsf1,
     d' genmatn2.  3')
      end if                                                            12d5s20
      if(iwpb1.lt.0)then                                                12d5s20
       nusedi=ibc(-iwpb1)/2                                             12d5s20
       if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
       icmp1=-iwpb1+1                                                   12d5s20
       icmp2=icmp1+ncsf1                                                12d5s20
       icmp3=icmp2+nusedi                                               12d5s20
       call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmp3),       12d5s20
     $       ncsf1,ncsf1,xout,ncsfb,ncsfb,mcol,ff,fmul)                 3d25s22
      else                                                              12d5s20
       call dgemm('n','n',ncsfb,mcol,ncsf1,fmul,bc(iwpb1),ncsfb,        2d23s21
     $     bc(itmp3),ncsf1,ff,xout,ncsfb,                               12d4s20
     d' genmatn2.  4')
      end if                                                            12d5s20
      else                                                              10d14s22
       if(iwpb1.lt.0)then                                               2d14s23
        nn=ibc(-iwpb1)                                                  2d14s23
        n123=nn*ncsfm                                                   2d14s23
        n124=nn*ncsf2
        n125=nn*mcsf                                                    2d14s23
        n126=nn*mcol
        if(iwpk1.lt.0)then                                              2d14s23
         n123=n123/2                                                    2d14s23
        end if                                                          2d14s23
       else                                                             2d14s23
        if(iwpk1.lt.0)then                                              2d14s23
         nn=ibc(-iwpk1)                                                 2d14s23
         n123=ncsfb*nn                                                  2d14s23
        else                                                            2d14s23
         n123=ncsfb*ncsf1*ncsfm                                          2d14s23
        end if                                                          2d14s23
        n124=ncsfb*ncsf1*ncsf2
        n125=ncsfb*ncsf1*mcsf                                           2d14s23
        n126=ncsfb*ncsf1*mcol                                            9d30s22
       end if
       if(iwpk1.lt.0)then
        nn=ibc(-iwpk1)
        n234=nn*ncsf2
        n235=nn*mcsf                                                    2d14s23
        n236=nn*mcol
       else
        n234=ncsf1*ncsfm*ncsf2                                           9d30s22
        n235=ncsf1*ncsfm*mcsf                                           2d14s23
        n236=ncsf1*ncsfm*mcol                                            9d30s22
       end if
       if(iwpb2.lt.0)then                                               2d14s23
        nn=ibc(-iwpb2)                                                  2d14s23
        n134=ncsfb*nn                                                   2d14s23
        n345=nn*ncsfm                                                   2d14s23
       else                                                             2d14s23
        n134=ncsfb*ncsfm*ncsf2                                           9d30s22
        n345=ncsfm*ncsf2*mcsf                                           2d14s23
       end if                                                           2d14s23
       n145=ncsfb*ncsf2*mcsf                                            2d14s23
       n245=ncsf1*ncsf2*mcsf                                            2d14s23
       n456=ncsf2*mcsf*mcol                                             2d14s23
       n136=ncsfb*ncsfm*mcol                                            9d30s22
       n135=ncsfb*ncsfm*mcsf                                            2d14s23
       n146=ncsfb*ncsf2*mcol                                            9d30s22
       n156=ncsfb*ncsfk*mcol                                            9d30s22
       n246=ncsf1*ncsf2*mcol                                            9d30s22
       n256=ncsf1*ncsfk*mcol                                            9d30s22
       n346=ncsfm*ncsf2*mcol                                            9d30s22
       n356=ncsfm*ncsfk*mcol                                            9d30s22
cccc
       ntry1=n145+n156                                                  9d30s22
       ntry2=n146+n456                                                  9d30s22
       if(ntry1.le.ntry2)then                                           9d30s22
        nuse1=ntry1                                                     9d30s22
        icase1=1                                                        9d30s22
       else                                                             9d30s22
        nuse1=ntry2                                                     9d30s22
        icase1=2                                                        9d30s22
       end if                                                           9d30s22
       ntry1=n345+n356                                                  9d30s22
       ntry2=n346+n456                                                  9d30s22
       if(ntry1.le.ntry2)then                                           9d30s22
        nuse2=ntry1                                                     9d30s22
        icase2=3                                                        9d30s22
       else                                                             9d30s22
        nuse2=ntry2                                                     9d30s22
        icase2=5                                                        9d30s22
       end if                                                           9d30s22
       ntry1=n134+nuse1                                                 9d30s22
       ntry2=n136+nuse2                                                 9d30s22
       if(ntry1.le.ntry2)then                                           9d30s22
        nuse12=ntry1                                                    9d30s22
        icase12=icase1                                                  9d30s22
       else                                                             9d30s22
        nuse12=ntry2                                                    9d30s22
        icase12=icase2                                                  9d30s22
       end if                                                           9d30s22
       ntry=n135+n156+n345                                              9d30s22
       if(ntry.le.nuse12)then                                           9d30s22
        nuse12=ntry                                                     9d30s22
        icase12=4                                                       9d30s22
       end if                                                           9d30s22
       nuse12=nuse12+n123                                               9d30s22
       ntry3=n245+n256                                                  9d30s22
       ntry4=n246+n456                                                  9d30s22
       if(ntry3.le.ntry4)then                                           9d30s22
        nuse3=ntry3                                                     9d30s22
        icase3=9                                                        9d30s22
       else                                                             9d30s22
        nuse3=ntry4                                                     9d30s22
        icase3=10                                                       9d30s22
       end if                                                           9d30s22
       ntry3=n345+n356                                                  9d30s22
       ntry4=n346+n456                                                  9d30s22
       if(ntry3.le.ntry4)then                                           9d30s22
        nuse4=ntry3                                                     9d30s22
        icase4=11                                                       9d30s22
       else                                                             9d30s22
        nuse4=ntry4                                                     9d30s22
        icase4=14                                                       9d30s22
       end if                                                           9d30s22
       ntry3=n234+nuse3                                                 9d30s22
       ntry4=n236+nuse4                                                 9d30s22
       if(ntry3.le.ntry4)then                                           9d30s22
        nuse34=ntry3                                                    9d30s22
        icase34=icase3                                                  9d30s22
       else                                                             9d30s22
        nuse34=ntry4                                                    9d30s22
        icase34=icase4                                                  9d30s22
       end if                                                           9d30s22
       ntry=n235+n256+n345                                              9d30s22
       if(ntry.le.nuse34)then                                           9d30s22
        icase34=13                                                      9d30s22
        nuse34=ntry                                                     9d30s22
       end if                                                           9d30s22
       nuse34=nuse34+n126                                               9d30s22
       if(nuse12.le.nuse34)then                                         9d30s22
        nuse1234=nuse12                                                 9d30s22
        icase1234=icase12                                               9d30s22
       else                                                             9d30s22
        nuse1234=nuse34                                                 9d30s22
        icase1234=icase34                                               9d30s22
       end if                                                           9d30s22
       ntry5=n146+n456                                                  9d30s22
       ntry6=n145+n156                                                  9d30s22
       if(ntry5.le.ntry6)then                                           9d30s22
        nuse5=ntry5                                                     9d30s22
        icase5=6                                                        9d30s22
       else                                                             9d30s22
        nuse5=ntry6                                                     9d30s22
        icase5=7                                                        9d30s22
       end if                                                           9d30s22
       nuse5=nuse5+n124+n234                                            9d30s22
       if(nuse5.le.nuse1234)then                                        9d30s22
        nuse1234=nuse5                                                  9d30s22
        icase1234=icase5                                                9d30s22
       end if                                                           9d30s22
       ntry5=n234+n245                                                  9d30s22
       ntry6=n235+n345                                                  9d30s22
       if(ntry5.le.ntry6)then                                           9d30s22
        nuse6=ntry5                                                     9d30s22
        icase6=8                                                        9d30s22
       else                                                             9d30s22
        nuse6=ntry6                                                     9d30s22
        icase6=12                                                       9d30s22
       end if                                                           9d30s22
       nuse6=nuse6+n125+n156                                            9d30s22
       if(nuse6.le.nuse1234)then                                        9d30s22
        icase1234=icase6                                                9d30s22
       end if                                                           9d30s22
       if(icase1234.eq.3.or.icase1234.eq.4.or.icase1234.eq.11.or.       9d30s22
     $    icase1234.eq.13.or.icase1234.eq.12)then                       9d30s22
        itmpcd=ibcoff                                                   9d30s22
        ibcoff=itmpcd+ncsfm*mcsf                                        10d6s22
        call enough('genmatn2.  2',bc,ibc)
        if(iwpk2.lt.0)then                                                12d5s20
         nusedi=ibc(-iwpk2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk2+1                                                   12d5s20
         icmp2=icmp1+ncsfk                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         itrial=ibcoff                                                  9d30s22
         ibcoff=itrial+ncsf2*ncsfk                                      9d30s22
         call enough('genmatn2.  3',bc,ibc)
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsf2,ncsfk)                                              9d30s22
        else                                                            9d30s22
         itrial=iwpk2                                                   9d30s22
        end if                                                          9d30s22
        jtrial=itrial+ncsf2*(icsf0-1)                                   10d7s22
        if(iwpb2.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpb2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb2+1                                                   12d5s20
         icmp2=icmp1+ncsf2                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(jtrial),    9d30s22
     $       ncsf2,ncsf2,bc(itmpcd),ncsfm,ncsfm,mcsf,0d0,1d0)           10d7s22
        else                                                            9d30s22
         call dgemm('n','n',ncsfm,mcsf,ncsf2,1d0,bc(iwpb2),ncsfm,       10d6s22
     $        bc(jtrial),ncsf2,0d0,bc(itmpcd),ncsfm,                    10d7s22
     d' genmatn2.  1')
        end if                                                          9d30s22
        if(iwpk2.lt.0)ibcoff=itrial                                     9d30s22
       end if                                                           9d30s22
       if(icase1234.le.5.or.icase1234.eq.7)then                         10s3s22
        itmpab=ibcoff                                                   9d30s22
        ibcoff=itmpab+ncsfb*ncsfm                                       9d30s22
        call enough('genmatn2.  4',bc,ibc)
        call prodn(iwpb1,iwpk1,ncsfb,ncsfm,ncsf1,bc(itmpab),bc,ibc,1d0, 2d14s23
     $       0d0)                                                       2d14s23
       end if                                                           9d30s22
       if(icase1234.eq.3.or.icase1234.eq.4.or.icase1234.eq.11.or.       10s3s22
     $      icase1234.eq.13)then                                        10s3s22
        itmpcde=ibcoff                                                  9d30s22
        ibcoff=itmpcde+ncsfm*mcol                                       9d30s22
        call enough('genmatn2.  6',bc,ibc)
        call dgemm('n','n',ncsfm,mcol,mcsf,1d0,bc(itmpcd),ncsfm,        10d7s22
     $       vec,ndim,0d0,bc(itmpcde),ncsfm,                            9d30s22
     d' genmatn2.  3')
       end if                                                           9d30s22
       if(icase1234.eq.12.or.icase1234.eq.13)then                       9d30s22
        itmpbcd=ibcoff                                                  9d30s22
        ibcoff=itmpbcd+ncsf1*mcsf                                       10d7s22
        call enough('genmatn2.  7',bc,ibc)
        if(iwpk1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpcd),    9d30s22
     $       ncsfm,ncsfm,bc(itmpbcd),ncsf1,ncsf1,mcsf,0d0,1d0)          10d7s22
        else                                                            9d30s22
         call dgemm('n','n',ncsf1,mcsf,ncsfm,1d0,bc(iwpk1),ncsf1,       10d7s22
     $        bc(itmpcd),ncsfm,0d0,bc(itmpbcd),ncsf1,                   9d30s22
     d' genmatn2.  4')
        end if                                                          9d30s22
       end if                                                           9d30s22
       if(icase1234.le.2)then                                           10d14s22
        itmpabc=ibcoff                                                  9d30s22
        ibcoff=itmpabc+ncsfb*ncsf2                                      9d30s22
        call enough('genmatn2.  8',bc,ibc)
        mrow=ncsfb
        ncol=ncsf2
        mwant=itmpabc
        mgot=itmpabc                                                    2d14s23
        call prodn(itmpab,iwpb2,ncsfb,ncsf2,ncsfm,bc(mgot),bc,ibc,      2d14s23
     $       1d0,0d0)                                                   2d14s23
       end if                                                           9d30s22
       if(icase1234.ge.6.and.icase1234.le.10)then                       9d30s22
        itmpbc=ibcoff                                                   9d30s22
        ibcoff=itmpbc+ncsf1*ncsf2                                       9d30s22
        call enough('genmatn2. 10',bc,ibc)
        mrow=ncsf1
        ncol=ncsf2
        mwant=itmpbc
        mgot=itmpbc                                                     2d14s23
        call prodn(iwpk1,iwpb2,ncsf1,ncsf2,ncsfm,bc(mgot),bc,ibc,       2d14s23
     $       1d0,0d0)                                                   2d14s23
        if(icase1234.eq.8.or.icase1234.eq.9)then                        9d30s22
         itmpbcd=ibcoff                                                 9d30s22
         ibcoff=itmpbcd+ncsf1*mcsf                                      10d6s22
         call enough('genmatn2. 12',bc,ibc)
         if(iwpk2.lt.0)then                                             9d30s22
          itrial=ibcoff                                                 9d30s22
          ibcoff=itrial+ncsf2*ncsfk                                     9d30s22
          call enough('genmatn2. 13',bc,ibc)
          nusedi=ibc(-iwpk2)/2                                             12d5s20
          if(2*nusedi.ne.ibc(-iwpk2))nusedi=nusedi+1                       12d5s20
          icmp1=-iwpk2+1                                                   12d5s20
          icmp2=icmp1+ncsfk                                                12d5s20
          icmp3=icmp2+nusedi                                               12d5s20
          call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsf2,ncsfk)                                              9d30s22
         else                                                           9d30s22
          itrial=iwpk2                                                  9d30s22
         end if                                                         9d30s22
         jtrial=itrial+ncsf2*(icsf0-1)                                  10d6s22
         call dgemm('n','n',ncsf1,mcsf,ncsf2,1d0,bc(itmpbc),ncsf1,      10d6s22
     $        bc(jtrial),ncsf2,0d0,bc(itmpbcd),ncsf1,                   10d7s22
     d' genmatn2.  7')
         if(iwpk2.lt.0)ibcoff=itrial                                    9d30s22
        end if                                                          9d30s22
        if(icase1234.eq.6.or.icase1234.eq.7)then                        10d14s22
         itmpabc=ibcoff                                                 10d14s22
         ibcoff=itmpabc+ncsfb*ncsf2                                     10d14s22
         call enough('genmatn2. 14',bc,ibc)
         if(iwpb1.lt.0)then                                                12d5s20
          nusedi=ibc(-iwpb1)/2                                             12d5s20
          if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
          icmp1=-iwpb1+1                                                   12d5s20
          icmp2=icmp1+ncsf1                                                12d5s20
          icmp3=icmp2+nusedi                                               12d5s20
          call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpbc),   10d14s22
     $         ncsf1,ncsf1,bc(itmpabc),ncsfb,ncsfb,ncsf2,0d0,1d0)       10d14s22
         else                                                              12d5s20
          call dgemm('n','n',ncsfb,ncsf2,ncsf1,1d0,bc(iwpb1),ncsfb,
     $     bc(itmpbc),ncsf1,0d0,bc(itmpabc),ncsfb,                      10d14s22
     d' genmatn2.  8')
         end if                                                            12d5s20
        end if                                                          10d14s22
       end if                                                           9d30s22
       if(icase1234.eq.2.or.icase1234.eq.5.or.icase1234.eq.10.or.       10s3s22
     $      icase1234.eq.14.or.icase1234.eq.6)then                      10s3s22
        itmpde=ibcoff                                                   10s3s22
        ibcoff=itmpde+ncsf2*mcol                                        10s3s22
        call enough('genmatn2. 15',bc,ibc)
        if(iwpk2.lt.0)then                                              10s3s22
         nusedi=ibc(-iwpk2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk2+1                                                   12d5s20
         icmp2=icmp1+ncsfk                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),vec,             12d5s20
     $       ndim,ncsfk,bc(itmpde),ncsf2,ncsf2,mcol,0d0,icsf0,mcsf,bc,  11d14s22
     $        ibc)                                                      11d14s22
        else                                                            10s3s22
         jwpk2=iwpk2+ncsf2*(icsf0-1)                                    10d6s22
         call dgemm('n','n',ncsf2,mcol,mcsf,1d0,bc(jwpk2),ncsf2,        10d6s22
     $        vec,ndim,0d0,bc(itmpde),ncsf2,                            10s3s22
     d' genmatn2.  9')
        end if                                                          10s3s22
       end if                                                           10s3s22
       if(icase1234.eq.5.or.icase1234.eq.14)then                        10s3s22
        itmpccde=ibcoff                                                 10s3s22
        ibcoff=itmpccde+ncsfm*mcol                                      10s3s22
        call enough('genmatn2. 16',bc,ibc)
        if(iwpb2.lt.0)then                                              10s3s22
         nusedi=ibc(-iwpb2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb2+1                                                   12d5s20
         icmp2=icmp1+ncsf2                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpde),    10s3s22
     $       ncsf2,ncsf2,bc(itmpccde),ncsfm,ncsfm,mcol,0d0,1d0)         10s3s22
        else                                                            10s3s22
         call dgemm('n','n',ncsfm,mcol,ncsf2,1d0,bc(iwpb2),ncsfm,       10s3s22
     $        bc(itmpde),ncsf2,0d0,bc(itmpccde),ncsfm,                  10s3s22
     d' genmatn2. 10')
        end if                                                          10s3s22
       end if                                                           10s3s22
       if(icase1234.eq.1)then                                           9d30s22
c (((ab)c)d)e
        itmpad=ibcoff                                                   9d30s22
        ibcoff=itmpad+ncsfb*mcsf                                        10d6s22
        call enough('genmatn2. 17',bc,ibc)
        if(iwpk2.lt.0)then                                              9d30s22
         itrial=ibcoff                                                  9d30s22
         ibcoff=itrial+ncsf2*ncsfk                                      9d30s22
         call enough('genmatn2. 18',bc,ibc)
         nusedi=ibc(-iwpk2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk2+1                                                   12d5s20
         icmp2=icmp1+ncsfk                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsf2,ncsfk)                                              9d30s22
        else                                                            9d30s22
         itrial=iwpk2                                                   9d30s22
        end if                                                          9d30s22
        jtrial=itrial+ncsf2*(icsf0-1)                                   10d6s22
        call dgemm('n','n',ncsfb,mcsf,ncsf2,1d0,bc(itmpabc),ncsfb,      10d6s22
     $       bc(jtrial),ncsf2,0d0,bc(itmpad),ncsfb,                     9d30s22
     d' genmatn2. 11')
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpad),ncsfb,       10d6s22
     $       vec,ndim,ff,xout,ncsfb,                                    10s3s22
     d' genmatn2. 12')
        if(iwpk2.lt.0)ibcoff=itrial                                     9d30s22
       else if(icase1234.eq.2)then                                           9d30s22
c ((ab)c)(de)
        call dgemm('n','n',ncsfb,mcol,ncsf2,fmul,bc(itmpabc),ncsfb,     10d7s22
     $       bc(itmpde),ncsf2,ff,xout,ncsfb,                            10s3s22
     d' genmatn2. 13')
       else if(icase1234.eq.3)then                                           9d30s22
c (ab)((cd)e)
        call dgemm('n','n',ncsfb,mcol,ncsfm,fmul,bc(itmpab),ncsfb,      10d7s22
     $       bc(itmpcde),ncsfm,ff,xout,ncsfb,                           10s3s22
     d' genmatn2. 14')
       else if(icase1234.eq.4)then                                      10s3s22
c ((ab)(cd))e                                                           10s3s22
        call dgemm('n','n',ncsfb,mcol,ncsfm,fmul,bc(itmpab),ncsfb,      10d7s22
     $       bc(itmpcde),ncsfm,ff,xout,ncsfb,                           10s3s22
     d' genmatn2. 15')
       else if(icase1234.eq.5)then                                           9d30s22
c (ab)(c(de))
        call dgemm('n','n',ncsfb,mcol,ncsfm,fmul,bc(itmpab),ncsfb,      10d7s22
     $       bc(itmpccde),ncsfm,ff,xout,ncsfb,                          10s3s22
     d' genmatn2. 16')
       else if(icase1234.eq.6)then                                           9d30s22
c (a(bc))(de)
        call dgemm('n','n',ncsfb,mcol,ncsf2,fmul,bc(itmpabc),ncsfb,     10d14s22
     $       bc(itmpde),ncsf2,ff,xout,ncsfb,                            10d14s22
     d' genmatn2. 17')
       else if(icase1234.eq.7)then                                           9d30s22
c ((a(bc))d)e
        itmpx=ibcoff                                                    10s3s22
        ibcoff=itmpx+ncsfb*mcsf                                         10d6s22
        call enough('genmatn2. 19',bc,ibc)
        if(iwpk2.lt.0)then                                              10s3s22
         itrial=ibcoff                                                  10s3s22
         ibcoff=itrial+ncsf2*ncsfk                                      10s3s22
         call enough('genmatn2. 20',bc,ibc)
         nusedi=ibc(-iwpk2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk2+1                                                   12d5s20
         icmp2=icmp1+ncsfk                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsf2,ncsfk)                                              9d30s22
        else                                                            10s3s22
         itrial=iwpk2                                                   10s3s22
        end if                                                          10s3s22
        jtrial=itrial+ncsf2*(icsf0-1)                                   10d6s22
        call dgemm('n','n',ncsfb,mcsf,ncsf2,1d0,bc(itmpabc),ncsfb,      10d6s22
     $       bc(jtrial),ncsf2,0d0,bc(itmpx),ncsfb,                      10d6s22
     d' genmatn2. 18')
        if(iwpk2.lt.0)ibcoff=itrial                                     10s3s22
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpx),ncsfb,        10d7s22
     $       vec,ndim,ff,xout,ncsfb,                                    10s3s22
     d' genmatn2. 19')
       else if(icase1234.eq.8)then                                           9d30s22
c (a((bc)d))e
        itmpad=ibcoff                                                   9d30s22
        ibcoff=itmpad+ncsfb*mcsf                                        10d7s22
        call enough('genmatn2. 21',bc,ibc)
        if(iwpb1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpbcd),   9d30s22
     $       ncsf1,ncsf1,bc(itmpad),ncsfb,ncsfb,mcsf,0d0,1d0)           10d7s22
        else                                                            9d30s22
         call dgemm('n','n',ncsfb,mcsf,ncsf1,1d0,bc(iwpb1),ncsfb,       10d7s22
     $        bc(itmpbcd),ncsf1,0d0,bc(itmpad),ncsfb,                   9d30s22
     d' genmatn2. 20')
        end if                                                          9d30s22
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpad),ncsfb,       10d7s22
     $       vec,ndim,ff,xout,ncsfb,                                    10s3s22
     d' genmatn2. 21')
       else if(icase1234.eq.9)then                                           9d30s22
c  a(((bc)d)e)
        itmpx=ibcoff                                                    10d14s22
        ibcoff=itmpx+ncsf1*mcol                                         10d14s22
        call enough('genmatn2. 22',bc,ibc)
        call dgemm('n','n',ncsf1,mcol,mcsf,1d0,bc(itmpbcd),ncsf1,       10d14s22
     $       vec,ndim,0d0,bc(itmpx),ncsf1,                              10d14s22
     d' genmatn2. 22')
        call prodn(iwpb1,itmpx,ncsfb,mcol,ncsf1,xout,bc,ibc,fmul,       2d14s23
     $       ff)                                                        2d14s23
       else if(icase1234.eq.10)then                                           9d30s22
c a((bc)*(de))
        itmpx=ibcoff                                                    10d14s22
        ibcoff=itmpx+ncsf1*mcol                                         10d14s22
        call enough('genmatn2. 24',bc,ibc)
        call dgemm('n','n',ncsf1,mcol,ncsf2,1d0,bc(itmpbc),ncsf1,       10d14s22
     $       bc(itmpde),ncsf2,0d0,bc(itmpx),ncsf1,                      10d14s22
     d' genmatn2. 24')
        call prodn(iwpb1,itmpx,ncsfb,mcol,ncsf1,xout,bc,ibc,fmul,       2d14s23
     $       ff)                                                        2d14s23
       else if(icase1234.eq.11)then                                           9d30s22
c a(b((cd)e))
        itmpbe=ibcoff                                                   9d30s22
        ibcoff=itmpbe+ncsfm*mcol                                        9d30s22
        call enough('genmatn2. 26',bc,ibc)
        if(iwpk1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpcde),   9d30s22
     $       ncsfm,ncsfm,bc(itmpbe),ncsf1,ncsf1,mcol,0d0,1d0)            3d25s22
        else                                                            9d30s22
         call dgemm('n','n',ncsf1,mcol,ncsfm,1d0,bc(iwpk1),ncsf1,       9d30s22
     $        bc(itmpcde),ncsfm,0d0,bc(itmpbe),ncsf1,                   9d30s22
     d' genmatn2. 26')
        end if                                                          9d30s22
        if(iwpb1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpbe),    9d30s22
     $       ncsf1,ncsf1,xout,ncsfb,ncsfb,mcol,ff,1d0)                  10s3s22
        else                                                            9d30s22
         call dgemm('n','n',ncsfb,mcol,ncsf1,fmul,bc(iwpb1),ncsfb,      10d7s22
     $        bc(itmpbe),ncsf1,ff,xout,ncsfb,                           10s3s22
     d' genmatn2. 27')
        end if                                                          9d30s22
       else if(icase1234.eq.12)then                                           9d30s22
c (a(b(cd)))e
        itmpad=ibcoff                                                   9d30s22
        ibcoff=itmpad+mcsf*ncsfb                                        10d7s22
        call enough('genmatn2. 27',bc,ibc)
        if(iwpb1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpbcd),   9d30s22
     $       ncsf1,ncsf1,bc(itmpad),ncsfb,ncsfb,mcsf,0d0,1d0)           10d7s22
        else                                                            9d30s22
         call dgemm('n','n',ncsfb,mcsf,ncsf1,1d0,bc(iwpb1),ncsfb,       10d7s22
     $        bc(itmpbcd),ncsf1,0d0,bc(itmpad),ncsfb,                   9d30s22
     d' genmatn2. 28')
        end if                                                          9d30s22
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpad),ncsfb,vec,   10d7s22
     $       ndim,ff,xout,ncsfb,                                        10s3s22
     d' genmatn2. 29')
       else if(icase1234.eq.13)then                                           9d30s22
c a((b(cd))e)
        itmpx=ibcoff                                                    10d14s22
        ibcoff=itmpx+ncsf1*mcol                                         10d14s22
        call enough('genmatn2. 28',bc,ibc)
        call dgemm('n','n',ncsf1,mcol,mcsf,1d0,bc(itmpbcd),ncsf1,       10d14s22
     $       vec,ndim,0d0,bc(itmpx),ncsf1,                              10d14s22
     d' genmatn2. 30')
        call prodn(iwpb1,itmpx,ncsfb,mcol,ncsf1,xout,bc,ibc,fmul,       2d14s23
     $       ff)                                                        2d14s23
       else if(icase1234.eq.14)then                                           9d30s22
c a(b(c(de)))
        itmpx=ibcoff                                                    10s3s22
        ibcoff=itmpx+ncsf1*mcol                                         10s3s22
        call enough('genmatn2. 30',bc,ibc)
        if(iwpk1.lt.0)then                                              10s3s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpccde),  10s3s22
     $       ncsfm,ncsfm,bc(itmpx),ncsf1,ncsf1,mcol,0d0,1d0)            10s3s22
        else                                                            10s3s22
         call dgemm('n','n',ncsf1,mcol,ncsfm,1d0,bc(iwpk1),ncsf1,       10s3s22
     $        bc(itmpccde),ncsfm,0d0,bc(itmpx),ncsf1,                   10s3s22
     d' genmatn2. 32')
        end if                                                          10s3s22
        if(iwpb1.lt.0)then                                              10s3s22
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpx),     10s3s22
     $       ncsf1,ncsf1,xout,ncsfb,ncsfb,mcol,ff,1d0)                  10s3s22
        else                                                            10s3s22
         call dgemm('n','n',ncsfb,mcol,ncsf1,fmul,bc(iwpb1),ncsfb,      10d7s22
     $        bc(itmpx),ncsf1,ff,xout,ncsfb,                            10s3s22
     d' genmatn2. 33')
        end if                                                          10s3s22
       else                                                             9d30s22
        write(6,*)('unknown case1234 '),icase1234
        stop 'case1234'
       end if                                                           9d30s22
      end if                                                            9d30s22
      ibcoff=ibcoffo                                                    10s3s22
      return
      end
