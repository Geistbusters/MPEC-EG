c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genmatn3(ncsfb,ncsfm,ncsfk,ncsf1,iwpb1,iwpk1,ncsf2,    2d23s21
c                           1     3     5    2    vects          4
     $     iwpb2,iwpk2,vec,ndim,mcol,xout,ff,icsf0,mcsf,jcsf0,kcsf,fmul,11d14s22
     $     bc,ibc)                                                      11d14s22
c          iprod
      implicit real*8 (a-h,o-z)                                         11d13s20
c
c     like genmatn2 when jcsf0=1 and kcsf=ncsfb.                        11d2s22
c     multiply iwpb1*iwpk1*iwpb2*iwpk2*vec and add to xout              11d26s20
c     iwpb1: ncsfb,kcsf
c     iwpk1: ncsf1,ncsfm
c     iwpb2: ncsfm,ncsf2
c     iwpk2: ncsf2,mcsf
c     vec: ndim,mcol
c                                                                       11d13s20
      dimension vec(ndim,*),xout(ncsfb,*)                               11d4s22
      include "common.store"                                            11d13s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
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
      call enough('genmatn3.  1',bc,ibc)
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
     d' genmatn3.  1')
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
       jtmp3=itmp3+jcsf0-1                                               11d3s22
      else                                                              12d5s20
       jwpk1=iwpk1+jcsf0-1                                              11d2s22
       call dgemm('n','n',kcsf,mcol,ncsfm,1d0,bc(jwpk1),ncsf1,          11d2s22
     $     bc(itmp2),ncsfm,0d0,bc(itmp3),ncsf1,                         11d4s22
     d' genmatn3.  3')
       jtmp3=itmp3
      end if                                                            12d5s20
      icpy=ibcoff
      ibcoff=icpy+ncsfb*mcol
      jcpy=icpy-1
      do i=1,mcol
       do j=1,ncsfb
        bc(jcpy+j)=xout(j,i)
       end do
      end do
      call dgemm('n','n',ncsfb,mcol,kcsf,fmul,bc(iwpb1),ncsfb,          11d3s22
     $     bc(jtmp3),ncsf1,ff,xout,ncsfb,                                11d2s22
     d' genmatn3.  4')
      else                                                              10d14s22
       n123=ncsfb*kcsf*ncsfm
       n124=ncsfb*kcsf*ncsf2
       n125=ncsfb*kcsf*mcsf                                             10d6s22
       n126=ncsfb*kcsf*mcol                                             9d30s22
       n134=ncsfb*ncsfm*ncsf2                                           9d30s22
       n136=ncsfb*ncsfm*mcol                                            9d30s22
       n135=ncsfb*ncsfm*mcsf                                            10d6s22
       n145=ncsfb*ncsf2*mcsf                                            10d6s22
       n146=ncsfb*ncsf2*mcol                                            9d30s22
       n156=ncsfb*mcsf*mcol                                             10d6s22
       n234=kcsf*ncsfm*ncsf2                                            9d30s22
       n235=kcsf*ncsfm*mcsf                                             10d6s22
       n236=kcsf*ncsfm*mcol                                             9d30s22
       n245=kcsf*ncsf2*mcsf                                             10d6s22
       n246=kcsf*ncsf2*mcol                                             9d30s22
       n256=kcsf*mcsf*mcol                                              10d6s22
       n345=ncsfm*ncsf2*mcsf                                            10d6s22
       n346=ncsfm*ncsf2*mcol                                            9d30s22
       n356=ncsfm*mcsf*mcol                                             10d6s22
       n456=ncsf2*mcsf*mcol                                             10d6s22
       if(iwpk2.lt.0)then                                               10s3s22
        n345=n345*3                                                     10s3s22
        n456=n456/2                                                     10s3s22
       end if                                                           10s3s22
       if(iwpb2.lt.0)then                                               10s3s22
        n234=n234*3                                                     10s3s22
        n345=n345/2                                                     10s3s22
        n346=n346/2                                                     10s3s22
       end if                                                           10s3s22
       if(iwpk1.lt.0)then                                               10s3s22
        n123=n123*3                                                     10s3s22
        n234=n234/2                                                     10s3s22
        n235=n235/2                                                     10s3s22
        n236=n236/2                                                     10s3s22
       end if                                                           10s3s22
       if(iwpb1.lt.0)then                                               10s3s22
        n123=n123/2                                                     10s3s22
        n124=n124/2                                                     10s3s22
        n125=n125/2                                                     10s3s22
        n126=n126/2                                                     10s3s22
       end if                                                           10s3s22
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
        call enough('genmatn3.  2',bc,ibc)
        if(iwpk2.lt.0)then                                                12d5s20
         nusedi=ibc(-iwpk2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk2+1                                                   12d5s20
         icmp2=icmp1+ncsfk                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         itrial=ibcoff                                                  9d30s22
         ibcoff=itrial+ncsf2*ncsfk                                      9d30s22
         call enough('genmatn3.  3',bc,ibc)
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
     d' genmatn3.  1')
        end if                                                          9d30s22
        if(iwpk2.lt.0)ibcoff=itrial                                     9d30s22
       end if                                                           9d30s22
       if(icase1234.le.5.or.icase1234.eq.7)then                         10s3s22
        itmpab=ibcoff                                                   9d30s22
        ibcoff=itmpab+ncsfb*ncsfm                                       9d30s22
        call enough('genmatn3.  4',bc,ibc)
        if(iwpk1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         itrial=ibcoff                                                  9d30s22
         ibcoff=itrial+ncsf1*ncsfm                                      9d30s22
         call enough('genmatn3.  5',bc,ibc)
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsf1,ncsfm)                                              9d30s22
        else                                                            9d30s22
         itrial=iwpk1                                                   9d30s22
        end if                                                          9d30s22
        jtrial=itrial+jcsf0-1                                           11d3s22
        if(iwpb1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),    9d30s22
     $       ncsf1,ncsf1,bc(itmpab),ncsfb,ncsfb,ncsfm,0d0,1d0)          9d30s22
        else                                                            9d30s22
         call dgemm('n','n',ncsfb,ncsfm,kcsf,1d0,bc(iwpb1),ncsfb,       11d3s22
     $        bc(jtrial),ncsf1,0d0,bc(itmpab),ncsfb,                    9d30s22
     d' genmatn3.  2')
        end if                                                          9d30s22
        if(iwpk1.lt.0)ibcoff=itrial                                     9d30s22
       end if                                                           9d30s22
       if(icase1234.eq.3.or.icase1234.eq.4.or.icase1234.eq.11.or.       10s3s22
     $      icase1234.eq.13)then                                        10s3s22
        itmpcde=ibcoff                                                  9d30s22
        ibcoff=itmpcde+ncsfm*mcol                                       9d30s22
        call enough('genmatn3.  6',bc,ibc)
        call dgemm('n','n',ncsfm,mcol,mcsf,1d0,bc(itmpcd),ncsfm,        10d7s22
     $       vec,ndim,0d0,bc(itmpcde),ncsfm,                            9d30s22
     d' genmatn3.  3')
       end if                                                           9d30s22
       if(icase1234.eq.12.or.icase1234.eq.13)then                       9d30s22
        itmpbcd=ibcoff                                                  9d30s22
        ibcoff=itmpbcd+ncsf1*mcsf                                       10d7s22
        call enough('genmatn3.  7',bc,ibc)
        if(iwpk1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpcd),    9d30s22
     $       ncsfm,ncsfm,bc(itmpbcd),ncsf1,ncsf1,mcsf,0d0,1d0)          10d7s22
         jtmpbcd=itmpbcd+jcsf0-1                                        11d4s22
        else                                                            9d30s22
         jwpk1=iwpk1+jcsf0-1
         call dgemm('n','n',kcsf,mcsf,ncsfm,1d0,bc(jwpk1),ncsf1,        11d4s22
     $        bc(itmpcd),ncsfm,0d0,bc(itmpbcd),ncsf1,                   9d30s22
     d' genmatn3.  4')
         jtmpbcd=itmpbcd
        end if                                                          9d30s22
       end if                                                           9d30s22
       if(icase1234.le.2)then                                           10d14s22
        itmpabc=ibcoff                                                  9d30s22
        ibcoff=itmpabc+ncsfb*ncsf2                                      9d30s22
        call enough('genmatn3.  8',bc,ibc)
        if(iwpb2.lt.0)then                                              9d30s22
         itrial=ibcoff                                                  9d30s22
         ibcoff=itrial+ncsfm*ncsf2                                      9d30s22
         call enough('genmatn3.  9',bc,ibc)
         nusedi=ibc(-iwpb2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb2+1                                                   12d5s20
         icmp2=icmp1+ncsf2                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsfm,ncsf2)                                              9d30s22
        else                                                            9d30s22
         itrial=iwpb2                                                   9d30s22
        end if                                                          9d30s22
        call dgemm('n','n',ncsfb,ncsf2,ncsfm,1d0,bc(itmpab),ncsfb,      9d30s22
     $       bc(itrial),ncsfm,0d0,bc(itmpabc),ncsfb,                    9d30s22
     d' genmatn3.  5')
        if(iwpb2.lt.0)ibcoff=itrial                                     9d30s22
       end if                                                           9d30s22
       if(icase1234.ge.6.and.icase1234.le.10)then                       9d30s22
        itmpbc=ibcoff                                                   9d30s22
        ibcoff=itmpbc+ncsf1*ncsf2                                       9d30s22
        call enough('genmatn3. 10',bc,ibc)
        if(iwpb2.lt.0)then                                              9d30s22
         itrial=ibcoff                                                  9d30s22
         ibcoff=itrial+ncsfm*ncsf2                                      9d30s22
         call enough('genmatn3. 11',bc,ibc)
         nusedi=ibc(-iwpb2)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb2))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb2+1                                                   12d5s20
         icmp2=icmp1+ncsf2                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsfm,ncsf2)                                              9d30s22
        else                                                            9d30s22
         itrial=iwpb2                                                   9d30s22
        end if                                                          9d30s22
        if(iwpk1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),    9d30s22
     $       ncsfm,ncsfm,bc(itmpbc),ncsf1,ncsf1,ncsf2,0d0,1d0)          9d30s22
         jtmpbc=itmpbc+jcsf0-1                                          11d4s22
        else                                                            9d30s22
         jwpk1=iwpk1+jcsf0-1                                            11d4s22
         call dgemm('n','n',kcsf,ncsf2,ncsfm,1d0,bc(jwpk1),ncsf1,       11d4s22
     $        bc(itrial),ncsfm,0d0,bc(itmpbc),ncsf1,                    9d30s22
     d' genmatn3.  6')
         jtmpbc=itmpbc                                                  11d4s22
        end if                                                          9d30s22
        if(iwpb2.lt.0)ibcoff=itrial                                     9d30s22
        if(icase1234.eq.8.or.icase1234.eq.9)then                        9d30s22
         itmpbcd=ibcoff                                                 9d30s22
         ibcoff=itmpbcd+kcsf*mcsf                                       11d4s22
         call enough('genmatn3. 12',bc,ibc)
         if(iwpk2.lt.0)then                                             9d30s22
          itrial=ibcoff                                                 9d30s22
          ibcoff=itrial+ncsf2*ncsfk                                     9d30s22
          call enough('genmatn3. 13',bc,ibc)
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
         call dgemm('n','n',kcsf,mcsf,ncsf2,1d0,bc(jtmpbc),ncsf1,       11d4s22
     $        bc(jtrial),ncsf2,0d0,bc(itmpbcd),kcsf,                    11d4s22
     d' genmatn3.  7')
         if(iwpk2.lt.0)ibcoff=itrial                                    9d30s22
        end if                                                          9d30s22
        if(icase1234.eq.6.or.icase1234.eq.7)then                        10d14s22
         itmpabc=ibcoff                                                 10d14s22
         ibcoff=itmpabc+ncsfb*ncsf2                                     10d14s22
         call enough('genmatn3. 14',bc,ibc)
         if(iwpb1.lt.0)then                                                12d5s20
          nusedi=ibc(-iwpb1)/2                                             12d5s20
          if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
          icmp1=-iwpb1+1                                                   12d5s20
          icmp2=icmp1+ncsf1                                                12d5s20
          icmp3=icmp2+nusedi                                               12d5s20
          call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(jtmpbc),   10d14s22
     $         ncsf1,kcsf,bc(itmpabc),ncsfb,ncsfb,ncsf2,0d0,1d0)        11d4s22
         else                                                              12d5s20
          call dgemm('n','n',ncsfb,ncsf2,kcsf,1d0,bc(iwpb1),ncsfb,      11d4s22
     $     bc(jtmpbc),ncsf1,0d0,bc(itmpabc),ncsfb,                      10d14s22
     d' genmatn3.  8')
         end if                                                            12d5s20
        end if                                                          10d14s22
       end if                                                           9d30s22
       if(icase1234.eq.2.or.icase1234.eq.5.or.icase1234.eq.10.or.       10s3s22
     $      icase1234.eq.14.or.icase1234.eq.6)then                      10s3s22
        itmpde=ibcoff                                                   10s3s22
        ibcoff=itmpde+ncsf2*mcol                                        10s3s22
        call enough('genmatn3. 15',bc,ibc)
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
     d' genmatn3.  9')
        end if                                                          10s3s22
       end if                                                           10s3s22
       if(icase1234.eq.5.or.icase1234.eq.14)then                        10s3s22
        itmpccde=ibcoff                                                 10s3s22
        ibcoff=itmpccde+ncsfm*mcol                                      10s3s22
        call enough('genmatn3. 16',bc,ibc)
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
     d' genmatn3. 10')
        end if                                                          10s3s22
       end if                                                           10s3s22
       if(icase1234.eq.1)then                                           9d30s22
c (((ab)c)d)e
        itmpad=ibcoff                                                   9d30s22
        ibcoff=itmpad+ncsfb*mcsf                                        10d6s22
        call enough('genmatn3. 17',bc,ibc)
        if(iwpk2.lt.0)then                                              9d30s22
         itrial=ibcoff                                                  9d30s22
         ibcoff=itrial+ncsf2*ncsfk                                      9d30s22
         call enough('genmatn3. 18',bc,ibc)
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
     d' genmatn3. 11')
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpad),ncsfb,       10d6s22
     $       vec,ndim,ff,xout,ncsfb,                                    10s3s22
     d' genmatn3. 12')
        if(iwpk2.lt.0)ibcoff=itrial                                     9d30s22
       else if(icase1234.eq.2)then                                           9d30s22
c ((ab)c)(de)
        call dgemm('n','n',ncsfb,mcol,ncsf2,fmul,bc(itmpabc),ncsfb,     10d7s22
     $       bc(itmpde),ncsf2,ff,xout,ncsfb,                            10s3s22
     d' genmatn3. 13')
       else if(icase1234.eq.3)then                                           9d30s22
c (ab)((cd)e)
        call dgemm('n','n',ncsfb,mcol,ncsfm,fmul,bc(itmpab),ncsfb,      10d7s22
     $       bc(itmpcde),ncsfm,ff,xout,ncsfb,                           10s3s22
     d' genmatn3. 14')
       else if(icase1234.eq.4)then                                      10s3s22
c ((ab)(cd))e                                                           10s3s22
        call dgemm('n','n',ncsfb,mcol,ncsfm,fmul,bc(itmpab),ncsfb,      10d7s22
     $       bc(itmpcde),ncsfm,ff,xout,ncsfb,                           10s3s22
     d' genmatn3. 15')
       else if(icase1234.eq.5)then                                           9d30s22
c (ab)(c(de))
        call dgemm('n','n',ncsfb,mcol,ncsfm,fmul,bc(itmpab),ncsfb,      10d7s22
     $       bc(itmpccde),ncsfm,ff,xout,ncsfb,                          10s3s22
     d' genmatn3. 16')
       else if(icase1234.eq.6)then                                           9d30s22
c (a(bc))(de)
        call dgemm('n','n',ncsfb,mcol,ncsf2,fmul,bc(itmpabc),ncsfb,     10d14s22
     $       bc(itmpde),ncsf2,ff,xout,ncsfb,                            10d14s22
     d' genmatn3. 17')
       else if(icase1234.eq.7)then                                           9d30s22
c ((a(bc))d)e
        itmpx=ibcoff                                                    10s3s22
        ibcoff=itmpx+ncsfb*mcsf                                         10d6s22
        call enough('genmatn3. 19',bc,ibc)
        if(iwpk2.lt.0)then                                              10s3s22
         itrial=ibcoff                                                  10s3s22
         ibcoff=itrial+ncsf2*ncsfk                                      10s3s22
         call enough('genmatn3. 20',bc,ibc)
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
     d' genmatn3. 18')
        if(iwpk2.lt.0)ibcoff=itrial                                     10s3s22
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpx),ncsfb,        10d7s22
     $       vec,ndim,ff,xout,ncsfb,                                    10s3s22
     d' genmatn3. 19')
       else if(icase1234.eq.8)then                                           9d30s22
c (a((bc)d))e
        itmpad=ibcoff                                                   9d30s22
        ibcoff=itmpad+ncsfb*mcsf                                        10d7s22
        call enough('genmatn3. 21',bc,ibc)
        call dgemm('n','n',ncsfb,mcsf,kcsf,1d0,bc(iwpb1),ncsfb,         11d4s22
     $        bc(itmpbcd),kcsf,0d0,bc(itmpad),ncsfb,                    11d4s22
     d' genmatn3. 20')
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpad),ncsfb,       10d7s22
     $       vec,ndim,ff,xout,ncsfb,                                    10s3s22
     d' genmatn3. 21')
       else if(icase1234.eq.9)then                                           9d30s22
c  a(((bc)d)e)
        itmpx=ibcoff                                                    10d14s22
        ibcoff=itmpx+kcsf*mcol                                          11d4s22
        call enough('genmatn3. 22',bc,ibc)
        call dgemm('n','n',kcsf,mcol,mcsf,1d0,bc(itmpbcd),kcsf,         11d4s22
     $       vec,ndim,0d0,bc(itmpx),kcsf,                               11d4s22
     d' genmatn3. 22')
        call dgemm('n','n',ncsfb,mcol,kcsf,fmul,bc(iwpb1),ncsfb,        11d4s22
     $       bc(itmpx),kcsf,ff,xout,ncsfb,                              11d4s22
     d' genmatn3. 23')
       else if(icase1234.eq.10)then                                           9d30s22
c a((bc)*(de))
        itmpx=ibcoff                                                    10d14s22
        ibcoff=itmpx+kcsf*mcol                                          11d4s22
        call enough('genmatn3. 24',bc,ibc)
        call dgemm('n','n',kcsf,mcol,ncsf2,1d0,bc(jtmpbc),ncsf1,        11d4s22
     $       bc(itmpde),ncsf2,0d0,bc(itmpx),kcsf,                       11d4s22
     d' genmatn3. 24')
        if(iwpb1.lt.0)then                                                12d5s20
         itrial=ibcoff                                                  10d14s22
         ibcoff=itrial+ncsfb*ncsf1                                      10d14s22
         call enough('genmatn3. 25',bc,ibc)
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsfb,ncsf1)                                              10d14s22
        else                                                              12d5s20
         itrial=iwpb1                                                   10d14s22
        end if                                                          10d14s22
        call dgemm('n','n',ncsfb,mcol,kcsf,fmul,bc(itrial),ncsfb,       11d4s22
     $       bc(itmpx),kcsf,ff,xout,ncsfb,                              11d4s22
     d' genmatn3. 25')
       else if(icase1234.eq.11)then                                           9d30s22
c a(b((cd)e))
        itmpbe=ibcoff                                                   9d30s22
        ibcoff=itmpbe+ncsfm*mcol                                        9d30s22
        call enough('genmatn3. 26',bc,ibc)
        if(iwpk1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpcde),   9d30s22
     $       ncsfm,ncsfm,bc(itmpbe),ncsf1,ncsf1,mcol,0d0,1d0)            3d25s22
         jtmpbe=itmpbe+jcsf0-1                                          11d4s22
        else                                                            9d30s22
         jwpk1=iwpk1+jcsf0-1                                            11d4s22
         call dgemm('n','n',kcsf,mcol,ncsfm,1d0,bc(jwpk1),ncsf1,        11d4s22
     $        bc(itmpcde),ncsfm,0d0,bc(itmpbe),ncsf1,                   9d30s22
     d' genmatn3. 26')
         jtmpbe=itmpbe                                                  11d4s22
        end if                                                          9d30s22
        call dgemm('n','n',ncsfb,mcol,kcsf,fmul,bc(iwpb1),ncsfb,        11d4s22
     $        bc(jtmpbe),ncsf1,ff,xout,ncsfb,                           11d4s22
     d' genmatn3. 27')
       else if(icase1234.eq.12)then                                           9d30s22
c (a(b(cd)))e
        itmpad=ibcoff                                                   9d30s22
        ibcoff=itmpad+mcsf*ncsfb                                        10d7s22
        call enough('genmatn3. 27',bc,ibc)
        if(iwpb1.lt.0)then                                              9d30s22
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(jtmpbcd),   11d4s22
     $       ncsf1,kcsf,bc(itmpad),ncsfb,ncsfb,mcsf,0d0,1d0)            11d4s22
        else                                                            9d30s22
         call dgemm('n','n',ncsfb,mcsf,kcsf,1d0,bc(iwpb1),ncsfb,        11d4s22
     $        bc(jtmpbcd),ncsf1,0d0,bc(itmpad),ncsfb,                   9d30s22
     d' genmatn3. 28')
        end if                                                          9d30s22
        call dgemm('n','n',ncsfb,mcol,mcsf,fmul,bc(itmpad),ncsfb,vec,   10d7s22
     $       ndim,ff,xout,ncsfb,                                        10s3s22
     d' genmatn3. 29')
       else if(icase1234.eq.13)then                                           9d30s22
c a((b(cd))e)
        itmpx=ibcoff                                                    10d14s22
        ibcoff=itmpx+kcsf*mcol                                          11d4s22
        call enough('genmatn3. 28',bc,ibc)
        call dgemm('n','n',kcsf,mcol,mcsf,1d0,bc(jtmpbcd),ncsf1,        11d4s22
     $       vec,ndim,0d0,bc(itmpx),kcsf,                               11d4s22
     d' genmatn3. 30')
        if(iwpb1.lt.0)then                                                12d5s20
         itrial=ibcoff                                                  10d14s22
         ibcoff=itrial+ncsfb*ncsf1                                      10d14s22
         call enough('genmatn3. 29',bc,ibc)
         nusedi=ibc(-iwpb1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpb1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpb1+1                                                   12d5s20
         icmp2=icmp1+ncsf1                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),      9d30s22
     $        ncsfb,ncsf1)                                              10d14s22
        else                                                              12d5s20
         itrial=iwpb1                                                   10d14s22
        end if                                                          10d14s22
        call dgemm('n','n',ncsfb,mcol,kcsf,fmul,bc(itrial),ncsfb,       11d4s22
     $       bc(itmpx),kcsf,ff,xout,ncsfb,                              11d4s22
     d' genmatn3. 31')
       else if(icase1234.eq.14)then                                           9d30s22
c a(b(c(de)))
        itmpx=ibcoff                                                    10s3s22
        ibcoff=itmpx+ncsf1*mcol                                         10s3s22
        call enough('genmatn3. 30',bc,ibc)
        if(iwpk1.lt.0)then                                              10s3s22
         nusedi=ibc(-iwpk1)/2                                             12d5s20
         if(2*nusedi.ne.ibc(-iwpk1))nusedi=nusedi+1                       12d5s20
         icmp1=-iwpk1+1                                                   12d5s20
         icmp2=icmp1+ncsfm                                                12d5s20
         icmp3=icmp2+nusedi                                               12d5s20
         call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpccde),  10s3s22
     $       ncsfm,ncsfm,bc(itmpx),ncsf1,ncsf1,mcol,0d0,1d0)            10s3s22
         jtmpx=itmpx+jcsf0-1                                            11d4s22
        else                                                            10s3s22
         jwpk1=iwpk1+jcsf0-1                                            11d4s22
         call dgemm('n','n',kcsf,mcol,ncsfm,1d0,bc(jwpk1),ncsf1,        11d4s22
     $        bc(itmpccde),ncsfm,0d0,bc(itmpx),ncsf1,                   10s3s22
     d' genmatn3. 32')
         jtmpx=itmpx                                                    11d4s22
        end if                                                          10s3s22
        call dgemm('n','n',ncsfb,mcol,kcsf,fmul,bc(iwpb1),ncsfb,        11d4s22
     $        bc(jtmpx),ncsf1,ff,xout,ncsfb,                            10s3s22
     d' genmatn3. 33')
       else                                                             9d30s22
        write(6,*)('unknown case1234 '),icase1234
        stop 'case1234'
       end if                                                           9d30s22
      end if                                                            9d30s22
      ibcoff=ibcoffo                                                    10s3s22
      return
      end
