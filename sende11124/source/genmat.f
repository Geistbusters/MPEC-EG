c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genmat(ncsfb,ncsfm,ncsfk,ncsf1,iwpb1,iwpk1,ncsf2,iwpb2,11d13s20
     $     iwpk2,iprod,bc,ibc)                                          11d10s22
      implicit real*8 (a-h,o-z)                                         11d13s20
      integer*8 nwork1,nwork2,nwork3,nwork4,nwork5,nbest                2d17s21
c
c     multiply iwpb1*iwpk1*iwpb2*iwpk2 and return in iprod              11d13s20
c                                                                       11d13s20
      include "common.store"                                            11d13s20
      data icall/0/
      icall=icall+1
      iprod=ibcoff                                                      10d27s22
      ibcoff=iprod+ncsfb*ncsfk                                          10d27s22
      call enough('genmat.0',bc,ibc)                                           11d15s22
      ibcoffo=ibcoff                                                    10s3s22
      if(iwpb1.lt.0)then                                                2d10s23
       nn=ibc(-iwpb1)                                                   2d10s23
       n123=nn*ncsfm                                                    2d10s23
       if(iwpk1.lt.0)n123=n123/2                                        2d10s23
       n124=nn*ncsf2
       n125=nn*ncsfk
      else                                                              2d10s23
       if(iwpk1.lt.0)then                                               2d10s23
        nn=ibc(-iwpk1)                                                  2d10s23
        n123=nn*ncsfb                                                   2d10s23
       else                                                             2d10s23
        n123=ncsfb*ncsf1*ncsfm                                            10d6s22
       end if                                                           2d10s23
       n124=ncsfb*ncsf1*ncsf2
       n125=ncsfb*ncsf1*ncsfk
      end if                                                            2d10s23
      if(iwpb2.lt.0)then                                                2d10s23
       nn=ibc(-iwpb2)                                                   2d10s23
       n134=ncsfb*nn                                                    2d10s23
       n234=ncsfm*nn
       if(iwpk1.lt.0)n234=n234/2
       n345=nn*ncsfk                                                    2d10s23
      else                                                              2d10s23
       n134=ncsfb*ncsfm*ncsf2
       n345=ncsfm*ncsf2*ncsfk                                           2d10s23
       if(iwpk1.lt.0)then
        nn=ibc(-iwpk1)
        n234=nn*ncsf2
       else                                                             2d10s23
        n234=ncsf1*ncsfm*ncsf2
       end if                                                           2d10s23
      end if                                                            2d10s23
      n135=ncsfb*ncsfm*ncsfk                                            10d6s22
      if(iwpk2.lt.0)then                                                2d10s23
       nn=ibc(-iwpk2)
       n145=ncsfb*nn
       n245=ncsf1*nn
      else                                                              2d10s23
       n145=ncsfb*ncsf2*ncsfk                                            10d6s22
       n245=ncsf1*ncsf2*ncsfk                                            10d6s22
      end if                                                            2d10s23
      if(iwpk1.lt.0)then                                                2d10s23
       nn=ibc(-iwpk1)                                                   2d10s23
       n235=nn*ncsfk
      else
       n235=ncsf1*ncsfm*ncsfk                                            10d6s22
      end if
      ntry1=n134+n145                                                   10d6s22
      ntry2=n345+n135                                                   10d6s22
      if(ntry1.le.ntry2)then                                            10d6s22
       nuse=ntry1                                                       10d6s22
       icase1=1                                                         10d6s22
      else                                                              10d6s22
       nuse=ntry2                                                       10d6s22
       icase1=2                                                         10d6s22
      end if
      nsofar=n123+nuse                                                  10d6s22
      ntry1=n124+n145                                                   10d6s22
      ntry2=n245+n125                                                   10d6s22
      if(ntry1.le.ntry2)then                                            10d6s22
       nuse=ntry1                                                       10d6s22
       icase2=3                                                         10d6s22
      else                                                              10d6s22
       nuse=ntry2                                                       10d6s22
       icase2=4                                                         10d6s22
      end if                                                            10d6s22
      ntry=nuse+n234                                                    10d6s22
      if(ntry.le.nsofar)then                                            10d6s22
       icase1=icase2                                                    10d6s22
       nsofar=ntry                                                      10d6s22
      end if                                                            10d6s22
      ntry=n345+n235+n125                                               10d6s22
      if(ntry.le.nsofar)then                                            10d6s22
       icase1=5                                                         10d6s22
      end if                                                            10d6s22
      if(icase1.le.2)then                                               10d6s22
       itmpab=ibcoff                                                    10d6s22
       ibcoff=itmpab+ncsfb*ncsfm                                        10d6s22
       call enough('genmat.  8',bc,ibc)
       call prodn(iwpb1,iwpk1,ncsfb,ncsfm,ncsf1,bc(itmpab),bc,ibc,      2d13s23
     $      1d0,0d0)                                                    2d13s23
      end if                                                            10d6s22
      if(icase1.ge.3.and.icase1.le.4)then                               10d6s22
       itmpbc=ibcoff                                                    10d6s22
       ibcoff=itmpbc+ncsf1*ncsf2                                        10d6s22
       call enough('genmat. 10',bc,ibc)
       call prodn(iwpk1,iwpb2,ncsf1,ncsf2,ncsfm,bc(itmpbc),bc,ibc,      2d13s23
     $      1d0,0d0)                                                    2d13s23
      end if                                                            10d6s22
      if(icase1.eq.2.or.icase1.eq.5)then                                10d6s22
       itmpcd=ibcoff                                                    10d6s22
       ibcoff=itmpcd+ncsfm*ncsfk                                        10d6s22
       call prodn(iwpb2,iwpk2,ncsfm,ncsfk,ncsf2,bc(itmpcd),bc,ibc,      2d13s23
     $      1d0,0d0)                                                    2d13s23
      end if                                                            10d6s22
      if(icase1.eq.1)then                                               10d6s22
       itmpabc=ibcoff                                                   10d6s22
       ibcoff=itmpabc+ncsfb*ncsf2                                       10d6s22
       call enough('genmat. 14',bc,ibc)
       call prodn(itmpab,iwpb2,ncsfb,ncsf2,ncsfm,bc(itmpabc),bc,ibc,    2d13s23
     $      1d0,0d0)                                                    2d13s23
       call prodn(itmpabc,iwpk2,ncsfb,ncsfk,ncsf2,bc(iprod),bc,ibc,     2d13s23
     $      1d0,0d0)                                                    2d13s23
      else if(icase1.eq.2)then                                          10d6s22
       call dgemm('n','n',ncsfb,ncsfk,ncsfm,1d0,bc(itmpab),ncsfb,       10d6s22
     $      bc(itmpcd),ncsfm,0d0,bc(iprod),ncsfb,                       10d27s22
     d' genmat.  6')
      else if(icase1.eq.3)then                                          10d6s22
       itmpabc=ibcoff                                                   10d6s22
       ibcoff=itmpabc+ncsfb*ncsf2                                       2d10s23
       call enough('genmat. 17',bc,ibc)
       call prodn(iwpb1,itmpbc,ncsfb,ncsf2,ncsf1,bc(itmpabc),bc,ibc,    2d13s23
     $      1d0,0d0)                                                    2d13s23
       call prodn(itmpabc,iwpk2,ncsfb,ncsfk,ncsf2,bc(iprod),bc,ibc,     2d13s23
     $      1d0,0d0)                                                    2d13s23
      else if(icase1.eq.4)then                                          10d6s22
       itmpbd=ibcoff                                                    2d10s23
       ibcoff=itmpbd+ncsf1*ncsfk                                        2d10s23
       call enough('genmat {bc}d ',bc,ibc)                              2d10s23
       call prodn(itmpbc,iwpk2,ncsf1,ncsfk,ncsf2,bc(itmpbd),bc,ibc,     2d13s23
     $      1d0,0d0)                                                    2d13s23
       call prodn(iwpb1,itmpbd,ncsfb,ncsfk,ncsf1,bc(iprod),bc,ibc,      2d13s23
     $      1d0,0d0)                                                    2d13s23
       ibcoff=itmpbd                                                    2d10s23
      else if(icase1.eq.5)then                                          10d6s22
       itmpbd=ibcoff                                                    10d6s22
       ibcoff=itmpbd+ncsf1*ncsfk                                        10d6s22
       call enough('genmat. 19',bc,ibc)
       call prodn(iwpk1,itmpcd,ncsf1,ncsfk,ncsfm,bc(itmpbd),bc,ibc,     2d13s23
     $      1d0,0d0)                                                    2d13s23
       call prodn(iwpb1,itmpbd,ncsfb,ncsfk,ncsf1,bc(iprod),bc,ibc,      2d13s23
     $      1d0,0d0)                                                    2d13s23
      else                                                              10d6s22
       write(6,*)('case '),icase1
       stop 'genmat'
      end if                                                            10d6s22
      ibcoff=ibcoffo                                                    10d6s22
      return
      end
