c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genhsoa2(iwaveb,iwavek,nspc,mlb2,mlk2,msb2,msk2,aout,  5d18s21
     $     nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,iri,itypex,    8d30s21
     $     nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,
     $     isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,iifmx,    11d15s21
     $     ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,  12d20s20
     $         isstor,isym,ascale,idorel,iorb,l2e,l2es,lri,bc,ibc,n4vso,3d6s24
     $     irori,irev)                                                  3d27s24
      implicit real*8 (a-h,o-z)                                         5d14s21
      integer*1 ipack1(4)                                               5d14s21
      integer*8 ipack8                                                  5d18s21
      equivalence (npack4,ipack1),(ipack8,ipack4)                       5d14s21
      dimension iwaveb(nspc,*),ipack4(2),multh(8,8),irefo(*),           5d18s21
     $     iwavek(nspc,*),ih0a(8,2),i4so2(8,8,8,2),irel(*),ism(*),      1d19s22
     $     aout(*),itypex(4),nbasp(*),nbaspc(*)                         12d20s20
      include "common.store"                                            5d14s21
      data icall/0/
      icall=icall+1
      npack4=iwaveb(6,1)                                                5d18s21
      if(ipack1(2).eq.6)then                                            5d27s21
       nllb=2*ipack1(3)+1                                                 5d14s21
      else if(ipack1(2).eq.2)then                                       5d27s21
       if(ipack1(3).eq.0)then                                           5d27s21
        nllb=1                                                          5d27s21
       else                                                             5d27s21
        nllb=2                                                          5d27s21
       end if                                                           5d27s21
      else                                                              5d27s21
       nllb=1                                                           5d27s21
      end if                                                            5d27s21
      llb=ipack1(3)                                                     5d20s21
      mlba=iabs(mlb2)/2                                                 5d14s21
      mlka=iabs(mlk2)/2                                                 5d14s21
      do i=1,nllb                                                        5d14s21
       npack4=iwaveb(6,i)                                               5d18s21
       if(ipack1(4).eq.mlba)ibr=i                                       5d24s21
       if(ipack1(4).eq.-mlba)ibi=i                                      5d24s21
      end do                                                            5d14s21
      npack4=iwavek(6,1)                                                5d18s21
      if(ipack1(2).eq.6)then                                            5d27s21
       nllk=2*ipack1(3)+1                                                 5d14s21
      else if(ipack1(2).eq.2)then                                       5d27s21
       if(ipack1(3).eq.0)then                                           5d27s21
        nllk=1                                                          5d27s21
       else                                                             5d27s21
        nllk=2                                                          5d27s21
       end if                                                           5d27s21
      else                                                              5d27s21
       nllk=1                                                           5d27s21
      end if                                                            5d27s21
      llk=ipack1(3)                                                     5d20s21
      do i=1,nllk                                                        5d14s21
       npack4=iwavek(6,i)                                               5d18s21
       if(ipack1(4).eq.mlka)ikr=i                                       5d24s21
       if(ipack1(4).eq.-mlka)iki=i                                      5d24s21
      end do                                                            5d14s21
      if(ipack1(2).eq.6)then                                            5d27s21
c
c     atomic case
c
c     use phase convention L -Ml = (-1)^Ml L Ml* for ordinary states
c     and L -Ml = -(-1)^Ml L Ml* for extraordinary states
c     we have stored data for L Ml.
c
       iordb=0                                                           5d20s21
       if(iwaveb(15,1).ne.0)iordb=1                                      5d20s21
       if(mod(llb,2).ne.0)iordb=iordb+1                                     5d19s21
       iordb=mod(iordb,2)                                                  5d19s21
       pbr=1d0                                                           5d14s21
       pbi=1d0                                                           5d14s21
       if(mlb2.lt.0)then                                                 5d14s21
        if(mod(mlba+iordb,2).eq.0)then                                   5d20s21
         pfb=1d0                                                          5d14s21
        else                                                              5d14s21
         pfb=-1d0                                                         5d14s21
        end if                                                           5d14s21
        pbr=pfb                                                          5d14s21
        pbi=-pfb                                                         5d14s21
       end if                                                            5d14s21
       iordk=0                                                           5d20s21
       if(iwavek(15,1).ne.0)iordk=1                                      5d20s21
       if(mod(llk,2).ne.0)iordk=iordk+1                                  5d20s21
       iordk=mod(iordk,2)                                                5d20s21
       pkr=1d0                                                           5d14s21
       pki=1d0                                                           5d14s21
       if(mlk2.lt.0)then                                                 5d14s21
        if(mod(mlka+iordk,2).eq.0)then                                   5d20s21
         pfk=1d0                                                          5d14s21
        else                                                              5d14s21
         pfk=-1d0                                                         5d14s21
        end if                                                           5d14s21
        pkr=pfk                                                          5d14s21
        pki=-pfk                                                         5d14s21
       end if                                                            5d14s21
      else if(ipack1(2).eq.2)then                                       5d27s21
       iordb=0                                                           5d20s21
       pbr=1d0                                                          5d27s21
       pbi=1d0                                                          5d27s21
       if(llb.gt.0)then                                                 3d9s22
        if(mod(llb,2).eq.0)then                                          3d9s22
         if(mlb2.gt.0)pbi=-1d0                                          3d9s22
        else                                                            3d9s22
         pbi=-1d0                                                       3d9s22
         if(mlb2.lt.0)pbr=-1d0                                          3d9s22
        end if                                                          3d9s22
       end if                                                           3d9s22
       if(llb.eq.0.and.(iwaveb(2,1).eq.8.or.iwaveb(2,1).eq.4))then       5d27s21
        iordb=1                                                         5d27s21
       end if                                                           5d27s21
       iordk=0                                                           5d20s21
       pkr=1d0                                                          5d27s21
       pki=1d0                                                          5d27s21
       if(llk.gt.0)then                                                 3d9s22
        if(mod(llk,2).eq.0)then                                          3d9s22
         if(mlk2.gt.0)pki=-1d0                                          3d9s22
        else                                                            3d9s22
         pki=-1d0                                                       3d9s22
         if(mlk2.lt.0)pkr=-1d0                                          3d9s22
        end if                                                          3d9s22
       end if                                                           3d9s22
       if(llk.eq.0.and.(iwavek(2,1).eq.8.or.iwavek(2,1).eq.4))then       5d27s21
        iordk=1                                                         5d27s21
       end if                                                           5d27s21
      else
c
c     general case                                                      5d27s21
c
       npack4=iwaveb(6,1)                                               3d28s22
       iordb=ipack1(3)                                                  3d28s22
       npack4=iwavek(6,1)                                               3d28s22
       iordk=ipack1(3)                                                  3d28s22
       pbr=1d0                                                          5d27s21
       pkr=1d0                                                          5d27s21
       pbi=1d0                                                          5d27s21
       pki=1d0                                                          5d27s21
       ibr=1                                                            3d24s22
       ibi=1                                                            3d24s22
       ikr=1                                                            3d24s22
       iki=1                                                            3d24s22
      end if                                                            5d27s21
      mdoopb=((nec-iwaveb(1,1)+1)/2)+1                                  5d18s21
      nffbr=iwaveb(9,ibr)+iwaveb(4,ibr)                                 5d18s21
      iffbr=iwaveb(10,ibr)+iwaveb(4,ibr)                                5d18s21
      iwbr=iwaveb(13,ibr)+iwaveb(4,ibr)                                 5d18s21
      ipack8=ibc(iwbr)                                                  5d14s21
      ncsftbr=ipack4(1)                                                 5d14s21
      iwbr=iwbr+1                                                       5d14s21
      nffbi=iwaveb(9,ibi)+iwaveb(4,ibi)                                 5d18s21
      iffbi=iwaveb(10,ibi)+iwaveb(4,ibi)                                5d18s21
      iwbi=iwaveb(13,ibi)+iwaveb(4,ibi)                                 5d18s21
      ipack8=ibc(iwbi)                                                  5d14s21
      ncsftbi=ipack4(1)                                                 5d14s21
      iwbi=iwbi+1                                                       5d14s21
      mdoopk=((nec-iwavek(1,1)+1)/2)+1                                  5d18s21
      nffkr=iwavek(9,ikr)+iwavek(4,ikr)                                 5d18s21
      iffkr=iwavek(10,ikr)+iwavek(4,ikr)                                5d18s21
      iwkr=iwavek(13,ikr)+iwavek(4,ikr)                                 5d18s21
      ipack8=ibc(iwkr)                                                  5d14s21
      ncsftkr=ipack4(1)                                                 5d14s21
      iwkr=iwkr+1                                                       5d14s21
      nffki=iwavek(9,iki)+iwavek(4,iki)                                 5d18s21
      iffki=iwavek(10,iki)+iwavek(4,iki)                                5d18s21
      iwki=iwavek(13,iki)+iwavek(4,iki)                                 5d18s21
      ipack8=ibc(iwki)                                                  5d14s21
      ncsftki=ipack4(1)                                                 5d14s21
      iwki=iwki+1                                                       5d14s21
      isymc=multh(iwaveb(2,ibr),iwavek(2,ikr))                          5d18s21
      nhdimb=iwaveb(3,ibr)                                              5d18s21
      nhdimk=iwavek(3,ikr)                                              5d18s21
      ihsor=ibcoff                                                      5d18s21
      ihsoi=ihsor+nhdimb*nhdimk                                           5d14s21
      ibcoff=ihsoi+nhdimb*nhdimk                                        5d18s21
      ihfullr=ibcoff                                                    5d14s21
      ihfulli=ihfullr+nhdimb*nhdimk                                     3d6s24
      ibcoff=ihfulli+nhdimb*nhdimk                                      3d6s24
      call enough('genhsoa2.  1',bc,ibc)
      do i=ihsor,ibcoff-1                                               5d14s21
       bc(i)=0d0                                                        5d14s21
      end do                                                            5d14s21
      nn=nhdimb*nhdimk*2                                                5d18s21
      iboff=1                                                           5d14s21
      ikoff=1                                                           5d14s21
      shift=0d0                                                         5d14s21
      nbb=iwaveb(1,ibr)-1                                               5d18s21
      nkk=iwavek(1,ikr)-1                                               5d18s21
      iptmp=ibcoff                                                      9d2s21
      ibcoff=iptmp+nn                                                   9d2s21
      call enough('genhsoa2.  2',bc,ibc)
      if((ibr.ne.ibi.or.iordb.eq.0).and.(ikr.ne.iki.or.iordk.eq.0))then 5d20s21
       if(l2e.eq.0.or.isymc.eq.l2es)then                                2d17s22
        ibcoff0=ibcoff                                                   10d20s21
        call psisopsi(iwaveb(1,ibr),iwaveb(3,ibr),msb2,isymc,            8d30s21
     $     iwavek(1,ikr),iwavek(3,ikr),msk2,ih0a(1,1),                  1d19s22
     $     i4so2,bc(iptmp),mdon,ism,irel,irefo,norb,nvirt,maxbx,maxbxd, 3d19s24
     $       srh,sr2,multh,nsymb,irori,irw0,irw1,irw2,                  3d19s24
     $     ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr, 11d15s21
     $     iifmx,ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,   12d20s20
     $      ibstor,isstor,isym,ascale,idorel,iorb,l2e,bc,ibc,n4vso)     3d29s24
       pr=pbr*pkr                                                        5d14s21
        if(l2e.eq.0.or.lri.eq.0)then                                    2d25s22
         do i=0,nhdimb*nhdimk-1                                              5d14s21
          bc(ihfullr+i)=bc(iptmp+i)*pr                                    10d20s21
         end do                                                           10d20s21
        else                                                            2d25s22
         do i=0,nhdimb*nhdimk-1                                         2d25s22
          bc(ihfulli+i)=bc(iptmp+i)*pr                                  2d25s22
         end do                                                         2d25s22
        end if                                                          2d25s22
        ibcoff=ibcoff0                                                   10d20s21
       end if                                                           2d17s22
       do i=0,nn-1                                                       5d18s21
        bc(ihsor+i)=0d0                                                  5d14s21
       end do                                                            5d14s21
      end if                                                            5d20s21
      if((ibr.ne.ibi.or.iordb.eq.0).and.(ikr.ne.iki.or.iordk.eq.1))then 5d20s21
       isymc=multh(iwaveb(2,ibr),iwavek(2,iki))                         5d18s21
       if(l2e.eq.0.or.isymc.eq.l2es)then                                2d17s22
        ibcoff0=ibcoff                                                   10d20s21
        call psisopsi(iwaveb(1,ibr),iwaveb(3,ibr),msb2,isymc,            8d30s21
     $     iwavek(1,iki),iwavek(3,ikr),msk2,ih0a(1,2),                  1d19s22
     $     i4so2,bc(iptmp),mdon,ism,irel,irefo,                         10d22s21
     $     norb,nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,2,irw0,irw1,irw2,10d8s21
     $      ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,11d15s21
     $      iifmx,ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,  12d20s20
     $      ibstor,isstor,isym,ascale,idorel,iorb,l2e,bc,ibc,n4vso)     3d29s24
        pr=pbr*pki                                                        5d14s21
        if(l2e.eq.0.or.lri.eq.1)then                                    2d25s22
         do i=0,nhdimb*nhdimk-1                                              5d14s21
          bc(ihfullr+i)=bc(ihfullr+i)-bc(iptmp+i)*pr                      10d20s21
         end do                                                           10d20s21
        else                                                            2d25s22
         do i=0,nhdimb*nhdimk-1                                         2d25s22
          bc(ihfulli+i)=bc(ihfulli+i)+bc(iptmp+i)*pr                    2d25s22
         end do                                                         2d25s22
        end if                                                          2d25s22
        ibcoff=ibcoff0                                                   10d20s21
       end if                                                           2d17s22
       do i=0,nn-1                                                      5d18s21
        bc(ihsor+i)=0d0                                                  5d14s21
       end do                                                            5d14s21
      end if                                                            5d14s21
      if((ibr.ne.ibi.or.iordb.eq.1).and.(ikr.ne.iki.or.iordk.eq.0))then 5d20s21
       isymc=multh(iwaveb(2,ibi),iwavek(2,ikr))                         5d18s21
       if(l2e.eq.0.or.isymc.eq.l2es)then                                2d17s22
        ibcoff0=ibcoff                                                   10d20s21
        call psisopsi(iwaveb(1,ibi),iwaveb(3,ibr),msb2,isymc,            8d30s21
     $     iwavek(1,ikr),iwavek(3,ikr),msk2,ih0a(1,2),                  1d19s22
     $     i4so2,bc(iptmp),mdon,ism,irel,irefo,                         10d22s21
     $     norb,nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,2,irw0,irw1,irw2,10d8s21
     $      ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,11d15s21
     $      iifmx,ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,  12d20s20
     $      ibstor,isstor,isym,ascale,idorel,iorb,l2e,bc,ibc,n4vso)     3d29s24
        pr=pbi*pkr                                                        5d14s21
        if(l2e.eq.0.or.lri.eq.1)then                                    2d25s22
         do i=0,nhdimb*nhdimk-1                                           10d20s21
          bc(ihfullr+i)=bc(ihfullr+i)+bc(iptmp+i)*pr                      10d20s21
         end do                                                           10d20s21
        else                                                            2d25s22
         do i=0,nhdimb*nhdimk-1                                           10d20s21
          bc(ihfulli+i)=bc(ihfulli+i)-bc(iptmp+i)*pr                    2d25s22
         end do                                                           10d20s21
        end if                                                          2d25s22
        ibcoff=ibcoff0                                                   10d20s21
       end if                                                           2d17s22
       do i=0,nn-1                                                      5d18s21
        bc(ihsor+i)=0d0                                                  5d14s21
       end do                                                            5d14s21
      end if                                                            5d14s21
      if((ibr.ne.ibi.or.iordb.eq.1).and.(ikr.ne.iki.or.iordk.eq.1))then 5d20s21
       isymc=multh(iwaveb(2,ibi),iwavek(2,iki))                         5d18s21
       if(l2e.eq.0.or.isymc.eq.l2es)then                                2d17s22
        ibcoff0=ibcoff                                                   10d20s21
        call psisopsi(iwaveb(1,ibi),iwaveb(3,ibr),msb2,isymc,            8d30s21
     $     iwavek(1,iki),iwavek(3,ikr),msk2,ih0a(1,1),                  1d19s22
     $     i4so2,bc(iptmp),mdon,ism,irel,irefo,                         10d22s21
     $     norb,nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,1,irw0,irw1,irw2,10d8s21
     $      ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,11d15s21
     $      iifmx,ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,  12d20s20
     $      ibstor,isstor,isym,ascale,idorel,iorb,l2e,bc,ibc,n4vso)     3d29s24
        pr=pbi*pki                                                        5d14s21
        if(l2e.eq.0.or.lri.eq.0)then                                    2d25s22
         do i=0,nhdimb*nhdimk-1                                           10d20s21
          bc(ihfullr+i)=bc(ihfullr+i)+bc(iptmp+i)*pr                      10d20s21
         end do                                                           10d20s21
        else                                                            2d25s22
         do i=0,nhdimb*nhdimk-1                                           10d20s21
          bc(ihfulli+i)=bc(ihfulli+i)+bc(iptmp+i)*pr                      10d20s21
         end do                                                           10d20s21
        end if                                                          2d25s22
        ibcoff=ibcoff0                                                   10d20s21
       end if                                                           2d17s22
      end if                                                            5d14s21
      ff=1d0
      if(ibr.ne.ibi)ff=srh                                              5d18s21
      if(ikr.ne.iki)ff=ff*srh                                           5d18s21
      sz=0d0                                                            5d18s21
      szr=0d0                                                           5d19s21
      do i=0,nhdimb*nhdimk-1                                            5d18s21
       ip=i+1                                                           5d18s21
       aout(ip)=bc(ihfullr+i)*ff                                        5d18s21
       sz=sz+bc(ihfulli+i)**2                                           5d18s21
       szr=szr+bc(ihfullr+i)**2                                           5d18s21
      end do                                                            5d18s21
      if(l2e.ne.0)then                                                  2d25s22
       do i=0,nhdimb*nhdimk-1                                            5d18s21
        ip=i+1+nhdimb*nhdimk                                            2d25s22
        aout(ip)=bc(ihfulli+i)*ff                                        5d18s21
       end do                                                           2d25s22
      end if                                                            2d25s22
      ibcoff=ihsor                                                      3d15s22
      sz=sqrt(sz/dfloat(nhdimb*nhdimk))                                 5d18s21
      szr=sqrt(szr/dfloat(nhdimb*nhdimk))                                 5d18s21
      call dws_bcast(sz,1)                                              5d18s21
      iri=0                                                             5d19s21
      if(sz.gt.1d-10)iri=1                                              3d6s24
      return                                                            5d14s21
      end                                                               5d14s21
