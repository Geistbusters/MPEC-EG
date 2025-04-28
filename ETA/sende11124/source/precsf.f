c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine precsf(nopx,nsymb,irel,ism,norb,                       12d28s19
     $     multh,ibasisp,icsfp,nfcnp,nctp,ixw1p,ixw2p,iptrbitp,ipoint2, 11d9s22
     $     bc,ibc)                                                      11d9s22
      implicit real*8 (a-h,o-z)
      integer*8 icsfxv,ipack8                                           7d30s22
      integer*4 ipack4(2)                                               7d30s22
      integer*2 ipack2(4)                                               8d2s22
      equivalence (ipack8,ipack4),(ipack4,ipack2)                       7d30s22
      include "common.cas"                                              12d28s19
      dimension irel(*),ism(*),ibasis(8),multh(8,8)                     12d28s19
     $     ,nfcn(8),nct(8),nfill(4),nhole(3)                            4d19s23
      include "common.store"                                            12d28s19
      ii=0                                                              12d28s19
      do isb=1,nsymb                                                    12d28s19
       do i=1,iacto(isb)                                                12d28s19
        ii=ii+1                                                         12d28s19
        irel(ii)=i                                                      12d28s19
        ism(ii)=isb                                                     12d28s19
       end do                                                           12d28s19
      end do                                                            12d28s19
      mdoo=0                                                            12d28s19
      mdon=norb+1                                                       12d28s19
      nodd=0                                                            12d28s19
      nevn=0d0                                                          12d28s19
      ispino=0                                                          12d28s19
      ispine=0                                                          12d28s19
      nopex=0
      nopox=0                                                           12d28s19
      nec=isinfo(3,1)                                                   12d28s19
      ismult=isinfo(2,1)                                                12d28s19
      spin=0.5d0*(dfloat(ismult)-1d0)                                   2d19s19
      iuniq=ibcoff                                                      7d30s22
      ipoint2=iuniq+nstate                                              7d30s22
      ibcoff=ipoint2+nstate                                             7d30s22
      call enough('precsf.  1',bc,ibc)
      nuniq=0                                                           7d30s22
      ierr=0                                                            12d28s19
      do i=1,nstate
       im=i-1                                                           7d30s22
       do j=0,nuniq-1                                                   7d30s22
        ipack8=ibc(iuniq+j)                                             7d30s22
        if(isinfo(3,i).eq.ipack2(1).and.isinfo(2,i).eq.ipack2(2))then   7d30s22
         ibc(ipoint2+im)=j                                              7d30s22
         go to 57                                                       7d30s22
        end if                                                          7d30s22
       end do                                                           7d30s22
       ipack2(1)=isinfo(3,i)                                            7d30s22
       ipack2(2)=isinfo(2,i)                                            7d30s22
       ibc(iuniq+nuniq)=ipack8                                          7d30s22
       ibc(ipoint2+im)=nuniq                                            8d2s22
       nuniq=nuniq+1                                                    7d30s22
   57  continue                                                         7d30s22
      end do                                                            7d30s22
      icsfp=ibcoff                                                      8d2s22
      nfcnp=icsfp+nuniq                                                 8d2s22
      ibasisp=nfcnp+nuniq*nsymb                                         8d2s22
      iptrbitp=ibasisp+nuniq*nsymb                                         8d2s22
      ixw1p=iptrbitp+nuniq                                              8d2s22
      ixw2p=ixw1p+nuniq                                                 8d2s22
      nctp=ixw2p+nuniq                                                  8d2s22
      ibcoff=nctp+nuniq*nsymb                                           8d2s22
      call enough('precsf.  2',bc,ibc)
      do iq=0,nuniq-1                                                   8d2s22
       ipack8=ibc(iuniq+iq)                                             8d2s22
       nec=ipack2(1)                                                    8d2s22
       ismult=ipack2(2)                                                 8d2s22
       spin=0.5d0*(dfloat(ismult)-1d0)                                   2d19s19
       icsf=ibcoff                                                       6d5s19
       idet=icsf+norb                                                    11d19s20
       ibcoff=idet+norb                                                  11d19s20
       call enough('precsf.  3',bc,ibc)
       call gencsf1(nec,spin,norb,nos,nod,mdon,mdoo,ibc(icsf),nopx,      12d28s19
     $     ixsoo,2,ibc(idet),0)                                         4d13s21
       ipack2(3)=mdoo                                                   8d2s22
       ipack2(4)=mdon                                                   8d2s22
       ibc(iuniq+iq)=ipack8                                             8d2s22
       mdoop=mdoo+1                                                      10d28s20
       ibcoff=icsf+mdoop-mdon                                            4d12s21
       icsf2=ibcoff                                                      10d28s20
       ibcoff=icsf2+4*(mdoop-mdon)                                       4d12s21
       xms=spin                                                          12d28s19
       if(mod(ixsoo,2).eq.0)then                                         2d14s20
        ixsoo=max(ixsoo,2)                                               2d14s20
       else                                                              2d14s20
        ixsoo=max(ixsoo,3)                                               2d14s20
       end if                                                            2d14s20
       nfcnpx=ibcoff                                                     9d20s20
       ibcoff=nfcnpx+4*nsymb*mdoop                                       10d28s20
       idorb=ibcoff
       nod8=nsymb*nod/8                                                  7d15s19
       if(nod8*8.lt.nod*nsymb)nod8=nod8+1                                7d16s19
       isorb=idorb+nod8                                                  6d12s19
       nos=nsymb*2*nos*8                                                 7d15s19
       nos8=nos/8                                                        6d12s19
       if(nos8*8.lt.nos)nos8=nos8+1                                      6d12s19
       iptr=isorb+nos8                                                   6d12s19
       ibasis(1)=iptr+mdoop*4*nsymb                                      4d12s21
       do i=1,4                                                         4d19s23
        nfill(i)=0                                                       12d28s19
       end do                                                            12d28s19
       do i=1,3                                                          12d28s19
        nhole(i)=0                                                       12d28s19
       end do                                                            12d28s19
       call enough('precsf.  4',bc,ibc)
       call gencsf2(nec,spin,norb,ibc(idorb),ibc(isorb),ibc(iptr),       6d4s19
     $     nsymb,mdon,mdoo,ibc(ibasis(1)),multh,ism,ibc(icsf),nfcn,nct, 7d15s19
     $     nod,nos,nopx,nhole,ihole,nfill,ifill,1,0,0,                  12d28s19
     $     ismult,idum,idum,idum,idum,idum,idum,idum,idum,0,            12d28s19
     $     nodu,nosu,nsymb,idum,idum,idum,idum,idum,idum,ibc(nfcnpx),   10d28s20
     $     iiiidum,idum,ibc(icsf2),1,bc,ibc)                            11d14s22
       need=0                                                            7d15s19
       jbasis=ibasis(1)                                                  7d15s19
       infcnp=nfcnp+iq*nsymb-1                                          8d2s22
       jbasisp=ibasisp+iq*nsymb-1                                       8d2s22
       inctp=nctp+iq*nsymb-1                                            8d2s22
       do isb=1,nsymb                                                    7d15s19
        ibc(infcnp+isb)=nfcn(isb)                                       8d2s22
        ibc(inctp+isb)=nct(isb)                                         8d2s22
        nhere=(nfcn(isb)*3)+mod(nfcn(isb),2)                             7d15s19
        need=need+nhere                                                  7d15s19
        ibasis(isb)=jbasis                                               7d15s19
        ibc(jbasisp+isb)=jbasis                                         8d2s22
        jbasis=jbasis+(nhere/2)                                          7d15s19
       end do                                                            7d15s19
       need8=need/2                                                      7d15s19
       ibcoff=ibasis(1)+need8                                            7d15s19
       iptrbit=ibcoff                                                    5d7s20
       ibcoff=iptrbit+nsymb*mdoop                                        5d7s20
       ibc(iptrbitp+iq)=iptrbit                                         8d2s22
       icsfpd=ibcoff                                                     6d11s19
       ibcoff=icsfpd+mdoop-mdon                                          4d12s21
       icsf2=ibcoff                                                      8d2s19
       ibcoff=icsf2+4*(mdoop-mdon)                                       4d12s21
       call enough('precsf.  5',bc,ibc)
       call int1tobit(ibc(iptr),ibc(iptrbit),ibc(idorb),ibc(isorb),      5d7s20
     $     idorbbit,isorbbit,mdoop,nsymb,nec,bc,ibc)                    11d10s22
       ibcb4=ibcoff                                                      4d12s21
       call stepwisecsf(ixsoo,spin,xms,dum,idum,idum,1,idata,bc,ibc)    11d9s22
       call enough('precsf.  6',bc,ibc)
       call gencsf3(nec,spin,mdon,mdoo,norb,ibc(icsfpd),ismult,          8d2s19
     $     ibc(icsf2),.false.,ixw1,ixw2,ibc(icsf),1,ibcb4,icsfxv,bc,ibc)11d11s22
       ibc(icsfp+iq)=icsf                                               8d2s22
       ibc(ixw1p+iq)=ixw1                                               8d2s22
       ibc(ixw2p+iq)=ixw2                                               8d2s22
      end do                                                            8d2s22
      return                                                            12d28s19
      end                                                               12d28s19
