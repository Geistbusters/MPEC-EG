c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cont2ecsf(idorb,isorb,mdon,mdoo,ibasis,iptr,ncsf,      8d14s19
     $     ncsf2,nfcn,                                                  8d14s19
     $     nct,nec,icsfpd,multh,ivintref,irwt,nsymb,thrs,               11d19s20
     $     nct2,nct2c,lprint,iff2,iptrto,nff2,nfdat,ndetc,nff22,icsfxv, 4d14s21
     $     ipaddout,bc,ibc,nder,idvi,idvmt,idervcv)                     8d7s24
      implicit real*8 (a-h,o-z)
      external second
      integer*1 idorb(*),isorb(*)                                       11d18s20
      integer*8 icsfxv,idvmt(2,*)                                       9d1s23
      logical ldebug,lprint                                             12d12s19
      dimension ibasis(*),ncsf(*),nfcn(*),icsfpd(*),multh(8,8),
     $     ivintref(*),irwt(*),nct2(*),nct(*),                          11d19s20
     $     ivintr(8),nok(3,2),ncsf2(4,*),nff22(mdoo+1,2,nsymb),         12d18s20
     $     nct2c(4,8),ivcc(4),nff2(mdoo+1,nsymb,8),joff(4),joff0(4),    1d12s23
     $     nfdat(5,4,8),ndetc(*),idvi(*),idervcv(2,8)                   8d7s24
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.mrci"
      include "common.store"
      include "common.print"                                            1d5s20
      if(iprtr(9).eq.0)then                                             1d5s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d5s20
       ldebug=.true.                                                    1d5s20
      end if                                                            1d5s20
      if(ldebug)write(6,*)('hi, my name is cont2ecsf ')                 12d12s19
      joff=0                                                            11d16s20
      koff=1                                                            11d16s20
      do isb=1,nsymb                                                    11d16s20
       joff0=joff                                                       11d16s20
       do l=1,4                                                         11d19s20
        nct2c(l,isb)=0                                                  11d19s20
       end do                                                           11d19s20
       do nclo=mdon,mdoo                                                11d16s20
        nclop=nclo+1                                                    11d16s20
        iarg=nclop-mdon                                                 11d16s20
c     this is offset in vector
        do l=1,4                                                        11d16s20
         nff2(nclop,isb,2+l)=joff(l)-joff0(l)                                    11d16s20
         joff(l)=joff(l)+ncsf2(l,iarg)*nff2(nclop,isb,1)                11d16s20
         nct2c(l,isb)=nct2c(l,isb)+ncsf2(l,iarg)*nff2(nclop,isb,1)      11d19s20
        end do                                                          11d16s20
c     this is offset in list
        nff2(nclop,isb,2)=koff                                          11d16s20
        koff=koff+nff2(nclop,isb,1)                                     11d16s20
       end do                                                           11d16s20
      end do                                                            11d16s20
      ismultm=ismult-1                                                  9d28s20
c     for everything
      js3j=ismultm                                                      9d28s20
      msj3j=ismultm                                                     9d28s20
c     for singlet coupled pairs
      ipair3j=0                                                         9d28s20
      is3j=ismultm                                                      9d28s20
      igamma3j=0                                                        9d28s20
c     igam+msi=msj, mis=msj-igam
      msi3j=msj3j-igamma3j                                              9d28s20
      iphs3j=is3j-ipair3j-msj3j                                         9d28s20
      iphs3j=iabs(iphs3j)/2                                             10d26s20
      fact3j=f3j(ipair3j,is3j,js3j,igamma3j,msi3j,-msj3j,1)             9d28s20
      fact3j=fact3j*sqrt(dfloat(js3j+1))                                9d28s20
      if(mod(iphs3j,2).ne.0)fact3j=-fact3j                              9d28s20
      fact3js=fact3j                                                    9d28s20
c     for same spin triplet coupled pairs
      ipair3j=2                                                         9d28s20
      iphs3j=is3j-ipair3j-msj3j                                         9d28s20
      iphs3j=iabs(iphs3j)/2                                             10d26s20
      fact3j=f3j(ipair3j,is3j,js3j,igamma3j,msi3j,-msj3j,1)             9d28s20
      fact3j=fact3j*sqrt(dfloat(js3j+1))                                9d28s20
      if(mod(iphs3j,2).ne.0)fact3j=-fact3j                              9d28s20
      fact3jst=fact3j                                                   9d28s20
c     for high spin triplet coupled pairs
      ipair3j=2                                                         9d28s20
      is3j=ismultm+2                                                    9d28s20
      igamma3j=0                                                        9d28s20
      msi3j=msj3j-igamma3j                                              9d28s20
      iphs3j=is3j-ipair3j-msj3j                                         9d28s20
      iphs3j=iabs(iphs3j)/2                                             10d26s20
      fact3j=f3j(ipair3j,is3j,js3j,igamma3j,msi3j,-msj3j,1)             9d28s20
      fact3j=fact3j*sqrt(dfloat(js3j+1))                                9d28s20
      if(mod(iphs3j,2).ne.0)fact3j=-fact3j                              9d28s20
      fact3jht=fact3j                                                   9d28s20
c     now for ms one less
      igamma3j=+2                                                        9d28s20
c     for high spin triplet coupled pairs
      ipair3j=2                                                         9d28s20
      is3j=ismultm+2                                                    9d28s20
      msi3j=msj3j-igamma3j                                              9d28s20
      iphs3j=is3j-ipair3j-msj3j                                         9d28s20
      iphs3j=iabs(iphs3j)/2                                             10d26s20
      fact3j=f3j(ipair3j,is3j,js3j,igamma3j,msi3j,-msj3j,1)             9d28s20
      fact3j=fact3j*sqrt(dfloat(js3j+1))                                9d28s20
      if(mod(iphs3j,2).ne.0)fact3j=-fact3j                              9d28s20
      fact3jht2=fact3j                                                   9d28s20
c     for same spin triplet coupled pairs
      ipair3j=2                                                         9d28s20
      is3j=ismultm                                                      9d28s20
      msi3j=msj3j-igamma3j                                              9d28s20
      iphs3j=is3j-ipair3j-msj3j                                         9d28s20
      iphs3j=iabs(iphs3j)/2                                             10d26s20
      fact3j=f3j(ipair3j,is3j,js3j,igamma3j,msi3j,-msj3j,1)             9d28s20
      fact3j=fact3j*sqrt(dfloat(js3j+1))                                9d28s20
      if(mod(iphs3j,2).ne.0)fact3j=-fact3j                              9d28s20
      fact3jst2=fact3j                                                   9d28s20
c     for low spin triplet coupled pairs
      ipair3j=2                                                         9d28s20
      is3j=ismultm-2                                                    9d28s20
      msi3j=is3j
      msi3j=msj3j-igamma3j                                              9d28s20
      iphs3j=is3j-ipair3j-msj3j                                         9d28s20
      iphs3j=iabs(iphs3j)/2                                             10d26s20
      fact3j=f3j(ipair3j,is3j,js3j,igamma3j,msi3j,-msj3j,1)             9d28s20
      fact3j=fact3j*sqrt(dfloat(js3j+1))                                9d28s20
      if(mod(iphs3j,2).ne.0)fact3j=-fact3j                              9d28s20
      fact3jlt=fact3j                                                   9d28s20
      if(mynowprog.eq.0)then                                            9d29s20
       write(6,21)fact3js                                               2d2s21
   21  format('3j for sp: ',f12.6)                                      2d2s21
       write(6,22)fact3jst,fact3jht                                     2d2s21
   22  format('3j for tp, same ms: ',2f12.6)                            2d2s21
       write(6,23)fact3jlt,fact3jst2,fact3jht2                          2d2s21
   23  format('3j for tp, ms-1: ',3f12.6)                               2d2s21
      end if                                                            9d29s20
      if(abs(fact3jst).gt.1d-12)fact3jst=1d0/fact3jst                   9d28s20
      if(abs(fact3jht).gt.1d-12)fact3jht=1d0/fact3jht                   9d28s20
      if(abs(fact3jlt).gt.1d-12)fact3jlt=1d0/fact3jlt                   9d28s20
      if(abs(fact3jst2).gt.1d-12)fact3jst2=1d0/fact3jst2                9d28s20
      if(abs(fact3jht2).gt.1d-12)fact3jht2=1d0/fact3jht2                9d28s20
      srh=sqrt(0.5d0)                                                   9d28s20
      if(ldebug)then                                                    12d12s19
       write(6,*)('nsymb = '),nsymb
       write(6,*)('norb = '),norb
      end if                                                            12d12s19
      nrtot=0                                                           8d1s19
      nderp=nder+1                                                      8d30s23
      do isb=1,nsymb
       idviu=idvi(isb)                                                  8d29s23
       if(ldebug)then                                                   12d12s19
        write(6,*)('for isb = '),isb
        write(6,*)('nfcn,nct2 '),nfcn(isb),nct(isb),
     $      nct2(isb)
        write(6,*)('irwt, irxinfo: '),irwt(isb),irxinfo(isb)
        write(6,*)('ivintref: '),ivintref(isb)
       end if                                                           12d12s19
       if(isb.eq.isymmrci)then                                          8d1s19
        nrootu=max(irxinfo(isb),nroot)                                  8d1s19
       else                                                             8d1s19
        nrootu=irxinfo(isb)                                             8d1s19
       end if                                                           8d1s19
       nrtot=nrtot+nrootu                                               8d1s19
       if(nrootu.gt.0)then                                              8d1s19
        ivintr(isb)=ibcoff                                              8d1s19
        ibcoff=ivintr(isb)+nct(isb)*nrootu*nderp                        8d30s23
        call enough('cont2ecsf.  1',bc,ibc)
        do j=0,nrootu-1                                                 8d1s19
         iad1=ivintr(isb)+nct(isb)*nderp*j                              8d30s23
         iad2=ivintref(isb)+nct(isb)*j                                  8d1s19
         do k=0,nct(isb)-1                                              8d1s19
          bc(iad1+k)=bc(iad2+k)*bc(irwt(isb)+j)                         8d1s19
         end do                                                         8d1s19
        end do                                                          8d1s19
        ladd=nrootu                                                     8d29s23
        do l=1,nder                                                     8d29s23
         do j=0,nrootu-1                                                8d29s23
          jp=l+nderp*j                                                  8d30s23
          iad1=ivintr(isb)+nct(isb)*jp                                  8d29s23
          iad2=ivintref(isb)+nct(isb)*j
          iad3=idviu+nrootu+nct(isb)*j                                  8d29s23
          do k=0,nct(isb)-1                                             8d29s23
           bc(iad1+k)=bc(iad2+k)*bc(idviu+j)+bc(iad3+k)*bc(irwt(isb)+j) 8d31s23
          end do                                                        8d29s23
         end do                                                         8d29s23
         ladd=ladd+nrootu                                               8d29s23
         idviu=idviu+nrootu*(1+nct(isb))                                8d29s23
        end do                                                          8d29s23
       end if                                                           8d1s19
      end do
      ibcoffo=ibcoff                                                    11d21s19
      if(ldebug)write(6,*)('nrtot = '),nrtot
      nrow=(norb*(norb+1))/2                                            8d1s19
      call ilimts(nrow,1,mynprocg,mynowprog,ilrow,ihrow,i1s,i1e,i2s,i2e)10d16s20
      ihrow=ihrow-1                                                     10d17s20
      ilrow=ilrow-1                                                     10d17s20
      nrowh=max(1,ihrow-ilrow+1)                                        10d17s20
      do isb=1,nsymb                                                    7d26s19
       nsum=nct2c(1,isb)                                                8d14s19
       do i=2,4                                                         8d14s19
        nsum=nsum+nct2c(i,isb)                                          8d14s19
       end do                                                           8d14s19
       if(ldebug)
     $ write(6,*)('nct2c vs nct2: '),(nct2c(i,isb),i=1,4),nsum,nct2(isb)12d12s19
       do i=1,4                                                         8d14s19
        ivcc(i)=ibcoff                                                   8d14s19
        ibcoff=ivcc(i)+nct2c(i,isb)*nrowh*nrtot*nderp                   8d30s23
        if(ldebug)write(6,*)('ivcc: '),i,ivcc(i),nct2c(i,isb),nrow,nrtot12d12s19
     $       ,nct2c(i,isb)*nrowh*nrtot
       end do                                                           8d14s19
       call enough('cont2ecsf.  2',bc,ibc)
       do i=ivcc(1),ibcoff-1                                            8d14s19
        bc(i)=0d0                                                       8d14s19
       end do                                                           8d14s19
       ivc=ibcoff                                                       8d1s19
       ibcoff=ivc+nct2(isb)*nrowh*nrtot*nderp                           8d30s23
       if(ldebug)write(6,*)('space for ivc '),nct2(isb)*nrowh*nrtot*    8d30s23
     $      nderp                                                       8d30s23
       call enough('cont2ecsf.  3',bc,ibc)
       do i=0,nct2(isb)*nrowh*nrtot*nderp-1                             8d30s23
        bc(ivc+i)=0d0                                                   8d1s19
       end do                                                           8d1s19
       jvc=ivc                                                          8d1s19
       irooto=0                                                         2d18s20
       jptrto=1                                                         11d16s20
       do jsb=1,nsymb                                                   7d26s19
        if(jsb.eq.isymmrci.or.irxinfo(jsb).gt.0)then                    5d10s21
         if(ldebug)write(6,*)('starting from vectors of symmetry '),jsb,12d12s19
     $       ivintref(jsb)
         if(jsb.ne.isymmrci)then                                        7d26s19
          nrootu=irxinfo(jsb)                                           7d26s19
         else                                                           7d26s19
          nrootu=max(nroot,irxinfo(jsb))                                7d26s19
         end if                                                         7d26s19
         if(ldebug)then                                                 12d12s19
          call prntm2(bc(ivintref(jsb)),nct(jsb),nrootu,                8d30s23
     $        nct(jsb))                                                 8d30s23
          write(6,*)('weighted '),ivintr(jsb)
          call prntm2(bc(ivintr(jsb)),nct(jsb),nrootu*nderp,nct(jsb))   8d30s23
         end if                                                         12d12s19
         jptr=iptr+(mdoo+1)*2*(jsb-1)
         call cont2ecsfa(ibc(ibasis(jsb)),ibc(jptr),idorb,isorb,        7d26s19
     $        nfcn(jsb),nec,nrootu*nderp,icsfpd,                        8d30s23
     $        bc(ivintr(jsb)),nct(jsb),mdon,ncsf,ncsf2,isb,             8d14s19
     $        jsb,bc(jvc),nct2(isb),nrowh,ivcc,nct2c(1,isb),            10d17s20
     $       multh,irooto*(1+nder),fact3jst,fact3jht,fact3jlt,fact3jst2,8d30s23
     $        fact3jht2,srh,ilrow,ihrow,ibc(iff2),ibc(iptrto),jptrto,   11d16s20
     $        nff2,mdoo,nsymb,ndetc,icsfxv,nrtot*(1+nder),bc,ibc)       8d30s23
         irooto=irooto+nrootu                                           8d30s23
         jvc=jvc+nct2(isb)*nrowh*nrootu*nderp                           8d30s23
        end if                                                          7d26s19
       end do                                                           7d26s19
       call cont2ecsfb(bc(ivc),nct2(isb),nrtot,nrowh,                   11d19s20
     $      mdon,ncsf,icsfpd,nec,thrs,nok,iscp,ivcc,                    11d19s20
     $      nct2c(1,isb),ncsf2,lprint,ilrow,                            11d19s20
     $      ihrow,ibc(iff2),nff2,mdoo,isb,nsymb,nfdat(1,1,isb),needcv,  12d18s20
     $      nff22(1,1,isb),bc,ibc,nder,idvmt(1,isb),ntackon)            8d7s24
c
c     consolidate storage.
c
       if(ldebug)write(6,*)('first word '),ibc(nfdat(5,1,isb))          1d25s21
       noffset=ibcoffo-nfdat(5,1,isb)                                   11d18s20
       do i=0,needcv-1                                                  11d18s20
        bc(ibcoffo+i)=bc(nfdat(5,1,isb)+i)                              11d18s20
       end do                                                           11d18s20
       ibcoffo=ibcoffo+needcv                                           11d18s20
       lastadd=0                                                        8d7s24
       do l=1,4                                                         11d18s20
        nfdat(5,l,isb)=nfdat(5,l,isb)+noffset                           3d1s21
        if(nfdat(2,l,isb).gt.0)then                                     11d18s20
         ivec=nfdat(4,l,isb)                                            11d18s20
         do i=0,nfdat(2,l,isb)*nfdat(3,l,isb)-1                         9d5s23
          bc(ibcoffo+i)=bc(ivec+i)                                      11d18s20
         end do                                                         11d18s20
         lastadd=ivec+nfdat(2,l,isb)*nfdat(3,l,isb)                     8d7s24
         nfdat(4,l,isb)=ibcoffo                                         11d18s20
         ibcoffo=ibcoffo+nfdat(2,l,isb)*nfdat(3,l,isb)                  9d5s23
        end if                                                          11d18s20
       end do                                                           11d18s20
       idervcv(1,isb)=ibcoffo                                           8d7s24
       idervcv(2,isb)=ntackon                                           8d7s24
       if(ntackon.gt.0)then                                             8d7s24
        do i=0,ntackon-1                                                8d7s24
         bc(ibcoffo+i)=bc(lastadd+i)                                    8d7s24
        end do                                                          8d7s24
        ibcoffo=ibcoffo+ntackon                                         8d7s24
       end if                                                           8d7s24
       ibcoff=ibcoffo                                                   11d18s20
       if(ldebug)then                                                   1d25s21
        write(6,*)('after move, first word '),ibc(nfdat(5,1,isb)),
     $      nfdat(5,1,isb)
       end if                                                           1d25s21
      end do                                                            7d26s19
      return                                                            7d26s19
      end                                                               7d26s19
