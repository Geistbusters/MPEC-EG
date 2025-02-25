c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genfcnp(nec,mdon,mdoo,multh,nff,iff,nffp,iffp,nsymb,   11d12s20
     $     iskeep,nct,ncsf,ntot,bc,ibc)                                 11d11s22
      implicit real*8 (a-h,o-z)
      integer*1 iorb
      integer*8 iff(*)                                                  11d12s20
c
c     add one electron
c
      dimension multh(8,8),nff(mdoo+1,nsymb),                           11d12s20
     $     nffp(mdoo+1,nsymb,3),iorb(64),nct(8),ncsf(*)                 12d11s20
      include "common.store"
      include "common.mrci"                                             9d10s19
      include "common.print"                                            1d5s20
      ismultm=ismult-1                                                  11d12s20
      iffp=ibcoff
      mffp=0                                                            11d12s20
      do isb=1,nsymb
       do ii=1,mdoo+1
        nffp(ii,isb,1)=0
       end do
      end do
      ioff=0                                                            11d12s20
      do isb=1,nsymb                                                    11d12s20
       do ii=1,mdoo+1                                                   11d12s20
        if(nff(ii,isb).gt.0)then                                        11d12s20
         do i=1,nff(ii,isb)
          nclo=0
          nopen=0
          do j=1,norb
           jp=j+32                                                      11d18s20
           if(btest(iff(ioff+i),j))then
            nclo=nclo+1
            iorb(nclo)=j
           end if
           if(btest(iff(ioff+i),jp))then
            nopen=nopen+1
            iorb(norb+nopen)=j
           end if
          end do
          nech=2*nclo+nopen+1                                           11d12s20
          necx=nec-nech                                                 11d12s20
          nopenx=nopen+necx                                             11d13s20
          nok=0
          do j=1,norb                                                   11d12s20
           if(.not.btest(iff(ioff+i),j))then                            11d12s20
            nok=nok+1                                                   11d12s20
            iorb(nok)=j                                                 11d12s20
           end if                                                       11d12s20
          end do                                                        11d12s20
          ibcoff=iffp+mffp+nok                                          11d12s20
          call enough('genfcnp.  1',bc,ibc)
          do j=1,nok                                                    11d12s20
           jffp=iffp+mffp                                               11d12s20
           jp=iorb(j)+32                                                11d18s20
           isbh=multh(isb,ism(iorb(j)))                                 11d12s20
           nclop=ii                                                     11d12s20
           if(btest(iff(ioff+i),jp))then                                11d12s20
            ibc(jffp)=ibset(iff(ioff+i),iorb(j))                        11d12s20
            ibc(jffp)=ibclr(ibc(jffp),jp)                               11d12s20
            nclop=ii+1                                                  11d12s20
            nopenh=nopenx-1                                             11d13s20
           else                                                         11d12s20
            ibc(jffp)=ibset(iff(ioff+i),jp)                             11d12s20
            nopenh=nopenx+1                                             11d13s20
           end if                                                       11d12s20
           if(nopenh.ge.ismultm)then                                    11d13s20
            nffp(nclop,isbh,1)=nffp(nclop,isbh,1)+1                     11d13s20
            mffp=mffp+1                                                 11d13s20
           end if                                                       11d13s20
          end do                                                        11d12s20
         end do
         ioff=ioff+nff(ii,isb)                                          11d12s20
        end if                                                          11d12s20
       end do                                                           11d12s20
      end do                                                            11d12s20
      nsum=0                                                            11d12s20
      iffc=ibcoff                                                       11d12s20
      iffc0=iffc                                                        11d12s20
      do isb=1,nsymb                                                    11d12s20
       do ii=1,mdoo+1                                                   11d12s20
        nffp(ii,isb,2)=iffc                                             11d12s20
        iffc=iffc+nffp(ii,isb,1)                                        11d12s20
        nsum=nsum+nffp(ii,isb,1)                                        11d12s20
       end do                                                           11d12s20
      end do                                                            11d12s20
      ibcoff=iffc                                                       11d12s20
      do i=0,mffp-1                                                     11d12s20
       jffp=iffp+i                                                      11d12s20
       isbh=1                                                           11d12s20
       nclo=0                                                           11d12s20
       do j=1,norb                                                      11d12s20
        jp=j+32                                                         11d18s20
        if(btest(ibc(jffp),j))nclo=nclo+1                               11d12s20
        if(btest(ibc(jffp),jp))isbh=multh(isbh,ism(j))                  11d12s20
       end do                                                           11d12s20
       nclop=nclo+1                                                     11d12s20
       ibc(nffp(nclop,isbh,2))=ibc(jffp)                                11d12s20
       nffp(nclop,isbh,2)=nffp(nclop,isbh,2)+1                          11d12s20
      end do                                                            11d12s20
      iffc=iffc0                                                        11d12s20
      jffp=iffp                                                         11d12s20
      ntot=0                                                            11d12s20
      do isb=1,nsymb                                                    11d12s20
       if(iprtr(20).ne.0)write(6,*)('for symmetry '),isb
       nsumh=0                                                          11d12s20
       nsumcsf=0                                                        11d19s20
       do ii=1,mdoo+1                                                   11d12s20
        nufcn=0                                                         11d12s20
        if(nffp(ii,isb,1).gt.0.and.(iskeep.eq.0.or.iskeep.eq.isb))then  11d12s20
         isort=ibcoff                                                   11d12s20
         ibcoff=isort+nffp(ii,isb,1)                                    11d12s20
         call idsortdws(ibc(iffc),ibc(isort),nffp(ii,isb,1))            1d18s23
         ibcoff=isort                                                   11d12s20
         ireff=0
    1    continue                                                       11d12s20
         kref=iffc+ireff                                                11d12s20
         do i=ireff,nffp(ii,isb,1)-1                                    11d12s20
          if(ibc(kref).ne.ibc(iffc+i))then
           ireff=i                                                      11d12s20
           go to 2                                                      11d12s20
          end if                                                        11d12s20
         end do                                                         11d12s20
         ireff=nffp(ii,isb,1)                                           11d12s20
    2    continue
         ibc(jffp)=ibc(kref)                                            11d12s20
         jffp=jffp+1                                                    11d12s20
         nufcn=nufcn+1
         nsumh=nsumh+1                                                  11d12s20
         iarg=ii-mdon                                                   11d19s20
         nsumcsf=nsumcsf+ncsf(iarg)                                     11d19s20
         nclo=0                                                         11d12s20
         nopen=0                                                        11d12s20
         do j=1,norb                                                    11d12s20
          if(btest(ibc(kref),j))then                                    11d12s20
           nclo=nclo+1                                                  11d12s20
           iorb(nclo)=j                                                 11d12s20
          end if                                                        11d12s20
         end do                                                         11d12s20
         jjo=nclo+1                                                     11d12s20
         jj=jjo                                                         11d12s20
         do j=1,norb                                                    11d12s20
          jp=j+32                                                       11d18s20
          if(btest(ibc(kref),jp))then                                   11d12s20
           nopen=nopen+1                                                11d12s20
           iorb(jj)=j                                                   11d12s20
           jj=jj+1                                                      11d12s20
          end if                                                        11d12s20
         end do                                                         11d12s20
         if(iprtr(20).ne.0)write(6,*)nsumh,ibc(kref),('c:o'),
     $         (iorb(j),j=1,nclo),(':'),(iorb(j),j=jjo,jj-1)
         if(ireff.lt.nffp(ii,isb,1))go to 1                             11d12s20
        end if                                                          11d12s20
        iffc=iffc+nffp(ii,isb,1)                                        11d12s20
        nffp(ii,isb,1)=nufcn                                            11d12s20
       end do                                                           11d12s20
       if(iprtr(20).ne.0)                                               1d26s21
     $      write(6,*)('total no. of fcns of this symmetry = '),nsumh   1d26s21
       if(iskeep.eq.0)then                                              11d18s20
        nct(isb)=nsumcsf                                                11d19s20
       else if(iskeep.eq.isb)then                                       11d18s20
        nct(1)=nsumcsf                                                  11d19s20
       end if                                                           11d18s20
       ntot=ntot+nsumh
      end do                                                            11d12s20
      if(iskeep.eq.0)then                                               11d26s20
       ioff=1                                                           11d26s20
       do isb=1,nsymb                                                   11d26s20
        do ii=1,mdoo+1                                                  11d26s20
         nffp(ii,isb,2)=ioff                                            11d26s20
         ioff=ioff+nffp(ii,isb,1)*2                                     11d26s20
        end do                                                          11d26s20
       end do                                                           11d26s20
      else                                                              11d26s20
       do ii=1,mdoo+1
        nffp(ii,1,1)=nffp(ii,iskeep,1)                                  11d26s20
       end do                                                           11d26s20
       ioff=1                                                           11d26s20
       ioffc=1                                                          12d11s20
       do ii=mdon+1,mdoo+1                                              12d11s20
c                                                                       1d23s23
c     one would think indices should be ii,1,2 and ii,1,3, but          1d23s23
c     if we are here, we are going after internal part which will       1d23s23
c     have no symmetry index.                                           1d23s23
c                                                                       1d23s23
        if(nsymb.eq.1)then                                              2d10s23
         nffp(ii,1,2)=ioff                                               1d23s23
         nffp(ii,1,3)=ioffc                                              1d23s23
        else if(nsymb.eq.2)then                                         2d10s23
         nffp(ii,2,1)=ioff                                               1d23s23
         nffp(ii,1,2)=ioffc                                              1d23s23
        else                                                            2d10s23
         nffp(ii,2,1)=ioff                                               1d23s23
         nffp(ii,3,1)=ioffc                                              1d23s23
        end if                                                          2d10s23
        ioff=ioff+nffp(ii,1,1)*2                                          11d26s20
        iarg=ii-mdon                                                    12d11s20
        ioffc=ioffc+ncsf(iarg)*nffp(ii,1,1)                             12d11s20
       end do                                                           11d26s20
      end if                                                            11d26s20
      ibcoff=jffp
      return
      end
