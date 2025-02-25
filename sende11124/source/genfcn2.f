c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genfcn2(ibasis,nbasis,iptr,idorb,isorb,nfcn,nec,       11d12s20
     $     mdon,mdoo,multh,isbg,nsbg,ncount,nsymb,iffn,ipto2,nct2,ncsf, 11d11s22
     $     bc,ibc)                                                      11d11s22
      implicit real*8 (a-h,o-z)
      integer*1 idorb(*),isorb(*),iorb(64)                              11d12s20
      integer*8 itesta                                                  11d18s20
      logical ldebug                                                    12d12s19
      dimension ibasis(*),iptr(4,mdoo+1,*),multh(8,8),nfcn(*),isbg(*),  1d26s23
     $     nbasis(*),ncount(mdoo+1,nsymb,2),itesta4(2),nct2(*),ncsf(*)  11d19s20
      equivalence (itesta,itesta4)                                      11d18s20
c
c     take reference set of functions, and generate N-2 set of functions
c
      include "common.store"
      include "common.mrci"                                             9d10s19
      include "common.print"                                            1d5s20
       ldebug=.false.                                                    12d12s19
      ismultm=max(0,ismult-3)                                           11d12s20
      if(ldebug)then                                                    12d12s19
       write(6,*)('hi, my name is genfcn2')
      end if                                                            12d12s19
      do isb=1,nsymb                                                    11d12s20
       do i=1,mdoo+1                                                    11d12s20
        ncount(i,isb,1)=0                                               11d12s20
       end do                                                           11d12s20
      end do                                                            11d12s20
      iffn=ibcoff                                                       11d12s20
      necmm=nec-2                                                       11d12s20
      nffn=0                                                            11d12s20
      do ib=1,nsbg                                                      11d12s20
       isb=isbg(ib)                                                     11d12s20
       ioff=1                                                            7d26s19
       if(ldebug)write(6,*)('for symmetry '),isb
       do if=1,nfcn(isb)                                                11d12s20
        iad=1+nbasis(isb)+3*(if-1)                                      1d26s23
        nclo=ibasis(iad)
        nclom=nclo-1                                                     9d10s19
        nclomm=nclom-1                                                   9d10s19
        nclop=nclo+1
        nopen=nec-2*nclo
        nopenpp=nopen+2                                                  9d10s19
        nopenmm=nopen-2                                                  9d10s19
        ic=iptr(2,nclop,isb)+nclo*(ibasis(iad+1)-1)                     1d26s23
        io=iptr(4,nclop,isb)+nopen*(ibasis(iad+2)-1)
        iarg=nclop-mdon
        if(ldebug)then                                                   12d12s19
         write(6,*)('for function ')
         write(6,*)if,(idorb(ic+j),j=0,nclo-1),(':'),
     $      (isorb(io+j),j=0,nopen-1)
        end if                                                           12d12s19
c
c     there are 4 scenarios:
c     1: excite both electrons from a single cs orbital
c     2: excite electrons from two distinct cs orbitals
c     3: excite one electron from a cs orbital and one from an os orbital
c     4: excite electrons from two distinct os orbitals.
c
        if(nclo.gt.0)then                                                11d12s20
         isbh=1                                                         11d12s20
         do j=0,nopen-1                                                 11d12s20
          isbh=multh(isbh,ism(isorb(io+j)))                             11d12s20
         end do                                                         11d12s20
         do i=1,nclo                                                     9d10s19
          im=i-1                                                         9d10s19
          jj=0                                                           9d10s19
          jffn=iffn+nffn                                                 11d12s20
          nffn=nffn+1                                                    11d12s20
          ibcoff=iffn+nffn                                               11d12s20
          ncount(nclo,isbh,1)=ncount(nclo,isbh,1)+1                     11d12s20
          call enough('genfcn2.  1',bc,ibc)
          ibc(jffn)=0                                                    11d12s20
          do j=1,nclo                                                    9d10s19
           if(j.ne.i)then                                                9d10s19
            jm=j-1                                                       9d10s19
            jj=jj+1                                                       9d10s19
            iorb(jj)=idorb(ic+jm)                                         9d10s19
            ibc(jffn)=ibset(ibc(jffn),idorb(ic+jm))                      11d12s20
           end if                                                        9d10s19
          end do                                                         9d10s19
          do j=0,nopen-1                                                 11d12s20
           ibc(jffn)=ibset(ibc(jffn),32+isorb(io+j))                    11d18s20
          end do                                                         11d12s20
          if(ldebug)write(6,*)('looking for '),(iorb(ij),ij=1,nclom),   12d12s19
     $         (':'),(iorb(ij),ij=nclo,nclo+nopen),(' nffn = '),nffn-1             5d10s21
          icol=((idorb(ic+im)*(idorb(ic+im)+1))/2)-1                     9d12s19
         end do                                                          9d10s19
        end if                                                           9d10s19
        if(nclo.gt.1)then                                                9d10s19
         do i=1,nclo                                                     9d10s19
          im=i-1                                                         9d10s19
          do j=1,i-1                                                     9d10s19
           jm=j-1                                                        9d10s19
           if(idorb(ic+im).lt.idorb(ic+jm))then                          9d10s19
            in=idorb(ic+im)                                              9d10s19
            ix=idorb(ic+jm)                                              9d10s19
           else                                                          9d10s19
            in=idorb(ic+jm)                                              9d10s19
            ix=idorb(ic+im)                                              9d10s19
           end if                                                        9d10s19
           if(ldebug)write(6,*)('scenario 2')
           if(ldebug)write(6,*)('in,ix: '),in,ix                         12d12s19
           icol=((ix*(ix-1))/2)+in-1                                     9d12s19
           if(ldebug)write(6,*)('store to col '),icol
           jj=0                                                          9d10s19
           do k=1,nclo                                                   9d10s19
            if(k.ne.i.and.k.ne.j)then                                    9d10s19
             km=k-1                                                      11d12s20
             jj=jj+1
             iorb(jj)=idorb(ic+km)                                       9d10s19
            end if                                                       9d10s19
           end do                                                         9d10s19
           iput=0                                                        9d10s19
           jjo=jj+1                                                      9d10s19
           do k=0,nopen-1                                                9d10s19
            if(isorb(io+k).lt.in)then                                    9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=isorb(io+k)                                        9d10s19
            else if(isorb(io+k).lt.ix)then                               9d10s19
             if(iput.eq.0)then                                           9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=in                                                9d10s19
              iput=1                                                     9d10s19
             end if                                                      9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=isorb(io+k)                                        9d10s19
            else if(isorb(io+k).gt.ix)then                               9d10s19
             if(iput.eq.0)then                                           9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=in                                                9d10s19
              iput=1                                                     9d10s19
             end if                                                      9d10s19
             if(iput.eq.1)then                                           9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=ix                                                9d10s19
              iput=2                                                     9d10s19
             end if                                                      9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=isorb(io+k)                                        9d10s19
            end if                                                       9d10s19
           end do                                                        9d10s19
           if(iput.eq.0)then                                             9d10s19
            jj=jj+1                                                      9d10s19
            iorb(jj)=in                                                  9d10s19
            iput=1                                                       9d10s19
           end if                                                        9d10s19
           if(iput.eq.1)then                                             9d10s19
            jj=jj+1                                                      9d10s19
            iorb(jj)=ix                                                  9d10s19
           end if                                                        9d10s19
           if(ldebug)write(6,*)('looking for '),(iorb(ij),ij=1,nclomm),  12d12s19
     $         (':'),(iorb(ij),ij=jjo,jj),(' nffn = '),nffn             5d10s21
           jffn=iffn+nffn                                                 11d12s20
           nffn=nffn+1                                                    11d12s20
           ibcoff=iffn+nffn                                               11d12s20
           call enough('genfcn2.  2',bc,ibc)
           ibc(jffn)=0                                                    11d12s20
           do ij=1,nclomm                                                11d12s20
            ibc(jffn)=ibset(ibc(jffn),iorb(ij))                          11d12s20
           end do                                                        11d12s20
           isbh=1                                                       11d12s20
           do ij=jjo,jj                                                  11d12s20
            ibc(jffn)=ibset(ibc(jffn),32+iorb(ij))                      11d18s20
            isbh=multh(isbh,ism(iorb(ij)))                              11d12s20
           end do                                                        11d12s20
           ncount(nclom,isbh,1)=ncount(nclom,isbh,1)+1                  11d12s20
          end do                                                         9d10s19
         end do                                                          9d10s19
        end if                                                           9d10s19
        if(nclo.gt.0.and.nopen.gt.0)then                                 9d10s19
         do i=0,nclo-1                                                   9d10s19
          do j=0,nopen-1                                                 9d10s19
           if(ldebug)                                                    12d12s19
     $        write(6,*)('try exciting out of '),idorb(ic+i),isorb(io+j)12d12s19
           if(ldebug)write(6,*)('scenario 3')                             12d12s19
           if(idorb(ic+i).lt.isorb(io+j))then                            9d12s19
            in=idorb(ic+i)                                               9d12s19
            ix=isorb(io+j)                                               9d12s19
           else                                                          9d12s19
            ix=idorb(ic+i)                                               9d12s19
            in=isorb(io+j)                                               9d12s19
           end if                                                        9d12s19
           icol=((ix*(ix-1))/2)+in-1                                     9d12s19
           jj=0                                                          9d10s19
           do k=0,nclo-1                                                 9d10s19
            if(k.ne.i)then                                               9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=idorb(ic+k)                                        9d10s19
            end if                                                       9d10s19
           end do                                                        9d10s19
           jjo=jj+1                                                      9d10s19
           iput=0                                                        9d10s19
           do k=0,nopen-1                                                9d10s19
            if(k.ne.j)then                                               9d10s19
             if(isorb(io+k).lt.idorb(ic+i))then                          9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=isorb(io+k)                                       9d10s19
             else                                                        9d10s19
              if(iput.eq.0)then                                          9d10s19
               jj=jj+1                                                   9d10s19
               iorb(jj)=idorb(ic+i)                                      9d10s19
               iput=1                                                    9d10s19
              end if                                                     9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=isorb(io+k)                                       9d10s19
             end if                                                      9d10s19
            else if(iput.eq.0)then                                       11d18s19
             if(idorb(ic+i).lt.isorb(io+k))then
              jj=jj+1                                                    11d18s19
              iorb(jj)=idorb(ic+i)                                       11d18s19
              iput=1                                                     11d18s19
             end if                                                      11d18s19
            end if                                                       9d10s19
           end do                                                        9d10s19
           if(iput.eq.0)then                                             9d10s19
            jj=jj+1                                                      9d10s19
            iorb(jj)=idorb(ic+i)                                         9d10s19
           end if                                                        9d10s19
           if(idorb(ic+i).lt.isorb(io+j))then                            11d18s19
            iextra=1                                                     9d16s19
           else                                                          9d16s19
            iextra=0                                                     9d16s19
           end if                                                        9d16s19
           if(ldebug)write(6,*)('looking for '),(iorb(ij),ij=1,nclom),   12d12s19
     $         (':'),(iorb(ij),ij=jjo,jj),(' nffn = '),nffn                               12d12s19
           jffn=iffn+nffn                                                 11d12s20
           nffn=nffn+1                                                    11d12s20
           ibcoff=iffn+nffn                                               11d12s20
           call enough('genfcn2.  3',bc,ibc)
           ibc(jffn)=0                                                    11d12s20
           do ij=1,nclom                                                 11d12s20
            ibc(jffn)=ibset(ibc(jffn),iorb(ij))                          11d12s20
           end do                                                        11d12s20
           isbh=1                                                       11d12s20
           do ij=jjo,jj                                                  11d12s20
            ibc(jffn)=ibset(ibc(jffn),32+iorb(ij))                      11d18s20
            isbh=multh(isbh,ism(iorb(ij)))                              11d12s20
           end do                                                        11d12s20
           ncount(nclo,isbh,1)=ncount(nclo,isbh,1)+1                    11d12s20
          end do
         end do                                                          9d10s19
        end if                                                           9d10s19
        if(nopenmm.ge.ismultm)then                                      11d12s20
         do i=0,nopen-1                                                  9d10s19
          do j=0,i-1                                                     9d10s19
           if(ldebug)write(6,*)('scenario 4')                             12d12s19
           icol=((isorb(io+i)*(isorb(io+i)-1))/2)+isorb(io+j)-1          9d12s19
           jj=0                                                          9d10s19
           do k=0,nopen-1                                                9d10s19
            if(k.ne.i.and.k.ne.j)then                                    9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=isorb(io+k)                                        9d10s19
            end if                                                       9d10s19
           end do                                                        9d10s19
           jffn=iffn+nffn                                                 11d12s20
           if(ldebug)write(6,*)('looking for '),(iorb(ij),ij=1,nclom),   12d12s19
     $         (':'),(iorb(ij),ij=jjo,jj),(' nffn = '),nffn                               12d12s19
           nffn=nffn+1                                                    11d12s20
           ibcoff=iffn+nffn                                               11d12s20
           call enough('genfcn2.  4',bc,ibc)
           ibc(jffn)=0                                                    11d12s20
           do ij=0,nclo-1                                                11d12s20
            ibc(jffn)=ibset(ibc(jffn),idorb(ic+ij))                           11d12s20
           end do                                                        11d12s20
           isbh=1                                                       11d12s20
           do ij=1,nopenmm                                               11d12s20
            ibc(jffn)=ibset(ibc(jffn),32+iorb(ij))                      11d18s20
            isbh=multh(isbh,ism(iorb(ij)))                              11d12s20
           end do                                                        11d12s20
           ncount(nclop,isbh,1)=ncount(nclop,isbh,1)+1                  11d12s20
          end do                                                         9d10s19
         end do                                                          9d10s19
        end if                                                           9d10s19
       end do
      end do
      nsum=0
      iffnc0=ibcoff                                                     11d12s20
      iffnc=ibcoff                                                      11d12s20
      do isb=1,nsymb                                                    11d12s20
       nsumh=0                                                          11d12s20
       do ii=1,mdoo+1
        ncount(ii,isb,2)=iffnc                                          11d12s20
        iffnc=iffnc+ncount(ii,isb,1)*2                                  11d16s20
        nsumh=nsumh+ncount(ii,isb,1)                                    11d12s20
       end do
       nsum=nsum+nsumh                                                  11d12s20
      end do                                                            11d12s20
      ibcoff=iffnc                                                      11d12s20
      ipto2=ibcoff                                                      11d16s20
      ibcoff=ipto2+nffn                                                 11d16s20
      call enough('genfcn2.  5',bc,ibc)
      do i=0,nffn-1                                                     11d16s20
       ibc(ipto2+i)=0                                                   11d16s20
      end do                                                            11d16s20
      do i=0,nffn-1                                                     11d12s20
       jffn=iffn+i                                                      11d12s20
       isbh=1                                                           11d12s20
       nopen=0                                                          11d12s20
       do j=33,32+norb                                                  11d18s20
        if(btest(ibc(jffn),j))then                                      11d12s20
         nopen=nopen+1                                                  11d12s20
         isbh=multh(isbh,ism(j-32))                                     11d18s20
        end if                                                          11d12s20
       end do                                                           11d12s20
       nclo=(necmm-nopen)/2                                             11d12s20
       nclop=nclo+1                                                     11d12s20
       ibc(ncount(nclop,isbh,2))=ibc(jffn)                              11d12s20
       ncp=ncount(nclop,isbh,2)+ncount(nclop,isbh,1)                    11d16s20
       ibc(ncp)=jffn+1-iffn                                             11d16s20
       ncount(nclop,isbh,2)=ncount(nclop,isbh,2)+1                      11d12s20
      end do                                                            11d12s20
      iffnc=iffnc0                                                      11d12s20
      jffn=iffn                                                         11d12s20
      ntot=0                                                            11d12s20
      do isb=1,nsymb                                                    11d12s20
       if(iprtr(20).ne.0)write(6,*)('for symmetry '),isb                1d26s21
       nsumh=0                                                          11d12s20
       nsumcsf=0                                                        11d19s20
       do ii=1,mdoo+1                                                   11d12s20
        nufcn=0                                                         11d12s20
        if(ncount(ii,isb,1).gt.0)then                                   11d12s20
         isort=ibcoff                                                   11d12s20
         ibcoff=isort+ncount(ii,isb,1)                                  11d12s20
         call idsortdws(ibc(iffnc),ibc(isort),ncount(ii,isb,1))         1d18s23
         ibcoff=isort                                                   11d12s20
         ireff=0
    1    continue                                                       11d12s20
         kref=iffnc+ireff                                                11d12s20
         do i=ireff,ncount(ii,isb,1)-1                                   11d12s20
          kp=ibc(isort+i)-1+iffnc+ncount(ii,isb,1)                      11d16s20
          if(ibc(kref).ne.ibc(iffnc+i))then
           ireff=i                                                      11d12s20
           go to 2                                                      11d12s20
          else                                                          11d16s20
           jpto2=ipto2+ibc(kp)-1                                        11d16s20
           ibc(jpto2)=jffn+1-iffn
          end if                                                        11d12s20
         end do                                                         11d12s20
         ireff=ncount(ii,isb,1)                                         11d12s20
    2    continue
         ibc(jffn)=ibc(kref)                                            11d12s20
         jffn=jffn+1                                                    11d12s20
         nufcn=nufcn+1
         iarg=ii-mdon                                                   11d19s20
         nsumh=nsumh+1                                                  11d12s20
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
         if(iprtr(20).ne.0)write(6,*)nsumh,ibc(kref),('c:o'),                             11d26s20
     $         (iorb(j),j=1,nclo),(':'),(iorb(j),j=jjo,jj-1)            11d12s20
         if(ireff.lt.ncount(ii,isb,1))go to 1                           11d12s20
        end if                                                          11d12s20
        iffnc=iffnc+ncount(ii,isb,1)*2                                  11d16s20
        ncount(ii,isb,1)=nufcn                                          11d12s20
       end do                                                           11d12s20
       if(iprtr(20).ne.0)
     $      write(6,*)('total no. of fcns of this symmetry = '),nsumh        11d12s20
       nct2(isb)=nsumcsf                                                11d19s20
       ntot=ntot+nsumh
      end do                                                            11d12s20
      jpto2=jffn                                                        11d16s20
      do i=0,nffn-1                                                     11d16s20
       ibc(jpto2+i)=ibc(ipto2+i)                                        11d16s20
      end do                                                            11d16s20
      ipto2=jpto2                                                       11d16s20
      ibcoff=ipto2+nffn                                                 11d16s20
      iff=1                                                             6d22s21
      do isb=1,nsymb                                                    6d22s21
       do ii=mdon+1,mdoo+1                                              6d22s21
        ncount(ii,isb,2)=iff                                            6d22s21
        iff=iff+ncount(ii,isb,1)*2                                      6d22s21
       end do                                                           6d22s21
      end do                                                            6d22s21
      return
      end
