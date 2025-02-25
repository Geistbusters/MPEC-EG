c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine make1dm(nsymb,noc4,nbasdws,lprint,maxddi,mjmat,mkmat,
     $     mhmat,jdenpt,isblk,isblkk,isblkh,nsdlk,nsdlkk,nsdlkh,itab,
     $     idbk,isblk1,nsdlk1,isblkd,nsdlkd,ifull4x)                    7d5s21
      logical lprint
      dimension noc4(nsymb),nbasdws(nsymb),jdenpt(*),isblk(4,*),        7d13s23
     $     isblkk(4,*),mjmat(*),mkmat(*),mhmat(*),isblkd(4,*),          7d13s23
     $     isblkh(4,*),itab(8,8),isblk1(4,*)                            7d13s23
      nsdlk=0                                                           6d28s04
      nsdlkk=0
      nsdlk1=0                                                          3d9s12
      nsdlkd=0                                                          3d1s16
      do isb=1,nsymb
       n1=(noc4(isb)*(noc4(isb)+1))/2                                   6d28s04
       n2=nbasdws(isb)*nbasdws(isb)                                     6d28s04
       n3=noc4(isb)*noc4(isb)                                           7d12s04
c
c     cuz we also index 3x as onex and we may have unoccupied symmetries8d24s16
c
       if(n2.gt.0)then                                                  8d24s16
        nsdlk1=nsdlk1+1                                                 3d9s12
        do j=1,4                                                        3d9s12
         isblk1(j,nsdlk1)=isb                                           3d9s12
        end do                                                          8d24s16
       end if                                                           8d24s16
       if(n1.gt.0)then                                                  6d28s04
        nsdlk=nsdlk+1                                                   6d28s04
        iarg1=n1
        iarg2=n2                                                        1d25s05
  536   format('jmat ',i3,' dimensions ',3i8,' type ',4i1)
        jdenpt(nsdlk)=0
        nsdlkk=nsdlkk+1
        iarg1=n3
        iarg2=n2
  535   format('kmat ',i3,' dimensions ',2i8)
        maxddi=max(maxddi,mkmat(nsdlkk))                                9d17s04
        isblk(1,nsdlk)=isb                                              6d28s04
        isblk(2,nsdlk)=isb                                              6d28s04
        isblk(3,nsdlk)=isb                                              6d28s04
        isblk(4,nsdlk)=isb                                              6d28s04
        isblkk(1,nsdlkk)=isb                                              6d28s04
        isblkk(2,nsdlkk)=isb                                              6d28s04
        isblkk(3,nsdlkk)=isb                                              6d28s04
        isblkk(4,nsdlkk)=isb                                              6d28s04
       end if                                                           6d28s04
      end do
      do isa=1,nsymb
       do isb=1,nsymb                                                   6d28s04
        if(isa.ne.isb)then                                              6d28s04
        n1=(noc4(isa)*(noc4(isa)+1))/2                                  6d28s04
        n2=nbasdws(isb)*nbasdws(isb)                                    6d28s04
        n3=noc4(isa)*noc4(isa)                                          7d12s04
c
c     cuz we also index 3x as onex and we may have unoccupied symmetries8d24s16
c
       if(n2.ne.0.and.nbasdws(isa)*nbasdws(isa).ne.0)then               8d24s16
         nsdlk1=nsdlk1+1                                                3d9s12
         isblk1(1,nsdlk1)=isa                                           3d9s12
         isblk1(2,nsdlk1)=isa                                           3d9s12
         isblk1(3,nsdlk1)=isb                                           3d9s12
         isblk1(4,nsdlk1)=isb                                           3d9s12
        end if                                                          8d24s16
        if(n1.gt.0.and.n2.gt.0)then                                     8d7s07
         nsdlk=nsdlk+1                                                  6d28s04
         jdenpt(nsdlk)=0
         iarg1=n1                                                       1d25s05
         iarg2=n2                                                       1d25s05
        end if                                                          7d9s07
        if(n3.gt.0.and.n2.gt.0)then                                     8d7s07
         nsdlkk=nsdlkk+1
         iarg1=n3                                                       1d25s05
         iarg2=n2                                                       1d25s05
         maxddi=max(maxddi,mkmat(nsdlkk))                               9d17s04
         isblk(1,nsdlk)=isa                                             6d28s04
         isblk(2,nsdlk)=isa                                             6d28s04
         isblk(3,nsdlk)=isb                                             6d28s04
         isblk(4,nsdlk)=isb                                             6d28s04
         isblkk(1,nsdlkk)=isa                                             6d28s04
         isblkk(2,nsdlkk)=isa                                             6d28s04
         isblkk(3,nsdlkk)=isb                                             6d28s04
         isblkk(4,nsdlkk)=isb                                             6d28s04
        end if                                                          6d28s04
        end if                                                          6d28s04
       end do
      end do
      do isa=1,nsymb
       do isb=1,isa-1
        n1=noc4(isa)*noc4(isb)
        n2=nbasdws(isa)*nbasdws(isb)
         nsdlk=nsdlk+1                                                    6d28s04
        jdenpt(nsdlk)=0
c
c     cuz we also index 3x as onex and we may have unoccupied symmetries8d24s16
c
        if(n2.ne.0.and.noc4(isa).ne.0)then                              8d24s16
         nsdlk1=nsdlk1+1                                                3d9s12
         isblk1(1,nsdlk1)=isa                                           3d9s12
         isblk1(2,nsdlk1)=isb                                           3d9s12
         isblk1(3,nsdlk1)=isa                                           3d9s12
         isblk1(4,nsdlk1)=isb                                           3d9s12
        end if                                                          8d24s16
        if(n2.ne.0.and.noc4(isb).ne.0)then                              8d24s16
         nsdlk1=nsdlk1+1                                                3d9s12
         isblk1(1,nsdlk1)=isa                                           3d9s12
         isblk1(2,nsdlk1)=isb                                           3d9s12
         isblk1(3,nsdlk1)=isb                                           3d9s12
         isblk1(4,nsdlk1)=isa                                           3d9s12
        end if                                                          8d24s16
        if(n1.gt.0.and.n2.gt.0)then                                     8d7s07
         iarg1=n1                                                       1d25s05
         iarg2=n2                                                       1d25s05
        else
        end if                                                          6d28s04
         nsdlkk=nsdlkk+1
        if(n1.gt.0.and.n2.gt.0)then                                     8d7s07
         iarg1=n1                                                       1d25s05
         iarg2=n2                                                       1d25s05
         maxddi=max(maxddi,mkmat(nsdlkk))                               9d17s04
        else
        end if                                                          6d28s04
         isblk(1,nsdlk)=isa
         isblk(2,nsdlk)=isb
         isblk(3,nsdlk)=isa
         isblk(4,nsdlk)=isb
         isblkk(1,nsdlkk)=isa
         isblkk(2,nsdlkk)=isb
         isblkk(3,nsdlkk)=isa
         isblkk(4,nsdlkk)=isb
         nsdlkk=nsdlkk+1
        if(n1.gt.0.and.n2.gt.0)then                                     8d7s07
         iarg1=n1                                                       1d25s05
         iarg2=n2                                                       1d25s05
         maxddi=max(maxddi,mkmat(nsdlkk))                                9d17s04
        else
        end if                                                          6d28s04
         isblkk(1,nsdlkk)=isb
         isblkk(2,nsdlkk)=isa
         isblkk(3,nsdlkk)=isa
         isblkk(4,nsdlkk)=isb
       end do
      end do
      do isa=1,nsymb
       do isb=1,isa
        iab=itab(isb,isa)
        do isc=1,isa
         if(isc.eq.isa)then
          itop=isb
         else
          itop=isc
         end if
         do isd=1,itop
          if(isa.eq.isb.and.isc.eq.isd)go to 2199
          if(isa.eq.isc.and.isb.eq.isd)go to 2199
          icd=itab(isd,isc)
          iabcd=itab(icd,iab)
          if(iabcd.eq.1)then
 1199      format(4i3,2x,2i3,2x,i3)
           n1=noc4(isa)*noc4(isb)
           n2=nbasdws(isc)*nbasdws(isd)
           nvirtd=nbasdws(isd)-noc4(isd)                                9d30s14
c
c     cuz we also index 3x as onex and we may have unoccupied symmetries8d24s16
c
           if(noc4(isc)*nbasdws(isa)*nbasdws(isb)*nvirtd.gt.0)then      8d24s16
            nsdlk1=nsdlk1+1                                             3d9s12
            isblk1(1,nsdlk1)=isa                                        3d9s12
            isblk1(2,nsdlk1)=isb                                        3d9s12
            isblk1(3,nsdlk1)=isc                                        3d9s12
            isblk1(4,nsdlk1)=isd                                        3d9s12
           end if                                                       3d9s12
           nvirtc=nbasdws(isc)-noc4(isc)                                9d30s14
           if(noc4(isd)*nbasdws(isa)*nbasdws(isb)*nvirtc.gt.0)then      8d24s16
            nsdlk1=nsdlk1+1                                             3d9s12
            isblk1(1,nsdlk1)=isa                                        3d9s12
            isblk1(2,nsdlk1)=isb                                        3d9s12
            isblk1(3,nsdlk1)=isd                                        3d9s12
            isblk1(4,nsdlk1)=isc                                        3d9s12
           end if                                                       3d9s12
           nvirtb=nbasdws(isb)-noc4(isb)                                9d30s14
           if(nbasdws(isc)*noc4(isa)*nbasdws(isd)*nvirtb.gt.0)then      8d24s16
            nsdlk1=nsdlk1+1                                             3d9s12
            isblk1(1,nsdlk1)=isc                                        3d9s12
            isblk1(2,nsdlk1)=isd                                        3d9s12
            isblk1(3,nsdlk1)=isa                                        3d9s12
            isblk1(4,nsdlk1)=isb                                        3d9s12
           end if                                                       3d9s12
           nvirta=nbasdws(isa)-noc4(isa)                                9d30s14
           if(nbasdws(isd)*noc4(isb)*nbasdws(isc)*nvirta.gt.0)then      8d24s16
            nsdlk1=nsdlk1+1                                             3d9s12
            isblk1(1,nsdlk1)=isc                                        3d9s12
            isblk1(2,nsdlk1)=isd                                        3d9s12
            isblk1(3,nsdlk1)=isb                                        3d9s12
            isblk1(4,nsdlk1)=isa                                        3d9s12
           end if                                                       3d9s12
            nsdlk=nsdlk+1                                               6d28s04
           jdenpt(nsdlk)=0
            if(nsdlk.gt.idbk)stop 'idbk'                                6d28s04
           if(n1*n2.gt.0)then                                           6d28s04
            iarg1=n1                                                    1d25s05
            iarg2=n2                                                    1d25s05
           else                                                         8d11s04
           end if                                                       6d28s04
            nsdlkk=nsdlkk+1
           if(n1*n2.gt.0)then                                           6d28s04
            iarg1=n1                                                    1d25s05
            iarg2=n2                                                    1d25s05
            maxddi=max(maxddi,mkmat(nsdlkk))                            9d17s04
           else
           end if                                                       6d28s04
            isblk(1,nsdlk)=isa
            isblk(2,nsdlk)=isb
            isblk(3,nsdlk)=isc
            isblk(4,nsdlk)=isd
            isblkk(1,nsdlkk)=isa
            isblkk(2,nsdlkk)=isb
            isblkk(3,nsdlkk)=isc
            isblkk(4,nsdlkk)=isd
            nsdlkk=nsdlkk+1
           if(n1*n2.gt.0)then                                           6d28s04
            iarg1=n1                                                    1d25s05
            iarg2=n2                                                    1d25s05
            maxddi=max(maxddi,mkmat(nsdlkk))                            9d17s04
           else                                                         8d11s04
           end if                                                       6d28s04
            isblkk(1,nsdlkk)=isb                                        7d28s04
            isblkk(2,nsdlkk)=isa                                        7d28s04
            isblkk(3,nsdlkk)=isc                                        7d28s04
            isblkk(4,nsdlkk)=isd                                        7d28s04
           n1=noc4(isc)*noc4(isd)
           n2=nbasdws(isa)*nbasdws(isb)
            nsdlk=nsdlk+1                                               6d28s04
            jdenpt(nsdlk)=0
            if(nsdlk.gt.idbk)stop 'idbk'                                6d28s04
           if(n1*n2.gt.0)then                                           6d28s04
            iarg1=n1
            iarg2=n2
           else                                                         8d11s04
           end if                                                       6d28s04
            nsdlkk=nsdlkk+1
           if(n1*n2.gt.0)then                                           6d28s04
            iarg1=n1
            iarg2=n2
            maxddi=max(maxddi,mkmat(nsdlkk))                            9d17s04
           else                                                         8d11s04
           end if                                                       6d28s04
            isblk(1,nsdlk)=isc
            isblk(2,nsdlk)=isd
            isblk(3,nsdlk)=isa
            isblk(4,nsdlk)=isb
            isblkk(1,nsdlkk)=isc
            isblkk(2,nsdlkk)=isd
            isblkk(3,nsdlkk)=isa
            isblkk(4,nsdlkk)=isb
            nsdlkk=nsdlkk+1
           if(n1*n2.gt.0)then                                           6d28s04
            iarg1=n1
            iarg2=n2
            maxddi=max(maxddi,mkmat(nsdlkk))                            9d17s04
           else                                                         8d11s04
           end if                                                       6d28s04
            isblkk(1,nsdlkk)=isd                                        7d28s04
            isblkk(2,nsdlkk)=isc                                        7d28s04
            isblkk(3,nsdlkk)=isa                                        7d28s04
            isblkk(4,nsdlkk)=isb                                        7d28s04
          end if
 2199     continue
         end do
        end do
       end do
      end do
      nsdlkh=0                                                          8d11s04
      do isb=1,nsymb
       if(ifull4x.eq.0)then                                             7d5s21
        n4=noc4(isb)*nbasdws(isb)                                        7d28s04
       else                                                             7d5s21
        n4=nbasdws(isb)**2                                              7d5s21
       end if                                                           7d5s21
       n5=(nbasdws(isb)*(nbasdws(isb)+1))/2                             7d28s04
       if(min(n4,n5).gt.0)then                                          7d13s23
        iarg1=n4                                                        1d25s05
        iarg2=n5                                                        1d25s05
        nsdlkh=nsdlkh+1
        maxddi=max(maxddi,mhmat(nsdlkh))                                9d17s04
        isblkh(1,nsdlkh)=isb                                            8d11s04
        isblkh(2,nsdlkh)=isb                                            8d11s04
        isblkh(3,nsdlkh)=isb                                            8d11s04
        isblkh(4,nsdlkh)=isb                                            8d11s04
        nsdlkd=nsdlkd+1                                                 3d1s16
        isblkd(1,nsdlkd)=isb                                            3d1s16
        isblkd(2,nsdlkd)=isb                                            3d1s16
        isblkd(3,nsdlkd)=isb                                            3d1s16
        isblkd(4,nsdlkd)=isb                                            3d1s16
       end if                                                           6d28s04
      end do
      do isa=1,nsymb
       do isb=1,nsymb                                                   6d28s04
        if(isa.ne.isb)then                                              6d28s04
         if(ifull4x.eq.0)then                                           7d5s21
          n4=noc4(isa)*nbasdws(isa)                                      7d28s04
         else                                                           7d5s21
          n4=nbasdws(isa)**2                                            7d5s21
         end if                                                         7d5s21
         n5=(nbasdws(isb)*(nbasdws(isb)+1))/2                           7d28s04
         if(min(n4,n5).gt.0)then                                        7d13s23
          iarg1=n4
          iarg2=n5                                                       1d25s05
          nsdlkh=nsdlkh+1                                               12d6s05
          maxddi=max(maxddi,mhmat(nsdlkh))                               9d17s04
          isblkh(1,nsdlkh)=isa                                           8d11s04
          isblkh(2,nsdlkh)=isa                                           8d11s04
          isblkh(3,nsdlkh)=isb                                           8d11s04
          isblkh(4,nsdlkh)=isb                                           8d11s04
          nsdlkd=nsdlkd+1                                               3d1s16
          isblkd(1,nsdlkd)=isa                                          3d1s16
          isblkd(2,nsdlkd)=isa                                          3d1s16
          isblkd(3,nsdlkd)=isb                                          3d1s16
          isblkd(4,nsdlkd)=isb                                          3d1s16
         end if                                                          6d28s04
        end if                                                          6d28s04
       end do
      end do
      do isa=1,nsymb
       do isb=1,isa-1
        if(ifull4x.eq.0)then                                            7d5s21
         n3=noc4(isa)*nbasdws(isb)                                       8d11s04
        else                                                            7d5s21
         n3=nbasdws(isa)*nbasdws(isb)                                   7d5s21
        end if                                                          7d5s21
        if(n3.gt.0)then                                                 8d11s04
         nsdlkh=nsdlkh+1                                                8d11s04
         iarg1=n3                                                       1d25s05
         iarg2=nbasdws(isa)*nbasdws(isb)                                3d22s06
         maxddi=max(maxddi,mhmat(nsdlkh))                               9d17s04
         isblkh(1,nsdlkh)=isa                                           8d11s04
         isblkh(2,nsdlkh)=isb                                           8d11s04
         isblkh(3,nsdlkh)=isa                                           8d11s04
         isblkh(4,nsdlkh)=isb                                           8d11s04
         nsdlkd=nsdlkd+1                                                3d1s16
         isblkd(1,nsdlkd)=isa                                           3d1s16
         isblkd(2,nsdlkd)=isb                                           3d1s16
         isblkd(3,nsdlkd)=isa                                           3d1s16
         isblkd(4,nsdlkd)=isb                                           3d1s16
         nsdlkd=nsdlkd+1                                                3d1s16
         isblkd(1,nsdlkd)=isa                                           3d1s16
         isblkd(2,nsdlkd)=isb                                           3d1s16
         isblkd(3,nsdlkd)=isb                                           3d1s16
         isblkd(4,nsdlkd)=isa                                           3d1s16
        end if                                                          8d11s04
        if(ifull4x.eq.0)then                                            7d5s21
         n4=noc4(isb)*nbasdws(isa)                                       8d11s04
        else                                                            7d5s21
         n4=nbasdws(isb)*nbasdws(isa)                                   7d5s21
        end if                                                          7d5s21
        if(n4.gt.0)then                                                 8d11s04
         nsdlkh=nsdlkh+1                                                8d11s04
         iarg1=n4                                                       1d25s05
         iarg2=nbasdws(isa)*nbasdws(isb)                                3d22s06
         maxddi=max(maxddi,mhmat(nsdlkh))                               9d17s04
         isblkh(1,nsdlkh)=isb                                           8d11s04
         isblkh(2,nsdlkh)=isa                                           8d11s04
         isblkh(3,nsdlkh)=isa                                           8d11s04
         isblkh(4,nsdlkh)=isb                                           8d11s04
         nsdlkd=nsdlkd+1                                                3d1s16
         isblkd(1,nsdlkd)=isb                                           3d1s16
         isblkd(2,nsdlkd)=isa                                           3d1s16
         isblkd(3,nsdlkd)=isa                                           3d1s16
         isblkd(4,nsdlkd)=isb                                           3d1s16
         nsdlkd=nsdlkd+1                                                3d1s16
         isblkd(1,nsdlkd)=isb                                           3d1s16
         isblkd(2,nsdlkd)=isa                                           3d1s16
         isblkd(3,nsdlkd)=isb                                           3d1s16
         isblkd(4,nsdlkd)=isa                                           3d1s16
        end if                                                          8d11s04
       end do
      end do
      do isa=1,nsymb
       do isb=1,isa
        iab=itab(isb,isa)
        do isc=1,isa
         if(isc.eq.isa)then
          itop=isb
         else
          itop=isc
         end if
         do isd=1,itop
          if(isa.eq.isb.and.isc.eq.isd)go to 12199
          if(isa.eq.isc.and.isb.eq.isd)go to 12199
          icd=itab(isd,isc)
          iabcd=itab(icd,iab)
          if(iabcd.eq.1)then
           if(ifull4x.eq.0)then                                         7d5s21
            n3=noc4(isa)*nbasdws(isb)                                    8d11s04
           else                                                         7d5s21
            n3=nbasdws(isa)*nbasdws(isb)                                7d5s21
           end if                                                       7d5s21
           n2=nbasdws(isc)*nbasdws(isd)
           if(min(n3,n2).gt.0)then                                      7d13s23
            nsdlkh=nsdlkh+1                                             8d11s04
            iarg1=n3                                                    1d25s05
            iarg2=n2                                                    1d25s05
            maxddi=max(maxddi,mhmat(nsdlkh))                            9d17s04
            isblkh(1,nsdlkh)=isa                                        8d11s04
            isblkh(2,nsdlkh)=isb                                        8d11s04
            isblkh(3,nsdlkh)=isc                                        8d11s04
            isblkh(4,nsdlkh)=isd                                        8d11s04
            nsdlkd=nsdlkd+1                                             3d1s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isc                                        4d18s16
            isblkd(2,nsdlkd)=isd                                        4d18s16
            isblkd(3,nsdlkd)=isb                                        4d18s16
            isblkd(4,nsdlkd)=isa                                        4d18s16
            nsdlkd=nsdlkd+1                                             4d18s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isd                                        4d18s16
            isblkd(2,nsdlkd)=isc                                        4d18s16
            isblkd(3,nsdlkd)=isb                                        4d18s16
            isblkd(4,nsdlkd)=isa                                        4d18s16
           end if                                                       8d11s04
           if(ifull4x.eq.0)then                                         7d5s21
            n4=noc4(isb)*nbasdws(isa)                                    8d11s04
           else                                                         7d5s21
            n4=nbasdws(isb)*nbasdws(isa)                                7d5s21
           end if                                                       7d5s21
           if(min(n4,n2).gt.0)then                                      7d13s23
            nsdlkh=nsdlkh+1                                             8d11s04
            iarg1=n4                                                    1d25s05
            iarg2=n2                                                    1d25s05
            maxddi=max(maxddi,mhmat(nsdlkh))                            9d17s04
            isblkh(1,nsdlkh)=isb                                        8d11s04
            isblkh(2,nsdlkh)=isa                                        8d11s04
            isblkh(3,nsdlkh)=isc                                        8d11s04
            isblkh(4,nsdlkh)=isd                                        8d11s04
            nsdlkd=nsdlkd+1                                                3d1s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isc                                        4d18s16
            isblkd(2,nsdlkd)=isd                                        4d18s16
            isblkd(3,nsdlkd)=isa                                        4d18s16
            isblkd(4,nsdlkd)=isb                                        4d18s16
            nsdlkd=nsdlkd+1                                             4d18s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isd                                        4d18s16
            isblkd(2,nsdlkd)=isc                                        4d18s16
            isblkd(3,nsdlkd)=isa                                        4d18s16
            isblkd(4,nsdlkd)=isb                                        4d18s16
           end if                                                       8d11s04
           n2=nbasdws(isa)*nbasdws(isb)
           if(ifull4x.eq.0)then                                         7d5s21
            n3=noc4(isc)*nbasdws(isd)                                    8d11s04
           else                                                         7d5s21
            n3=nbasdws(isc)*nbasdws(isd)                                7d5s21
           end if                                                       7d5s21
           if(min(n3,n2).gt.0)then                                      7d13s23
            nsdlkh=nsdlkh+1                                             8d11s04
            iarg1=n3
            iarg2=n2
            maxddi=max(maxddi,mhmat(nsdlkh))                            9d17s04
            isblkh(1,nsdlkh)=isc                                        8d11s04
            isblkh(2,nsdlkh)=isd                                        8d11s04
            isblkh(3,nsdlkh)=isa                                        8d11s04
            isblkh(4,nsdlkh)=isb                                        8d11s04
            nsdlkd=nsdlkd+1                                                3d1s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isa                                        4d18s16
            isblkd(2,nsdlkd)=isb                                        4d18s16
            isblkd(3,nsdlkd)=isd                                        4d18s16
            isblkd(4,nsdlkd)=isc                                        4d18s16
            nsdlkd=nsdlkd+1                                             4d18s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isb                                        4d18s16
            isblkd(2,nsdlkd)=isa                                        4d18s16
            isblkd(3,nsdlkd)=isd                                        4d18s16
            isblkd(4,nsdlkd)=isc                                        4d18s16
           end if                                                       8d11s04
           if(ifull4x.eq.0)then                                         7d5s21
            n4=noc4(isd)*nbasdws(isc)                                    7d28s04
           else                                                         7d5s21
            n4=nbasdws(isd)*nbasdws(isc)                                7d5s21
           end if                                                       7d5s21
           if(min(n4,n2).gt.0)then                                      7d13s23
            nsdlkh=nsdlkh+1                                             8d11s04
             iarg1=n4
             iarg2=n2
            maxddi=max(maxddi,mhmat(nsdlkh))                            9d17s04
            isblkh(1,nsdlkh)=isd                                        8d11s04
            isblkh(2,nsdlkh)=isc                                        8d11s04
            isblkh(3,nsdlkh)=isa                                        8d11s04
            isblkh(4,nsdlkh)=isb                                        8d11s04
            nsdlkd=nsdlkd+1                                                3d1s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isa                                        4d18s16
            isblkd(2,nsdlkd)=isb                                        4d18s16
            isblkd(3,nsdlkd)=isc                                        4d18s16
            isblkd(4,nsdlkd)=isd                                        4d18s16
            nsdlkd=nsdlkd+1                                             4d18s16
            if(nsdlkd.gt.idbk)stop 'idbk'                               4d18s16
            isblkd(1,nsdlkd)=isb                                        4d18s16
            isblkd(2,nsdlkd)=isa                                        4d18s16
            isblkd(3,nsdlkd)=isc                                        4d18s16
            isblkd(4,nsdlkd)=isd                                        4d18s16
           end if                                                       8d11s04
          end if
12199     continue
         end do
        end do
       end do
      end do
      return
      end
