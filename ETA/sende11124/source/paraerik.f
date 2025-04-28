c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paraerik(natom,ngaus,ibdat,nbasis,ihmat,ipair,         11d27s18
     $     nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,iter,idwsdeb,   8d24s15
     $     idorel,ascale,nbasisp,iaddr,nfcn,nbasispc,ncd,iorb,hcdoub,   12d13s18
     $     ndoub,nrootu,sr2,srh,bc,ibc)                                 11d9s22
c
c     form 4v Hdd*Vd calculation in so basis.                           11d21s18
c     input
c
      implicit real*8 (a-h,o-z)
      include "common.hf"
      include "common.store"
      include "common.basis"
      include "common.mrci"                                             12d12s18
      integer*8 ibstor,isstor                                           5d10s10
      dimension isym(3,1),iapair(3,1),ibstor(1),isstor(1),              11d23s18
     $     carta(3),cartb(3),cartc(3),cartd(3),ieraw(8,8),itmpab(8),    12d7s18
     $     multh(8,8),iptoh(8,8),nbasispc(8),ntmpab(8),                 12d7s18
     $     nbasisp(8),iaddr(36,2,4),nfcn(36,3),ncd(3,8,2),iorb(8),      11d10s20
     $     hcdoub(*)                                                    12d13s18
      dimension nvprt(8),nsz(2)                                         12d13s18
      common/timerocm/tovr,telapo(15)                                   4d26s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      loop=0
      ntmpmat=2                                                         10d28s20
      if(idorel.eq.0.or.idorel.eq.1)then                                10d28s20
       nerimat=1                                                        10d28s20
       ntmpmat=1                                                        10d28s20
      else if(idorel.eq.2.or.idorel.eq.-2.or.idorel.eq.-4)then          10d28s20
       nerimat=3                                                        10d28s20
      else                                                              10d28s20
       nerimat=4                                                        10d28s20
      end if                                                            10d28s20
      xnan=-1d0
      ibcoffo=ibcoff                                                    2d19s10
      if(idorel.eq.0)then                                               8d24s15
       ncomp=1                                                          8d24s15
      else                                                              8d24s15
       ncomp=2                                                          8d24s15
      end if                                                            8d24s15
      nbasallc=nbasall*ncomp                                            8d26s15
      if(idwsdeb.gt.10)then
      write(6,*)('in paraerik '),ibcoff,iter,idorel
      write(6,*)('natom etc '),natom,ngaus,ibdat,nbasis,ihmat,
     $     ipair,nhcolt,isym,iter,idwsdeb,idorel,ascale
      end if
c
c     build shell pair order
c
      nn=(ngaus*(ngaus+1))/2                                            2d19s10
      nn2=nn*2                                                          2d19s10
      ipair=ibcoff                                                      2d19s10
      if(idwsdeb.gt.10)then                                             7d27s16
       write(6,*)('nn,nn2,ibcoff '),nn,nn2,ibcoff
      end if                                                            7d27s16
      ibcoff=ipair+nn*2                                                 3d4s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      memneed(1:3)='i12'
      call enough('paraerik.  1',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus3=ngaus*3                                                    5d12s10
      ngaus7=ngaus*7                                                    5d3s10
c
c     what's under ibdat
c     for each gaussian...
c     1: l
c     2: exponential parameter
c     3: normalization coefficient
c     4: offset in total basis function list
c     5: x of center
c     6: y of center
c     7: z of center
c     8: nhere -> address in iapair of possible symmetry equivalent
c        center
c
c     what's under iapair:
c     1: if nonzero, iabs is symmetry equivalent basis fcn
c     2: group operation relating symmetry equivalent basis fcns
c     3: phase?
c
c
c     estimate work for shell pair
c
      do i1=1,ngaus                                                     2d19s10
       do i2=i1,ngaus                                                   12d10s18
        ibc(j12)=-((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                  2d19s10
        if(iapair(1,ibc(jbdat+i1+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        if(iapair(1,ibc(jbdat+i2+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        ibc(j12+nn)=i1                                                  2d19s10
        ibc(j12+nn2)=i2                                                 2d19s10
        j12=j12+1                                                       2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nn*3                                                     2d19s10
      nn8=nn                                                            2d19s10
c
c     sort, so we can assign most expensive blocks first
c
      if(idwsdeb.gt.10)then
       write(6,*)('calling ihpsort next '),nn8,nn
      end if
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
c
c     order shell pair indices
c
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       ibc(jpair)=ibc(ii1+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
       ibc(jpair)=ibc(ii2+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
      end do                                                            2d19s10
      needt=0                                                           12d12s18
      ibptoh=ibcoff                                                     12d12s18
      do isb=1,nsymb                                                    12d12s18
       isbv12=multh(isb,isymmrci)                                       12d12s18
       do isbv1=1,nsymb                                                 12d12s18
        isbv2=multh(isbv1,isbv12)                                       12d12s18
         iptoh(isbv1,isbv2)=ibcoff                                      12d12s18
         need=nvirt(isbv1)*nbasisp(isbv2)                               11d2s20
         need=need*nrootu*                                               12d13s18
     $        (ncd(1,isb,1)+ncd(1,isb,2)+ncd(2,isb,2)+ncd(3,isb,2))     12d13s18
         ibcoff=ibcoff+need*ntmpmat                                     11d2s20
         needt=needt+need*ntmpmat                                       11d5s20
       end do                                                           12d12s18
      end do                                                            12d12s18
      call enough('paraerik.  2',bc,ibc)
      do i=0,needt-1                                                    12d12s18
       bc(ibptoh+i)=0d0                                                 12d12s18
      end do                                                            12d12s18
      idoit=10
      ndoit=0
      nxtot=0
      do i=1+mynowprog,nn,mynprocg                                      3d4s10
       jpair=ipair+2*(i-1)                                              2d19s10
       if(idwsdeb.gt.20)then
        write(6,2)i,ibc(jpair),ibc(jpair+1),ibc(jpair+2),tcur            3d14s12
        write(6,*)('this is i = '),i
       end if
    2  format('I am going to do ',i5,3i3,f10.4)                               2d19s10
       idoit=idoit+1
       ja=jbdat+ibc(jpair)                                              2d19s10
       ja2=ja+ngaus                                                     2d19s10
       ja3=ja2+ngaus                                                    2d19s10
       ja4=ja3+ngaus                                                    2d19s10
       ja5=ja4+ngaus                                                    2d19s10
       ja6=ja5+ngaus                                                    2d19s10
       ja7=ja6+ngaus                                                    2d19s10
       ja8=ja7+ngaus                                                    5d11s10
       nsza=2*ibc(ja)+1                                                 5d11s10
       nsza0=nsza                                                       11d26s18
       nsza=nsza*ncomp                                                  8d24s15
       nszaa=nsza                                                       5d12s10
       nszaa0=nsza0                                                     10d31s20
       if(iapair(1,ibc(ja8)).gt.0)then                                  5d11s10
        nszaa=nszaa*2                                                   5d12s10
        nszaa0=nszaa0*2                                                 10d31s20
       end if                                                           5d11s10
       jb=jbdat+ibc(jpair+1)                                            2d19s10
       jb2=jb+ngaus                                                     2d19s10
       jb3=jb2+ngaus                                                    2d19s10
       jb4=jb3+ngaus                                                    2d19s10
       jb5=jb4+ngaus                                                    2d19s10
       jb6=jb5+ngaus                                                    2d19s10
       jb7=jb6+ngaus                                                    2d19s10
       jb8=jb7+ngaus                                                    5d11s10
       nszb=2*ibc(jb)+1                                                 5d11s10
       nszb0=nszb                                                       11d26s18
       nszb=nszb*ncomp                                                  8d24s15
       nszba=nszb                                                       5d12s10
       nszba0=nszb0                                                     11d2s20
 3352  format(6i5)
       if(iapair(1,ibc(jb8)).gt.0)then                                  5d11s10
        nszba=nszba*2                                                   5d12s10
        nszba0=nszba0*2                                                 10d31s20
       end if                                                           5d11s10
       ncsz=nsza*nszb                                                   5d11s10
       ncsz0=nsza0*nszb0                                                10d31s20
       ncsza=nszaa*nszba                                                1d10s11
       ncsza0=nszaa0*nszba0                                             10d31s20
c
c     index of of itmpab is symmetry of contracted n-2 electron fcns
c
       do isd=1,nsymb                                                   12d7s18
        itmpab(isd)=ibcoff                                              12d7s18
        is12=multh(isd,isymmrci)                                        12d7s18
        nhere=0                                                         12d7s18
        do isc=1,nsymb                                                  12d7s18
         ise=multh(is12,isc)                                            12d7s18
         if(ise.ge.isc)then                                             12d7s18
          iiii=((ise*(ise-1))/2)+isc                                     12d7s18
          if(nfcn(iiii,3).eq.isd)nhere=max(nhere,                       12d7s18
     $         nfcn(iiii,1)+nfcn(iiii,2))                               12d7s18
         end if                                                         12d7s18
        end do                                                          12d7s18
        ibcoff=itmpab(isd)+nhere*ncsza0*nerimat                         11d2s20
        ntmpab(isd)=nhere                                               12d7s18
        if(ibcoff.lt.0d0)write(6,*)isd,itmpab(isd),nhere,ncsza0,
     $       nerimat,ibcoffo,nszaa0,nszba0,ibc(ja),ja,ibc(jb),jb
        call enough('paraerik.  3',bc,ibc)
        do iz=itmpab(isd),ibcoff-1                                      10d31s20
         bc(iz)=0d0                                                     10d31s20
        end do                                                          12d7s18
       end do                                                           5d12s10
       do isd=1,nsymb                                                   5d12s10
        do isc=1,nsymb                                                  5d12s10
         ieraw(isc,isd)=ibcoff                                          5d12s10
         ibcoff=ieraw(isc,isd)+ncsz0*nbasisp(isc)*nbasisp(isd)*nerimat  10d31s20
         call enough('paraerik.  4',bc,ibc)
         do iz=ieraw(isc,isd),ibcoff-1                                  10d31s20
          bc(iz)=xnan                                                   10d31s20
         end do
        end do                                                          5d12s10
       end do                                                           5d12s10
       memneed(1:5)='tmpab'
       call enough('paraerik.  5',bc,ibc)
       do ipassb=1,nszba/nszb                                           5d12s10
        if(ipassb.eq.1)then
         cartb(1)=bc(jb5)                                               5d12s10
         cartb(2)=bc(jb6)                                               5d12s10
         cartb(3)=bc(jb7)                                               5d12s10
         iob=0                                                          5d12s10
        else                                                            5d12s10
         cartb(1)=bc(jb5)*dfloat(isym(1,iapair(2,ibc(jb8))))            5d11s10
         cartb(2)=bc(jb6)*dfloat(isym(2,iapair(2,ibc(jb8))))            5d11s10
         cartb(3)=bc(jb7)*dfloat(isym(3,iapair(2,ibc(jb8))))            5d11s10
         iob=nszb0                                                      11d2s20
        end if                                                          5d12s10
        do ipassa=1,nszaa/nsza                                          5d12s10
         if(ipassa.eq.1)then                                            5d12s10
          carta(1)=bc(ja5)                                              5d12s10
          carta(2)=bc(ja6)                                              8d15s11
          carta(3)=bc(ja7)                                              8d15s11
          ioa=0                                                         5d12s10
         else                                                           5d12s10
          carta(1)=bc(ja5)*dfloat(isym(1,iapair(2,ibc(ja8))))           5d12s10
          carta(2)=bc(ja6)*dfloat(isym(2,iapair(2,ibc(ja8))))           5d12s10
          carta(3)=bc(ja7)*dfloat(isym(3,iapair(2,ibc(ja8))))           5d12s10
          ioa=nsza0                                                     11d2s20
         end if                                                         5d12s10
         do jjc=1,ngaus                                                 11d28s18
          jc=jbdat+jjc
          jc2=jc+ngaus                                                    2d19s10
          jc3=jc2+ngaus                                                   2d19s10
          jc4=jc3+ngaus                                                   2d19s10
          jc5=jc4+ngaus                                                   2d19s10
          jc6=jc5+ngaus                                                   2d19s10
          jc7=jc6+ngaus                                                   2d19s10
          jc8=jc7+ngaus                                                   5d11s10
          do jjd=1,ngaus                                                11d28s18
           jd=jbdat+jjd                                                 11d23s18
           jd2=jd+ngaus                                                    2d19s10
           jd3=jd2+ngaus                                                   2d19s10
           jd4=jd3+ngaus                                                   2d19s10
           jd5=jd4+ngaus                                                   2d19s10
           jd6=jd5+ngaus                                                   2d19s10
           jd7=jd6+ngaus                                                   2d19s10
           jd8=jd7+ngaus                                                   5d11s10
           nszc=2*ibc(jc)+1                                              2d7s12
           nszc0=nszc                                                   10d28s20
           nszc=nszc*ncomp                                               8d24s15
           nszca=nszc                                                    2d7s12
           nszca0=nszc0                                                 10d31s20
           nszd=2*ibc(jd)+1                                              2d7s12
           nszd0=nszd                                                   10d28s20
           nszd=nszd*ncomp                                               8d24s15
           nszda=nszd                                                    2d7s12
           nszda0=nszd0                                                 10d31s20
           if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
            nszca=nszca*2                                                 5d12s10
            nszca0=nszca0*2                                             10d31s20
           end if                                                         5d12s10
           if(iapair(1,ibc(jd8)).gt.0)then                                5d12s10
            nszda=nszda*2                                                 5d12s10
            nszda0=nszda0*2                                             10d31s20
           end if                                                         5d12s10
           itmp=ibcoff                                                    5d12s10
           jtmp=itmp-1                                                  11d26s18
           ibcoff=itmp+ncsz0*nszca0*nszda0*nerimat                      10d31s20
           memneed(1:4)='itmp'
           call enough('paraerik.  6',bc,ibc)
           do iz=itmp,ibcoff-1                                          10d31s20
            bc(iz)=0d0                                                  10d31s20
           end do                                                       12d10s18
           jdoit=10
c     (ac|bd)=(ab||cd), inner loop over shells is c,d
           if(idorel.ne.0)then                                          10d28s20
           call restin(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3), 10d28s20
     $         ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),         10d28s20
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       10d28s20
     $         ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),         10d28s20
     $       ieri,jdoit,1d0,idorel,ascale,nmat,bc,ibc)                  11d9s22
           else                                                         10d28s20
           call erir(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),    10d4s16
     $         ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),         10d4s16
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       10d4s16
     $         ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),         10d4s16
     $         ieri,jdoit,1d0,bc,ibc)                                   11d14s22
           end if                                                       10d28s20
           if(idoit.eq.1.or.idwsdeb.gt.10)then
            write(6,*)('output from eri: '),ieri                        10d31s20
            call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,nszd0*nszc0) 10d31s20
           end if
c     for 2 component calculation, and optionally compute additional
c     2 electron integrals.
c                                           acbd
c     return nszd*nszc*nszb*nsza block of   LLLL ints,
c     and if idorel ne 1,
c     next the nszd*nszc*nszb*nsza block of LLSS ints,
c     next the nszd*nszc*nszb*nsza block of SSLL ints,
c     and if idorel indicates so,
c     next the nszd*nszc*nszb*nsza block of SSSS ints.
           do in=0,nerimat-1                                            10d31s20
            do id=0,nszd0-1                                               12d10s18
             do ic=0,nszc0-1                                              12d10s18
              do ib=0,nszb0-1                                             12d10s18
               do ia=0,nsza0-1                                            12d10s18
c     ig a b c d in
c     0  L L L L  0
c     1  L S L S  1
c     2  S L S L  2
c     3  S S S S  3
                iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+nsza0*in)))  2d16s21
                iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszca0*(id+nszda0*in))) 10d31s20
                bc(iad2)=bc(iad1)                                        12d10s18
               end do                                                    12d10s18
c     thus tmp is 1 LLLL, 2 LSLS, 3 SLSL, 4 SSSS
              end do                                                     12d10s18
             end do                                                      12d10s18
            end do                                                       12d10s18
           end do                                                       11d2s20
           if(idoit.eq.1.or.idwsdeb.gt.10)then
            write(6,*)('transposed ')
            call prntm2(bc(itmp),ncsz0,nszca0*nszda0*nerimat,ncsz0)      10d31s20
           end if
           if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
            cartc(1)=bc(jc5)*dfloat(isym(1,iapair(2,ibc(jc8))))           5d12s10
            cartc(2)=bc(jc6)*dfloat(isym(2,iapair(2,ibc(jc8))))           5d12s10
            cartc(3)=bc(jc7)*dfloat(isym(3,iapair(2,ibc(jc8))))           5d12s10
            if(idorel.eq.0)then                                         10d28s20
             call erir(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),  4d23s18
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),         4d23s18
     $         ieri,jdoit,1d0,bc,ibc)                                   11d14s22
            else                                                        10d28s20
             call restin(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),
     $         carta(3),ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),         4d23s18
     $         ieri,jdoit,1d0,idorel,ascale,nmat,bc,ibc)                11d9s22
            end if                                                      10d28s20
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('outputb from eri: '),ieri
             call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,nszd0*nszc0)
            end if
            do in=0,nerimat-1                                           10d31s20
             do id=0,nszd0-1                                               12d10s18
              do ic=0,nszc0-1                                              12d10s18
               do ib=0,nszb0-1                                             12d10s18
                do ia=0,nsza0-1                                            12d10s18
                 iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+nsza0*in))) 2d16s21
                 iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszc0+nszca0*(id      10d31s20
     $                +nszda0*in)))                                     10d31s20
                 bc(iad2)=bc(iad1)                                        12d10s18
                end do                                                    12d10s18
               end do                                                     12d10s18
              end do                                                      12d10s18
             end do                                                       12d10s18
            end do                                                      10d31s20
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('transposedb ')
             call prntm2(bc(itmp),ncsz0,nszca0*nszda0*nerimat,ncsz0)      10d31s20
            end if
           end if                                                         5d12s10
           if(iapair(1,ibc(jd8)).gt.0)then                              1d10s11
            cartd(1)=bc(jd5)*dfloat(isym(1,iapair(2,ibc(jd8))))           5d12s10
            cartd(2)=bc(jd6)*dfloat(isym(2,iapair(2,ibc(jd8))))           5d12s10
            cartd(3)=bc(jd7)*dfloat(isym(3,iapair(2,ibc(jd8))))           5d12s10
            if(idorel.eq.0)then                                         10d28s20
             call erir(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),  4d23s18
     $         ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         4d23s18
     $         ieri,jdoit,1d0,bc,ibc)                                   11d14s22
            else                                                        10d28s20
             call restin                                                10d28s20
     $            (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),     10d28s20
     $         ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         4d23s18
     $         ieri,jdoit,1d0,idorel,ascale,nmat,bc,ibc)                11d9s22
            end if                                                      10d28s20
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('outputc from eri: '),ieri,nszd,nszc,ncsz
             call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,nszd0*nszc0)10d31s20
            end if
            do in=0,nerimat-1                                           10d31s20
             do id=0,nszd0-1                                               12d10s18
              do ic=0,nszc0-1                                              12d10s18
               do ib=0,nszb0-1                                             12d10s18
                do ia=0,nsza0-1                                            12d10s18
                 iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+nsza0*in))) 2d16s21
                 iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszca0*(id+nszd0      10d31s20
     $                +nszda0*in)))                                     10d31s20
                 bc(iad2)=bc(iad1)                                        12d10s18
                end do                                                  10d31s20
               end do                                                    12d10s18
              end do                                                     12d10s18
             end do                                                      12d10s18
            end do                                                       12d10s18
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('transposedc ')
             call prntm2(bc(itmp),ncsz0,nszca0*nszda0*nerimat,ncsz0)      10d31s20
            end if
            if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
             if(idorel.eq.0)then                                        10d28s20
              call erir(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3), 4d23s18
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         4d23s18
     $         ieri,jdoit,1d0,bc,ibc)                                   11d14s22
             else                                                       10d28s20
              call restin                                               10d28s20
     $             (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),    10d28s20
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         4d23s18
     $         ieri,jdoit,1d0,idorel,ascale,nmat,bc,ibc)                11d9s22
             end if                                                     10d28s20
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('outputd from eri: '),ieri,nszd,nszc,ncsz
              call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,
     $             nszd0*nszc0)
             end if
             do in=0,nerimat-1                                          10d31s20
              do id=0,nszd0-1                                               12d10s18
               do ic=0,nszc0-1                                              12d10s18
                do ib=0,nszb0-1                                             12d10s18
                 do ia=0,nsza0-1                                            12d10s18
                  iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+nsza0*in)))2d16s21
                  iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszc0+nszca0*(id      10d31s20
     $                 +nszd0+nszda0*in)))                              10d31s20
                  bc(iad2)=bc(iad1)                                        12d10s18
                 end do                                                 11d2s20
                end do                                                    12d10s18
               end do                                                     12d10s18
              end do                                                      12d10s18
             end do                                                       12d10s18
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('transposedd ')
              call prntm2(bc(itmp),ncsz0,nszca0*nszda0*nerimat,ncsz0)      10d31s20
             end if
            end if                                                         5d12s10
           end if                                                         5d12s10
           nszabc0=ncsz0*nszca0                                             5d12s10
           if(nszd.ne.nszda)then                                          5d12s10
            do in=0,nerimat-1                                           10d31s20
             jtmp=itmp+nszabc0*nszda0*in                                10d31s20
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('before d sum/dif: matrix '),in
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)              10d31s20
             end if
             do id=0,nszd0-1                                            10d31s20
              it1=jtmp+nszabc0*id                                       10d31s20
              it2=jtmp+nszabc0*(id+nszd0)                               10d31s20
              do iabc=0,nszabc0-1                                       1d21s21
               sum=srh*(bc(it1+iabc)+bc(it2+iabc))                         5d12s10
               diff=srh*(-bc(it1+iabc)+bc(it2+iabc))                       5d12s10
               bc(it1+iabc)=sum                                            5d12s10
               bc(it2+iabc)=diff                                           5d12s10
              end do                                                       5d12s10
             end do                                                        5d12s10
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('after sum/dif: '),jtmp
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)
             end if
            end do                                                      10d31s20
           end if                                                         5d12s10
           if(nszc.ne.nszca)then                                          5d12s10
            do in=0,nerimat-1                                           10d31s20
             jtmp=itmp+nszabc0*nszda0*in                                2d11s21
             do id=0,nszda0-1                                           10d31s20
              do ic=0,nszc0-1                                           10d31s20
               it1=jtmp+ncsz0*(ic+nszca0*id)                            10d31s20
               it2=jtmp+ncsz0*(ic+nszc0+nszca0*id)                      10d31s20
               do iab=0,ncsz0-1                                         10d31s20
                sum=srh*(bc(it1+iab)+bc(it2+iab))                          5d12s10
                diff=srh*(-bc(it1+iab)+bc(it2+iab))                        5d12s10
                bc(it1+iab)=sum                                            5d12s10
                bc(it2+iab)=diff                                           5d12s10
               end do                                                      5d12s10
              end do                                                       5d12s10
             end do                                                        5d12s10
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('after c sum/dif: for matrix '),in,jtmp
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)
             end if
            end do                                                      10d31s20
           end if                                                         5d12s10
c
c     c&d correspond to b&d from my notes
c
           do id=1,nszda0                                               8d26s15
            idm=id-1                                                    10d31s20
            idd=id+ibc(jd4)                                               5d12s10
            isd=isstor(idd)                                               5d12s10
            iddd=ibstor(idd)-1                                          10d31s20
            do ic=1,nszca0                                              8d26s15
             icm=ic-1                                                   10d31s20
             icc=ic+ibc(jc4)                                              5d12s10
             isc=isstor(icc)                                              5d12s10
             iccc=ibstor(icc)-1                                         10d31s20
             do in=0,nerimat-1                                          10d31s20
c     in a b c d
c     0  L L L L
c     1  L S L S
c     2  S L S L
c     3  S S S S
c     thus eraw is 1 LLLL, 2 LSLS, 3 SLSL, 4 SSSS
              iad1=ieraw(isc,isd)+ncsz0*(iccc+nbasisp(isc)*(iddd        10d31s20
     $             +nbasisp(isd)*in))                                   10d31s20
              iad2=itmp+ncsz0*(icm+nszca0*(idm+nszda0*in))               10d31s20
              do iab=0,ncsz0-1                                          11d2s20
               bc(iad1+iab)=bc(iad2+iab)                                   5d12s10
              end do                                                       5d12s10
             end do                                                     10d31s20
            end do                                                        5d12s10
           end do                                                         5d12s10
           ibcoff=itmp                                                    5d12s10
          end do                                                        11d26s18
c
c     end of jjc loop
c
         end do                                                           3d4s10
c
c     at this point, we have all of cd (=db in my notes) ints,          11d27s18
c     and furthermore they have been symmetrized (turned into so's) and 11d27s18
c     stored under eraw.                                                11d27s18
c     so do the b<d restriction and then multiply by ci vector in ao    11d27s18
c     basis.
c                                                                       1d10s11
         do isd=1,nsymb                                                 11d27s18
          do isc=1,isd                                                  11d27s18
           iiii=((isd*(isd-1))/2)+isc                                   11d27s18
           if(isc.eq.isd)then                                           11d27s18
            nvv=(nbasispc(isc)*(nbasispc(isc)+1))/2                     11d27s18
            nvv0=(nbasisp(isc)*(nbasisp(isc)+1))/2                      10d31s20
           else                                                         11d27s18
            nvv=nbasispc(isc)*nbasispc(isd)                             11d27s18
            nvv0=nbasisp(isc)*nbasisp(isd)                              10d31s20
           end if                                                       11d27s18
           nvvx=nbasisp(isc)*nbasisp(isd)                               2d11s21
           mcd=ntmpab(nfcn(iiii,3))                                     12d7s18
           if(min(mcd,nvv0,ncsz0).gt.0)then                             10d31s20
            itmp=ibcoff                                                  11d27s18
            itmp2=itmp+nvvx*ncsz0                                       2d11s21
            ibcoff=itmp2+ncsz0*max(nfcn(iiii,1),nfcn(iiii,2))           10d31s20
            call enough('paraerik.  7',bc,ibc)
            do in=0,nerimat-1                                           10d31s20
             inp=in+1                                                   11d10s20
             iioff=0                                                     5d15s19
             do ist=1,2                                                   11d27s18
              if(nfcn(iiii,ist).gt.0)then                                11d28s18
               if(ist.eq.1)then                                            11d27s18
                fact=1d0                                                   11d27s18
               else                                                        11d27s18
                fact=-1d0                                                  11d27s18
               end if                                                      11d27s18
               do ii=0,nvvx*ncsz0-1
                bc(itmp+ii)=xnan
               end do
               if(isc.eq.isd)then                                       2d12s21
                if(ist.eq.1)then                                        2d12s21
                 nvv=(nbasispc(isc)*(nbasispc(isc)+1))/2                2d12s21
                 nvv0=(nbasisp(isc)*(nbasisp(isc)+1))/2                 2d12s21
                else                                                    2d12s21
                 nvv=(nbasispc(isc)*(nbasispc(isc)-1))/2                 2d11s21
                 nvv0=(nbasisp(isc)*(nbasisp(isc)-1))/2                  2d11s21
                end if                                                  2d12s21
               end if                                                   2d12s21
               ididit=0                                                 2d11s21
c     thus eraw is 1 LLLL, 2 LSLS, 3 SLSL, 4 SSSS
               if(in.eq.0.or.in.eq.3)then                               2d15s21
c     for in=0 or 3, get both ints with associated ci vector
                if(isc.eq.isd)then                                      2d11s21
                 jtmp=itmp                                                  11d27s18
                 do id=0,nbasisp(isd)-1                                    11d27s18
                  do ic=0,id-1                                              11d27s18
                   iadcd=ieraw(isc,isc)+ncsz0*(ic+nbasisp(isc)*(id       10d31s20
     $                 +nbasisp(isd)*in))                               10d31s20
                   iaddc=ieraw(isc,isc)+ncsz0*(id+nbasisp(isc)*(ic       10d31s20
     $                 +nbasisp(isc)*in))                               10d31s20
                   do ii=0,ncsz0-1                                            11d27s18
                    bc(jtmp+ii)=bc(iadcd+ii)+fact*bc(iaddc+ii)                 11d27s18
                   end do                                                   11d27s18
                   jtmp=jtmp+ncsz0                                       10d31s20
                  end do                                                    11d27s18
                  if(ist.eq.1)then                                          11d27s18
                   iaddd=ieraw(isc,isc)+ncsz0*(id+nbasisp(isc)*(id       10d31s20
     $                +nbasisp(isc)*in))                                11d2s20
                   do ii=0,ncsz0-1                                       10d31s20
                    bc(jtmp+ii)=bc(iaddd+ii)                                  11d27s18
                   end do                                                   11d27s18
                   jtmp=jtmp+ncsz0                                       10d31s20
                  end if                                                    11d27s18
                 end do                                                     11d27s18
                else                                                        11d27s18
                 do id=0,nbasisp(isd)-1                                    11d27s18
                  do ic=0,nbasisp(isc)-1                                   11d27s18
                   iadcd=ieraw(isc,isd)+ncsz0*(ic+nbasisp(isc)*(id       10d31s20
     $                 +nbasisp(isd)*in))                               10d31s20
                   iaddc=ieraw(isd,isc)+ncsz0*(id+nbasisp(isd)*(ic       10d31s20
     $                 +nbasisp(isc)*in))                               10d31s20
                   jtmp=itmp+ncsz0*(ic+nbasisp(isc)*id)                  10d31s20
                   do ii=0,ncsz0-1                                            11d27s18
                    bc(jtmp+ii)=bc(iadcd+ii)+fact*bc(iaddc+ii)                 11d27s18
                   end do                                                   11d27s18
                  end do                                                    11d27s18
                 end do                                                     11d27s18
                end if                                                      11d27s18
                if(min(ncsz,nfcn(iiii,ist),nvv0).gt.0)then               11d2s20
                 ididit=1                                               2d11s21
                 if(in.eq.0)then                                        2d11s21
                  jaddr=iaddr(iiii,ist,1)                               2d11s21
                 else                                                   2d11s21
                  jaddr=iaddr(iiii,ist,4)                               2d11s21
                 end if                                                 2d11s21
                 call dgemm('n','n',ncsz0,nfcn(iiii,ist),nvv0,1d0,       11d2s20
     $               bc(itmp),ncsz0,bc(jaddr),nvv0,0d0,bc(itmp2),ncsz0, 11d2s20
     d' paraerik.  1')
                end if                                                  2d11s21
               else
c     thus eraw is 1 LLLL, 2 LSLS, 3 SLSL, 4 SSSS
                jaddr=iaddr(iiii,ist,inp)                                 2d11s21
                iabcd=ieraw(isc,isd)+ncsz0*nbasisp(isc)*nbasisp(isd)*in 2d11s21
                iabdc=ieraw(isd,isc)+ncsz0*nbasisp(isd)*nbasisp(isc)*in 2d16s21
                if(min(ncsz,nfcn(iiii,ist),nvvx).gt.0)then              2d11s21
                 ididit=1                                               2d11s21
                 call dgemm('n','n',ncsz0,nfcn(iiii,ist),nvvx,1d0,      2d11s21
     $               bc(iabcd),ncsz0,bc(jaddr),nvvx,0d0,bc(itmp2),ncsz0,2d17s21
     d' paraerik.  2')
                 itransp=ibcoff                                         2d16s21
                 ibcoff=itransp+ncsz0*nvvx                              2d16s21
                 call enough('paraerik.  8',bc,ibc)
                 jtransp=itransp                                        2d16s21
                 do id=0,nbasisp(isd)-1                                 2d16s21
                  do ic=0,nbasisp(isc)-1                                2d16s21
                   jabdc=iabdc+ncsz0*(id+nbasisp(isd)*ic)               2d16s21
                   do j=0,ncsz0-1                                       2d16s21
                    bc(jtransp+j)=bc(jabdc+j)                           2d16s21
                   end do                                               2d16s21
                   jtransp=jtransp+ncsz0                                2d16s21
                  end do                                                2d16s21
                 end do                                                 2d16s21
                 if(in.eq.1)then                                        2d16s21
                  kaddr=iaddr(iiii,ist,3)                               2d16s21
                 else                                                   2d16s21
                  kaddr=iaddr(iiii,ist,2)                               2d16s21
                 end if                                                 2d16s21
                 call dgemm('n','n',ncsz0,nfcn(iiii,ist),nvvx,fact,     2d16s21
     $                bc(itransp),ncsz0,bc(kaddr),nvvx,1d0,bc(itmp2),   2d16s21
     $                ncsz0,                                            2d16s21
     d' paraerik.  3')
                 ibcoff=itransp                                         2d16s21
                end if                                                  2d11s21
               end if                                                   2d11s21
               if(ididit.ne.0)then                                      2d11s21
c     itmp2  is 1 LL, 2 LS, 3 SL, 4 SS thus
c     itmpab is 1 LL, 2 LS, 3 SL, 4 SS
                do ii=0,nfcn(iiii,ist)-1                                      11d27s18
                 iip=ii+iioff
                 do ib=0,nszb0-1                                          5d15s19
                  ibp=ib+iob                                             5d15s19
                  do ia=0,nsza0-1                                         5d15s19
                   iap=ia+ioa                                            5d15s19
                   iad1=itmpab(nfcn(iiii,3))+iip+mcd*(iap+nszaa0*(ibp   11d2s20
     $                  +nszba0*in))                                    11d2s20
                   iad2=itmp2+ia+nsza0*(ib+nszb0*ii)                     11d2s20
                   bc(iad1)=bc(iad1)+bc(iad2)                              12d7s18
                  end do                                                 5d15s19
                 end do                                                      11d27s18
                end do                                                       11d27s18
               end if                                                    5d11s19
              end if                                                     11d28s18
              iioff=iioff+nfcn(iiii,ist)                                 5d15s19
             end do                                                       11d27s18
            end do                                                      11d2s20
            ibcoff=itmp                                                  11d27s18
           end if                                                       11d28s18
          end do                                                        11d27s18
         end do                                                         11d27s18
c     end of ipassa loop
        end do                                                          5d12s10
c     end of ipassb loop
       end do                                                           5d12s10
c
c     complete ao to so transformation for ab                           1d10s11
c
       if(iapair(1,ibc(jb8)).gt.0)then                                  5d12s10
        do isd=1,nsymb                                                  5d12s10
         mcd=ntmpab(isd)                                                12d7s18
         do in=0,nerimat-1                                              11d2s20
          do ib=0,nszb0-1                                               11d2s20
           do ia=0,nszaa0-1                                             11d2s20
            it1=itmpab(isd)+mcd*(ia+nszaa0*(ib+nszba0*in))              11d2s20
            it2=itmpab(isd)+mcd*(ia+nszaa0*(ib+nszb0+nszba0*in))        11d2s20
            do icd=0,mcd-1                                              11d2s20
             sum=srh*(bc(it1+icd)+bc(it2+icd))                           12d7s18
             diff=srh*(-bc(it1+icd)+bc(it2+icd))                         12d7s18
             bc(it1+icd)=sum                                             12d7s18
             bc(it2+icd)=diff                                            12d7s18
            end do                                                       12d7s18
           end do                                                        12d7s18
          end do                                                         12d7s18
         end do                                                          5d12s10
        end do                                                          11d2s20
       end if                                                           5d12s10
       if(iapair(1,ibc(ja8)).gt.0)then                                  5d12s10
        do isd=1,nsymb                                                  5d12s10
         mcd=ntmpab(isd)                                                12d11s18
         do in=0,nerimat-1                                              11d2s20
          do ia=0,nsza0-1                                               11d2s20
           do ib=0,nszba0-1                                             11d2s20
            it1=itmpab(isd)+mcd*(ia+nszaa0*(ib+nszba0*in))              11d2s20
            it2=itmpab(isd)+mcd*(ia+nsza0+nszaa0*(ib+nszba0*in))        11d2s20
            do icd=0,mcd-1                                              11d2s20
             sum=srh*(bc(it1+icd)+bc(it2+icd))                           12d11s18
             diff=srh*(-bc(it1+icd)+bc(it2+icd))                         12d11s18
             bc(it1+icd)=sum                                             12d11s18
             bc(it2+icd)=diff                                            12d11s18
            end do                                                       12d11s18
           end do                                                        12d7s18
          end do                                                         5d12s10
         end do                                                          5d12s10
        end do                                                          11d2s20
       end if                                                           5d12s10
c
c     final storage                                                     1d10s11
c
       nmei=0                                                           12d11s18
       do ia=1,nszaa0                                                   8d26s15
        iaa=ia+ibc(ja4)                                                 2d3s12
        isa=isstor(iaa)                                                 2d3s12
        iaaa=ibstor(iaa)                                                12d11s18
        do ib=1,nszba0                                                  2d16s21
         ibb=ib+ibc(jb4)                                                2d3s12
         isb=isstor(ibb)                                                2d3s12
         ibbb=ibstor(ibb)                                               12d12s18
         ix=max(isa,isb)                                                12d11s18
         in=min(isa,isb)                                                12d11s18
         iiii=((ix*(ix-1))/2)+in                                        12d11s18
         mcd=ntmpab(nfcn(iiii,3))                                       12d7s18
c     itmpab is 1 LL, 2 LS, 3 SL, 4 SS
c     we want to transform v1 from so to mo basis. We require ibb ge
c     iaa, and isbv2 ge isbv1, so isbv2=max(isa,isb) and
c     isbv1=min(isa,isb).
c
         if(mcd.gt.0.and.ibb.ge.iaa)then                                12d12s18
           factab=1d0                                                   2d17s21
           if(iaa.eq.ibb)factab=0.5d0                                   2d17s21
c     isbv1 is a
           ii=iptoh(isa,isb)                                            2d17s21
           jorb=iorb(isa)+iaaa-1+nbasispc(isa)*(irefo(isa)+idoubo(isa)) 2d16s21
           iadd0=ii+mcd*nvirt(isa)*(ibbb-1)                             2d16s21
           do inn=0,nerimat-1                                           2d16s21
            if(inn.eq.2)jorb=jorb+nbasisp(isa)                          2d16s21
            jadd0=iadd0                                                 2d16s21
            if(inn.eq.1.or.inn.eq.3)jadd0=iadd0                         2d16s21
     $           +mcd*nvirt(isa)*nbasisp(isb)                           2d16s21
            jj=itmpab(nfcn(iiii,3))+mcd*(ia-1+nszaa0*(ib-1+nszba0*inn)) 2d16s21
            korb=jorb                                                   2d16s21
            do iv1=0,nvirt(isa)-1                                       2d16s21
             orb=bc(korb)*factab                                        2d17s21
             do j=0,mcd-1                                               2d16s21
              bc(jadd0+j)=bc(jadd0+j)+bc(jj+j)*orb                      2d17s21
             end do                                                     2d16s21
             korb=korb+nbasispc(isa)                                    2d16s21
             jadd0=jadd0+mcd                                            2d16s21
            end do                                                      2d16s21
           end do                                                       2d16s21
         end if                                                         12d12s18
        end do                                                          5d12s10
       end do                                                           5d12s10
       ibcoff=itmpab(1)                                                 12d7s18
c
c     end of do i loop                                                  11d27s18
c
      end do                                                            2d19s10
      call dws_gsumf(bc(ibptoh),needt)                                  12d12s18
      ihcdoub=1                                                         12d13s18
      do isb=1,nsymb                                                    12d12s18
       isbv12=multh(isb,isymmrci)                                       12d12s18
       do isbv1=1,nsymb                                                 12d12s18
        isbv2=multh(isbv1,isbv12)                                       12d12s18
        if(isbv2.ge.isbv1)then                                          12d12s18
         nvvh=nvirt(isbv1)*nbasisp(isbv2)                                11d2s20
         mcd=nrootu                                                     12d13s18
     $        *(ncd(1,isb,1)+ncd(1,isb,2)+ncd(2,isb,2)+ncd(3,isb,2))    12d13s18
         if(mcd*min(nvirt(isbv2),nvvh).gt.0)then                        11d5s20
          itmpa=ibcoff                                                  12d12s18
          itmpb=itmpa+mcd*nvirt(isbv1)*nvirt(isbv2)                     2d17s21
          ibcoff=itmpb+mcd*nvirt(isbv1)*nvirt(isbv2)                    2d17s21
          call enough('paraerik.  9',bc,ibc)
          nrow=mcd*nvirt(isbv1)                                         12d12s18
          fact=0d0                                                      2d17s21
          iorig=iptoh(isbv1,isbv2)                                      2d17s21
          jorb=iorb(isbv2)+nbasispc(isbv2)*(idoubo(isbv2)+irefo(isbv2)) 2d17s21
          do in=0,ntmpmat-1                                             2d17s21
           call dgemm('n','n',nrow,nvirt(isbv2),nbasisp(isbv2),1d0,     2d17s21
     $          bc(iorig),nrow,bc(jorb),nbasispc(isbv2),fact,           2d17s21
     $          bc(itmpa),nrow,                                         2d17s21
     d' paraerik.  4')
           fact=1d0                                                     2d17s21
           iorig=iorig+nrow*nbasisp(isbv2)                              2d17s21
           jorb=jorb+nbasisp(isbv2)                                     2d17s21
          end do                                                        2d17s21
          if(isbv1.eq.isbv2)then                                        12d12s18
           nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                       12d12s18
           nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       12d12s18
           jtmpa=itmpa                                                  2d17s21
           do ist=1,2                                                    12d13s18
            if(ist.eq.1)then                                             12d13s18
             nx=1                                                        12d13s18
             nvv=nvvs
            else                                                         12d13s18
             nx=3                                                        12d13s18
             nvv=nvvt
            end if                                                       12d13s18
            do is=1,nx                                                   12d13s18
             do iv2=0,nvirt(isbv2)-1                                    2d17s21
              if(ist.eq.1)then                                          2d17s21
               do iv1=0,iv2-1                                           2d17s21
                iv12=jtmpa+mcd*(iv1+nvirt(isbv1)*iv2)                   2d17s21
                iv21=jtmpa+mcd*(iv2+nvirt(isbv2)*iv1)                   2d17s21
                do j=0,ncd(is,isb,ist)-1                                2d17s21
                 do ir=0,nrootu-1                                       2d17s21
                  jhcdoub=ihcdoub+ndoub*ir                              2d17s21
                  jr=ir+nrootu*j                                        2d17s21
                  orig=hcdoub(jhcdoub)
                  hcdoub(jhcdoub)=hcdoub(jhcdoub)+bc(iv12+jr)           2d17s21
     $                 +bc(iv21+jr)                                     2d17s21
                 end do                                                 2d17s21
                 ihcdoub=ihcdoub+1                                      2d17s21
                end do                                                  2d17s21
               end do                                                   2d17s21
               ivv=jtmpa+mcd*(iv2+nvirt(isbv2)*iv2)                     2d17s21
               do j=0,ncd(is,isb,ist)-1                                 2d17s21
                do ir=0,nrootu-1                                        2d17s21
                 jhcdoub=ihcdoub+ndoub*ir                               2d17s21
                 jr=ir+nrootu*j                                         2d17s21
                 orig=hcdoub(jhcdoub)
                 hcdoub(jhcdoub)=hcdoub(jhcdoub)+bc(ivv+jr)*sr2         2d17s21
                end do                                                  2d17s21
                ihcdoub=ihcdoub+1                                       2d17s21
               end do                                                   2d17s21
              else                                                      2d17s21
               do iv1=0,iv2-1                                           2d17s21
                iv12=jtmpa+mcd*(iv1+nvirt(isbv1)*iv2)                   2d17s21
                iv21=jtmpa+mcd*(iv2+nvirt(isbv2)*iv1)                   2d17s21
                do j=0,ncd(is,isb,ist)-1                                2d17s21
                 do ir=0,nrootu-1                                       2d17s21
                  jhcdoub=ihcdoub+ndoub*ir                              2d17s21
                  jr=ir+nrootu*j                                        2d17s21
                  orig=hcdoub(jhcdoub)
                  hcdoub(jhcdoub)=hcdoub(jhcdoub)+bc(iv12+jr)           2d17s21
     $                 -bc(iv21+jr)                                     2d17s21
                 end do                                                 2d17s21
                 ihcdoub=ihcdoub+1                                      2d17s21
                end do                                                  2d17s21
               end do                                                   2d17s21
              end if                                                    2d17s21
             end do                                                     2d17s21
             jtmpa=jtmpa+nrootu*ncd(is,isb,ist)                         2d17s21
            end do                                                      2d17s21
           end do                                                       2d17s21
          else if(isbv1.ne.isbv2)then                                        2d17s21
           nrow=mcd*nvirt(isbv2)                                         12d12s18
           fact=0d0                                                      2d17s21
           iorig=iptoh(isbv2,isbv1)                                      2d17s21
           jorb=iorb(isbv1)+nbasispc(isbv1)*(idoubo(isbv1)+irefo(isbv1)) 2d17s21
           do in=0,ntmpmat-1                                             2d17s21
            call dgemm('n','n',nrow,nvirt(isbv1),nbasisp(isbv1),1d0,     2d17s21
     $          bc(iorig),nrow,bc(jorb),nbasispc(isbv1),fact,           2d17s21
     $          bc(itmpb),nrow,                                         2d17s21
     d' paraerik.  5')
            fact=1d0                                                     2d17s21
            iorig=iorig+nrow*nbasisp(isbv1)                             2d17s21
            jorb=jorb+nbasisp(isbv1)                                     2d17s21
           end do                                                        2d17s21
           jtmpa=itmpa                                                  2d17s21
           jtmpb=itmpb                                                  2d17s21
           do ist=1,2                                                    12d13s18
            if(ist.eq.1)then                                             12d13s18
             nx=1                                                        12d13s18
             fact=1d0                                                   2d17s21
            else                                                         12d13s18
             nx=3                                                        12d13s18
             fact=-1d0                                                  2d17s21
            end if                                                       12d13s18
            do is=1,nx                                                   12d13s18
             do iv2=0,nvirt(isbv2)-1                                      2d17s21
              do iv1=0,nvirt(isbv1)-1                                     2d17s21
               iv12=jtmpa+mcd*(iv1+nvirt(isbv1)*iv2)                      2d17s21
               iv21=jtmpb+mcd*(iv2+nvirt(isbv2)*iv1)                    2d17s21
               do j=0,ncd(is,isb,ist)-1                                 2d17s21
                do ir=0,nrootu-1                                        2d17s21
                 jhcdoub=ihcdoub+ndoub*ir                               2d17s21
                 jr=ir+nrootu*j                                         2d17s21
                 orig=hcdoub(jhcdoub)
                 hcdoub(jhcdoub)=hcdoub(jhcdoub)+bc(iv12+jr)            2d17s21
     $                +fact*bc(iv21+jr)                                 2d17s21
                end do                                                  2d17s21
                ihcdoub=ihcdoub+1                                       2d17s21
               end do                                                   2d17s21
              end do                                                    2d17s21
             end do                                                     2d17s21
             jtmpa=jtmpa+nrootu*ncd(is,isb,ist)                         2d17s21
             jtmpb=jtmpb+nrootu*ncd(is,isb,ist)                         2d17s21
            end do                                                      2d17s21
           end do                                                       2d17s21
          end if                                                        2d17s21
          ibcoff=itmpa                                                  2d17s21
         end if                                                         12d12s18
        end if                                                          12d12s18
       end do                                                           12d12s18
      end do                                                            12d12s18
      return
      end                                                               2d19s10
