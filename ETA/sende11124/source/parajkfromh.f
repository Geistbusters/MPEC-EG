c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parajkfromh(ipair,ngaus,ibdat,nhcolt,ihmat,nocc,icol,  6d22s10
     $     iorbx,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,     1d5s12
     $     ioooo,ionex,irefo,iter,idwsdeb,ncomp,nvirtc,i3xi,iflag,ih0,  5d31s18
     $     idoub,shift,nbasisp,idorel,ifull4x,i4x,i4xb,ibcb4,bc,ibc,    5d4s23
     $     isend,irecv,ipt,ipr,nder,iff)                                7d11s23
c
c     if iflag = 0, ih0, idoub, and shift are dummies.                  5d31s18
c     and irefo(i)=nocc(i).
c     if iflag = 1, then we are doing these transformation for mrci,    5d29s18
c     thus we exclude positron states and fold doubly occ orbs into     5d29s18
c     shift and h0 and 2e ints will contain no doubly occ orbs.         5d29s18
c     nocc(i) will include doubly occupied orbitals while irefo will    5d31s18
c     only be orbitals in active space.                                 5d31s18
c     ifull4x=1: cd is noc only.                                        8d18s23
c     ifull4x=2: d is noc, c is all. (normal operation)                 8d18s23
c     ifull4x=3: cd are both all.                                       8d18s23
c                                                                       5d29s18
      implicit real*8 (a-h,o-z)
      external second
      integer*8 isstor(1),ibstor(1)                                     4d14s21
      integer isend(*),irecv(*),ipt(*),ipr(*),idoub(*),irefo(*)         5d4s23
      logical ltest                                                     3d29s13
      include "common.store"
      include "common.hf"
      include "common.print"                                            1d6s20
      dimension nocc(1),iptoh(8,8,8),iapair(3,*),jmats(*),kmats(*),     6d11s21
     $     i4x(*),i4xb(*),nduse(8),iff(*)                               7d11s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension multh(8,8),iorbn(8),ihj(idbk),ilj(idbk),ihk(idbk),      12d13s11
     $     ilk(idbk),ijhalf(idbk),ikhalf(idbk),i3xi(*),                 7d26s16
     $     j1s(idbk),iph0(8),nbasisp(8),nduse2(8),                      8d18s23
     $     j1e(idbk),j2s(idbk),j2e(idbk),k1s(idbk),k1e(idbk),k2s(idbk), 6d23s10
     $     k2e(idbk),ioffab(8,8),joffab(8,8),iorbx(8),ioooo(1),         1d5s12
     $     ilo(idbk),iho(idbk),io1s(idbk),io1e(idbk),io2s(idbk),        1d5s12
     $     io2e(idbk),ionex(idbk),ix1s(idbk),ix1e(idbk),ix2s(idbk),     1d5s12
     $     ix2e(idbk),ilx(idbk),ihx(idbk),myh(idbk),nvirtc(8)           5d31s18
      common/ilimcm/ilimcode                                            5d7s19
      common/fnd4xcm/inv4x(2,8,8,8)                                     7d4s21
      if(ifull4x.eq.1)then                                              6d11s21
       do isb=1,nsymb                                                   6d11s21
        nduse(isb)=nocc(isb)                                            6d11s21
        nduse2(isb)=nocc(isb)                                           8d18s23
       end do                                                           6d11s21
      else if(ifull4x.eq.2)then                                         8d18s23
       do isb=1,nsymb                                                   6d11s21
        nduse(isb)=nocc(isb)                                            8d18s23
        nduse2(isb)=nbasdws(isb)                                        8d18s23
       end do                                                           6d11s21
      else if(ifull4x.eq.3)then                                         8d18s23
       do isb=1,nsymb                                                   6d11s21
        nduse(isb)=nbasdws(isb)                                         6d11s21
        nduse2(isb)=nbasdws(isb)                                        8d18s23
       end do                                                           6d11s21
      else                                                              8d18s23
       write(6,*)('unknown ifull4x in parajkfromh '),ifull4x            8d18s23
       stop 'parajkfromh'                                               8d18s23
      end if                                                            6d11s21
      ilimsave=ilimcode                                                 5d7s19
      ilimcode=1                                                        5d7s19
      call second(time1)                                                11d27s12
      ntmpmat=2                                                         10d28s20
      if(idorel.eq.0.or.idorel.eq.1)then                                10d28s20
       ntmpmat=1                                                        10d28s20
      end if                                                            10d30s20
c
c     in this version, do transpose of hmat, then transform full hmat
c     to mo basis, then store appropriate parts of the result.
c
c
      do i=1,idbk
       jmats(i)=0
       kmats(i)=0
       ionex(i)=0
       i3xi(i)=0                                                        7d26s16
      end do
      ibcoffo=ibcoff                                                    6d22s10
      iix=0
      nbasallp=0                                                        1d30s18
      do isb=1,nsymb
       nbasallp=nbasallp+nbasisp(isb)                                   1d30s18
       do isc=1,nsymb
        isbc=multh(isc,isb)
        do isd=1,nsymb
         iix=max(iix,iptoh(isd,isc,isb))
        end do
       end do
      end do
c
c     on pass 1, compute number of words going to each processor,
c     on pass 2, store by blocks for all2allvb and move
c
      ngaus7=ngaus*7                                                    5d13s10
      ngaus3=ngaus*3                                                    5d13s10
      nn=(ngaus*(ngaus+1))/2
      jbdat=ibdat-1
      do i=1,idbk                                                       11d30s12
       myh(i)=0                                                         11d30s12
      end do                                                            11d30s12
      do ip=1,mynprocg                                                  11d30s12
       isend(ip)=0                                                      11d30s12
       irecv(ip)=0                                                      11d30s12
      end do                                                            11d28s12
      nsend=0                                                           11d30s12
      nrecv=0                                                           11d30s12
      do ipass=1,2                                                      11d28s12
c
c     for sending
c
       ihmax=0
       do isb=1,nsymb                                                   11d30s12
        do isc=1,nsymb                                                  11d30s12
         isbc=multh(isc,isb)                                            11d30s12
         do isd=1,nsymb                                                 11d30s12
          isa=multh(isd,isbc)                                           11d30s12
          ii=iptoh(isd,isc,isb)                                         11d30s12
           ncd=nduse2(isc)*nduse(isd)                                   8d18s23
          if(ii.gt.0.and.ncd.gt.0)then                                  11d30s12
  323     format(4i3,4x,i5,es14.6,3i8)
           if(nhcol(ii).gt.0)then                                       11d30s12
            nhere=ncd/mynprocg                                            11d30s12
            nleft=ncd-nhere*mynprocg                                      11d30s12
            mcol=nhcol(ii)/ncd                                          11d30s12
            if(nleft.gt.0)nhere=nhere+1                                 11d30s12
            if(ipass.eq.1)then
             do ip=1,mynprocg                                           11d30s12
              isend(ip)=isend(ip)+mcol*nhere                            11d30s12
              if(ip.eq.nleft)nhere=nhere-1                              11d30s12
             end do                                                      11d30s12
            else
             ioff=ihcol(ii)-1-ncd
             ipt0=ipt(1)
             do ip=1,mynprocg                                           11d30s12
              do i=1,mcol
               do j=1,nhere
                iad=ioff+j+ncd*i
                bc(ipt(ip))=bc(iad)
                ihmax=max(ihmax,iad)
                ipt(ip)=ipt(ip)+1
               end do
              end do
              ioff=ioff+nhere
              if(ip.eq.nleft)nhere=nhere-1                              11d30s12
             end do
            end if                                                       11d30s12
           end if                                                        11d30s12
          end if                                                        11d30s12
         end do                                                         11d30s12
        end do                                                          11d30s12
       end do                                                           11d30s12
       if(ipass.eq.2)then                                               11d30s12
        ipt(1)=0                                                        11d30s12
        ipr(1)=0                                                        11d30s12
        do ip=2,mynprocg                                                11d30s12
         ipm=ip-1                                                       11d30s12
         ipt(ip)=ipt(ipm)+isend(ipm)                                    11d30s12
         ipr(ip)=ipr(ipm)+irecv(ipm)                                    11d30s12
        end do                                                          11d30s12
c      if(bc(132).ne.-132d0)then
c       call dws_sync
c       call dws_finalize
c       stop
c      end if
        call dws_all2allvb(bc(ibufs),isend,ipt,bc(ihmat),irecv,ipr)     11d30s12
        do ip=1,mynprocg                                                11d30s12
         ipr(ip)=ipr(ip)+ihmat                                          11d30s12
        end do                                                          11d30s12
        nsofar=myh(1)                                                   11d30s12
        myh(1)=ibufs                                                    11d30s12
        do i=2,idbk                                                     11d30s12
         nsum=nsofar+myh(i)                                             11d30s12
         myh(i)=ibufs+nsofar                                            11d30s12
         nsofar=nsum                                                    11d30s12
        end do                                                          11d30s12
        valnz=-1d0
        do i=0,nsofar-1
         bc(ibufs+i)=valnz
        end do
       end if
c
c     for receiving
c
       do ip=0,mynprocg-1                                               2d11s13
c`        call flushit
        do ii=1,1
         do isb=1,nsymb                                                  2d11s13
          nbasdwsb=nbasisp(isb)                                         10d30s20
          do isc=1,nsymb                                                 2d11s13
           isbc=multh(isc,isb)                                           2d11s13
           do isd=1,nsymb                                                2d11s13
            isa=multh(isd,isbc)                                          2d11s13
            nbasdwsa=nbasisp(isa)                                       10d30s20
            ntotrun=0
            ncd=nduse2(isc)*nduse(isd)                                  8d18s23
            iib=iptoh(isd,isc,isb)                                       2d12s13
 3323       format(4i2,4i8)
            iia=0
            if(max(iia,iib).gt.0.and.ncd.gt.0)then                       2d12s13
             nhere=ncd/mynprocg                                          2d11s13
             nleft=ncd-nhere*mynprocg                                    2d11s13
             if(mynowprog.lt.nleft)nhere=nhere+1                         3d5s13
             nsumx=0
             nsumx2=0
             nsumx3=0
             nsumx4=0
             do i=1+ip,nn,mynprocg                                       2d11s13
              jpair=ipair+2*(i-1)
              la=ibc(jbdat+ibc(jpair))
              na=(2*la+1)*ncomp                                          8d31s15
              lb=ibc(jbdat+ibc(jpair+1))
              nb=(2*lb+1)*ncomp                                          8d31s15
              ncsz=0                                                          5d13s10
              naa=na                                                          5d13s10
              nba=nb                                                          5d13s10
              if(iapair(1,ibc(jbdat+ibc(jpair)+ngaus7)).gt.0)then             5d13s10
               naa=naa*2                                                      5d13s10
              end if                                                          5d13s10
              if(iapair(1,ibc(jbdat+ibc(jpair+1)+ngaus7)).gt.0)then           5d13s10
               nba=nba*2                                                      5d13s10
              end if                                                          5d13s10
              naa0=naa/ncomp                                             8d31s15
              nba0=nba/ncomp                                             8d31s15
              do ia=1,naa0                                               8d31s15
               iaa=ia+ibc(jbdat+ibc(jpair)+ngaus3)                           6d4s10
               issa=isstor(iaa)                                               5d14s10
               iaaa=ibstor(iaa)                                              11d30s12
               do ib=1,nba0                                              8d31s15
                ibb=ib+ibc(jbdat+ibc(jpair+1)+ngaus3)                          6d4s10
                issb=isstor(ibb)                                                5d13s10
                ibbb=ibstor(ibb)                                               11d30s12
                if(ibb.le.iaa)then                                       2d12s13
                 if(issb.eq.isb.and.issa.eq.isa.and.iib.gt.0)then        2d12s13
                  if(issb.eq.issa)then                                   2d12s13
                   in=min(ibbb,iaaa)                                        11d30s12
                   ix=max(ibbb,iaaa)                                        11d30s12
                   ix=((ix*(ix-1))/2)+in-1                               2d26s13
                   nxn=(nbasdwsb*(nbasdwsb+1))/2                         8d31s15
                  else
                   nxn=nbasdwsb*nbasdwsa                                 8d31s15
                   ix=iaaa-1+nbasdwsa*(ibbb-1)                           8d31s15
                  end if                                                 2d11s13
                  if(ipass.eq.1)then                                     2d11s13
                   irecv(ip+1)=irecv(ip+1)+nhere*ntmpmat                10d30s20
                   myh(iib)=myh(iib)+nhere*ntmpmat                      10d30s20
                   nsumx=nsumx+1
                  else                                                   2d11s13
                   do in=0,ntmpmat-1                                    10d30s20
                    do j=0,nhere-1                                         2d11s13
                     iad=myh(iib)+ix+nxn*(j+nhere*in)                   10d30s20
                     bc(iad)=bc(ipr(ip+1))                                 2d11s13
                     ntotrun=ntotrun+1
                     ipr(ip+1)=ipr(ip+1)+1                                 2d11s13
                    end do                                              10d30s20
                   end do                                                 2d11s13
                  end if                                                 2d11s13
                 else if(issa.eq.isb.and.issb.eq.isa.and.                2d12s13
     $                iib.gt.0)then                                     2d12s13
                  nxn=nbasdwsb*nbasdwsa                                  8d31s15
                  ix=ibbb-1+nbasdwsa*(iaaa-1)                            8d31s15
                  if(ipass.eq.1)then                                     2d11s13
                   irecv(ip+1)=irecv(ip+1)+nhere*ntmpmat                10d30s20
                   myh(iib)=myh(iib)+nhere*ntmpmat                      10d30s20
                   nsumx=nsumx+1
                  else                                                   2d11s13
                   do in=0,ntmpmat-1                                    10d30s20
                    do j=0,nhere-1                                         2d11s13
                     iad=myh(iib)+ix+nxn*(j+nhere*in)                   10d30s20
                     bc(iad)=bc(ipr(ip+1))                                 2d11s13
                     ntotrun=ntotrun+1
                     ipr(ip+1)=ipr(ip+1)+1                                 2d11s13
                    end do                                                 2d11s13
                   end do                                               10d30s20
                  end if                                                 2d12s13
                 end if                                                  2d11s13
                end if                                                   2d12s13
               end do
              end do            !ia
             end do    !i
            end if                                                      9d2s15
           end do      !isd                                                 2d11s13
          end do       !isc                                                  2d11s13
         end do        !isb                                                  2d11s13
        end do         !ii                                                  2d11s13
       end do                   !ip                                                  2d15s13
       if(ipass.eq.1)then                                               11d28s12
        nsend=0                                                         11d28s12
        nrecv=0                                                         11d28s12
        do ip=1,mynprocg                                                11d30s12
         nsend=nsend+isend(ip)                                          11d30s12
         nrecv=nrecv+irecv(ip)                                          11d30s12
        end do                                                          11d28s12
        ngot=ibcoff+1-ihmat                                              11d30s12
        if(nrecv.gt.ngot)then                                           11d30s12
         more=nrecv-ngot                                                11d30s12
         ibcoff=ihmat+nrecv                                             12d3s12
         call enough('parajkfromh.  1',bc,ibc)
        end if                                                          11d30s12
        ibufs=ibcoff                                                    11d28s12
        ibcoff=ibufs+max(nsend,nrecv)                                   11d30s12
        call enough('parajkfromh.  2',bc,ibc)
        ipt(1)=ibufs                                                    11d30s12
        do ip=2,mynprocg                                                11d30s12
         ipm=ip-1                                                       11d28s12
         ipt(ip)=ipt(ipm)+isend(ipm)                                    11d30s12
        end do                                                          11d28s12
       end if                                                           11d28s12
       if(ipass.eq.1)then
        bc(ibcoff)=dfloat(nsend)
        bc(ibcoff+1)=dfloat(nrecv)
        call dws_gsumf(bc(ibcoff),2)
        nsendt=nint(bc(ibcoff))
        nrecvt=nint(bc(ibcoff+1))
       end if
        if(nsend.ne.nhcolt)then
         write(6,*)('send doesn''t match nhcolt: '),nsend,nhcolt
         call dws_sync
         call dws_finalize
         stop
        end if
        if(nsendt.ne.nrecvt)then                                        8d11s14
        write(6,*)('nsendt '),nsendt,(' vs. nrecvt '),nrecvt
         call dws_sync                                                  8d11s14
         call dws_finalize                                              8d11s14
         stop                                                           8d11s14
        end if                                                          8d11s14
      end do                                                            11d28s12
c
c     ihmat - raw ints: nduse2(c)*nduse(d);nn                           8d18s23
c     ibufs - ordered ints to send: nduse2(c)*nduse(d);nn               8d18s23
c     ihmat - ints that have been recieved: nduse2(c)*nduse(d);nn       8d18s23
c     ibufs - re-ordered ints: nduse2(c)*nduse(d);nn                    8d18s23
c     ihmat - transformed ints
c
      call second(time2)
      ipdeb=100
      ihlast=ihmat                                                      3d5s13
      do isb=1,nsymb                                                    11d30s12
       nbasdwsb=nbasisp(isb)                                            10d30s20
       do isc=1,nsymb                                                   11d30s12
        isbc=multh(isc,isb)                                             11d30s12
        do isd=1,nsymb                                                  11d30s12
         isa=multh(isd,isbc)                                            11d30s12
         iis=iptoh(isd,isc,isb)                                         12d21s12
         nbasdwsa=nbasisp(isa)                                          10d30s20
         ncd=nduse2(isc)*nduse(isd)                                     8d18s23
         nbasdwsf=nbasdws(isa)                                          5d29s18
         if(ncd.gt.0.and.iis.gt.0)then                                  12d21s12
          nhere=ncd/mynprocg                                            11d30s12
          nleft=ncd-nhere*mynprocg                                      11d30s12
          if(mynowprog.lt.nleft)nhere=nhere+1                           11d30s12
          if(nhere.gt.0)then                                            11d30s12
           if(isa.eq.isb)then                                           11d30s12
            itmp=ibcoff                                                 11d30s12
            itmp2=itmp+nbasdwsa*nbasdwsb                                10d30s20
            itmp3=itmp2+nbasdwsa*nbasdwsb                               10d30s20
            ibcoff=itmp3+nbasdwsa*nbasdwsb                              10d30s20
            call enough('parajkfromh.  3',bc,ibc)
            nxn=(nbasdwsa*(nbasdwsa+1))/2                               8d31s15
             nxnf=(nbasdws(isa)*(nbasdws(isa)+1))/2                     5d30s18
            ihcol(iis)=ihlast                                            3d5s13
            ihlast=ihlast+nxnf*nhere                                    5d30s18
            if(ihlast.gt.ibufs)then                                      3d5s13
             write(6,*)('ihlast now exceeds ibufs!!! '),ihlast,ibufs     3d5s13
             write(6,*)ihlast-ihmat
             call dws_sync                                               3d5s13
             call dws_finalize                                           3d5s13
             stop                                                        3d5s13
            end if                                                       3d5s13
            jtmp=itmp-1                                                 11d30s12
            jtmp3=itmp3-1                                               10d30s20
c
c     in myh we have triangle of fcns, nhere, then LL or SS
            do icol=1,nhere                                             11d30s12
             factr=0d0                                                  10d30s20
             do in=0,ntmpmat-1                                          10d30s20
              ij=0
              nan=0
              do i=1,nbasdwsa                                            8d31s15
               do j=1,i                                                  11d30s12
                iad0=myh(iis)+ij+nxn*(icol-1+nhere*in)                  10d30s20
                ij=ij+1                                                  2d12s13
                iad1=jtmp+j+nbasdwsa*(i-1)                               8d31s15
                iad2=jtmp+i+nbasdwsa*(j-1)                               8d31s15
                if(bc(iad0).ne.bc(iad0))nan=nan+1
                bc(iad1)=bc(iad0)                                        11d30s12
                bc(iad2)=bc(iad0)                                        11d30s12
               end do                                                    11d30s12
              end do                                                     11d30s12
              ipdeb=ipdeb+1
              if(in.eq.0)then                                           10d30s20
               jorb=iorbx(isa)                                          10d30s20
              else                                                      10d30s20
               jorb=iorbx(isa)+nbasdwsa                                 10d30s20
              end if                                                    10d30s20
              if(ipdeb.le.10)then
               write(6,*)('multiply from the right ')
               call prntm2(bc(itmp),nbasdwsa,nbasdwsa,
     $              nbasdwsa)
               if(icol.eq.1)then
                write(6,*)('by  ')
                call prntm2(bc(jorb),nbasdwsa,nbasdws(isa),             10d30s20
     $             nbasdwsa*ncomp)                                      10d30s20
               end if
              end if
              call dgemm('n','n',nbasdwsa,nbasdwsf,nbasdwsa,             5d29s18
     $            1d0,bc(itmp),nbasdwsa,bc(jorb),nbasdwsa*ncomp,        10d30s20
     $            0d0,bc(itmp2),nbasdwsa,                               5d29s18
     d' parajkfromh.  1')
              do i=1,nbasdwsf                                            5d29s18
               do j=1,nbasdwsa                                           8d31s15
                iad1=itmp2+j-1+nbasdwsa*(i-1)                            8d31s15
                iad2=itmp+i-1+nbasdwsf*(j-1)                             5d29s18
                bc(iad2)=bc(iad1)                                        11d30s12
               end do                                                    11d30s12
              end do                                                     11d30s12
              call dgemm('n','n',nbasdwsf,nbasdwsf,nbasdwsa,             5d29s18
     $            1d0,bc(itmp),nbasdwsf,bc(jorb),nbasdwsa*ncomp,        10d30s20
     $            factr,bc(itmp3),nbasdwsf,                             10d30s20
     d' parajkfromh.  2')
              factr=1d0                                                 10d30s20
              if(ipdeb.le.10)then
               write(6,*)('to yield '),iis,ihcol(iis),nxn
               call prntm2(bc(itmp3),nbasdwsf,nbasdwsf,
     $            nbasdwsf)
              end if
             end do                                                     10d30s20
             do i=1,nbasdwsf                                            5d29s18
              do j=1,i                                                  11d30s12
               iad0=ihcol(iis)+((i*(i-1))/2)+j-1+nxnf*(icol-1)          5d30s18
               iad1=jtmp3+j+nbasdwsf*(i-1)                              10d30s20
               bc(iad0)=bc(iad1)                                        11d30s12
               if(abs(bc(iad1)-2.8810023347952813d0).lt.1d-10)
     $              write(6,*)('4vinta '),iad1,bc(iad1)
              end do                                                    11d30s12
             end do                                                     11d30s12
            end do                                                      11d30s12
            ibcoff=itmp                                                 11d30s12
           else                                                         11d30s12
            nxn=nbasdwsa*nbasdwsb                                       8d31s15
             nxnf=nbasdws(isa)*nbasdws(isb)                             5d30s18
            ihcol(iis)=ihlast                                            3d5s13
            ihlast=ihlast+nxnf*nhere                                    5d30s18
            if(ihlast.gt.ibufs)then                                      3d5s13
             write(6,*)('ihlast now exceeds ibufs!!! '),ihlast,ibufs     3d5s13
             call dws_sync                                               3d5s13
             call dws_finalize                                           3d5s13
             stop                                                        3d5s13
            end if                                                       3d5s13
            itmp=ibcoff                                                 11d30s12
            itmp2=itmp+nxn                                              12d21s12
            itmp3=itmp2+nxn                                             10d30s20
            ibcoff=itmp3+nxn                                            10d30s20
            call enough('parajkfromh.  4',bc,ibc)
             nbasdwsaf=nbasdws(isa)                                     5d29s18
             nbasdwsbf=nbasdws(isb)                                     5d29s18
            do icol=1,nhere                                             11d30s12
             iads=ihcol(iis)+nxnf*(icol-1)                              5d30s18
             factr=0d0                                                  10d30s20
             jorbb=iorbx(isb)                                           10d30s20
             jorba=iorbx(isa)                                           10d30s20
             do in=0,ntmpmat-1                                          10d30s20
              iad=myh(iis)+nxn*(icol-1+nhere*in)                        10d30s20
              call dgemm('n','n',nbasdwsa,nbasdwsbf,nbasdwsb,            5d29s18
     $            1d0,bc(iad),nbasdwsa,bc(jorbb),nbasdwsb*ncomp,        10d30s20
     $            0d0,bc(itmp),nbasdwsa,                                8d31s15
     d' parajkfromh.  3')
              do i=1,nbasdwsa                                            8d31s15
               do j=1,nbasdwsbf                                          5d29s18
                iad1=iad+j-1+nbasdwsbf*(i-1)                             5d29s18
                iad2=itmp+i-1+nbasdwsa*(j-1)                             8d31s15
                bc(iad1)=bc(iad2)                                        11d30s12
               end do                                                    11d30s12
              end do                                                     11d30s12
              call dgemm('n','n',nbasdwsbf,nbasdwsaf,nbasdwsa,           5d29s18
     $            1d0,bc(iad),nbasdwsbf,bc(jorba),nbasdwsa*ncomp,       10d30s20
     $            factr,bc(itmp3),nbasdwsbf,                               5d29s18
     d' parajkfromh.  4')
              factr=1d0                                                 10d30s20
              jorbb=jorbb+nbasdwsb                                      10d30s20
              jorba=jorba+nbasdwsa                                      10d30s20
             end do                                                      10d30s20
             do i=1,nbasdwsbf                                           5d29s18
              do j=1,nbasdwsaf                                          5d29s18
               iad1=iads+j-1+nbasdwsaf*(i-1)                            5d29s18
               iad2=itmp3+i-1+nbasdwsbf*(j-1)                           10d30s20
               if(abs(bc(iad2)-2.8810023347952813d0).lt.1d-10)
     $              write(6,*)('4vintb '),iad2,bc(iad2)
               bc(iad1)=bc(iad2)                                        11d30s12
              end do                                                    11d30s12
             end do                                                     11d30s12
            end do                                                      11d30s12
            ibcoff=itmp                                                 11d30s12
           end if                                                       11d30s12
          end if                                                        11d30s12
         end if                                                         11d30s12
        end do                                                          11d30s12
       end do                                                           11d30s12
      end do                                                            11d30s12
      call second(time3)
      if(iflag.ne.0)then                                                5d29s18
       if(iprtr(10).ne.0)write(6,*)('computing shift and folding ')     10d28s20
       shift1=0d0                                                       5d29s18
       shift2=0d0                                                       5d29s18
       shift2k=0d0                                                      5d10s24
       jh0=ih0                                                          5d29s18
       kh0=ih0                                                          5d30s18
       mh0=0                                                            5d31s18
       ndoub=0                                                          5d31s18
       scale=1d0/dfloat(mynprocg)                                       5d29s18
       do isb=1,nsymb                                                   5d29s18
        ndoub=ndoub+idoub(isb)                                          5d31s18
        nbasdwsf=nbasisp(isb)*ncomp                                     2d15s19
        if(iprtr(10).ne.0)then                                          10d28s20
         write(6,*)('for symmetry block '),isb,(' h0 is ')
         call prntm2(bc(jh0),nbasdws(isb),nbasdws(isb),nbasdws(isb))     7d6s18
         if(nbasdws(isb).gt.1)then                                      10d28s20
          sym=0d0                                                        10d28s20
          do i=1,nbasdws(isb)-1                                          10d28s20
           do j=0,i-1                                                    10d28s20
            ji=jh0+j+nbasdws(isb)*i                                      10d28s20
            ij=jh0+i+nbasdws(isb)*j                                      10d28s20
            diff=bc(ji)-bc(ij)                                           10d28s20
            sym=sym+diff**2                                              10d28s20
           end do                                                        10d28s20
          end do                                                         10d28s20
          sym=sqrt(sym/dfloat(nbasdws(isb)**2))                         10d28s20
          write(6,*)('symmetry test '),sym                              10d28s20
         end if                                                         10d28s20
        end if                                                          10d28s20
c
c     first divide h0 by no. of procs, for we will be doing global sum  5d29s18
c     afterwords                                                        5d29s18
c                                                                       5d29s18
        do i=0,nbasdws(isb)*nbasdws(isb)-1                              6d6s18
         bc(jh0+i)=bc(jh0+i)*scale                                      5d29s18
        end do                                                          5d29s18
        do i=0,idoub(isb)-1                                             5d29s18
         ii=jh0+(nbasdws(isb)+1)*i                                      6d6s18
         shift1=shift1+2d0*bc(ii)                                       5d29s18
        end do                                                          5d29s18
c                                                                       5d30s18
c     compress h0                                                       5d30s18
c                                                                       5d30s18
        if(iprtr(10).ne.0)then                                          10d28s20
         write(6,*)('before compression, we have ')                      5d30s18
         call prntm2(bc(jh0),nbasdws(isb),nbasdws(isb),nbasdws(isb))
        end if                                                          10d28s20
        nhere=nbasdws(isb)-idoub(isb)                                   5d30s18
        do i=idoub(isb)+1,nbasdws(isb)                                  5d30s18
         im=i-idoub(isb)-1                                              5d30s18
         do j=idoub(isb)+1,nbasdws(isb)                                 5d30s18
          jm=j-idoub(isb)-1                                             5d30s18
          iad1=kh0+jm+nhere*im                                          5d30s18
          iad2=jh0+j-1+nbasdws(isb)*(i-1)                               6d6s18
          bc(iad1)=bc(iad2)                                             5d30s18
         end do                                                         5d30s18
        end do                                                          5d30s18
        iph0(isb)=kh0                                                   5d30s18
        if(iprtr(10).ne.0)then                                          10d28s20
         write(6,*)('after compression, we have '),kh0
         call prntm2(bc(kh0),nhere,nhere,nhere)                          5d30s18
         if(nhere.gt.1)then                                             10d28s20
          sym=0d0                                                        10d28s20
          do i=1,nhere-1                                                10d28s20
           do j=0,i-1                                                    10d28s20
            ji=kh0+j+nhere*i                                            10d28s20
            ij=kh0+i+nhere*j                                            10d28s20
            diff=bc(ji)-bc(ij)                                           10d28s20
            sym=sym+diff**2                                              10d28s20
           end do                                                        10d28s20
          end do                                                         10d28s20
          sym=sqrt(sym/dfloat(nhere**2))                                10d28s20
          write(6,*)('symmetry test '),sym                              10d28s20
         end if                                                         10d28s20
        end if                                                          10d28s20
        jh0=jh0+nbasdwsf*nbasdwsf                                       5d29s18
        kh0=kh0+nhere*nhere                                             5d30s18
        mh0=mh0+nhere*nhere                                             5d31s18
       end do                                                           5d29s18
       loopit=0
       do isb=1,nsymb                                                   5d30s18
        do isc=1,nsymb                                                  5d30s18
         isbc=multh(isc,isb)                                            5d30s18
         do isd=1,nsymb                                                 5d30s18
          loopit=loopit+1
          isa=multh(isd,isbc)                                           5d30s18
          iis=iptoh(isd,isc,isb)                                        5d30s18
          call ilimts(nduse(isd),nduse2(isc),mynprocg,mynowprog,ilh,    8d18s23
     $         ihh,i1s,i1e,i2s,i2e)                                     6d11s21
          ncd=ihh+1-ilh                                                 5d30s18
          if(iis.gt.0.and.ncd.gt.0)then                                 5d30s18
           if(isa.eq.isb)then                                           5d30s18
            nrowh=(nbasdws(isa)*(nbasdws(isa)+1))/2                     5d30s18
            iswitch=0                                                   5d30s18
           else                                                         5d30s18
            nrowh=nbasdws(isa)*nbasdws(isb)                             5d30s18
            iswitch=1                                                   5d30s18
           end if                                                       5d30s18
c
c     form J sum d (dd|vv')
c
           nh0=nbasdws(isa)-idoub(isa)                                  7d6s18
           if(isd.eq.isc)then                                           5d30s18
            jjs=ihcol(iis)-1                                            5d30s18
            i10=i1s                                                     5d30s18
            i1n=nduse(isd)                                              6d11s21
            do i2=i2s,i2e                                               5d30s18
             i2m=i2-1-idoub(isd)                                        5d31s18
             if(i2.eq.i2e)i1n=i1e                                       5d30s18
             do i1=i10,i1n                                              5d30s18
              if(i1.eq.i2.and.i2.le.idoub(isc))then                     5d30s18
               do i34=1,idoub(isa)                                      5d30s18
                iad=jjs+((i34*(i34+1))/2)
                shift2=shift2+bc(iad)*2d0                               5d30s18
               end do                                                   5d30s18
               do i4=0,nh0-1                                            7d6s18
                jh0=iph0(isa)+i4*nh0                                    7d6s18
                do i3=0,nh0-1                                           7d6s18
                 ix=max(i3,i4)+idoub(isa)                               7d6s18
                 in=min(i3,i4)+idoub(isa)                               7d6s18
                 iad=((ix*(ix+1))/2)+in+jjs+1                           7d6s18
                 orig=bc(jh0+i3)
                 bc(jh0+i3)=bc(jh0+i3)+2d0*bc(iad)                      7d6s18
                end do                                                  7d6s18
               end do                                                   7d6s18
              end if                                                    5d30s18
              jjs=jjs+nrowh                                             5d30s18
             end do                                                     5d30s18
             i10=1                                                      5d30s18
            end do                                                      5d30s18
           end if
c
c     for K: -sum d (dv|dv')
c
           if(isd.eq.isa)then                                           5d30s18
            nh0=nbasdws(isb)-idoub(isb)                                 7d6s18
c     ie (dv|dv')
            jjs=ihcol(iis)-1                                            5d30s18
            i10=i1s                                                     5d30s18
            i1n=nduse(isd)                                              6d11s21
            do i2=i2s,i2e                                               5d30s18
             i2m=i2-1                                                   5d31s18
             if(i2.eq.i2e)i1n=i1e                                       5d30s18
             do i1=i10,i1n                                              5d30s18
              if(i1.le.idoub(isd).and.i2.le.idoub(isc))then             5d30s18
               ix=max(i1,i2)                                            5d30s18
               in=min(i1,i2)                                            5d30s18
               ieq=((ix*(ix-1))/2)+in                                   5d30s18
               inot=i1+nbasdws(isa)*i2m                                 5d31s18
               iad=jjs+(inot-ieq)*iswitch+ieq                           5d30s18
               shift2k=shift2k-bc(iad)                                    5d30s18
              else if(i1.le.idoub(isd).and.i2.gt.idoub(isc))then        7d6s18
               jh0=iph0(isb)+nh0*(i2m-idoub(isc))                       7d6s18
               do i4=0,nh0-1                                            7d6s18
                i4p=i4+idoub(isb)                                       7d6s18
                i3=i1-1                                                 7d6s18
                inot=i3+nbasdws(isa)*i4p                                7d6s18
                ix=max(i4p,i3)                                          7d6s18
                in=min(i4p,i3)                                          7d6s18
                ieq=((ix*(ix+1))/2)+in                                  7d6s18
                iad=jjs+(inot-ieq)*iswitch+ieq+1                        7d6s18
                orig=bc(jh0+i4)                                         7d6s18
                bc(jh0+i4)=bc(jh0+i4)-bc(iad)                           7d6s18
               end do                                                   7d6s18
              end if                                                    5d30s18
              jjs=jjs+nrowh                                             5d30s18
             end do                                                     5d30s18
             i10=1                                                      5d30s18
            end do                                                      5d30s18
           else if(isd.eq.isb)then                                      5d30s18
c     ie (vd|dv')
            nh0=nbasdws(isa)-idoub(isa)                                 7d6s18
            jjs=ihcol(iis)-1                                            5d30s18
            i10=i1s                                                     5d30s18
            i1n=nduse(isd)                                              6d11s21
            do i2=i2s,i2e                                               5d30s18
             i2m=i2-1                                                   5d31s18
             if(i2.eq.i2e)i1n=i1e                                       5d30s18
             do i1=i10,i1n                                              5d30s18
              if(i1.le.idoub(isd).and.i2.le.idoub(isc))then             5d30s18
               iad=jjs+i2+nbasdws(isa)*(i1-1)                           5d31s18
               shift2k=shift2k-bc(iad)                                    5d30s18
              else if(i1.le.idoub(isd).and.i2.gt.idoub(isc))then        7d6s18
               i4=i1-1                                                  7d6s18
               jh0=iph0(isa)+nh0*(i2m-idoub(isc))                       7d6s18
               do i3=0,nh0-1                                            7d6s18
                inot=i3+idoub(isa)+nbasdws(isa)*i4                      7d6s18
                ix=max(i3+idoub(isa),i4)                                7d6s18
                in=min(i3+idoub(isa),i4)                                7d6s18
                ieq=((ix*(ix+1))/2)+in                                  7d6s18
                iad=jjs+(inot-ieq)*iswitch+ieq+1                        7d6s18
                orig=bc(jh0+i3)                                         7d6s18
                bc(jh0+i3)=bc(jh0+i3)-bc(iad)                           7d6s18
               end do                                                   7d6s18
              end if                                                    5d30s18
              jjs=jjs+nrowh                                             5d30s18
             end do                                                     5d30s18
             i10=1                                                      5d30s18
            end do                                                      5d30s18
           end if                                                       5d30s18
          end if                                                        5d30s18
         end do                                                         5d30s18
        end do                                                          5d30s18
       end do                                                           5d30s18
c
c     if there are no doubly occupied orbitals, then there will be no
c     shift.
c     if there are doubly occupied orbitals, then the new h0 will       5d31s18
c     be smaller than the original, so we can piggy back the shift      5d31s18
c     on the h0 global sum.                                             5d31s18
c                                                                       5d31s18
       if(ndoub.eq.0)then                                               5d31s18
        shift1=0d0                                                      5d31s18
        shift2=0d0                                                      5d31s18
        shift2k=0d0                                                     5d10s24
        call dws_gsumf(bc(ih0),mh0)                                     5d31s18
       else                                                             5d31s18
        bc(ih0+mh0)=shift1                                              5d31s18
        bc(ih0+mh0+1)=shift2+shift2k                                    5d10s24
        call dws_gsumf(bc(ih0),mh0+2)                                    5d31s18
        shift1=bc(ih0+mh0)                                              5d31s18
        shift2=bc(ih0+mh0+1)                                            5d31s18
       end if                                                           5d31s18
       if(iprtr(10).ne.0)then                                           6d12s24
        write(6,*)('at end of folding we have: ')
        do isb=1,nsymb
         nhere=nbasdws(isb)-idoub(isb)                                  5d30s18
         if(nhere.gt.0)then                                             6d12s24
          write(6,*)('for symmetry block '),isb,iph0(isb)               6d12s24
          call prntm2(bc(iph0(isb)),nhere,nhere,nhere)                  6d12s24
         end if                                                         6d12s24
        end do                                                          6d12s24
       end if                                                           6d12s24
       shift=shift1+shift2                                              5d31s18
      end if                                                            5d29s18
      do i=1,mynprocg                                                   3d4s13
       isend(i)=0                                                       3d4s13
       irecv(i)=0                                                       3d4s13
      end do                                                            3d4s13
c
c     if iflag = 0, then nvirtc is nvirtc. other wise it is nvirt, which5d31s18
c     is exactly what we what.
c
      ilimcode=ilimsave                                                 5d7s19
      do is=1,nsdlk1                                                    3d4s13
       nvc=nvirtc(isblk1(4,is))                                         8d31s15
       call ilimts(irefo(isblk1(3,is)),nvc,mynprocg,                    5d31s18
     $      mynowprog,ilx(is),ihx(is),ix1s(is),ix1e(is),ix2s(is),       3d4s13
     $      ix2e(is))                                                   3d4s13
      end do                                                            3d4s13
      do is=1,nsdlkk                                                    3d8s13
       nv3c=nvirtc(isblkk(3,is))                                        8d31s15
       nv4c=nvirtc(isblkk(4,is))                                        8d31s15
       call ilimts(nv3c,nv4c,mynprocg,                                  8d31s15
     $      mynowprog,ilk(is),ihk(is),k1s(is),k1e(is),k2s(is),k2e(is))
      end do                                                            3d8s13
      do is=1,nsdlk                                                     4d2s13
       nv3c=nvirtc(isblk(3,is))                                         8d31s15
       nv4c=nvirtc(isblk(4,is))                                         8d31s15
       call ilimts(nv3c,nv4c,mynprocg,                                  8d31s15
     $      mynowprog,ilj(is),ihj(is),j1s(is),j1e(is),j2s(is),j2e(is))  4d2s13
       call ilimts(irefo(isblk(3,is)),irefo(isblk(4,is)),mynprocg,      5d31s18
     $      mynowprog,ilo(is),iho(is),io1s(is),io1e(is),io2s(is),       4d2s13
     $      io2e(is))                                                   4d2s13
      end do                                                            4d2s13
      iso=ibcoff
      isj=iso+mynprocg*idbk
      isk=isj+mynprocg
      isx=isk+mynprocg
      iro=isx+mynprocg
      irj=iro+mynprocg*idbk
      irk=irj+mynprocg
      irx=irk+mynprocg
      ibcoff=irx+mynprocg
      is3=ibcoff                                                        7d27s16
      ir3=is3+mynprocg                                                  7d27s16
      ibcoff=ir3+mynprocg                                               7d27s16
      do i=iso,ibcoff-1                                                 6d1s18
       ibc(i)=0                                                         6d1s18
      end do                                                            6d1s18
      i4xsnd=ibcoff
      i4xrcv=i4xsnd+n4x
      ibcoff=i4xrcv+n4x
      i4ysnd=ibcoff
      i4yrcv=i4ysnd+nsdlk
      ibcoff=i4yrcv+nsdlk
      i4osnd=ibcoff-1
      i4orcv=i4osnd+nsdlk
      ijsnd=i4orcv+nsdlk
      ijrcv=ijsnd+nsdlk
      iksnd=ijrcv+nsdlkk
      ikrcv=iksnd+nsdlkk
      i3xsnd=ikrcv+nsdlkk
      i3xrcv=i3xsnd+nsdlk1
      i1xsnd=i3xrcv+nsdlk1
      i1xrcv=i1xsnd+nsdlk1
      ibcoff=i1xrcv+nsdlk1+1
      call enough('parajkfromh.  5',bc,ibc)
      do i=i4xsnd,ibcoff-1
       ibc(i)=0
      end do
      j4xsnd=i4xsnd-1
      j4xrcv=i4xrcv-1
      j4ysnd=i4ysnd-1
      j4yrcv=i4yrcv-1
      if(ifull4x.eq.3)then                                              8d18s23
       do i=1,n4x                                                        7d6s21
        i4xb(i)=1                                                        7d6s21
       end do                                                            7d6s21
      end if                                                            7d6s21
      do ipass=1,2                                                      3d4s13
       loopit=0
       do isb=1,nsymb                                                    1d10s11
         nbasdwsb=nbasdws(isb)                                          5d31s18
        do isc=1,nsymb                                                   1d10s11
         isbc=multh(isc,isb)                                             1d10s11
         jtmpl=ibcoff                                                   5d17s19
         jtmph=jtmpl+mynprocg                                           5d17s19
         ibcoff=jtmph+mynprocg                                          5d17s19
         call enough('parajkfromh.  6',bc,ibc)
         do ip=0,mynprocg-1                                             5d17s19
          call ilimts(nvirtc(isc),nvirtc(isb),mynprocg,ip,iilx,iihx,    5d17s19
     $         i1sx,i1ex,i2sx,i2ex)                                     5d17s19
          ibc(jtmpl+ip)=iilx-1                                          5d17s19
          ibc(jtmph+ip)=iihx-1                                          5d17s19
         end do                                                         5d17s19
         do isd=1,nsymb                                                  1d10s11
          itmpl=ibcoff                                                  5d15s19
          itmph=itmpl+mynprocg                                          5d15s19
          ibcoff=itmph+mynprocg                                         5d15s19
          call enough('parajkfromh.  7',bc,ibc)
          do ip=0,mynprocg-1                                            5d15s19
           call ilimts(irefo(isd),nvirt(isc),mynprocg,ip,iilx,iihx,     5d17s19
     $          i1sx,i1ex,i2sx,i2ex)                                    5d17s19
           ibc(itmpl+ip)=iilx-1                                          5d15s19
           ibc(itmph+ip)=iihx-1                                          5d15s19
          end do                                                        5d15s19
          isa=multh(isd,isbc)                                            1d10s11
          ktmpl=ibcoff                                                  5d17s19
          ktmph=ktmpl+mynprocg                                          5d17s19
          ibcoff=ktmph+mynprocg                                         5d17s19
          call enough('parajkfromh.  8',bc,ibc)
          do ip=0,mynprocg-1                                            5d15s19
           call ilimts(nvirtc(isc),nvirt(isa),mynprocg,ip,iilx,iihx,     5d17s19
     $          i1sx,i1ex,i2sx,i2ex)                                    5d17s19
           ibc(ktmpl+ip)=iilx-1                                         5d17s19
           ibc(ktmph+ip)=iihx-1                                         5d17s19
          end do                                                        5d15s19
           nbasdwsa=nbasdws(isa)                                        5d31s18
          ii=iptoh(isd,isc,isb)                                          1d10s11
          nbc=nduse2(isc)                                               8d18s23
          ilimcode=1                                                    5d7s19
          call ilimts(nduse(isd),nbc,mynprocg,mynowprog,ilh,            6d11s21
     $           ihh,i1s,i1e,i2s,i2e)                                   3d4s13
          ncd=ihh+1-ilh                                                 3d5s13
          loopit=loopit+1
          if(ii.gt.0.and.ncd.gt.0)then
c                                                                       11d17s17
c     nhcol is for uneven distribution by shells in paraeri. once we hav11d17s17
c     the data for final transformations, which we have already complete11d17s17
c     it is no longer relavant.                                         11d17s17
c                                                                       11d17s17
           if(nhcol(ii).ne.-1492)then                                        1d10s11
            if(isa.eq.isb)then
             nab=(nbasdwsa*(nbasdwsa+1))/2                              8d31s15
            else                                                         12d5s12
             nab=nbasdwsa*nbasdwsb                                      8d31s15
            end if                                                       12d5s12
            ntot=nduse(isd)*nvirtc(isc)                                 6d11s21
            nhere=ntot/mynprocg                                         3d5s13
            nleft=ntot-nhere*mynprocg                                   3d5s13
            nhere0=nhere                                                3d5s13
            if(nleft.gt.0)nhere=nhere+1                                 3d5s13
            ilimcode=ilimsave                                           5d7s19
            do is=1,nsdlk                                               4d2s13
             jso=iso
             if(isa.eq.isblk(3,is).and.isb.eq.isblk(4,is).and.isc.eq.    2d27s13
     $           isblk(1,is).and.isd.eq.isblk(2,is))then                2d27s13
              i10=i1s                                                    3d4s13
              i1x=nduse(isd)                                            6d11s21
              ioff=0
              do i2=i2s,min(i2e,nduse(isc))                             6d11s21
               if(i2.eq.i2e)i1x=i1e
               if(iflag.eq.0)then                                       5d31s18
                i2m=i2                                                  5d31s18
               else                                                     5d31s18
                i2m=i2-idoub(isc)                                       5d31s18
               end if                                                   5d31s18
               do i1=i10,i1x
                if(i1.le.nocc(isd).and.i2.le.nocc(isc))then             6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(isa.eq.isb)then
                  if(i1.ge.i2.and.i2m.gt.0)then                          5d31s18
                   if(ipass.eq.1)then
                    do ip=0,mynprocg-1
                     nva=nvirtc(isa)                                     8d31s15
                     if(ifull4x.gt.1)then                               8d18s23
                      call ilimts(nva,nva,mynprocg,ip,ilt,                8d31s15
     $                   iht,idum1,idum2,idum3,idum4)
                      isend(ip+1)=isend(ip+1)+iht+1-ilt
                      ibc(ijsnd+is)=ibc(ijsnd+is)+iht+1-ilt
                      ibc(isj+ip)=ibc(isj+ip)+iht+1-ilt
                     end if                                             8d18s23
                     call ilimts(irefo(isa),irefo(isa),mynprocg,ip,ilt,  5d31s18
     $                   iht,idum1,idum2,idum3,idum4)
                     isend(ip+1)=isend(ip+1)+iht+1-ilt
                     ibc(i4osnd+is)=ibc(i4osnd+is)+iht+1-ilt
                     ibc(jso+ip)=ibc(jso+ip)+iht+1-ilt
                    end do
                   else
                    do ip=0,mynprocg-1
                     nva=nvirtc(isa)                                     8d31s15
                     if(ifull4x.gt.1)then                               8d18s23
                      call ilimts(nva,nva,mynprocg,ip,ilt,                8d31s15
     $                   iht,i3s,i3e,i4s,i4e)
                      i30=i3s
                      i3x=nvirtc(isa)                                     8d31s15
                      do i4=i4s,i4e
                       ii4=i4+nocc(isa)
                       if(i4.eq.i4e)i3x=i3e
                       do i3=i30,i3x
                        ii3=i3+nocc(isa)
                        ix=max(ii3,ii4)
                        in=min(ii3,ii4)
                        iad=ihcol(ii)+((ix*(ix-1))/2)+in-1+ioff*nab
                        bc(ipt(ip+1))=bc(iad)
                        ipt(ip+1)=ipt(ip+1)+1
                       end do
                       i30=1
                      end do
                     end if                                             8d18s23
                     call ilimts(irefo(isa),irefo(isa),mynprocg,ip,ilt,  5d31s18
     $                   iht,i3s,i3e,i4s,i4e)
                     i30=i3s
                     i3x=irefo(isa)                                      5d31s18
                     if(iflag.eq.0)then                                  5d31s18
                      iaddd=0                                            5d31s18
                     else                                                5d31s18
                      iaddd=idoub(isa)                                   5d31s18
                     end if                                              5d31s18
                     do i4=i4s,i4e
                      if(i4.eq.i4e)i3x=i3e
                      do i3=i30,i3x
                       ix=max(i3,i4)+iaddd                               5d31s18
                       in=min(i3,i4)+iaddd                               5d31s18
                       iad=ihcol(ii)+((ix*(ix-1))/2)+in-1+ioff*nab
                       bc(ipt(ip+1))=bc(iad)
                       ipt(ip+1)=ipt(ip+1)+1
                      end do
                      i30=1
                     end do
                    end do
                   end if
                  end if
                 else if(min(i1m,i2m).gt.0)then                          5d31s18
                  if(ipass.eq.1)then
                   do ip=0,mynprocg-1
                    nva=nvirtc(isa)                                      8d31s15
                    nvb=nvirtc(isb)                                      8d31s15
                    if(ifull4x.gt.1)then                                8d18s23
                     call ilimts(nva,nvb,mynprocg,ip,ilt,                 8d31s15
     $                   iht,idum1,idum2,idum3,idum4)
                     isend(ip+1)=isend(ip+1)+iht+1-ilt
                     ibc(ijsnd+is)=ibc(ijsnd+is)+iht+1-ilt
                     ibc(isj+ip)=bc(isj+ip)+iht+1-ilt
                    end if                                              8d18s23
                    call ilimts(irefo(isa),irefo(isb),mynprocg,ip,ilt,   5d31s18
     $                   iht,idum1,idum2,idum3,idum4)
                    isend(ip+1)=isend(ip+1)+iht+1-ilt
                    ibc(i4osnd+is)=ibc(i4osnd+is)+iht+1-ilt
                    ibc(jso+ip)=ibc(jso+ip)+iht+1-ilt
                   end do
                  else
                   do ip=0,mynprocg-1
                    nva=nvirtc(isa)                                      8d31s15
                    nvb=nvirtc(isb)                                      8d31s15
                    if(ifull4x.gt.1)then                                8d18s23
                     call ilimts(nva,nvb,mynprocg,ip,ilt,                 8d31s15
     $                   iht,i3s,i3e,i4s,i4e)
                     i30=i3s
                     i3x=nvirtc(isa)                                      8d31s15
                     do i4=i4s,i4e
                      ii4=i4+nocc(isb)
                      if(i4.eq.i4e)i3x=i3e
                      do i3=i30,i3x
                       ii3=i3+nocc(isa)
                       iad=ihcol(ii)+ii3-1+nbasdwsa*(ii4-1)+ioff*nab      8d31s15
                       bc(ipt(ip+1))=bc(iad)
                       ipt(ip+1)=ipt(ip+1)+1
                      end do
                      i30=1
                     end do
                    end if                                              8d18s23
                    call ilimts(irefo(isa),irefo(isb),mynprocg,ip,ilt,   5d31s18
     $                  iht,i3s,i3e,i4s,i4e)
                    i30=i3s
                    i3x=irefo(isa)                                       5d31s18
                    if(iflag.eq.0)then                                   5d31s18
                     iadda=0                                             5d31s18
                     iaddb=0                                             5d31s18
                    else                                                 5d31s18
                     iadda=idoub(isa)                                    5d31s18
                     iaddb=idoub(isb)                                    5d31s18
                    end if                                               5d31s18
                    do i4=i4s,i4e
                     i4p=i4+iaddb-1                                      5d31s18
                     if(i4.eq.i4e)i3x=i3e
                     do i3=i30,i3x
                      i3p=i3+iadda-1                                     5d31s18
                      iad=ihcol(ii)+i3p+nbasdwsa*i4p+ioff*nab            5d31s18
                      bc(ipt(ip+1))=bc(iad)
                      ipt(ip+1)=ipt(ip+1)+1
                     end do
                     i30=1
                    end do
                   end do
                  end if
                 end if
                end if                                                  6d11s21
                ioff=ioff+1
               end do
               i10=1
              end do
             else if(isb.eq.isblk(3,is).and.isa.eq.isblk(4,is).and.      2d28s13
     $            isc.eq.isblk(1,is).and.isd.eq.isblk(2,is))then        2d28s13
              i10=i1s                                                    3d4s13
              i1x=nduse(isd)                                            6d11s21
              ioff=0
              do i2=i2s,min(i2e,nduse(isc))                             6d11s21
               if(iflag.eq.0)then                                       5d31s18
                i2m=i2                                                  5d31s18
               else                                                     5d31s18
                i2m=i2-idoub(isc)                                       5d31s18
               end if                                                   5d31s18
               if(i2.eq.i2e)i1x=i1e
               do i1=i10,i1x
                if(i1.le.nocc(isd).and.i2.le.nocc(isc))then             6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(ipass.eq.1.and.min(i1m,i2m).gt.0)then                5d31s18
                  do ip=0,mynprocg-1
                   if(ifull4x.gt.1)then                                 8d18s23
                    nva=nvirtc(isa)                                       8d31s15
                    nvb=nvirtc(isb)                                       8d31s15
                    call ilimts(nvb,nva,mynprocg,ip,ilt,                  8d31s15
     $                   iht,idum1,idum2,idum3,idum4)
                    ibc(ijsnd+is)=ibc(ijsnd+is)+iht+1-ilt
                    isend(ip+1)=isend(ip+1)+iht+1-ilt
                    ibc(isj+ip)=ibc(isj+ip)+iht+1-ilt
                   end if                                               8d18s23
                   call ilimts(irefo(isb),irefo(isa),mynprocg,ip,ilt,    5d31s18
     $                   iht,idum1,idum2,idum3,idum4)
                   isend(ip+1)=isend(ip+1)+iht+1-ilt
                   ibc(i4osnd+is)=ibc(i4osnd+is)+iht+1-ilt
                   ibc(jso+ip)=ibc(jso+ip)+iht+1-ilt
                  end do
                 else if(min(i1m,i2m).gt.0)then                          5d31s18
                  do ip=0,mynprocg-1
                   if(ifull4x.gt.1)then                                 8d18s23
                    nva=nvirtc(isa)                                       8d31s15
                    nvb=nvirtc(isb)                                       8d31s15
                    call ilimts(nvb,nva,mynprocg,ip,ilt,                  8d31s15
     $                  iht,i3s,i3e,i4s,i4e)
                    i30=i3s
                    i3x=nvirtc(isb)                                       8d31s15
                    do i4=i4s,i4e
                     ii4=i4+nocc(isa)
                     if(i4.eq.i4e)i3x=i3e
                     do i3=i30,i3x
                      ii3=i3+nocc(isb)
                      iad=ihcol(ii)+ii4-1+nbasdwsa*(ii3-1)+ioff*nab       8d31s15
                      bc(ipt(ip+1))=bc(iad)
                      ipt(ip+1)=ipt(ip+1)+1
                     end do
                     i30=1
                    end do
                   end if                                               8d18s23
                   call ilimts(irefo(isb),irefo(isa),mynprocg,ip,ilt,    5d31s18
     $                  iht,i3s,i3e,i4s,i4e)
                   i30=i3s
                   i3x=irefo(isb)                                        5d31s18
                   if(iflag.eq.0)then                                    5d31s18
                    iadda=0                                              5d31s18
                    iaddb=0                                              5d31s18
                   else                                                  5d31s18
                    iadda=idoub(isa)                                     5d31s18
                    iaddb=idoub(isb)                                     5d31s18
                   end if                                                5d31s18
                   do i4=i4s,i4e
                    i4p=i4-1+iadda                                       5d31s18
                    if(i4.eq.i4e)i3x=i3e
                    do i3=i30,i3x
                     i3p=i3-1+iaddb                                      5d31s18
                     iad=ihcol(ii)+i4p+nbasdwsa*i3p+ioff*nab             5d31s18
                     bc(ipt(ip+1))=bc(iad)
                     ipt(ip+1)=ipt(ip+1)+1
                    end do
                    i30=1
                   end do
                  end do
                 end if
                end if                                                  6d11s21
                ioff=ioff+1
               end do
               i10=1
              end do
             end if                                                      2d27s13
            end do                                                      4d2s13
            if(ifull4x.eq.3)then                                        8d18s23
             if(isc.ge.isd)then                                         7d4s21
              itricd=((isc*(isc-1))/2)+isd                              7d4s21
              iisx=max(isa,isb)                                         6d13s24
              isn=min(isa,isb)                                           7d4s21
              itriab=((iisx*(iisx-1))/2)+isn                            6d13s24
              if(itriab.le.itricd)then                                  7d4s21
               if(isc.eq.isd)then                                       7d5s21
                nrowcd=(nvirtc(isc)*(nvirtc(isc)+1))/2                  7d5s21
               else                                                     7d5s21
                nrowcd=nvirtc(isc)*nvirtc(isd)                          7d5s21
               end if                                                   7d5s21
               if(isa.eq.isb)then                                       7d5s21
                nrowab=(nvirtc(isa)*(nvirtc(isa)+1))/2                  7d5s21
                iswab=0                                                 7d5s21
               else                                                     7d5s21
                nrowab=nvirtc(isa)*nvirtc(isb)                          7d5s21
                iswab=1                                                 7d5s21
               end if                                                   7d5s21
               if(itriab.eq.itricd)then                                 7d5s21
                nabcd=(nrowab*(nrowab+1))/2                             7d5s21
                isw=0                                                   7d5s21
               else                                                     7d5s21
                nabcd=nrowab*nrowcd                                     7d5s21
                isw=1                                                   7d5s21
               end if                                                   7d5s21
               nhere=nabcd/mynprocg                                     7d5s21
               nleft=nabcd-nhere*mynprocg                               7d5s21
               nhere0=nhere                                               3d5s13
               if(nleft.gt.0)nhere=nhere+1                                3d5s13
c     abcd: i1+na*(i2+nb*(i3+nc*i4)
c     abab: ((ix*(ix-1))/2)+in, ix=i1+na*i2 ge in=i3+na*i4
c     aabb: ((i1*(i1-1))/2)+i2+naa*[((i3*(i3-1))/2)+i4]
c     aaaa: ((ix*(ix-1))/2)+in, ix=((i1*(i1-1))/2)+i2, in=((i3*(i3-1))/2+i4.
               i2eu=inv4x(1,isc,isd,iisx)                               6d13s24
               i10=i1s                                                  7d4s21
               i1x=nduse(isd)                                           7d4s21
               ioff=0                                                   7d5s21
               do ii2=i2s,i2e                                           7d4s21
                i2=ii2-nocc(isc)                                        11d8s23
                i2m=i2-1                                                7d6s21
                if(ii2.eq.i2e)i1x=i1e                                   7d4s21
                do ii1=i10,i1x                                          7d4s21
                 i1=ii1-nocc(isd)                                       11d8s23
                 i1m=i1-1                                               7d6s21
                 if(i1.gt.0.and.i2.gt.0.and.
     $                (isd.ne.isc.or.i2.ge.i1))then                     7d4s21
                  if(isc.eq.isd)then                                    7d4s21
                   itri12=((i2*(i2-1))/2)+i1-1                          7d5s21
                  else                                                  7d4s21
                   itri12=i2-1+nvirtc(isc)*(i1-1)                       7d5s21
                  end if                                                7d4s21
                  naabb=nrowab-1                                        7d5s21
                  if(itricd.eq.itriab)naabb=itri12                      7d4s21
                  jtriab=0
                  if(isa.ge.isb)then                                    7d5s21
                   do ib=0,nvirtc(isb)-1                                7d5s21
                    ibp=ib+nocc(isb)                                    11d8s23
c     jtriab+ia+1 le naabb, ia le naabb-1-jtriab
                    imxa=naabb-jtriab                                   7d5s21
                    itopa=min(imxa,ib+iswab*(nvirtc(isa)-1-ib))          7d5s21
                    do ia=0,itopa                                        7d5s21
                     iap=ia+nocc(isa)                                   11d8s23
                     ix=max(itri12,jtriab)                               7d5s21
                     in=min(itri12,jtriab)                               7d5s21
                     itri=((ix*(ix+1))/2)+in                             7d5s21
                     irec=jtriab+nrowab*itri12                           7d5s21
                     qq=get4int(i4xb,isc,isd,isa,isb,                   7d6s21
     $                     i2m,i1m,ia,ib,nvirtc,icol,bc,ibc)            11d10s22
                     icol=icol-1                                        7d6s21
                     jtriab=jtriab+1                                     7d5s21
                     ip=icol/nhere                                           3d5s13
                     if(ip.ge.nleft)then                                     3d5s13
                      ip=((icol-nleft*nhere)/nhere0)+nleft                   3d5s13
                     end if                                                  3d5s13
                     ipp=ip+1
                     if(ipass.eq.1)then                                  7d5s21
                      if(ipp.lt.1.or.ipp.gt.mynprocg)then
                       write(6,*)('ippa error !!! '),icol,ip,nabcd,
     $                      nhere,
     $                      nleft,nhere0
                      end if
                      isend(ipp)=isend(ipp)+1                            7d5s21
                      ibc(j4xsnd+i2eu)=ibc(j4xsnd+i2eu)+1
                     else                                                7d5s21
                      itrih=((ibp*(ibp+1))/2)+iap                        7d5s21
                      irech=iap+nbasdwsa*ibp                             7d5s21
                      iad=ihcol(ii)+itrih+iswab*(irech-itrih)+nab*ioff   7d5s21
                      bc(ipt(ipp))=bc(iad)                               7d5s21
                      ipt(ipp)=ipt(ipp)+1                                7d5s21
                     end if                                              7d5s21
                    end do                                               7d5s21
                   end do                                                7d5s21
                  else                                                  7d5s21
                   do ia=0,nvirtc(isa)-1                                7d5s21
                    iap=ia+nocc(isa)                                    11d8s23
                    imxb=naabb-jtriab                                   7d5s21
                    ibtop=min(nvirtc(isb)-1,imxb)                       7d5s21
                    do ib=0,ibtop
                     ibp=ib+nocc(isb)                                   11d8s23
                     ix=max(itri12,jtriab)                               7d5s21
                     in=min(itri12,jtriab)                               7d5s21
                     itri=((ix*(ix+1))/2)+in                             7d5s21
                     irec=jtriab+nrowab*itri12                           7d5s21
                     icol=itri+isw*(irec-itri)                          7d6s21
                     qq=get4int(i4xb,isc,isd,isb,isa,                   7d6s21
     $                     i2m,i1m,ib,ia,nvirtc,icol,bc,ibc)            11d10s22
                     icol=icol-1                                        7d6s21
                     jtriab=jtriab+1                                    7d5s21
                     ip=icol/nhere                                           3d5s13
                     if(ip.ge.nleft)then                                     3d5s13
                      ip=((icol-nleft*nhere)/nhere0)+nleft                   3d5s13
                     end if                                                  3d5s13
                     ipp=ip+1
                     if(ipass.eq.1)then                                  7d5s21
                      if(ipp.lt.1.or.ipp.gt.mynprocg)then
                       write(6,*)('ippb error !!! '),icol,ip,nabcd,
     $                      nhere,
     $                      nleft,nhere0
                      end if
                      isend(ipp)=isend(ipp)+1                            7d5s21
                      ibc(j4xsnd+i2eu)=ibc(j4xsnd+i2eu)+1
                     else                                                7d5s21
                      irech=iap+nbasdwsa*ibp                             7d5s21
                      iad=ihcol(ii)+irech+nab*ioff                      7d5s21
                      bc(ipt(ipp))=bc(iad)                               7d5s21
                      ipt(ipp)=ipt(ipp)+1                                7d5s21
                     end if                                              7d5s21
                    end do                                              7d5s21
                   end do                                               7d5s21
                  end if                                                7d5s21
                 end if                                                 7d4s21
                 ioff=ioff+1                                            7d5s21
                end do                                                  7d4s21
                i10=1                                                   7d4s21
               end do                                                   7d4s21
              end if                                                    7d4s21
             end if                                                     7d4s21
            end if                                                      7d4s21
            if(ifull4x.eq.1)then                                        8d21s23
             do is=1,nsdlk1                                             8d21s23
              if(isc.eq.isblk1(1,is).and.isd.eq.isblk1(2,is))then       8d21s23
               if(isc.eq.isd)then                                       8d21s23
                iswitch=0                                               8d21s23
               else                                                     8d21s23
                iswitch=1                                               8d21s23
               end if                                                   8d21s23
               if(isa.eq.isblk1(3,is))then                              8d21s23
                i10=i1s                                                  8d21s23
                i1n=nduse(isd)                                           8d21s23
                jhcol=ihcol(ii)                                         8d21s23
                do i2=i2s,i2e                                            8d21s23
                 if(i2.eq.i2e)i1n=i1e                                    8d21s23
                 i2m=i2-1                                                8d21s23
                 if(iflag.ne.0)i2m=i2m-idoub(isc)                       8d21s23
                 do i1=i10,i1n                                           8d21s23
                  i1m=i1-1                                               8d21s23
                  if(iflag.ne.0)i1m=i1m-idoub(isd)                      8d21s23
                  i1top=i2m+iswitch*(nduse(isd)-1-i2m)                   8d21s23
                  if(min(i2m,i1m).ge.0.and.i1m.le.i1top)then             8d21s23
                   if(ipass.eq.1)then                                    8d21s23
                    do ip=0,mynprocg-1                                  8d21s23
                     call ilimts(irefo(isblk1(3,is)),                   8d21s23
     $                    nvirt(isblk1(4,is)),mynprocg,ip,ll,lh,l1s,l1e,8d21s23
     $                    l2s,l2e)                                      8d21s23
                     nhere=lh+1-ll                                      8d21s23
                     isend(ip+1)=isend(ip+1)+nhere                      8d21s23
                    end do                                              8d21s23
                   else                                                  8d21s23
                    do ip=0,mynprocg-1                                  8d21s23
                     call ilimts(irefo(isblk1(3,is)),                   8d21s23
     $                    nvirt(isblk1(4,is)),mynprocg,ip,ll,lh,l1s,l1e,8d21s23
     $                    l2s,l2e)                                      8d21s23
                     j10=l1s                                            8d21s23
                     j1n=irefo(isblk1(3,is))                            8d21s23
                     do j2=l2s,l2e                                      8d21s23
                      if(j2.eq.l2e)j1n=l1e                              8d21s23
                      j2m=j2-1+nocc(isblk1(4,is))                       8d21s23
                      do j1=j10,j1n                                     8d21s23
                       j1m=j1-1                                         8d21s23
                       irec=j1m+nbasdwsa*j2m                            8d21s23
                       itri=((j2m*(j2m+1))/2)+j1m                       8d21s23
                       itri=jhcol+itri+iswitch*(irec-itri)              8d21s23
                       bc(ipt(ip+1))=bc(itri)                           8d21s23
                       ipt(ip+1)=ipt(ip+1)+1                            8d21s23
                      end do                                            8d21s23
                      j10=1                                             8d21s23
                     end do                                             8d21s23
                    end do                                              8d21s23
                   end if                                                8d21s23
                  end if                                                 8d21s23
                  jhcol=jhcol+nab                                       8d21s23
                 end do                                                  8d21s23
                 i10=1                                                   8d21s23
                end do                                                   8d21s23
               else if(isa.eq.isblk1(4,is))then                         8d21s23
                i10=i1s                                                  8d21s23
                i1n=nduse(isd)                                           8d21s23
                jhcol=ihcol(ii)                                         8d21s23
                do i2=i2s,i2e                                            8d21s23
                 if(i2.eq.i2e)i1n=i1e                                    8d21s23
                 i2m=i2-1                                                8d21s23
                 if(iflag.ne.0)i2m=i2m-idoub(isc)                       8d21s23
                 do i1=i10,i1n                                           8d21s23
                  i1m=i1-1                                               8d21s23
                  if(iflag.ne.0)i1m=i1m-idoub(isd)                      8d21s23
                  i1top=i2m+iswitch*(nduse(isd)-1-i2m)                   8d21s23
                  if(min(i2m,i1m).ge.0.and.i1m.le.i1top)then             8d21s23
                   if(ipass.eq.1)then                                    8d21s23
                    do ip=0,mynprocg-1                                  8d21s23
                     call ilimts(irefo(isblk1(3,is)),                   8d21s23
     $                    nvirt(isblk1(4,is)),mynprocg,ip,ll,lh,l1s,l1e,8d21s23
     $                    l2s,l2e)                                      8d21s23
                     nhere=lh+1-ll                                      8d21s23
                     isend(ip+1)=isend(ip+1)+nhere                      8d21s23
                    end do                                              8d21s23
                   else                                                  8d21s23
                    do ip=0,mynprocg-1                                  8d21s23
                     call ilimts(irefo(isblk1(3,is)),                   8d21s23
     $                    nvirt(isblk1(4,is)),mynprocg,ip,ll,lh,l1s,l1e,8d21s23
     $                    l2s,l2e)                                      8d21s23
                     j10=l1s                                            8d21s23
                     j1n=irefo(isblk1(3,is))                            8d21s23
                     do j2=l2s,l2e                                      8d21s23
                      if(j2.eq.l2e)j1n=l1e                              8d21s23
                      j2m=j2-1+nocc(isblk1(4,is))                       8d21s23
                      do j1=j10,j1n                                     8d21s23
                       j1m=j1-1                                         8d21s23
                       irec=jhcol+j2m+nbasdwsa*j1m                      8d21s23
                       bc(ipt(ip+1))=bc(irec)                           8d21s23
                       ipt(ip+1)=ipt(ip+1)+1                            8d21s23
                      end do                                            8d21s23
                      j10=1                                             8d21s23
                     end do                                             8d21s23
                    end do                                              8d21s23
                   end if                                                8d21s23
                  end if                                                 8d21s23
                  jhcol=jhcol+nab                                       8d21s23
                 end do                                                  8d21s23
                 i10=1                                                   8d21s23
                end do                                                   8d21s23
               end if                                                   8d21s23
              end if                                                    8d21s23
             end do                                                     8d21s23
            else                                                        8d21s23
            do is=1,nsdlk1                                               2d28s13
             ididit=0
             if(isa.eq.isblk1(1,is).and.isb.eq.isblk1(2,is).and.         2d28s13
     $            isd.eq.isblk1(3,is).and.isc.eq.isblk1(4,is))then         2d28s13
              ididit=1
              ntot=irefo(isd)*nvirtc(isc)                               5d31s18
              nhere=ntot/mynprocg                                        3d5s13
              nleft=ntot-nhere*mynprocg                                  3d5s13
              nhere0=nhere                                               3d5s13
              if(nleft.gt.0)nhere=nhere+1                                3d5s13
              i10=i1s                                                    3d4s13
              i1x=nduse(isd)                                            6d11s21
              ioff=0
              do ii2=i2s,i2e                                             3d5s13
               i2=ii2-nocc(isc)
               if(ii2.eq.i2e)i1x=i1e                                      3d4s13
               do i1=i10,i1x                                             3d4s13
                if(i1.le.nocc(isd))then                                 6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(i2.gt.0.and.i1m.gt.0)then                            5d31s18
                  itry=i1m-1+irefo(isd)*(i2-1)                           5d31s18
                  do iq=0,mynprocg-1                                     5d15s19
                   if(itry.ge.ibc(itmpl+iq).and.itry.le.ibc(itmph+iq))   5d15s19
     $                 then                                             5d15s19
                    go to 2020                                           5d15s19
                   end if                                                5d15s19
                  end do                                                 5d15s19
                  write(6,*)('iq failure ')
                  call dws_sync
                  call dws_finalize
                  stop
 2020             continue
                  ip=iq                                                  5d17s19
                  if(ip.lt.0.or.ip.ge.mynprocg)then                      6d1s18
                   write(6,*)('ip: '),ip,ntot,itry,nhere,nleft
                  end if
                  if(ipass.eq.1)then                                      3d5s13
                   if(isa.eq.isb)then                                      3d5s13
                    nxn=((irefo(isa)*(irefo(isa)+1))/2)                  5d31s18
                    isend(ip+1)=isend(ip+1)+nxn                           3d7s13
                    ibc(i1xsnd+is)=ibc(i1xsnd+is)+nxn
                    ibc(isx+ip)=ibc(isx+ip)+nxn                           3d7s13
                    if(ifull4x.gt.1)then                                8d18s23
                     nxn=(nvirtc(isa)*(nvirtc(isa)+1))/2                  7d26s16
                     isend(ip+1)=isend(ip+1)+nxn                          7d26s16
                     ibc(i3xsnd+is)=ibc(i3xsnd+is)+nxn
                     ibc(is3+ip)=ibc(is3+ip)+nxn                          7d27s16
                    end if                                              8d18s23
                   else                                                   3d5s13
                    nxn=irefo(isa)*irefo(isb)                            5d31s18
                    isend(ip+1)=isend(ip+1)+nxn                           3d7s13
                    ibc(i1xsnd+is)=ibc(i1xsnd+is)+nxn
                    ibc(isx+ip)=ibc(isx+ip)+nxn
                    if(ifull4x.gt.1)then                                8d18s23
                     nxn=nvirtc(isa)*nvirtc(isb)                          7d26s16
                     isend(ip+1)=isend(ip+1)+nxn                          7d26s16
                     ibc(i3xsnd+is)=ibc(i3xsnd+is)+nxn
                     ibc(is3+ip)=ibc(is3+ip)+nxn                          7d27s16
                    end if                                              8d18s23
                   end if                                                 3d5s13
                  else                                                    3d5s13
                   iad=ihcol(ii)+nab*ioff                                 3d5s13
                   lad=iad                                               7d26s16
                   if(isa.eq.isb)then                                     3d5s13
                    if(iflag.eq.0)then                                   5d31s18
                     iadda=0                                             5d31s18
                    else                                                 5d31s18
                     iadda=idoub(isa)                                    5d31s18
                    end if                                               5d31s18
                    do i3=1,irefo(isa)                                   5d31s18
                     i3p=i3+iadda                                        5d31s18
                     do i4=1,i3                                            3d5s13
                      i4p=i4+iadda                                       5d31s18
                      jad=iad+((i3p*(i3p-1))/2)+i4p-1                    5d31s18
                      bc(ipt(ip+1))=bc(jad)                              5d31s18
                      ipt(ip+1)=ipt(ip+1)+1                               3d5s13
                     end do                                                3d5s13
                    end do                                                 3d5s13
                    kad=lad-1                                            7d26s16
                    if(ifull4x.gt.1)then                                8d18s23
                     do i3=1,nvirtc(isa)                                  7d26s16
                      i3p=i3+nocc(isa)                                    5d31s18
                      do i4=1,i3                                            3d5s13
                       i4p=i4+nocc(isa)                                   5d31s18
                       jad=kad+((i3p*(i3p-1))/2)+i4p                      7d26s16
                       bc(ipt(ip+1))=bc(jad)                              7d26s16
                       ipt(ip+1)=ipt(ip+1)+1                               3d5s13
                      end do                                                3d5s13
                     end do                                                 3d5s13
                    end if                                              8d18s23
                   else                                                   3d5s13
                    if(iflag.eq.0)then                                   5d31s18
                     iadda=0                                             5d31s18
                     iaddb=0                                             5d31s18
                    else                                                 5d31s18
                     iadda=idoub(isa)                                    5d31s18
                     iaddb=idoub(isb)                                    5d31s18
                    end if                                               5d31s18
                    do i3=1,irefo(isa)                                   5d31s18
                     i3p=i3-1+iadda                                      5d31s18
                     do i4=1,irefo(isb)                                  5d31s18
                      i4p=i4-1+iaddb                                     5d31s18
                      iad2=i3p+nbasdwsa*i4p+iad                          5d31s18
                      bc(ipt(ip+1))=bc(iad2)                                3d5s13
                      ipt(ip+1)=ipt(ip+1)+1                               3d5s13
                     end do                                               3d5s13
                    end do                                                3d5s13
                    if(ifull4x.gt.1)then                                8d18s23
                     do i3=1,nvirtc(isa)                                  7d26s16
                      i3p=i3+nocc(isa)                                    5d31s18
                      do i4=1,nvirtc(isb)                                 7d26s16
                       i4p=i4+nocc(isb)                                   5d31s18
                       iad2=i3p-1+nbasdwsa*(i4p-1)+iad                    7d26s16
                       bc(ipt(ip+1))=bc(iad2)                             7d26s16
                       ipt(ip+1)=ipt(ip+1)+1                              7d26s16
                      end do                                              7d26s16
                     end do                                               7d26s16
                    end if                                              8d18s23
                   end if                                                 3d5s13
                  end if                                                  3d5s13
                 end if                                                    3d4s13
                end if                                                  6d11s21
                ioff=ioff+1
               end do                                                    3d4s13
               i10=1                                                     3d4s13
              end do                                                     3d4s13
             else if(isb.eq.isblk1(1,is).and.isa.eq.isblk1(2,is).and.         2d28s13
     $         isd.eq.isblk1(3,is).and.isc.eq.isblk1(4,is))then         2d28s13
              ididit=2
              ntot=irefo(isd)*nvirtc(isc)                               5d31s18
              nhere=ntot/mynprocg                                        3d5s13
              nleft=ntot-nhere*mynprocg                                  3d5s13
              nhere0=nhere                                               3d5s13
              if(nleft.gt.0)nhere=nhere+1                                3d5s13
              i10=i1s                                                    3d4s13
              i1x=nduse(isd)                                            6d11s21
              ioff=0
              do ii2=i2s,i2e                                             3d5s13
               i2=ii2-nocc(isc)
               if(ii2.eq.i2e)i1x=i1e                                      3d4s13
               do i1=i10,i1x                                             3d4s13
                if(i1.le.nocc(isd))then                                 6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(i2.gt.0.and.i1m.gt.0)then                            5d31s18
                  itry=i1m-1+irefo(isd)*(i2-1)                           6d1s18
                  do iq=0,mynprocg-1                                     5d15s19
                   if(itry.ge.ibc(itmpl+iq).and.itry.le.ibc(itmph+iq))   5d15s19
     $                 then                                             5d15s19
                    go to 2021                                           5d15s19
                   end if                                                5d15s19
                  end do                                                 5d15s19
                  write(6,*)('iq failure ')
                  call dws_sync
                  call dws_finalize
                  stop
 2021             continue
                  ip=iq                                                  5d17s19
                  if(ipass.eq.1)then                                      3d5s13
                   nxn=irefo(isa)*irefo(isb)                             5d31s18
                   isend(ip+1)=isend(ip+1)+nxn                           3d7s13
                   ibc(i1xsnd+is)=ibc(i1xsnd+is)+nxn
                   ibc(isx+ip)=ibc(isx+ip)+nxn
                   if(ifull4x.gt.1)then                                 8d18s23
                    nxn=nvirtc(isa)*nvirtc(isb)                           7d26s16
                    isend(ip+1)=isend(ip+1)+nxn                           7d26s16
                    ibc(i3xsnd+is)=ibc(i3xsnd+is)+nxn
                    ibc(is3+ip)=ibc(is3+ip)+nxn                           7d27s16
                   end if                                               8d18s23
                  else                                                    3d5s13
                   iad=ihcol(ii)+nab*ioff                                 3d5s13
                   if(iflag.eq.0)then                                    5d31s18
                    iadda=0                                              5d31s18
                    iaddb=0                                              5d31s18
                   else                                                  5d31s18
                    iadda=idoub(isa)                                     5d31s18
                    iaddb=idoub(isb)                                     5d31s18
                   end if                                                5d31s18
                   do i3=1,irefo(isa)                                    5d31s18
                    i3p=i3-1+iadda                                       5d31s18
                    do i4=1,irefo(isb)                                   5d31s18
                     i4p=i4-1+iaddb                                      5d31s18
                     iad2=i3p+nbasdwsa*i4p+iad                           5d31s18
                     bc(ipt(ip+1))=bc(iad2)                                3d5s13
                     ipt(ip+1)=ipt(ip+1)+1                               3d5s13
                    end do                                               3d5s13
                   end do                                                3d5s13
                   if(ifull4x.gt.1)then                                 8d18s23
                    do i3=1,nvirtc(isa)                                   7d26s16
                     i3p=i3+nocc(isa)                                     7d26s16
                     do i4=1,nvirtc(isb)                                  7d26s16
                      i4p=i4+nocc(isb)                                    7d26s16
                      iad2=i3p-1+nbasdwsa*(i4p-1)+iad                     7d26s16
                      bc(ipt(ip+1))=bc(iad2)                              7d26s16
                      ipt(ip+1)=ipt(ip+1)+1                               7d26s16
                     end do                                               7d26s16
                    end do                                                7d26s16
                   end if                                               8d18s23
                  end if                                                 3d5s13
                 end if                                                  3d5s13
                end if                                                  6d11s21
                ioff=ioff+1
               end do                                                    3d4s13
               i10=1                                                     3d4s13
              end do                                                     3d4s13
             end if                                                      2d28s13
            end do                                                       2d28s13
            end if                                                      8d21s23
            if(ifull4x.gt.1)then                                        8d18s23
             do is=1,nsdlkk                                               2d28s13
              if(isblkk(1,is).eq.isa.and.isblkk(2,is).eq.isd.and.         2d28s13
     $         isblkk(3,is).eq.isc.and.isblkk(4,is).eq.isb)then         2d28s13
               ntot=nvirtc(isb)*nvirtc(isc)                              8d31s15
               nhere=ntot/mynprocg                                        3d5s13
               nleft=ntot-nhere*mynprocg                                  3d5s13
               nhere0=nhere                                               3d5s13
               if(nleft.gt.0)nhere=nhere+1                                3d5s13
               i10=i1s                                                    3d4s13
               i1x=nduse(isd)                                            6d11s21
               ioff=0
               do ii2=i2s,i2e                                             3d5s13
                i2=ii2-nocc(isc)
                if(ii2.eq.i2e)i1x=i1e                                      3d4s13
                do i1=i10,i1x                                             3d4s13
                 if(i1.le.nocc(isd))then                                 6d11s21
                  if(iflag.eq.0)then                                      5d31s18
                   i1m=i1                                                 5d31s18
                  else                                                    5d31s18
                   i1m=i1-idoub(isd)                                      5d31s18
                  end if                                                  5d31s18
                  if(i2.gt.0.and.i1m.gt.0)then                            5d31s18
                   do ii4=nocc(isb)+1,nbasdwsb                            8d31s15
                    i4=ii4-nocc(isb)                                      3d8s13
                    itry=i2-1+nvirtc(isc)*(i4-1)                          8d31s15
                    do iq=0,mynprocg-1                                     5d15s19
                     if(itry.ge.ibc(jtmpl+iq).and.itry.le.ibc(jtmph+iq))   5d15s19
     $                 then                                             5d15s19
                      go to 2022                                          5d17s19
                     end if                                                5d15s19
                    end do                                                 5d15s19
                    write(6,*)('iq failure ')
                    call dws_sync
                    call dws_finalize
                    stop
 2022               continue
                    ip=iq                                                 5d17s19
                    if(ipass.eq.1)then                                      3d5s13
                     isend(ip+1)=isend(ip+1)+irefo(isa)                   5d31s18
                     ibc(iksnd+is)=ibc(iksnd+is)+irefo(isa)
                     ibc(isk+ip)=ibc(isk+ip)+irefo(isa)                   5d31s18
                    else                                                  3d8s13
                     if(isa.eq.isb)then                                   3d8s13
                      iad=ihcol(ii)+nab*ioff+((ii4*(ii4-1))/2)-1           3d8s13
                     else                                                 3d8s13
                      iad=ihcol(ii)+nab*ioff+nbasdwsa*(ii4-1)-1           8d31s15
                     end if                                               3d8s13
                     if(iflag.eq.0)then                                   5d31s18
                      iadda=0                                             5d31s18
                     else                                                 5d31s18
                      iadda=idoub(isa)                                    5d31s18
                     end if                                               5d31s18
                     do i3=1,irefo(isa)                                   5d31s18
                      i3p=i3+iadda                                        5d31s18
                      bc(ipt(ip+1))=bc(iad+i3p)                           5d31s18
                      ipt(ip+1)=ipt(ip+1)+1                               3d8s13
                     end do
                    end if
                   end do
                  end if
                 end if                                                  6d11s21
                 ioff=ioff+1
                end do
                i10=1
               end do
              else if(isblkk(1,is).eq.isb.and.isblkk(2,is).eq.isd.and.         2d28s13
     $         isblkk(3,is).eq.isc.and.isblkk(4,is).eq.isa)then         2d28s13
               nvc=nvirtc(isc)                                           8d31s15
               ntot=nvirtc(isa)*nvc                                      8d31s15
               nhere=ntot/mynprocg                                        3d5s13
               nleft=ntot-nhere*mynprocg                                  3d5s13
               nhere0=nhere                                               3d5s13
               if(nleft.gt.0)nhere=nhere+1                                3d5s13
               i10=i1s                                                    3d4s13
               i1x=nduse(isd)                                            6d11s21
               ioff=0
               do ii2=i2s,i2e                                             3d5s13
                i2=ii2-nocc(isc)
                if(ii2.eq.i2e)i1x=i1e                                      3d4s13
                do i1=i10,i1x                                             3d4s13
                 if(i1.le.nocc(isd))then                                 6d11s21
                  if(iflag.eq.0)then                                      5d31s18
                   i1m=i1                                                 5d31s18
                  else                                                    5d31s18
                   i1m=i1-idoub(isd)                                      5d31s18
                  end if                                                  5d31s18
                  if(i2.gt.0.and.i1m.gt.0)then                            5d31s18
                   do ii3=nocc(isa)+1,nbasdwsa                            8d31s15
                    i3=ii3-nocc(isa)                                      4d9s14
                    itry=i2-1+nvc*(i3-1)                                  8d31s15
                    do iq=0,mynprocg-1                                     5d15s19
                     if(itry.ge.ibc(ktmpl+iq).and.itry.le.ibc(ktmph+iq))   5d15s19
     $                 then                                             5d15s19
                      go to 2023                                          5d17s19
                     end if                                                5d15s19
                    end do                                                 5d15s19
                    write(6,*)('iq failure '),itry,ibc(ktmpl),
     $                 ibc(ktmph+mynprocg-1)
                    call dws_sync
                    call dws_finalize
                    stop
 2023               continue
                    ip=iq                                                 5d17s19
                    if(ipass.eq.1)then                                      3d5s13
                     isend(ip+1)=isend(ip+1)+irefo(isb)                   5d31s18
                     ibc(iksnd+is)=ibc(iksnd+is)+irefo(isb)
                     ibc(isk+ip)=ibc(isk+ip)+irefo(isb)                   5d31s18
                    else                                                  3d8s13
                     iad=ihcol(ii)+nab*ioff+ii3-1-nbasdwsa                8d31s15
                     if(iflag.eq.0)then                                   5d31s18
                      iaddb=0                                             5d31s18
                     else                                                 5d31s18
                      iaddb=idoub(isb)                                    5d31s18
                     end if                                               5d31s18
                     do i4=1,irefo(isb)                                   5d31s18
                      i4p=i4+iaddb                                        5d31s18
                      bc(ipt(ip+1))=bc(iad+i4p*nbasdwsa)                  5d31s18
                      ipt(ip+1)=ipt(ip+1)+1                               3d8s13
                     end do
                    end if
                   end do                                                 3d8s13
                  end if                                                  3d8s13
                 end if                                                  6d11s21
                 ioff=ioff+1                                             3d8s13
                end do                                                   3d8s13
                i10=1                                                    3d8s13
               end do                                                    3d8s13
              end if                                                     3d8s13
             end do                                                      3d8s13
            end if                                                      8d18s23
           end if                                                        1d10s11
          end if                                                         1d10s11
          ibcoff=itmpl                                                   5d17s19
         end do                                                          1d10s11
         ibcoff=jtmpl                                                    5d17s19
        end do                                                           1d10s11
       end do                                                            1d10s11
       if(ipass.eq.1)then                                               3d5s13
        ipt(1)=ibufs                                                    3d5s13
        do i=2,mynprocg                                                 3d5s13
         im=i-1                                                         3d5s13
         ipt(i)=ipt(im)+isend(im)                                       3d5s13
        end do                                                          3d5s13
        ibcoff=max(ibcoff,ipt(mynprocg)+isend(mynprocg))                9d11s17
       else                                                             3d5s13
        ipt(1)=0                                                        3d6s13
        do i=2,mynprocg                                                 3d5s13
         im=i-1                                                         3d5s13
         ipt(i)=ipt(im)+isend(im)                                       3d5s13
        end do                                                          3d5s13
        ibuf2=ibcoff                                                     3d5s13
        ipr(1)=0                                                        3d5s13
        do i=2,mynprocg                                                  3d5s13
         im=i-1                                                          3d5s13
         ipr(i)=ipr(im)+irecv(im)                                        3d5s13
        end do                                                           3d5s13
        ibcoff=ibuf2+ipr(mynprocg)+irecv(mynprocg)                      3d5s13
        call enough('parajkfromh.  9',bc,ibc)
        do i=0,ipr(mynprocg)+irecv(mynprocg)-1
         bc(ibuf2+i)=sqrt(3.141d0)
        end do
        istot=0
        irtot=0
        do i=1,mynprocg
         istot=istot+isend(i)
         irtot=irtot+irecv(i)
   33    format(13i8)
        end do
        if(nsendt.ne.nrecvt)then
         write(6,*)('send and recieve do not match ')
         call dws_sync
         call dws_finalize
         stop
        end if
        call dws_all2allvb(bc(ibufs),isend,ipt,bc(ibuf2),irecv,ipr)     3d6s13
        do i=1,mynprocg                                                 3d5s13
         ipr(i)=ipr(i)+ibuf2                                            3d5s13
        end do                                                          3d5s13
        isavs=ibcoff                                                    7d7s21
        do is=1,nsdlk
         jmats(is)=ibcoff
         if(isblk(1,is).eq.isblk(2,is))then
          nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2            5d31s18
         else
          nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                    5d31s18
         end if
         if(ifull4x.gt.1)then                                           8d21s23
          ibcoff=jmats(is)+nrow*(ihj(is)+1-ilj(is))                     8d18s23
         end if                                                         8d21s23
         ioooo(is)=ibcoff                                               8d18s23
         ibcoff=ioooo(is)+nrow*(iho(is)+1-ilo(is))                      4d2s13
         call enough('parajkfromh.ioooo',bc,ibc)
        end do
        if(ifull4x.eq.3)then                                            8d18s23
         do is=1,n4x                                                    7d4s21
          i4xb(is)=ibcoff                                               7d4s21
          if(isblk4x(1,is).eq.isblk4x(2,is))then                        7d4s21
           nrowab=(nvirtc(isblk4x(1,is))*(nvirtc(isblk4x(1,is))+1))/2   7d4s21
          else                                                          7d4s21
           nrowab=nvirt(isblk4x(1,is))*nvirt(isblk4x(2,is))             7d4s21
          end if                                                        7d4s21
          itriab=((isblk4x(1,is)*(isblk4x(1,is)-1))/2)+isblk4x(2,is)    7d4s21
          if(isblk4x(3,is).eq.isblk4x(4,is))then                        7d4s21
           nrowcd=(nvirtc(isblk4x(3,is))*(nvirtc(isblk4x(3,is))+1))/2   7d4s21
          else                                                          7d4s21
           nrowcd=nvirtc(isblk4x(3,is))*nvirtc(isblk4x(4,is))           7d4s21
          end if                                                        7d4s21
          itricd=((isblk4x(3,is)*(isblk4x(3,is)-1))/2)+isblk4x(4,is)    7d4s21
          if(itriab.eq.itricd)then                                      7d4s21
           ntot=(nrowab*(nrowab+1))/2                                   7d4s21
          else                                                          7d4s21
           ntot=nrowab*nrowcd                                           7d4s21
          end if                                                        7d4s21
          call ilimts(1,ntot,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  7d4s21
          ibcoff=i4xb(is)+ih+1-il                                       7d4s21
          call enough('parajkfromh. 10',bc,ibc)
          do i=i4xb(is),ibcoff-1                                        7d4s21
           bc(i)=valnz                                                  7d4s21
          end do                                                        7d4s21
         end do                                                         7d4s21
        end if                                                          7d4s21
        do is=1,nsdlk1                                                  3d6s13
         ionex(is)=ibcoff                                               3d6s13
         if(isblk1(1,is).eq.isblk1(2,is))then                           8d4s14
          nrow=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2          5d31s18
          nrow3=(nvirtc(isblk1(1,is))*(nvirtc(isblk1(1,is))+1))/2       7d26s16
         else                                                           8d4s14
          nrow=irefo(isblk1(1,is))*irefo(isblk1(2,is))                  5d31s18
          nrow3=nvirtc(isblk1(1,is))*nvirtc(isblk1(2,is))               7d26s16
         end if                                                         8d4s14
         ncol=ihx(is)+1-ilx(is)                                         7d26s16
         nwds=nrow*ncol                                                 7d26s16
         ibcoff=ionex(is)+nwds                                          8d4s14
         do j=0,nwds-1                                                  8d4s14
          bc(ionex(is)+j)=-2d0
         end do
         i3xi(is)=ibcoff                                                7d26s16
         if(ifull4x.gt.1)then                                           8d18s23
          nwds=nrow3*ncol                                                7d26s16
          ibcoff=i3xi(is)+nwds                                           7d26s16
          do j=0,nwds-1                                                  7d26s16
           bc(i3xi(is)+j)=valnz                                          7d26s16
          end do                                                         7d26s16
         end if                                                         8d18s23
        end do                                                          3d6s13
        if(ifull4x.gt.1)then                                            8d18s23
         do is=1,nsdlkk                                                  3d8s13
          kmats(is)=ibcoff                                               3d8s13
          ibcoff=kmats(is)+irefo(isblkk(1,is))*irefo(isblkk(2,is))       5d31s18
     $        *(ihk(is)+1-ilk(is))                                      3d8s13
          call enough('parajkfromh. 11',bc,ibc)
          do j=1,irefo(isblkk(1,is))*irefo(isblkk(2,is))                 5d31s18
     $        *(ihk(is)+1-ilk(is))                                      3d8s13
           bc(kmats(is)+j-1)=valnz                                       3d8s13
          end do
         end do                                                          3d8s13
        end if                                                          8d18s23
        isave=ibcoff                                                    7d7s21
        call enough('parajkfromh. 12',bc,ibc)
       end if                                                           3d5s13
       do ip=0,mynprocg-1                                               3d5s13
        ipp=ip+1
        do isb=1,nsymb                                                    1d10s11
          nbasdwsb=nbasdws(isb)                                         5d31s18
         do isc=1,nsymb                                                   1d10s11
          isbc=multh(isc,isb)                                             1d10s11
          do isd=1,nsymb                                                  1d10s11
           isa=multh(isd,isbc)                                            1d10s11
           nbasdwsa=nbasdws(isa)                                        5d31s18
           ii=iptoh(isd,isc,isb)                                          1d10s11
           nbc=nduse2(isc)                                              8d18s23
           call ilimts(nduse(isd),nbc,mynprocg,ip,ilh,                  6d11s21
     $           ihh,i1s,i1e,i2s,i2e)                                   3d4s13
           ncd=ihh+1-ilh                                                 3d5s13
           if(ii.gt.0.and.ncd.gt.0)then
            do is=1,nsdlk                                               4d2s13
             jro=iro
             if(isa.eq.isblk(3,is).and.isb.eq.isblk(4,is).and.isc.eq.   4d2s13
     $          isblk(1,is).and.isd.eq.isblk(2,is))then                 4d2s13
              i10=i1s                                                     3d4s13
              i1x=nduse(isd)                                            6d11s21
              if(isa.eq.isb)then
               nrow=(irefo(isc)*(irefo(isc)+1))/2                       5d31s18
               mrow=(nvirtc(isc)*(nvirtc(isc)+1))/2                     6d11s21
              else
               nrow=irefo(isc)*irefo(isd)                               5d31s18
               mrow=nvirtc(isc)*nvirtc(isd)                             6d11s21
              end if
              do i2=i2s,min(i2e,nduse(isc))                             6d11s21
               if(iflag.eq.0)then                                       5d31s18
                i2m=i2                                                  5d31s18
               else                                                     5d31s18
                i2m=i2-idoub(isc)                                       5d31s18
               end if                                                   5d31s18
               if(i2.eq.i2e)i1x=i1e
               do i1=i10,i1x
                if(i1.le.nocc(isd).and.i2.le.nocc(isc))then             6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(isa.eq.isb)then
                  if(i1.ge.i2.and.i2m.gt.0)then                          5d31s18
                   if(ipass.eq.1)then
                    if(ifull4x.gt.1)then                                8d18s23
                     irecv(ipp)=irecv(ipp)+ihj(is)+1-ilj(is)               4d2s13
                     ibc(ijrcv+is)=ibc(ijrcv+is)+ihj(is)+1-ilj(is)
                     ibc(irj+ipp-1)=ibc(irj+ipp-1)+ihj(is)+1-ilj(is)               4d2s13
                    end if                                              8d18s23
                    irecv(ipp)=irecv(ipp)+iho(is)+1-ilo(is)               4d2s13
                    ibc(i4orcv+is)=ibc(i4orcv+is)+iho(is)+1-ilo(is)
                    ibc(jro+ipp-1)=ibc(jro+ipp-1)+iho(is)+1-ilo(is)               4d2s13
                   else                                                   4d2s13
                    if(ifull4x.gt.1)then                                8d18s23
                     i12=((i1m*(i1m-1))/2)+jmats(is)+i2m-1                5d31s18
                     ncol=ihj(is)+1-ilj(is)
                     do i34=1,ncol
                      bc(i12)=bc(ipr(ipp))
                      i12=i12+nrow
                      ipr(ipp)=ipr(ipp)+1
                     end do
                    end if                                              8d18s23
                    i12=((i1m*(i1m-1))/2)+ioooo(is)+i2m-1                5d31s18
                    ncol=iho(is)+1-ilo(is)
                    do i34=1,ncol
                     bc(i12)=bc(ipr(ipp))
                     i12=i12+nrow
                     ipr(ipp)=ipr(ipp)+1
                    end do
                   end if
                  end if
                 else if(min(i1m,i2m).gt.0)then                          5d31s18
                  if(ipass.eq.1)then
                   if(ifull4x.gt.1)then                                 8d18s23
                    irecv(ipp)=irecv(ipp)+ihj(is)+1-ilj(is)               4d2s13
                    ibc(ijrcv+is)=ibc(ijrcv+is)+ihj(is)+1-ilj(is)
                    ibc(irj+ipp-1)=ibc(irj+ipp-1)+ihj(is)+1-ilj(is)               4d2s13
                   end if                                               8d18s23
                   irecv(ipp)=irecv(ipp)+iho(is)+1-ilo(is)               4d2s13
                   ibc(i4orcv+is)=ibc(i4orcv+is)+iho(is)+1-ilo(is)
                   ibc(jro+ipp-1)=ibc(jro+ipp-1)+iho(is)+1-ilo(is)               4d2s13
                  else
                   if(ifull4x.gt.1)then                                 8d18s23
                    i12=i2m+irefo(isc)*(i1m-1)+jmats(is)-1                5d31s18
                    ncol=ihj(is)+1-ilj(is)
                    do i34=1,ncol
                     bc(i12)=bc(ipr(ipp))
                     i12=i12+nrow
                     ipr(ipp)=ipr(ipp)+1
                    end do
                   end if                                               8d18s23
                   i12=i2m+irefo(isc)*(i1m-1)+ioooo(is)-1                5d31s18
                   ncol=iho(is)+1-ilo(is)
                   do i34=1,ncol
                    bc(i12)=bc(ipr(ipp))
                    i12=i12+nrow
                    ipr(ipp)=ipr(ipp)+1
                   end do
                  end if
                 end if
                end if                                                  6d11s21
               end do
               i10=1
              end do
             else if(isb.eq.isblk(3,is).and.isa.eq.isblk(4,is).and.     4d2s13
     $            isc.eq.isblk(1,is).and.isd.eq.isblk(2,is))then        4d2s13
              i10=i1s                                                     3d4s13
              i1x=nduse(isd)                                            6d11s21
              nrow=irefo(isc)*irefo(isd)                                5d31s18
              mrow=nvirtc(isc)*nvirtc(isd)                              6d11s21
              do i2=i2s,min(i2e,nduse(isc))                             6d11s21
               if(iflag.eq.0)then                                       5d31s18
                i2m=i2                                                  5d31s18
               else                                                     5d31s18
                i2m=i2-idoub(isc)                                       5d31s18
               end if                                                   5d31s18
               if(i2.eq.i2e)i1x=i1e
               do i1=i10,i1x
                if(i1.le.nocc(isd).and.i2.le.nocc(isc))then             6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(ipass.eq.1.and.min(i1m,i2m).gt.0)then                5d31s18
                  if(ifull4x.gt.1)then                                  8d18s23
                   irecv(ipp)=irecv(ipp)+ihj(is)+1-ilj(is)                4d2s13
                   ibc(ijrcv+is)=ibc(ijrcv+is)+ihj(is)+1-ilj(is)
                   ibc(irj+ipp-1)=ibc(irj+ipp-1)+ihj(is)+1-ilj(is)               4d2s13
                  end if                                                8d18s23
                  irecv(ipp)=irecv(ipp)+iho(is)+1-ilo(is)                4d2s13
                  ibc(i4orcv+is)=ibc(i4orcv+is)+iho(is)+1-ilo(is)
                  ibc(jro+ipp-1)=ibc(jro+ipp-1)+iho(is)+1-ilo(is)               4d2s13
                 else if(min(i1m,i2m).gt.0)then                          5d31s18
                  if(ifull4x.gt.1)then                                  8d18s23
                   i12=i2m+irefo(isc)*(i1m-1)+jmats(is)-1                 5d31s18
                   ncol=ihj(is)+1-ilj(is)
                   do i34=1,ncol
                    bc(i12)=bc(ipr(ipp))
                    i12=i12+nrow
                    ipr(ipp)=ipr(ipp)+1
                   end do
                  end if                                                8d18s23
                  i12=i2m+irefo(isc)*(i1m-1)+ioooo(is)-1                 5d31s18
                  ncol=iho(is)+1-ilo(is)
                  do i34=1,ncol
                   bc(i12)=bc(ipr(ipp))
                   i12=i12+nrow
                   ipr(ipp)=ipr(ipp)+1
                  end do
                 end if
                end if                                                  6d11s21
               end do
               i10=1
              end do
             end if                                                     4d2s13
            end do                                                      4d2s13
            if(ifull4x.eq.3)then                                        8d18s23
             if(isc.ge.isd)then                                         7d4s21
              itricd=((isc*(isc-1))/2)+isd                              7d4s21
              isx=max(isa,isb)                                           7d4s21
              isn=min(isa,isb)                                           7d4s21
              itriab=((isx*(isx-1))/2)+isn                              7d4s21
              if(itriab.le.itricd)then                                  7d4s21
               if(isc.eq.isd)then                                       7d5s21
                nrowcd=(nvirtc(isc)*(nvirtc(isc)+1))/2                  7d5s21
               else                                                     7d5s21
                nrowcd=nvirtc(isc)*nvirtc(isd)                          7d5s21
               end if                                                   7d5s21
               if(isa.eq.isb)then                                       7d5s21
                nrowab=(nvirtc(isa)*(nvirtc(isa)+1))/2                  7d5s21
                iswab=0                                                 7d5s21
               else                                                     7d5s21
                nrowab=nvirtc(isa)*nvirtc(isb)                          7d5s21
                iswab=1                                                 7d5s21
               end if                                                   7d5s21
               if(itriab.eq.itricd)then                                 7d5s21
                nabcd=(nrowab*(nrowab+1))/2                             7d5s21
                isw=0                                                   7d5s21
               else                                                     7d5s21
                nabcd=nrowab*nrowcd                                     7d5s21
                isw=1                                                   7d5s21
               end if                                                   7d5s21
               call ilimts(1,nabcd,mynprocg,mynowprog,il,ih,i1sz,i1ez,
     $              i2sz,i2ez)                                          7d7s21
c     abcd: i1+na*(i2+nb*(i3+nc*i4)
c     abab: ((ix*(ix-1))/2)+in, ix=i1+na*i2 ge in=i3+na*i4
c     aabb: ((i1*(i1-1))/2)+i2+naa*[((i3*(i3-1))/2)+i4]
c     aaaa: ((ix*(ix-1))/2)+in, ix=((i1*(i1-1))/2)+i2, in=((i3*(i3-1))/2+i4.
               i2eu=inv4x(1,isc,isd,isx)                                7d4s21
               icase=inv4x(2,isc,isd,isx)                                7d4s21
               i10=i1s                                                  7d4s21
               i1x=nduse(isd)                                           7d4s21
               do ii2=i2s,i2e                                           7d4s21
                i2=ii2-nocc(isc)                                        11d8s23
                i2m=i2-1                                                7d6s21
                if(ii2.eq.i2e)i1x=i1e                                   7d4s21
                do ii1=i10,i1x                                          7d4s21
                 i1=ii1-nocc(isd)                                       11d8s23
                 i1m=i1-1                                               7d6s21
                 if(i1.gt.0.and.i2.gt.0.and.
     $                (isd.ne.isc.or.i2.ge.i1))then                     7d4s21
                  if(isc.eq.isd)then                                    7d4s21
                   itri12=((i2*(i2-1))/2)+i1-1                          7d5s21
                  else                                                  7d4s21
                   itri12=i2+nvirtc(isc)*(i1-1)-1                       7d5s21
                  end if                                                7d4s21
                  naabb=nrowab-1                                        7d5s21
                  if(itricd.eq.itriab)naabb=itri12                      7d4s21
                  jtriab=0
                  if(isa.ge.isb)then                                    7d5s21
                   do ib=0,nvirtc(isb)-1                                  7d5s21
c     jtriab+ia+1 le naabb, ia le naabb-1-jtriab
                    imxa=naabb-jtriab                                   7d5s21
                    itopa=min(imxa,ib+iswab*(nvirtc(isa)-1-ib))           7d5s21
                    do ia=0,itopa                                         7d5s21
                     jtriab=jtriab+1                                      7d5s21
                     qq=get4int(i4xb,isc,isd,isa,isb,                   7d6s21
     $                    i2m,i1m,ia,ib,nvirtc,iadqq,bc,ibc)            11d10s22
                     icol=iadqq
                     if(icol.ge.il.and.icol.le.ih)then                    7d5s21
                      if(ipass.eq.1)then                                   7d5s21
                       irecv(ipp)=irecv(ipp)+1                             7d5s21
                       ibc(j4xrcv+i2eu)=ibc(j4xrcv+i2eu)+1
                      else                                                 7d5s21
                       iad=i4xb(i2eu)+icol-il                            7d5s21
                       bc(iad)=bc(ipr(ipp))                              7d5s21
                       ipr(ipp)=ipr(ipp)+1                                 7d5s21
                      end if                                               7d5s21
                     end if                                              7d5s21
                    end do                                                7d5s21
                   end do                                                 7d5s21
                  else                                                  7d5s21
                   do ia=0,nvirtc(isa)-1                                 7d5s21
                    imxb=naabb-jtriab                                   7d5s21
                    ibtop=min(nvirtc(isb)-1,imxb)                       7d5s21
                    do ib=0,ibtop                                       7d5s21
                     ix=max(itri12,jtriab)                                7d5s21
                     in=min(itri12,jtriab)                                7d5s21
                     itri=((ix*(ix+1))/2)+in                              7d5s21
                     irec=jtriab+nrowab*itri12                            7d5s21
                     icol=itri+isw*(irec-itri)+1                          7d5s21
                     qq=get4int(i4xb,isc,isd,isb,isa,                   7d6s21
     $                    i2m,i1m,ib,ia,nvirtc,iadqq,bc,ibc)            11d10s22
                     icol=iadqq
                     jtriab=jtriab+1                                     7d5s21
                     if(icol.ge.il.and.icol.le.ih)then                  7d5s21
                      if(ipass.eq.1)then                                   7d5s21
                       irecv(ipp)=irecv(ipp)+1                             7d5s21
                       ibc(j4xrcv+i2eu)=ibc(j4xrcv+i2eu)+1
                      else                                                 7d5s21
                       iad=i4xb(i2eu)+icol-il                            7d5s21
                       bc(iad)=bc(ipr(ipp))                              7d5s21
                       ipr(ipp)=ipr(ipp)+1                                 7d5s21
                      end if                                               7d5s21
                     end if                                             7d5s21
                    end do                                               7d5s21
                   end do                                                7d5s21
                  end if                                                7d5s21
                 end if                                                 7d4s21
                end do                                                  7d4s21
                i10=1                                                   7d4s21
               end do                                                   7d4s21
              end if                                                    7d4s21
             end if                                                     7d4s21
            end if                                                      7d4s21
            if(ifull4x.eq.1)then                                         8d21s23
             do is=1,nsdlk1                                             8d21s23
              if(isc.eq.isblk1(1,is).and.isd.eq.isblk1(2,is))then       8d21s23
               if(isc.eq.isd)then                                       8d21s23
                iswitch=0                                               8d21s23
                nrow1=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2   8d21s23
               else                                                     8d21s23
                iswitch=1                                               8d21s23
                nrow1=irefo(isblk1(1,is))*irefo(isblk1(2,is))           8d21s23
               end if                                                   8d21s23
               if(isa.eq.isblk1(3,is).or.isa.eq.isblk1(4,is))then       8d21s23
                i10=i1s                                                 8d21s23
                i1n=nduse(isd)                                          8d21s23
                do i2=i2s,i2e                                           8d21s23
                 if(i2.eq.i2e)i1n=i1e                                   8d21s23
                 i2m=i2-1                                               8d21s23
                 if(iflag.ne.0)i2m=i2m-idoub(isc)                       8d21s23
                 do i1=i10,i1n                                          8d21s23
                  i1m=i1-1                                              8d21s23
                  if(iflag.ne.0)i1m=i1m-idoub(isd)                      8d21s23
                  i1top=i2m+iswitch*(nduse(isd)-1-i2m)                  8d21s23
                  if(min(i1m,i2m).ge.0.and.i1m.le.i1top)then            8d21s23
                   call ilimts(irefo(isblk1(3,is)),                     8d21s23
     $                    nvirt(isblk1(4,is)),mynprocg,mynowprog,ll,lh, 8d21s23
     $                 l1s,l1e,l2s,l2e)                                 8d21s23
                   nhere=lh+1-ll                                        8d21s23
                   if(ipass.eq.1)then                                   8d21s23
                    irecv(ipp)=irecv(ipp)+nhere                         8d21s23
                    ibc(i1xrcv+is)=ibc(i1xrcv+is)+nhere                 8d21s23
                    ibc(irx+ipp-1)=ibc(irx+ipp-1)+nhere                 8d21s23
                   else                                                 8d21s23
                    irec=i2m+irefo(isc)*i1m                              8d21s23
                    ix=max(i1m,i2m)                                      8d21s23
                    in=min(i1m,i2m)                                      8d21s23
                    itri=((ix*(ix+1))/2)+in                              8d21s23
                    itri=ionex(is)+itri+iswitch*(irec-itri)             8d21s23
                    do j12=0,nhere-1                                    8d21s23
                     bc(itri)=bc(ipr(ipp))                              8d21s23
                     ipr(ipp)=ipr(ipp)+1                                8d21s23
                     itri=itri+nrow1                                    8d21s23
                    end do                                              8d21s23
                   end if                                               8d21s23
                  end if                                                8d21s23
                 end do                                                 8d21s23
                 i10=1                                                  8d21s23
                end do                                                  8d21s23
               end if                                                   8d21s23
              end if                                                    8d21s23
             end do                                                     8d21s23
            else                                                        8d21s23
            do is=1,nsdlk1                                                2d28s13
             if(isa.eq.isblk1(1,is).and.isb.eq.isblk1(2,is).and.          2d28s13
     $           isd.eq.isblk1(3,is).and.isc.eq.isblk1(4,is))then          2d28s13
              i10=i1s                                                   3d5s13
              i1x=nduse(isd)                                            6d11s21
              do ii2=i2s,i2e                                            3d5s13
               i2=ii2-nocc(isc)                                         3d5s13
               if(ii2.eq.i2e)i1x=i1e                                    3d5s13
               do i1=i10,i1x                                            3d5s13
                if(i1.le.nocc(isd))then                                 6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(i2.gt.0.and.i1m.gt.0)then                            5d31s18
                  itry=i1m+irefo(isd)*(i2-1)                             5d31s18
                  if(itry.ge.ilx(is).and.itry.le.ihx(is))then            3d7s13
                   if(ipass.eq.1)then                                    3d5s13
                    if(isa.eq.isb)then                                   3d5s13
                     nxn=((irefo(isa)*(irefo(isa)+1))/2)                 5d31s18
                     irecv(ipp)=irecv(ipp)+nxn                           3d7s13
                     ibc(i1xrcv+is)=ibc(i1xrcv+is)+nxn
                     ibc(irx+ipp-1)=ibc(irx+ipp-1)+nxn
                     if(ifull4x.gt.1)then                               8d18s23
                      nxn=(nvirtc(isa)*(nvirtc(isa)+1))/2                 7d26s16
                      irecv(ipp)=irecv(ipp)+nxn                           7d26s16
                      ibc(i3xrcv+is)=ibc(i3xrcv+is)+nxn
                      ibc(ir3+ipp-1)=ibc(ir3+ipp-1)+nxn                   7d27s16
                     end if                                             8d18s23
                    else                                                 3d5s13
                     nxn=irefo(isa)*irefo(isb)                           5d31s18
                     irecv(ipp)=irecv(ipp)+nxn                           3d7s13
                     ibc(i1xrcv+is)=ibc(i1xrcv+is)+nxn
                     ibc(irx+ipp-1)=ibc(irx+ipp-1)+nxn
                     if(ifull4x.gt.1)then                               8d18s23
                      nxn=nvirtc(isa)*nvirtc(isb)                         7d26s16
                      irecv(ipp)=irecv(ipp)+nxn                           7d26s16
                      ibc(i3xrcv+is)=ibc(i3xrcv+is)+nxn
                      ibc(ir3+ipp-1)=ibc(ir3+ipp-1)+nxn                   7d27s16
                     end if                                             8d18s23
                    end if                                               3d5s13
                   else                                                  3d5s13
                    if(isa.eq.isb)then                                   3d5s13
                     nrow=(irefo(isa)*(irefo(isa)+1))/2                  5d31s18
                     i34=0                                               8d4s14
                     do i3=1,irefo(isa)                                  5d31s18
                      do i4=1,i3                                         3d5s13
                       iad1=ionex(is)+i34+nrow*(itry-ilx(is))            8d4s14
                       i34=i34+1                                         8d4s14
                       bc(iad1)=bc(ipr(ipp))                             3d5s13
                       ipr(ipp)=ipr(ipp)+1                               3d5s13
                      end do                                             3d5s13
                     end do                                              3d5s13
                     if(ifull4x.gt.1)then                               8d18s23
                      nrow=(nvirtc(isa)*(nvirtc(isa)+1))/2                7d26s16
                      i34=0                                               7d26s16
                      do i3=1,nvirtc(isa)                                 7d26s16
                       do i4=1,i3                                         7d26s16
                        iad1=i3xi(is)+i34+nrow*(itry-ilx(is))             7d26s16
                        i34=i34+1                                         7d26s16
                        bc(iad1)=bc(ipr(ipp))                             7d26s16
                        ipr(ipp)=ipr(ipp)+1                               7d26s16
                       end do                                             7d26s16
                      end do                                              7d26s16
                     end if                                             8d18s23
                    else                                                 3d5s13
                     do i3=1,irefo(isb)                                  5d31s18
                      do i4=1,irefo(isa)                                 5d31s18
                       iad1=ionex(is)+i4-1+irefo(isa)*(i3-1+irefo(isb)   5d31s18
     $                     *(itry-ilx(is)))                             3d5s13
                       bc(iad1)=bc(ipr(ipp))                             3d5s13
                       ipr(ipp)=ipr(ipp)+1                               3d5s13
                      end do                                             3d5s13
                     end do                                              3d5s13
                     if(ifull4x.gt.1)then                               8d18s23
                      do i3=1,nvirtc(isb)                                 7d26s16
                       do i4=1,nvirtc(isa)                                7d26s16
                        iad1=i3xi(is)+i4-1+nvirtc(isa)*(i3-1            8d18s23
     $                       +nvirtc(isb)*(itry-ilx(is)))               8d18s23
                        bc(iad1)=bc(ipr(ipp))                             7d26s16
                        ipr(ipp)=ipr(ipp)+1                               7d26s16
                       end do                                             7d26s16
                      end do                                              7d26s16
                     end if                                             8d18s23
                    end if                                               3d5s13
                   end if                                                3d5s13
                  end if                                                 3d5s13
                 end if                                                  3d5s13
                end if                                                  6d11s21
               end do                                                   3d5s13
               i10=1                                                    3d5s13
              end do
             else if(isb.eq.isblk1(1,is).and.isa.eq.isblk1(2,is).and.   3d8s13
     $        isd.eq.isblk1(3,is).and.isc.eq.isblk1(4,is))then          2d28s13
              i10=i1s                                                   3d5s13
              i1x=nduse(isd)                                            6d11s21
              do ii2=i2s,i2e                                            3d5s13
               i2=ii2-nocc(isc)                                         3d5s13
               if(ii2.eq.i2e)i1x=i1e                                    3d5s13
               do i1=i10,i1x                                            3d5s13
                if(i1.le.nocc(isd))then                                 6d11s21
                 if(iflag.eq.0)then                                      5d31s18
                  i1m=i1                                                 5d31s18
                 else                                                    5d31s18
                  i1m=i1-idoub(isd)                                      5d31s18
                 end if                                                  5d31s18
                 if(i2.gt.0.and.i1m.gt.0)then                            5d31s18
                  itry=i1m+irefo(isd)*(i2-1)                             5d31s18
                  if(itry.ge.ilx(is).and.itry.le.ihx(is))then            3d7s13
                   if(ipass.eq.1)then                                    3d5s13
                    nxn=irefo(isa)*irefo(isb)                            5d31s18
                    irecv(ipp)=irecv(ipp)+nxn                            3d7s13
                    ibc(i1xrcv+is)=ibc(i1xrcv+is)+nxn
                    ibc(irx+ipp-1)=ibc(irx+ipp-1)+nxn
c???
                    if(ifull4x.gt.1)then                                8d18s23
                     nxn=nvirtc(isa)*nvirtc(isb)                          7d26s16
                     irecv(ipp)=irecv(ipp)+nxn                            7d26s16
                     ibc(ir3+ipp-1)=ibc(ir3+ipp-1)+nxn                    7d27s16
                     ibc(i3xrcv+is)=ibc(i3xrcv+is)+nxn
                    end if                                              8d18s23
                   else                                                  3d8s13
                    do i4=1,irefo(isa)                                   5d31s18
                     do i3=1,irefo(isb)                                  5d31s18
                      iad1=ionex(is)+i3-1+irefo(isb)*(i4-1+irefo(isa)    5d31s18
     $                    *(itry-ilx(is)))                              3d8s13
                      bc(iad1)=bc(ipr(ipp))                              3d8s13
                      ipr(ipp)=ipr(ipp)+1                                3d8s13
                     end do                                              3d8s13
                    end do                                               3d8s13
                    if(ifull4x.gt.1)then                                8d18s23
                     do i4=1,nvirtc(isa)                                  7d26s16
                      do i3=1,nvirtc(isb)                                 7d26s16
                       iad1=i3xi(is)+i3-1+nvirtc(isb)*(i4-1+nvirtc(isa)   7d26s16
     $                    *(itry-ilx(is)))                              7d26s16
                       bc(iad1)=bc(ipr(ipp))                              7d26s16
                       ipr(ipp)=ipr(ipp)+1                                3d8s13
                      end do                                              3d8s13
                     end do                                               3d8s13
                    end if                                              8d18s23
                   end if                                                3d8s13
                  end if                                                 3d8s13
                 end if                                                  3d8s13
                end if                                                  6d11s21
               end do                                                   3d8s13
               i10=1                                                    3d8s13
              end do                                                    3d8s13
             end if                                                     3d5s13
            end do                                                      3d5s13
            end if                                                      8d21s23
            if(ifull4x.gt.1)then                                        8d18s23
             do is=1,nsdlkk                                             3d8s13
              if(isblkk(1,is).eq.isa.and.isblkk(2,is).eq.isd.and.         2d28s13
     $        isblkk(3,is).eq.isc.and.isblkk(4,is).eq.isb)then          2d28s13
               i10=i1s                                                  3d5s13
               i1x=nduse(isd)                                           6d11s21
               ioff=0
               do ii2=i2s,i2e                                           3d5s13
                i2=ii2-nocc(isc)                                        3d5s13
                if(ii2.eq.i2e)i1x=i1e                                   3d5s13
                do i1=i10,i1x                                           3d5s13
                 if(i1.le.nocc(isd))then                                6d11s21
                  if(iflag.eq.0)then                                     5d31s18
                   i1m=i1                                                5d31s18
                  else                                                   5d31s18
                   i1m=i1-idoub(isd)                                     5d31s18
                  end if                                                 5d31s18
                  if(i2.gt.0.and.i1m.gt.0)then                           5d31s18
                   do ii4=nocc(isb)+1,nbasdwsb                           8d31s15
                    i4=ii4-nocc(isb)                                     3d8s13
                    nvc=nvirtc(isc)                                      8d31s15
                    itry=i2+nvc*(i4-1)                                   8d31s15
                    if(itry.ge.ilk(is).and.itry.le.ihk(is))then          3d8s13
                     if(ipass.eq.1)then                                  3d8s13
                      irecv(ipp)=irecv(ipp)+irefo(isa)                   5d31s18
                      ibc(ikrcv+is)=ibc(ikrcv+is)+irefo(isa)
                      ibc(irk+ipp-1)=ibc(irk+ipp-1)+irefo(isa)            5d31s18
                     else                                                3d8s13
c
c k_{ad}^{cb}=(ab|dc)
c
                      iad=kmats(is)-1+irefo(isa)*(i1m-1+irefo(isd)       5d31s18
     $                   *(itry-ilk(is)))                               3d8s13
                      do i3=1,irefo(isa)                                 5d31s18
                       bc(iad+i3)=bc(ipr(ipp))                           3d8s13
                       ipr(ipp)=ipr(ipp)+1                               3d8s13
                      end do                                             3d8s13
                     end if                                              3d8s13
                    end if                                               3d8s13
                   end do                                                3d8s13
                  end if                                                 3d8s13
                 end if                                                 6d11s21
                 ioff=ioff+1
                end do                                                  3d8s13
                i10=1                                                   3d8s13
               end do                                                   3d8s13
              else if(isblkk(1,is).eq.isb.and.isblkk(2,is).eq.isd.and.         2d28s13
     $         isblkk(3,is).eq.isc.and.isblkk(4,is).eq.isa)then         2d28s13
               i10=i1s                                                  3d5s13
               i1x=nduse(isd)                                           6d11s21
               ioff=0
               do ii2=i2s,i2e                                           3d5s13
                i2=ii2-nocc(isc)                                        3d5s13
                if(ii2.eq.i2e)i1x=i1e                                   3d5s13
                do i1=i10,i1x                                           3d5s13
                 if(i1.le.nocc(isd))then                                6d11s21
                  if(iflag.eq.0)then                                     5d31s18
                   i1m=i1                                                5d31s18
                  else                                                   5d31s18
                   i1m=i1-idoub(isd)                                     5d31s18
                  end if                                                 5d31s18
                  if(i2.gt.0.and.i1m.gt.0)then                           5d31s18
                   do ii3=nocc(isa)+1,nbasdwsa                           8d31s15
                    i3=ii3-nocc(isa)                                     3d8s13
                    nvc=nvirtc(isc)                                      8d31s15
                    itry=i2+nvc*(i3-1)                                   8d31s15
                    if(itry.ge.ilk(is).and.itry.le.ihk(is))then          3d8s13
                     if( ipass.eq.1)then                                  3d8s13
                      irecv(ipp)=irecv(ipp)+irefo(isb)                   5d31s18
                      ibc(ikrcv+is)=ibc(ikrcv+is)+irefo(isb)
                      ibc(irk+ipp-1)=ibc(irk+ipp-1)+irefo(isb)            5d31s18
                     else                                                3d8s13
                      iad=kmats(is)-1+irefo(isb)*(i1m-1+irefo(isd)       5d31s18
     $                     *(itry-ilk(is)))                              3d8s13
                      do i4=1,irefo(isb)                                 5d31s18
                       bc(iad+i4)=bc(ipr(ipp))                           4d1s13
                       ipr(ipp)=ipr(ipp)+1                               3d8s13
                      end do                                             3d8s13
                     end if                                              3d8s13
                    end if                                               3d8s13
                   end do                                                3d8s13
                  end if                                                 3d8s13
                 end if                                                 6d11s21
                 ioff=ioff+1
                end do                                                  3d8s13
                i10=1                                                   3d8s13
               end do                                                   3d8s13
              end if                                                    3d8s13
             end do                                                     3d8s13
            end if                                                      8d18s23
           end if                                                       3d5s13
          end do                                                        3d5s13
         end do                                                         3d5s13
        end do                                                          3d5s13
       end do                                                           3d5s13
      end do                                                            3d5s13
      ishfy=ibcb4-isavs
      nmve=isave-isavs                                                  7d7s21
      do i=0,nmve-1                                                     7d7s21
       bc(ibcb4+i)=bc(isavs+i)                                          7d7s21
      end do                                                            7d7s21
      ibcoff=isave+ishfy                                                7d11s22
      do is=1,nsdlk                                                     7d7s21
       ioooo(is)=ioooo(is)+ishfy                                        7d7s21
       if(ifull4x.gt.1)then                                             8d18s23
        jmats(is)=jmats(is)+ishfy                                        7d7s21
       end if                                                           8d18s23
      end do                                                            7d7s21
      if(ifull4x.gt.1)then                                              8d18s23
       do is=1,nsdlkk                                                    7d7s21
        kmats(is)=kmats(is)+ishfy                                        7d7s21
       end do                                                            7d7s21
      end if                                                            8d18s23
      do is=1,nsdlk1                                                    7d7s21
       ionex(is)=ionex(is)+ishfy                                        7d7s21
       if(ifull4x.gt.1)then                                             8d18s23
        i3xi(is)=i3xi(is)+ishfy                                          7d7s21
       end if                                                           8d18s23
      end do                                                            7d7s21
      if(ifull4x.eq.3)then                                              8d18s23
       do is=1,n4x                                                      7d7s21
        i4xb(is)=i4xb(is)+ishfy                                         7d7s21
       end do                                                           7d7s21
      end if                                                            7d7s21
      if(idwsdeb.gt.10.or.iprtr(10).ne.0)then                           1d6s20
       write(6,*)('for 4o')
   12  format('integral type ',4i2,5x,i5)                               3d4s13
       do is=1,nsdlk                                                     12d13s11
        nhere=iho(is)+1-ilo(is)                                          1d6s11
        if(isblk(1,is).eq.isblk(2,is))then                               12d13s11
         nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              5d31s18
        else                                                             12d13s11
         nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      5d31s18
        end if                                                           12d13s11
        if(nrow*nhere.gt.0)then
         write(6,12)(isblk(j,is),j=1,4),is                                   12d13s11
         call prntm2(bc(ioooo(is)),nrow,nhere,nrow)                       1d5s12
        end if
       end do                                                            12d13s11
       if(ifull4x.gt.1)then                                             8d18s23
        write(6,*)('for J'),nsdlk,bc(2110993)
        do is=1,nsdlk                                                     12d13s11
         nhere=ihj(is)+1-ilj(is)                                          12d13s11
         if(isblk(1,is).eq.isblk(2,is))then                               12d13s11
          nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              5d31s18
         else                                                             12d13s11
          nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      5d31s18
         end if                                                           12d13s11
         if(nrow.gt.0.and.nhere.gt.0)then                                 2d28s13
          write(6,12)(isblk(j,is),j=1,4)                                   12d13s11
          if(jmats(is).gt.0)then
           write(6,*)('no. '),is,jmats(is)
           call prntm2(bc(jmats(is)),nrow,nhere,nrow)                      12d13s11
          end if
         end if
        end do                                                            12d13s11
       end if                                                           8d18s23
       if(ifull4x.eq.3)then                                             8d18s23
        write(6,*)('for 4xb ')                                           7d5s21
        do is=1,n4x                                                     7d5s21
         if(isblk4x(1,is).eq.isblk4x(2,is))then                         7d4s21
          nrowab=(nvirtc(isblk4x(1,is))*(nvirtc(isblk4x(1,is))+1))/2    7d4s21
         else                                                           7d4s21
          nrowab=nvirt(isblk4x(1,is))*nvirt(isblk4x(2,is))              7d4s21
         end if                                                         7d4s21
         itriab=((isblk4x(1,is)*(isblk4x(1,is)-1))/2)+isblk4x(2,is)     7d4s21
         if(isblk4x(3,is).eq.isblk4x(4,is))then                         7d4s21
          nrowcd=(nvirtc(isblk4x(3,is))*(nvirtc(isblk4x(3,is))+1))/2    7d4s21
         else                                                           7d4s21
          nrowcd=nvirtc(isblk4x(3,is))*nvirtc(isblk4x(4,is))            7d4s21
         end if                                                         7d4s21
         itricd=((isblk4x(3,is)*(isblk4x(3,is)-1))/2)+isblk4x(4,is)     7d4s21
         if(itriab.eq.itricd)then                                       7d4s21
          ntot=(nrowab*(nrowab+1))/2                                    7d4s21
         else                                                           7d4s21
          ntot=nrowab*nrowcd                                            7d4s21
         end if                                                         7d4s21
         call ilimts(1,ntot,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   7d4s21
         nhere=ih+1-il                                                  7d4s21
          write(6,12)(isblk4x(j,is),j=1,4)                               7d5s21
          call prntm2(bc(i4xb(is)),1,nhere,1)                            7d4s21
        end do                                                          7d4s21
       end if                                                           6d11s21
       write(6,*)('for onex'),ionex(1)
       do is=1,nsdlk1                                                    3d9s12
        nhere=ihx(is)+1-ilx(is)                                          12d13s11
        if(isblk1(1,is).eq.isblk1(2,is))then                             8d4s14
         nrow=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2            5d31s18
        else                                                             8d4s14
         nrow=irefo(isblk1(1,is))*irefo(isblk1(2,is))                    5d31s18
        end if                                                           8d4s14
        if(ionex(is).gt.0.and.nhere.gt.0.and.nrow.gt.0)then              2d28s13
         write(6,*)('entry no. '),is
         write(6,12)(isblk1(j,is),j=1,4)                                  3d9s12
         call prntm2(bc(ionex(is)),nrow,nhere,nrow)                       1d5s12
        end if                                                           2d28s13
       end do                                                            12d13s11
       if(ifull4x.gt.1)then                                             8d18s23
        write(6,*)('for 3x'),i3xi(1)                                      7d26s16
        do is=1,nsdlk1                                                    7d26s16
         nhere=ihx(is)+1-ilx(is)                                          7d26s16
         if(isblk1(1,is).eq.isblk1(2,is))then                             7d26s16
          nrow=(nvirtc(isblk1(1,is))*(nvirtc(isblk1(1,is))+1))/2          7d26s16
         else                                                             7d26s16
          nrow=nvirtc(isblk1(1,is))*nvirtc(isblk1(2,is))                  7d26s16
         end if                                                           7d26s16
         if(i3xi(is).gt.0.and.nhere.gt.0.and.nrow.gt.0)then               7d26s16
          write(6,*)is,i3xi(is)
          write(6,12)(isblk1(j,is),j=1,4)                                 7d26s16
          call prntm2(bc(i3xi(is)),nrow,nhere,nrow)                       7d26s16
         end if                                                           7d26s16
        end do                                                            7d26s16
        write(6,*)('for K')
        do is=1,nsdlkk                                                     12d13s11
         nhere=ihk(is)+1-ilk(is)                                          12d13s11
         nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))                     5d31s18
         if(kmats(is).gt.0.and.nhere.gt.0.and.nrow.gt.0)then              2d28s13
          write(6,12)(isblkk(j,is),j=1,4)                                   12d13s11
          call prntm2(bc(kmats(is)),nrow,nhere,nrow)                      12d13s11
         end if                                                           2d28s13
        end do                                                            12d13s11
       end if
      end if                                                            8d18s23
      return
      end
