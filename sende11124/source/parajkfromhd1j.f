c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parajkfromhd1j(ipair,ngaus,ibdat,nhcolt,ihmat,nocc,icol9d20s16
     $     ,iorbx,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,     1d5s12
     $     ioooo,ionex,noc,iter,idwsdeb,ncomp,nvirtc,iprop,itrans,      3d21s16
     $     i4od,ionexd,kmatd,jmatd,i4od2b,ionexd2,der2e,nbasdwsc,        8d3s16
     $     isblkder,isblkxder,nsblkder,nsblkxder,isblkxder1,nsblkxder1, 9d16s16
     $     isblkder1,nsblkder1,isblkkder,nsblkkder,nbasisp,idorel,bc,   11d10s22
     $     ibc)                                                         11d10s22
      implicit real*8 (a-h,o-z)
c
c     like jkfromhd1, working with (ab|dc') instead of (a'b|dc)
c     also, both d and c are occ only.
c
      external second
      integer*8 isstor(1),ibstor(1)                                     5d13s10
      parameter (idp=512)                                                11d30s12
      integer isend(idp),irecv(idp),ipt(idp),ipr(idp),nbasdwsc(8)       3d24s16
      logical ltest                                                     3d29s13
      include "common.store"
      include "common.hf"
      dimension nocc(1),iptoh(8,8,8),iapair(3,1),jmats(1),kmats(1)      1d4s12
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension multh(8,8),iorbn(8),itrans(1),iorbx(8),nbasisp(*),      4d7s22
     $     ioooo(1),isblkder(4,idbk),isblkxder(4,idbk),                 8d24s16
     $     ionex(idbk),noc(8),myh(idbk),nvirtc(8),isblkkder(4,idbk),    9d16s16
     $     i4od(idbk),ionexd(idbk),kmatd(idbk),jmatd(idbk),i4od2b(idbk), 7d21s16
     $     ionexd2(idbk),isblkxder1(4,idbk),isblkder1(4,idbk)           9d16s16
      call second(time1)                                                11d27s12
c
c     in this version, do transpose of hmat, then transform full hmat
c     to mo basis, then store appropriate parts of the result.
c
c
      ntmpmat=2                                                         10d28s20
      if(idorel.eq.0.or.idorel.eq.1)then                                10d28s20
       ntmpmat=1                                                        10d28s20
      end if                                                            10d30s20
      ibcoffo=ibcoff                                                    6d22s10
      iix=0
      do isb=1,nsymb
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
      nn=(ngaus*(ngaus+1))/2                                            9d26s16
      jbdat=ibdat-1
      if(mynprocg.gt.idp)then                                           11d30s12
       write(6,*)('in parajkfromhd, mynprocg = '),mynprocg
       write(6,*)('greater than idp = '),idp
       call dws_sync
       call dws_finalize
       stop
      end if
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
         isbc=multh(isc,isb)                                            3d2s16
         do isd=1,nsymb                                                 11d30s12
          isad=multh(isd,isbc)                                           11d30s12
          isa=multh(isad,iprop)                                         9d20s16
          ii=iptoh(isd,isc,isb)                                         11d30s12
          ncd=nocc(isc)*nocc(isd)                                       9d20s16
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
        msend=isend(1)                                                  4d29s16
        mrev=irecv(1)                                                   4d29s16
        do ip=2,mynprocg                                                11d30s12
         ipm=ip-1                                                       11d30s12
         ipt(ip)=ipt(ipm)+isend(ipm)                                    11d30s12
         ipr(ip)=ipr(ipm)+irecv(ipm)                                    11d30s12
         msend=msend+isend(ip)                                          4d29s16
         mrev=mrev+irecv(ip)                                            4d29s16
        end do                                                          11d30s12
        rms=0d0
        do i=0,msend-1
         rms=rms+bc(ibufs+i)**2
        end do
        if(rms.gt.1d10)then
         call dws_sync
         call dws_finalize
         stop
        end if
        call dws_all2allvb(bc(ibufs),isend,ipt,bc(ihmat),irecv,ipr)     11d30s12
        rms=0d0
        do i=0,mrev-1
         rms=rms+bc(ihmat+i)**2
        end do
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
        do isb=1,nsymb                                                  2d11s13
         nbasdwsb=nbasisp(isb)                                          5d31s22
         do isc=1,nsymb                                                 9d27s16
          isbc=multh(isc,isb)                                           9d27s16
          do isd=1,nsymb                                                2d11s13
           isad=multh(isd,isbc)                                         9d27s16
           isa=multh(isad,iprop)                                        9d27s16
           nbasdwsa=nbasisp(isa)                                        4d7s22
           ntotrun=0
           ncd=nocc(isc)*nocc(isd)                                      9d27s16
           iib=iptoh(isd,isc,isb)                                       9d27s16
           if(iib.gt.0.and.ncd.gt.0)then                                3d1s16
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
              issa=isstor(iaa)                                          3d1s16
              iaaa=ibstor(iaa)                                              11d30s12
              do ib=1,nba0                                              8d31s15
               ibb=ib+ibc(jbdat+ibc(jpair+1)+ngaus3)                          6d4s10
               issb=isstor(ibb)                                                5d13s10
               ibbb=ibstor(ibb)                                               11d30s12
 5653          format(4i8)
               if(ibb.le.iaa.and.isa.le.isb)then                        9d27s16
                if(issb.eq.isb.and.issa.eq.isa)then                      9d20s16
                 if(ipass.eq.1)then                                      2d11s13
                 irecv(ip+1)=irecv(ip+1)+nhere*ntmpmat                  5d31s22
                  myh(iib)=myh(iib)+nhere*ntmpmat                       5d31s22
                  nsumx=nsumx+1
                 else                                                    2d11s13
                  if(issb.eq.issa)then
                   in=min(ibbb,iaaa)
                   ix=max(ibbb,iaaa)
                   ix=((ix*(ix-1))/2)+in-1
                   nxn=(nbasdwsb*(nbasdwsb+1))/2
                  else
                   ix=iaaa-1+nbasdwsa*(ibbb-1)                             3d1s16
                   nxn=nbasdwsb*nbasdwsa                                      3d1s16
                  end if
                  do in=0,ntmpmat-1                                     10d30s20
                   do j=0,nhere-1                                         3d2s16
                    iad=myh(iib)+ix+nxn*(j+nhere*in)                    5d31s22
                    bc(iad)=bc(ipr(ip+1))                                 2d11s13
                    ntotrun=ntotrun+1
                    ipr(ip+1)=ipr(ip+1)+1                                 2d11s13
                   end do                                                 2d11s13
                  end do                                                5d31s22
                 end if
                else if(issa.eq.isb.and.issb.eq.isa)then                9d20s16
                 if(ipass.eq.1)then
                  irecv(ip+1)=irecv(ip+1)+nhere*ntmpmat                 5d31s22
                  myh(iib)=myh(iib)+nhere*ntmpmat                       5d31s22
                  nsumx=nsumx+1
                 else                                                   9d20s16
                  ix=ibbb-1+nbasdwsa*(iaaa-1)                           9d20s16
                  nxn=nbasdwsb*nbasdwsa                                      3d1s16
                  do in=0,ntmpmat-1                                     5d31s22
                   do j=0,nhere-1
                    iad=myh(iib)+ix+nxn*(j+nhere*in)                    5d31s22
                    bc(iad)=bc(ipr(ip+1))                                9d20s16
                    ntotrun=ntotrun+1
                    ipr(ip+1)=ipr(ip+1)+1                                9d20s16
                   end do
                  end do                                                5d31s22
                 end if                                                 9d20s16
                end if                                                  2d11s13
               end if                                                   3d1s16
              end do
              end do            !ia
             end do             !i
            end if                                                      9d2s15
           end do      !isd                                                 2d11s13
          end do       !isc                                                  2d11s13
         end do        !isb                                                  2d11s13
       end do                   !ip                                                  2d15s13
       if(ipass.eq.1)then                                               11d28s12
        nsend=0                                                         11d28s12
        nrecv=0                                                         11d28s12
        do ip=1,mynprocg                                                11d30s12
         nsend=nsend+isend(ip)                                          11d30s12
         nrecv=nrecv+irecv(ip)                                          11d30s12
        end do                                                          11d28s12
        ngot=ibcoff+1-ihmat                                              11d30s12
        if(iter.eq.1.and.mynowprog.eq.0)then
         write(6,*)('no. of words to receive:  '),nrecv                   11d30s12
         write(6,*)('no. of words under ihmat: '),ngot,ibcoff                   11d30s12
        end if
        if(nrecv.gt.ngot)then                                           11d30s12
         more=nrecv-ngot                                                11d30s12
         if(iter.eq.1.and.mynowprog.eq.0)then
         write(6,*)('need to make some space for the receive ')         11d30s12
         write(6,*)('need '),more,(' additional words ')                11d30s12
         end if
         ibcoff=ihmat+nrecv                                             12d3s12
         call enough('parajkfromhd1j.  1',bc,ibc)
        end if                                                          11d30s12
        ibufs=ibcoff                                                    11d28s12
        ibcoff=ibufs+max(nsend,nrecv)                                   11d30s12
        call enough('parajkfromhd1j.  2',bc,ibc)
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
       if(iter.eq.1.and.ipass.eq.1.and.mynowprog.eq.0)then               3d3s17
        write(6,*)('send/recv: ')                                       11d28s12
        do ip=0,mynprocg-1                                              11d28s12
         write(6,335)ip,isend(1+ip),irecv(1+ip),ipt(1+ip)               11d30s12
  335    format(i5,4i8)                                                 11d28s12
        end do                                                          11d28s12
        write(6,*)('totals '),nsend,nrecv
        if(nsend.ne.nhcolt)then
         write(6,*)('send doesn''t match nhcolt: '),nsend,nhcolt
         call dws_sync
         call dws_finalize
         stop
        end if
        write(6,*)('nsendt '),nsendt,(' vs. nrecvt '),nrecvt
        if(nsendt.ne.nrecvt)then                                        8d11s14
         call dws_sync                                                  8d11s14
         call dws_finalize                                              8d11s14
         stop                                                           8d11s14
        end if                                                          8d11s14
       end if                                                            3d15s12
      end do                                                            11d28s12
c
c     ihmat - raw ints: nocc(c)*nocc(d);nn
c     ibufs - ordered ints to send: nocc(c)*nocc(d);nn
c     ihmat - ints that have been recieved: nocc(c)*nocc(d);nn
c     ibufs - re-ordered ints: nocc(c)*nocc(d);nn
c     ihmat - transformed ints
c
      call second(time2)
      ipdeb=10
      ihlast=ihmat                                                      3d5s13
      do isb=1,nsymb                                                    11d30s12
       nbasdwsb=nbasisp(isb)                                            5d31s22
       do isc=1,nsymb                                                   11d30s12
        isbc=multh(isc,isb)                                             3d4s16
        do isd=1,nsymb                                                  11d30s12
         isad=multh(isd,isbc)                                           3d4s16
         isa=multh(isad,iprop)                                          3d4s16
         nbasdwsa=nbasisp(isa)                                          5d31s22
         nxn=nbasdwsa*nbasdwsb                                          3d1s16
         nxnc=nbasdws(isa)*nbasdws(isb)                                 4d12s22
         iis=iptoh(isd,isc,isb)                                         12d21s12
         ncd=nocc(isc)*nocc(isd)                                        9d20s16
         if(ncd.gt.0.and.iis.gt.0)then                                  12d21s12
          nhere=ncd/mynprocg                                            11d30s12
          nleft=ncd-nhere*mynprocg                                      11d30s12
          if(mynowprog.lt.nleft)nhere=nhere+1                           11d30s12
          if(nhere.gt.0)then                                            11d30s12
           if(isa.eq.isb)then                                           9d20s16
            itmp=ibcoff                                                 9d20s16
            itmp2=itmp+nbasdwsa*nbasdwsb                                9d20s16
            itmp3=itmp2+nbasdwsa*nbasdwsb                               5d31s22
            ibcoff=itmp3+nbasdwsa*nbasdwsb                              9d20s16
            call enough('parajkfromhd1j.  3',bc,ibc)
            nxn=(nbasdwsa*(nbasdwsa+1))/2                               9d20s16
            nxnc=(nbasdws(isa)*(nbasdws(isa)+1))/2                      4d12s22
            ihcol(iis)=ihlast                                            3d5s13
            ihlast=ihlast+nxnc*nhere                                    4d12s22
            if(ihlast.gt.ibufs)then                                      3d5s13
             write(6,*)('ihlast now exceeds ibufs!!! '),ihlast,ibufs     3d5s13
             write(6,*)('no. of words '),ihlast+1-ihmat
             write(6,*)('nxn '),nxn,nxnc
             write(6,*)('nhere '),nhere
             write(6,*)('nbasdwsa '),nbasdwsa
             write(6,*)('nbasdwsb '),nbasdwsb
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
              call dgemm('n','n',nbasdwsa,nbasdws(isa),nbasdwsa,        2d21s23
     $            1d0,bc(itmp),nbasdwsa,bc(jorb),nbasdwsa*ncomp,        10d30s20
     $            0d0,bc(itmp2),nbasdwsa,                               5d29s18
     d' parajkfromh.  1')
              do i=1,nbasdws(isa)                                       2d21s23
               do j=1,nbasdwsa                                           8d31s15
                iad1=itmp2+j-1+nbasdwsa*(i-1)                            8d31s15
                iad2=itmp+i-1+nbasdws(isa)*(j-1)                        2d21s23
                bc(iad2)=bc(iad1)                                        11d30s12
               end do                                                    11d30s12
              end do                                                     11d30s12
              call dgemm('n','n',nbasdws(isa),nbasdws(isa),nbasdwsa,    2d21s23
     $            1d0,bc(itmp),nbasdws(isa),bc(jorb),nbasdwsa*ncomp,    2d21s23
     $            factr,bc(itmp3),nbasdws(isa),                         2d21s23
     d' parajkfromh.  2')
              factr=1d0                                                 10d30s20
              if(ipdeb.le.10)then
               write(6,*)('to yield '),iis,ihcol(iis),nxn
               call prntm2(bc(itmp3),nbasdws(isa),nbasdws(isa),         2d21s23
     $            nbasdws(isa))                                         2d21s23
              end if
             end do                                                     10d30s20
             do i=1,nbasdws(isa)                                        2d21s23
              do j=1,i                                                  11d30s12
               iad0=ihcol(iis)+((i*(i-1))/2)+j-1+nxnc*(icol-1)          5d30s18
               iad1=jtmp3+j+nbasdws(isa)*(i-1)                          2d21s23
               bc(iad0)=bc(iad1)                                        11d30s12
              end do                                                    11d30s12
             end do                                                     11d30s12
            end do                                                      11d30s12
            ibcoff=itmp                                                 9d20s16
           else                                                         9d20s16
            itmp=ibcoff                                                  11d30s12
            itmp2=itmp+nxn                                               3d2s16
            itmp3=itmp2+nxn                                             3d2s16
            ibcoff=itmp3+nxn                                             3d2s16
            call enough('parajkfromhd1j.  4',bc,ibc)
            nxn=nbasdwsa*nbasdwsb                                       8d31s15
            nxnc=nbasdws(isa)*nbasdws(isb)                              4d12s22
            ihcol(iis)=ihlast                                            3d5s13
            ihlast=ihlast+nxnc*nhere                                    4d12s22
            if(ihlast.gt.ibufs)then                                      3d5s13
             write(6,*)('ihlast now exceeds ibufs!!! '),ihlast,ibufs     3d5s13
             call dws_sync                                               3d5s13
             call dws_finalize                                           3d5s13
             stop                                                        3d5s13
            end if                                                       3d5s13
            do icol=1,nhere                                             11d30s12
             iads=ihcol(iis)+nxnc*(icol-1)                              5d30s18
             factr=0d0                                                  10d30s20
             jorbb=iorbx(isb)                                           10d30s20
             jorba=iorbx(isa)                                           10d30s20
             do in=0,ntmpmat-1                                          10d30s20
              iad=myh(iis)+nxn*(icol-1+nhere*in)                        10d30s20
              call dgemm('n','n',nbasdwsa,nbasdwsc(isb),nbasdwsb,       2d21s23
     $            1d0,bc(iad),nbasdwsa,bc(jorbb),nbasdwsb*ncomp,        10d30s20
     $            0d0,bc(itmp),nbasdwsa,                                8d31s15
     d' parajkfromh.  3')
              do i=1,nbasdwsa                                            8d31s15
               do j=1,nbasdwsc(isb)                                     2d21s23
                iad1=iad+j-1+nbasdwsc(isb)*(i-1)                        2d21s23
                iad2=itmp+i-1+nbasdwsa*(j-1)                             8d31s15
                bc(iad1)=bc(iad2)                                        11d30s12
               end do                                                    11d30s12
              end do                                                     11d30s12
              call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isa),nbasdwsa,  2d21s23
     $            1d0,bc(iad),nbasdwsc(isb),bc(jorba),nbasdwsa*ncomp,   2d21s23
     $            factr,bc(itmp3),nbasdwsc(isb),                        2d21s23
     d' parajkfromh.  4')
              factr=1d0                                                 10d30s20
              jorbb=jorbb+nbasdwsb                                      10d30s20
              jorba=jorba+nbasdwsa                                      10d30s20
             end do                                                      10d30s20
             do i=1,nbasdwsc(isb)                                       2d21s23
              do j=1,nbasdwsc(isa)                                      2d21s23
               iad1=iads+j-1+nbasdwsc(isa)*(i-1)                        2d21s23
               iad2=itmp3+i-1+nbasdwsc(isb)*(j-1)                       2d21s23
               bc(iad1)=bc(iad2)                                        11d30s12
              end do                                                    11d30s12
             end do                                                     11d30s12
            end do                                                      11d30s12
            ibcoff=itmp                                                  11d30s12
           end if                                                       9d20s16
          end if                                                        11d30s12
         end if                                                         11d30s12
        end do                                                          11d30s12
       end do                                                           11d30s12
      end do                                                            11d30s12
c
c     for j need
c     (o'o|vv), (oo'|vv), (oo|v'v) and (oo|vv')
c     recall we have  (vv|oo')
c                      a b dc  stored under iptoh(d,c,b),
c     thus presently we can get the first two
c     for mixed onex we need
c     (o'a|ov), (o'o|av), (o'o|oa)
c     (ao'|ov), (oo'|av), (oo'|ov)
c     (ao|o'v), (oa|o'v), (oo|o'a)
c     (ao|ov'), (oa|ov'), (oo|av')
c     we can now get ('oo|av) and (oo'|av).
c                                                                       9d16s16
      do isb=1,nsblkkder                                                9d16s16
       nrow=nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb))               8d22s16
       ncol=nvirtc(isblkkder(3,isb))*nvirtc(isblkkder(4,isb))           9d20s16
       nall=nrow*ncol                                                   8d4s16
       call ilimts(nvirtc(isblkkder(3,isb)),nvirtc(isblkkder(4,isb)),   8d22s16
     $      mynprocg,mynowprog,il,ih,i1s,i1e,is,i2e)                    8d22s16
       jmatda=ibcoff
       ibcoff=jmatda+nall
       call enough('parajkfromhd1j.  5',bc,ibc)
       do i=0,nall-1                                                    8d4s16
        bc(jmatda+i)=0d0                                                8d4s16
       end do                                                           8d4s16
c
c     (o'o|vv) part ok
c
       itry=iptoh(isblkkder(2,isb),isblkkder(1,isb),isblkkder(3,isb))
       itry2=iptoh(isblkkder(2,isb),isblkkder(1,isb),isblkkder(4,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(2,isb)),nocc(isblkkder(1,isb)),      9d20s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(2,isb))
        ii0=ihcol(itry)
        if(isblkkder(3,isb).eq.isblkkder(4,isb))then
         nrowh=(nbasdwsc(isblkkder(3,isb))                              9d20s16
     $        *(nbasdwsc(isblkkder(3,isb))+1))/2                        9d20s16
        else                                                            9d20s16
         nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        end if                                                          9d20s16
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(isblkkder(4,isb).eq.isblkkder(3,isb))then                  9d20s16
           do ia=0,nvirtc(isblkkder(4,isb))-1
            iap=ia+nocc(isblkkder(4,isb))
            do ib=0,nvirtc(isblkkder(3,isb))-1
             ibp=ib+nocc(isblkkder(3,isb))
             ix=max(iap,ibp)+1                                          9d20s16
             in=min(iap,ibp)                                            9d20s16
             ix=((ix*(ix-1))/2)+in
             iadh=ii0+ix                                                9d20s16
             iadj=jmatda+icm+nocc(isblkkder(1,isb))*(idm
     $         +nocc(isblkkder(2,isb))*(ib+nvirtc(isblkkder(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          else                                                          9d20s16
           do ia=0,nvirtc(isblkkder(4,isb))-1
            iap=ia+nocc(isblkkder(4,isb))
            do ib=0,nvirtc(isblkkder(3,isb))-1
             ibp=ib+nocc(isblkkder(3,isb))
             iadh=ii0+iap+nbasdwsc(isblkkder(4,isb))*ibp
             iadj=jmatda+icm+nocc(isblkkder(1,isb))*(idm
     $         +nocc(isblkkder(2,isb))*(ib+nvirtc(isblkkder(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d20s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then
        call ilimts(nocc(isblkkder(2,isb)),nocc(isblkkder(1,isb)),      9d20s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ia=0,nvirtc(isblkkder(4,isb))-1
           iap=ia+nocc(isblkkder(4,isb))
           do ib=0,nvirtc(isblkkder(3,isb))-1
            ibp=ib+nocc(isblkkder(3,isb))
            iadh=ii0+ibp+nbasdwsc(isblkkder(3,isb))*iap                 9d27s16
            iadj=jmatda+icm+nocc(isblkkder(1,isb))*(idm
     $         +nocc(isblkkder(2,isb))*(ib+nvirtc(isblkkder(3,isb))*ia))
            bc(iadj)=bc(iadj)+bc(iadh)                                  10d31s16
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|vv) ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|vv) part ok
c
       itry=iptoh(isblkkder(1,isb),isblkkder(2,isb),isblkkder(3,isb))
       itry2=iptoh(isblkkder(1,isb),isblkkder(2,isb),isblkkder(4,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(1,isb)),nocc(isblkkder(2,isb)),      9d20s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(1,isb))
        ii0=ihcol(itry)
        if(isblkkder(3,isb).eq.isblkkder(4,isb))then
         nrowh=(nbasdwsc(isblkkder(3,isb))                              9d20s16
     $        *(nbasdwsc(isblkkder(3,isb))+1))/2                        9d20s16
        else                                                            9d20s16
         nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        end if                                                          9d20s16
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(isblkkder(4,isb).eq.isblkkder(3,isb))then                  9d20s16
           do ia=0,nvirtc(isblkkder(4,isb))-1
            iap=ia+nocc(isblkkder(4,isb))
            do ib=0,nvirtc(isblkkder(3,isb))-1
             ibp=ib+nocc(isblkkder(3,isb))
             ix=max(iap,ibp)+1                                          9d20s16
             in=min(iap,ibp)                                            9d20s16
             ix=((ix*(ix-1))/2)+in
             iadh=ii0+ix                                                9d20s16
             iadj=jmatda+idm+nocc(isblkkder(1,isb))*(icm
     $         +nocc(isblkkder(2,isb))*(ib+nvirtc(isblkkder(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          else                                                          9d20s16
           do ia=0,nvirtc(isblkkder(4,isb))-1
            iap=ia+nocc(isblkkder(4,isb))
            do ib=0,nvirtc(isblkkder(3,isb))-1
             ibp=ib+nocc(isblkkder(3,isb))
             iadh=ii0+iap+nbasdwsc(isblkkder(4,isb))*ibp
             iadj=jmatda+idm+nocc(isblkkder(1,isb))*(icm
     $         +nocc(isblkkder(2,isb))*(ib+nvirtc(isblkkder(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d20s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then
        call ilimts(nocc(isblkkder(1,isb)),nocc(isblkkder(2,isb)),      9d20s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(1,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ia=0,nvirtc(isblkkder(4,isb))-1
           iap=ia+nocc(isblkkder(4,isb))
           do ib=0,nvirtc(isblkkder(3,isb))-1
            ibp=ib+nocc(isblkkder(3,isb))
            iadh=ii0+ibp+nbasdwsc(isblkkder(3,isb))*iap                 9d27s16
            iadj=jmatda+idm+nocc(isblkkder(1,isb))*(icm
     $         +nocc(isblkkder(2,isb))*(ib+nvirtc(isblkkder(3,isb))*ia))
            bc(iadj)=bc(iadj)+bc(iadh)                                  10d31s16
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|vv) ')
        call dws_sync
        call dws_finalize
        stop
       end if
       call dws_gsumf(bc(jmatda),nall)
       call ilimts(nvirtc(isblkkder(3,isb)),nvirtc(isblkkder(4,isb)),   8d22s16
     $      mynprocg,mynowprog,il,ih,i1s,i1e,is,i2e)                    8d22s16
       do i12=il,ih                                                     8d22s16
        i12m=i12-1                                                      8d22s16
        i12k=i12-il                                                     8d22s16
        iad1=jmatd(isb)+i12k*nrow                                       8d22s16
        iad2=jmatda+i12m*nrow                                           8d22s16
        do i34=0,nrow-1                                                 8d22s16
         bc(iad1+i34)=bc(iad1+i34)+bc(iad2+i34)                         9d16s16
        end do                                                          8d22s16
       end do
       ncolj=nvirtc(isblkkder(3,isb))*nvirtc(isblkkder(4,isb))
       nrowj=nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb))
       if(idwsdeb.gt.10)then
        write(6,*)('final jmat der ')
        write(6,12)(isblkkder(j,isb),j=1,4)
        ncolk=ih+1-il
        write(6,*)('locs '),jmatd(isb),jmatda
        call prntm2(bc(jmatd(isb)),nrowj,ih+1-il,nrowj)
        write(6,*)('jmatda '),jmatda
        call prntm2(bc(jmatda),nrowj,ncol,nrowj)
       end if
       ibcoff=jmatda                                                    8d22s16
      end do                                                            9d16s16
   12 format('integral type ',4i2,5x,i5)                                3d4s13\
      do isb=1,nsblkxder1                                                9d2s16
       nrow=nocc(isblkxder1(1,isb))*nocc(isblkxder1(2,isb))
       ncol=nocc(isblkxder1(3,isb))*nvirtc(isblkxder1(4,isb))
       ntot=nrow*ncol
       isw1=multh(isblkxder1(1,isb),iprop)
       isw2=multh(isblkxder1(2,isb),iprop)
       isw3=multh(isblkxder1(3,isb),iprop)
       isw4=multh(isblkxder1(4,isb),iprop)
c
c     (o'o|o`v): (o'o|av) part ok
c
c     recall we have (vv|oo')
c                     ab dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),isblkxder1(4,isb))
       itry2=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),isw3)
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(2,isb)),nocc(isblkxder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry)
        if(isblkxder1(4,isb).eq.isw3)then
         nrowh=(nbasdwsc(isw3)*(nbasdwsc(isw3)+1))/2
        else
         nrowh=nbasdwsc(isw3)*nbasdwsc(isblkxder1(4,isb))
        end if
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(isblkxder1(4,isb).eq.isw3)then
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nbasdwsc(isw3)-1
             ix=max(ia,ibp)
             in=min(ia,ibp)
             ixn=((ix*(ix+1))/2)+in
             iadh=ii0+ixn
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $      +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ib))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          else
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nbasdwsc(isw3)-1
             iadh=ii0+ia+nbasdwsc(isw3)*ibp
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $      +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ib))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then
        call ilimts(nocc(isblkxder1(2,isb)),nocc(isblkxder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw3)*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nbasdwsc(isw3)-1
           do ia=0,nvirtc(isblkxder1(4,isb))-1
            iap=ia+nocc(isblkxder1(4,isb))
            iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
            ff=bc(iadh)*2d0
            do m3=0,nocc(isblkxder1(3,isb))-1
             iad1=itrans(isw3)+ib+nbasdwsc(isw3)*m3
             iadk=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $      +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ia))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|av) for ionexd2 ')
        write(6,*)('want: '),isblkxder1(1,isb),isblkxder1(2,isb),
     $       isw4,isblkxder1(4,isb)
        write(6,*)('from ')
        do id=1,nsymb
         do ic=1,nsymb
          do ib=1,nsymb
           write(6,*)id,ic,ib,iptoh(id,ic,ib)
          end do
         end do
        end do
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|o`v): (oo'|av) part ok
c
c     recall we have (vv|oo')
c                     ab dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isblkxder1(4,isb))
       itry2=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isw3)
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nocc(isblkxder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry)
        if(isblkxder1(4,isb).eq.isw3)then
         nrowh=(nbasdwsc(isw3)*(nbasdwsc(isw3)+1))/2
        else
         nrowh=nbasdwsc(isw3)*nbasdwsc(isblkxder1(4,isb))
        end if
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(isblkxder1(4,isb).eq.isw3)then
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nbasdwsc(isw3)-1
             ix=max(ia,ibp)
             in=min(ia,ibp)
             ixn=((ix*(ix+1))/2)+in
             iadh=ii0+ixn
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $      +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ib))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          else
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nbasdwsc(isw3)-1
             iadh=ii0+ia+nbasdwsc(isw3)*ibp
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $      +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ib))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nocc(isblkxder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw3)*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nbasdwsc(isw3)-1
           do ia=0,nvirtc(isblkxder1(4,isb))-1
            iap=ia+nocc(isblkxder1(4,isb))
            iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
            ff=bc(iadh)*2d0
            do m3=0,nocc(isblkxder1(3,isb))-1
             iad1=itrans(isw3)+ib+nbasdwsc(isw3)*m3
             iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $      +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ia))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|av) for ionexd2 ')
        write(6,*)('want: '),isblkxder1(1,isb),isblkxder1(2,isb),
     $       isw4,isblkxder1(4,isb)
        write(6,*)('from ')
        do id=1,nsymb
         do ic=1,nsymb
          do ib=1,nsymb
           write(6,*)id,ic,ib,iptoh(id,ic,ib)
          end do
         end do
        end do
        call dws_sync
        call dws_finalize
        stop
       end if
       if(idwsdeb.gt.10.and.min(nrow,ncol).gt.0)then
       write(6,*)('final onexd2 ooox')
       write(6,12)(isblkxder1(j,isb),j=1,4)
       write(6,*)('not gsummed '),ionexd2(isb)
       call prntm2(bc(ionexd2(isb)),nrow,ncol,nrow)
       itmp=ibcoff
       ibcoff=itmp+nrow*ncol
       call enough('parajkfromhd1j.  6',bc,ibc)
       do i=0,nrow*ncol-1
        bc(itmp+i)=bc(ionexd2(isb)+i)
       end do
       call dws_gsumf(bc(itmp),nrow*ncol)
       write(6,*)('gsummed ')
       call prntm2(bc(itmp),nrow,ncol,nrow)
       if(nsymb.eq.1)call printa(bc(itmp),nocc,0,nocc,0,nocc,0,
     $      nvirtc,nocc,bc(ibcoff))
       ibcoff=itmp                                                      12d2s16
       end if
      end do
      return
      end
