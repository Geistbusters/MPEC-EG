c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parajkfromhd1(ipair,ngaus,ibdat,nhcolt,ihmat,nocc,icol,  6d22s10
     $     iorbx,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,     1d5s12
     $     ioooo,ionex,noc,iter,idwsdeb,ncomp,nvirtc,iprop,itrans,      3d21s16
     $     i4od,ionexd,kmatd,jmatd,i4od2b,ionexd2,der2e,nbasdwsc,        8d3s16
     $     isblkder,isblkxder,nsblkder,nsblkxder,isblkxder1,nsblkxder1, 9d16s16
     $     isblkder1,nsblkder1,isblkkder,nsblkkder,nbasisp,ndercode,    4d28s22
     $     igoal,idorel,bc,ibc)                                         11d10s22
      implicit real*8 (a-h,o-z)
c
c     like jkfromh, except for operators of symmetry iprop.
c     this version is for derivative contributions with the noc
c     contracted fcns differentiated once.
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
     $     ionexd2(idbk),isblkxder1(4,idbk),isblkder1(4,idbk),ihitit(4) 6d29s22
      call second(time1)                                                11d27s12
      if(idwsdeb.gt.10)write(6,*)('in parajkfromhd1')
c
c     in this version, do transpose of hmat, then transform full hmat
c     to mo basis, then store appropriate parts of the result.
c
c
      if(idwsdeb.gt.10)write(6,*)('xgoal at top '),bc(igoal)
      ntmpmat=2                                                         5d31s22
      if(idorel.eq.0.or.idorel.eq.1)then                                5d31s22
       ntmpmat=1                                                        5d31s22
      end if                                                            5d31s22
      kgoal=2248504
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
      nn=ngaus*ngaus                                                    3d2s16
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
          ii=iptoh(isd,isc,isb)                                         11d30s12
          ncd=nbasdwsc(isc)*nocc(isd)                                   3d24s16
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
       iqqq=ibcoff
       ibcoff=iqqq+nsymb**4
       call enough('parajkfromhd1.  1',bc,ibc)
       do iz=iqqq,ibcoff-1
        ibc(iz)=0
       end do
       do ip=0,mynprocg-1                                               2d11s13
        do isb=1,nsymb                                                  2d11s13
         nbasdwsb=nbasisp(isb)                                          5d31s22
         do isc=1,nsymb                                                 2d11s13
          isbc=multh(isc,isb)                                           3d2s16
          do isd=1,nsymb                                                2d11s13
           isad=multh(isd,isbc)                                          3d1s16
           isa=multh(isad,iprop)                                        3d1s16
           nbasdwsa=nbasisp(isa)                                        5d31s22
           ntotrun=0
           ncd=nbasdwsc(isc)*nocc(isd)                                  3d24s16
           iib=iptoh(isd,isc,isb)                                       2d12s13
           jz=iqqq+isa-1+nsymb*(isb-1+nsymb*(isc-1+nsymb*(isd-1)))
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
             nxn=nbasdwsb*nbasdwsa                                      3d1s16
             do ia=1,naa0                                               8d31s15
              iaa=ia+ibc(jbdat+ibc(jpair)+ngaus3)                           6d4s10
              issa=isstor(iaa)                                          3d1s16
              issad=multh(issa,iprop)                                    3d1s16
              iaaa=ibstor(iaa)                                              11d30s12
              do ib=1,nba0                                              8d31s15
               ibb=ib+ibc(jbdat+ibc(jpair+1)+ngaus3)                          6d4s10
               issb=isstor(ibb)                                                5d13s10
               ibbb=ibstor(ibb)                                               11d30s12
               if(issb.eq.isb.and.issa.eq.isa.and.iib.gt.0)then         5d31s22
                nxn=nbasdwsb*nbasdwsa                                   5d31s22
                ix=iaaa-1+nbasdwsa*(ibbb-1)                             5d31s22
                if(ipass.eq.1)then                                      2d11s13
                 irecv(ip+1)=irecv(ip+1)+nhere*ntmpmat                  10d30s20
                 myh(iib)=myh(iib)+nhere*ntmpmat                        10d30s20
                 ibc(jz)=ibc(jz)+nhere*ntmpmat
                 nsumx=nsumx+1
                else                                                    2d11s13
                 do in=0,ntmpmat-1                                      10d30s20
                  do j=0,nhere-1                                         2d11s13
                   iad=myh(iib)+ix+nxn*(j+nhere*in)                     10d30s20
                   bc(iad)=bc(ipr(ip+1))                                 2d11s13
                   ntotrun=ntotrun+1
                   ipr(ip+1)=ipr(ip+1)+1                                 2d11s13
                  end do                                                10d30s20
                 end do                                                 2d11s13
                end if                                                  2d11s13
               end if                                                   2d12s13
              end do            !ib
             end do            !ia
            end do             !i
           end if                                                       9d2s15
          end do      !isd                                                 2d11s13
         end do       !isc                                                  2d11s13
        end do                 !isb                                                  2d11s13
       end do                   !ip                                                  2d15s13
       if(ipass.eq.1)then                                               11d28s12
        nsend=0                                                         11d28s12
        nrecv=0                                                         11d28s12
        do ip=1,mynprocg                                                11d30s12
         nsend=nsend+isend(ip)                                          11d30s12
         nrecv=nrecv+irecv(ip)                                          11d30s12
        end do                                                          11d28s12
        do isb=1,nsymb
         do isc=1,nsymb
          isbc=multh(isc,isb)
          do isd=1,nsymb
           isad=multh(isd,isbc)
           jz=iqqq+isad-1+nsymb*(isb-1+nsymb*(isc-1+nsymb*(isd-1)))
          end do
         end do
        end do
        ngot=ibcoff+1-ihmat                                              11d30s12
        if(iter.eq.1)then
         write(6,*)('no. of words to receive:  '),nrecv                   11d30s12
         write(6,*)('no. of words under ihmat: '),ngot,ibcoff                   11d30s12
        end if
        if(nrecv.gt.ngot)then                                           11d30s12
         more=nrecv-ngot                                                11d30s12
         if(iter.eq.1.and.mynowprog.eq.0)then                            3d3s17
         write(6,*)('need to make some space for the receive ')         11d30s12
         write(6,*)('need '),more,(' additional words ')                11d30s12
         end if
         ibcoff=ihmat+nrecv                                             12d3s12
         call enough('parajkfromhd1.  2',bc,ibc)
        end if                                                          11d30s12
        ibufs=ibcoff                                                    11d28s12
        ibcoff=ibufs+max(nsend,nrecv)                                   11d30s12
        call enough('parajkfromhd1.  3',bc,ibc)
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
c     ihmat - raw ints: nbasdws(c)*nocc(d);nn
c     ibufs - ordered ints to send: nbasdws(c)*nocc(d);nn
c     ihmat - ints that have been recieved: nbasdws(c)*nocc(d);nn
c     ibufs - re-ordered ints: nbasdws(c)*nocc(d);nn
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
         ncd=nbasdwsc(isc)*nocc(isd)                                    3d24s16
         if(ncd.gt.0.and.iis.gt.0)then                                  12d21s12
          nhere=ncd/mynprocg                                            11d30s12
          nleft=ncd-nhere*mynprocg                                      11d30s12
          if(mynowprog.lt.nleft)nhere=nhere+1                           11d30s12
          if(nhere.gt.0)then                                            11d30s12
           ihcol(iis)=ihlast                                            3d5s13
           ihlast=ihlast+nxnc*nhere                                     4d12s22
           if(ihlast.gt.ibufs)then                                      3d5s13
            write(6,*)('ihlast now exceeds ibufs!!! '),ihlast,ibufs     3d5s13
            write(6,*)('no. of words '),ihlast+1-ihmat
            write(6,*)('nxn '),nxn
            write(6,*)('nhere '),nhere
            write(6,*)('nbasdwsa '),nbasdwsa
            write(6,*)('nbasdwsb '),nbasdwsb
            call dws_sync                                               3d5s13
            call dws_finalize                                           3d5s13
            stop                                                        3d5s13
           end if                                                       3d5s13
           itmp=ibcoff                                                  11d30s12
           itmp2=itmp+nxn                                               3d2s16
           itmp3=itmp2+nxn                                              5d31s22
           ibcoff=itmp3+nxn                                             5d31s22
           call enough('parajkfromhd1.  4',bc,ibc)
           nbasdwsaf=nbasdws(isa)                                       5d29s18
           nbasdwsbf=nbasdws(isb)                                       5d29s18
           if(min(nbasdws(isa),nbasdws(isb)).gt.0)then                  7d11s22
            do icol=1,nhere                                              11d30s12
             iads=ihcol(iis)+nxnc*(icol-1)                               5d30s18
             factr=0d0                                                   10d30s20
             jorbb=iorbx(isb)                                            10d30s20
             jorba=iorbx(isa)                                            10d30s20
             do in=0,ntmpmat-1                                           10d30s20
              iad=myh(iis)+nxn*(icol-1+nhere*in)                         10d30s20
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
              factr=1d0                                                  10d30s20
              jorbb=jorbb+nbasdwsb                                       10d30s20
              jorba=jorba+nbasdwsa                                       10d30s20
             end do                                                       10d30s20
             do i=1,nbasdwsbf                                            5d29s18
              do j=1,nbasdwsaf                                           5d29s18
               iad1=iads+j-1+nbasdwsaf*(i-1)                             5d29s18
               iad2=itmp3+i-1+nbasdwsbf*(j-1)                            10d30s20
               bc(iad1)=bc(iad2)                                         11d30s12
              end do                                                     11d30s12
             end do                                                      11d30s12
            end do                                                       11d30s12
           end if                                                       7d11s22
           ibcoff=itmp                                                  11d30s12
          end if                                                        11d30s12
         end if                                                         11d30s12
        end do                                                          11d30s12
       end do                                                           11d30s12
      end do                                                            11d30s12
      if(idwsdeb.gt.10)then
       write(6,*)('we have completed transformation ')
      end if                                                            5d19s22
      call second(time3)
      n4o3o=0                                                             7d26s16
      if(idwsdeb.gt.10)write(6,*)('first derivative part '),bc(igoal)   5d19s22
      do is=1,nsblkxder
       nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
       ncol=nocc(isblkxder(3,is))*nocc(isblkxder(4,is))
       nn=nrow*ncol
       if(nn.gt.0)then
        if(idwsdeb.gt.10)then
         write(6,12)(isblkxder(j,is),j=1,4)
         write(6,*)('starting i4od: ')
         if(mynprocg.eq.1)then
          call prntm2(bc(i4od(is)),nrow,ncol,nrow)
         else
          do i=0,nrow*ncol-1                                            5d31s22
           bc(ibcoff+i)=bc(i4od(is)+i)                                  5d31s22
          end do                                                        5d31s22
          call dws_gsumf(bc(ibcoff),nrow*ncol)
          call prntm2(bc(ibcoff),nrow,ncol,nrow)                        5d31s22
         end if
        end if
        n4o3o=n4o3o+nn                                                      7d26s16
c
c     under ihcol, der is wrt to first index, syma
c     under i4o, want der to be last index, sym4
c     since 4o, we can swap isd and isc, ie 2st two args of
c     iptoh.
c
        itry=iptoh(isblkxder(1,is),isblkxder(2,is),isblkxder(3,is))
        if(itry.gt.0)then
         isd=isblkxder(1,is)
         isc=isblkxder(2,is)
         isb=isblkxder(3,is)
         isa=isblkxder(4,is)
         nh=nbasdwsc(isc)
         call ilimts(nocc(isd),nh,mynprocg,mynowprog,ilh,ihh,i1s,i1e,
     $        i2s,i2e)
         i10=i1s
         i1n=nocc(isd)
         ii0=ihcol(itry)
         nha=nbasdwsc(isa)
         nrowh=nha*nbasdwsc(isb)                                        3d24s16
         do i2=i2s,i2e
          if(i2.eq.i2e)i1n=i1e
          i2m=i2-1
          do i1=i10,i1n
           i1m=i1-1
           if(i2.le.nocc(isc))then
            do i3=0,nocc(isb)-1
             do i4=0,nocc(isa)-1
              iad1=ii0+i4+nha*i3
              iad5=i4od(is)+i1m+nocc(isd)*(i2m+nocc(isc)
     $             *(i3+nocc(isb)*i4))
              bc(iad5)=bc(iad5)+bc(iad1)                                7d26s16
             end do
            end do
           end if
           ii0=ii0+nrowh
          end do
          i10=1
         end do
        else
         itry2=iptoh(isblkxder(2,is),isblkxder(1,is),isblkxder(3,is))
         if(itry2.gt.0)then
          isd=isblkxder(2,is)
          isc=isblkxder(1,is)
          isb=isblkxder(3,is)
          isa=isblkxder(4,is)
          nh=nbasdwsc(isc)
          call ilimts(nocc(isd),nh,mynprocg,mynowprog,ilh,ihh,i1s,i1e,
     $        i2s,i2e)
          i10=i1s
          i1n=nocc(isd)
          ii0=ihcol(itry)
          nha=nbasdwsc(isa)
          nrowh=nha*nbasdwsc(isb)                                        3d24s16
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           i2m=i2-1
           do i1=i10,i1n
            i1m=i1-1
            if(i2.le.nocc(isc))then
             do i3=0,nocc(isb)-1
              do i4=0,nocc(isa)-1
               iad1=ii0+i4+nha*i3
               iad5=i4od(is)+i2m+nocc(isc)*(i1m+nocc(isd)               4d18s16
     $              *(i3+nocc(isb)*i4))
               bc(iad5)=bc(iad5)+bc(iad1)                               7d26s16
              end do
             end do
            end if
            ii0=ii0+nrowh
           end do
           i10=1
          end do
         else
          write(6,*)('unable to find match for this type in hmat '),
     $        (isblkxder(j,is),j=1,4)
          write(6,*)('searching on ')
          do is1=1,nsymb                                                8d21s23
           do is2=1,nsymb                                               8d21s23
            do is3=1,nsymb                                              8d21s23
             if(iptoh(is3,is2,is1).gt.0)write(6,*)is3,is2,is1
            end do
           end do
          end do
          call dws_sync
          call dws_finalize
          stop
         end if
        end if
       end if
      end do
      ilookat=ionexd(1)
      if(ionexd(1).gt.0)then                                            8d17s24
       if(idwsdeb.gt.10)write(6,*)('ionexd part '),bc(igoal)
       do is=1,nsblkxder
        nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
        ncol=nocc(isblkxder(3,is))*nvirtc(isblkxder(4,is))
        nn=nrow*ncol
        if(nn.gt.0)then
         n4o3o=n4o3o+nn                                                      7d26s16
         if(idwsdeb.gt.10)then
          write(6,*)('starting ionexd: '),ionexd(is),bc(ilookat)
          write(6,12)(isblkxder(j,is),j=1,4)
          if(mynprocg.eq.1)then                                          5d31s22
           call prntm2(bc(ionexd(is)),nrow,ncol,nrow)
          else                                                           5d31s22
           do i=0,nrow*ncol-1                                            5d31s22
            bc(ibcoff+i)=bc(ionexd(is)+i)                                5d31s22
           end do                                                        5d31s22
           call dws_gsumf(bc(ibcoff),nrow*ncol)                          5d31s22
           call prntm2(bc(ibcoff),nrow,ncol,nrow)                        5d31s22
          end if                                                         5d31s22
         end if
c
c     under ihcol, der is wrt to first index, syma
c     for onex, need o'oox+oo'ox+ooo'x+ooox'
c
         itry=iptoh(isblkxder(3,is),isblkxder(4,is),isblkxder(2,is))
         if(itry.gt.0)then
          isd=isblkxder(3,is)
          isc=isblkxder(4,is)
          isb=isblkxder(2,is)
          isa=isblkxder(1,is)
          nh=nbasdwsc(isc)                                               3d24s16
          call ilimts(nocc(isd),nh,mynprocg,mynowprog,ilh,ihh,i1s,i1e,
     $        i2s,i2e)
          i10=i1s
          i1n=nocc(isd)
          ii0=ihcol(itry)
          nha=nbasdwsc(isa)                                              3d24s16
          nrowh=nha*nbasdwsc(isb)                                        3d24s16
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           i2m=i2-1-nocc(isc)
           do i1=i10,i1n
            i1m=i1-1
            if(i2.gt.nocc(isc))then
             do i3=0,nocc(isb)-1
              do i4=0,nocc(isa)-1
               iad1=ii0+i4+nha*i3
               iad5=ionexd(is)+i4+nocc(isa)*(i3+nocc(isb)*(i1m
     $             +nocc(isd)*i2m))
               bc(iad5)=bc(iad5)+bc(iad1)                                7d26s16
              end do
             end do
            end if
            ii0=ii0+nrowh
           end do
           i10=1
          end do
         else if(nocc(isblkxder(3,is)).ne.0)then                         4d18s16
          write(6,*)('unable to find match for this type in hmat ')
          write(6,*)('o''oox'),(isblkxder(j,is),j=1,4)
          call dws_sync
          call dws_finalize
          stop
         end if
         itry=iptoh(isblkxder(3,is),isblkxder(4,is),isblkxder(1,is))
         nprog=nocc(multh(isblkxder(2,is),iprop))                        4d19s16
         if(itry.gt.0)then
          isd=isblkxder(3,is)
          isc=isblkxder(4,is)
          isb=isblkxder(1,is)
          isa=isblkxder(2,is)
          nh=nbasdwsc(isc)                                               3d24s16
          call ilimts(nocc(isd),nh,mynprocg,mynowprog,ilh,ihh,i1s,i1e,
     $        i2s,i2e)
          i10=i1s
          i1n=nocc(isd)
          ii0=ihcol(itry)
          nha=nbasdwsc(isa)                                              3d24s16
          nrowh=nha*nbasdwsc(isb)                                        3d24s16
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           i2m=i2-1-nocc(isc)
           do i1=i10,i1n
            i1m=i1-1
            if(i2.gt.nocc(isc))then
             do i3=0,nocc(isb)-1
              do i4=0,nocc(isa)-1
               iad1=ii0+i4+nha*i3
               iad5=ionexd(is)+i3+nocc(isb)*(i4+nocc(isa)*(i1m
     $             +nocc(isd)*i2m))
               bc(iad5)=bc(iad5)+bc(iad1)
              end do
             end do
            end if
            ii0=ii0+nrowh
           end do
           i10=1
          end do
         else if(nprog.ne.0)then                                         4d19s16
          write(6,*)('unable to find match for this type in hmat ')
          write(6,*)('oo''ox'),(isblkxder(j,is),j=1,4)
          write(6,*)('nprog = '),nprog
          call dws_sync
          call dws_finalize
          stop
         end if
         itry=iptoh(isblkxder(2,is),isblkxder(1,is),isblkxder(4,is))
         nprog=nocc(multh(isblkxder(3,is),iprop))                        4d19s16
         if(itry.gt.0)then
          isd=isblkxder(2,is)
          isc=isblkxder(1,is)
          isb=isblkxder(4,is)
          isa=isblkxder(3,is)
          nh=nbasdwsc(isc)                                               3d24s16
          call ilimts(nocc(isd),nh,mynprocg,mynowprog,ilh,ihh,i1s,i1e,
     $        i2s,i2e)
          i10=i1s
          i1n=nocc(isd)
          ii0=ihcol(itry)
          nha=nbasdwsc(isa)                                              3d24s16
          nrowh=nha*nbasdwsc(isb)                                        3d24s16
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           i2m=i2-1
           do i1=i10,i1n
            i1m=i1-1
            if(i2.le.nocc(isc))then
             do i3=0,nvirtc(isb)-1
              i3p=i3+nocc(isb)
              do i4=0,nocc(isa)-1
               iad1=ii0+i4+nha*i3p
               iad5=ionexd(is)+i2m+nocc(isc)*(i1m+nocc(isd)*(i4
     $             +nocc(isa)*i3))
               bc(iad5)=bc(iad5)+bc(iad1)
              end do
             end do
            end if
            ii0=ii0+nrowh
           end do
           i10=1
          end do
         else if(nprog.gt.0)then                                         4d19s16
          write(6,*)('unable to find match for this type in hmat ')
          write(6,*)('ooo''x'),(isblkxder(j,is),j=1,4)                  4d19s16
          write(6,*)('nprog = '),nprog                                  4d19s16
          call dws_sync
          call dws_finalize
          stop
         end if
         itry1=iptoh(isblkxder(2,is),isblkxder(1,is),isblkxder(3,is))
         itry2=iptoh(isblkxder(2,is),isblkxder(1,is),isblkxder(4,is))
         if(itry1.gt.0)then
          isd=isblkxder(2,is)
          isc=isblkxder(1,is)
          isb=isblkxder(3,is)
          isa=isblkxder(4,is)
          nh=nbasdwsc(isc)                                               3d24s16
          call ilimts(nocc(isd),nh,mynprocg,mynowprog,ilh,ihh,i1s,i1e,
     $        i2s,i2e)
          i10=i1s
          i1n=nocc(isd)
          ii0=ihcol(itry1)                                               7d11s22
          nha=nbasdwsc(isa)                                              3d24s16
          nrowh=nha*nbasdwsc(isb)                                        3d24s16
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           i2m=i2-1
           do i1=i10,i1n
            i1m=i1-1
            if(i2.le.nocc(isc))then
             do i3=0,nocc(isb)-1
              do i4=0,nvirtc(isa)-1
               iad1=ii0+i4+nocc(isa)+nha*i3
               iad5=ionexd(is)+i2m+nocc(isc)*(i1m+nocc(isd)*(i3
     $              +nocc(isb)*i4))
               bc(iad5)=bc(iad5)+bc(iad1)
              end do
             end do
            end if
            ii0=ii0+nrowh
           end do
           i10=1
          end do
         else if(itry2.gt.0)then                                         7d11s22
          isd=isblkxder(2,is)
          isc=isblkxder(1,is)
          isb=isblkxder(3,is)                                            7d11s22
          isa=isblkxder(4,is)                                            7d11s22
          nh=nbasdwsc(isc)                                               3d24s16
          call ilimts(nocc(isd),nh,mynprocg,mynowprog,ilh,ihh,i1s,i1e,
     $        i2s,i2e)
          i10=i1s
          i1n=nocc(isd)
          ii0=ihcol(itry2)                                               7d11s22
          nhb=nbasdwsc(isb)                                              7d11s22
          nrowh=nhb*nbasdwsc(isa)                                        7d11s22
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           i2m=i2-1
           do i1=i10,i1n
            i1m=i1-1
            if(i2.le.nocc(isc))then
             do i3=0,nocc(isb)-1
              do i4=0,nvirtc(isa)-1
               iad1=ii0+i3+nocc(isb)+nhb*i4                              7d11s22
               iad5=ionexd(is)+i2m+nocc(isc)*(i1m+nocc(isd)*(i3
     $             +nocc(isb)*i4))
               bc(iad5)=bc(iad5)+bc(iad1)
              end do
             end do
            end if
            ii0=ii0+nrowh
           end do
           i10=1
          end do
         else
          write(6,*)('unable to find match for this type in hmat ')
          write(6,*)('ooox'' '),(isblkxder(j,is),j=1,4)                  4d19s16
          itry1=iptoh(isblkxder(2,is),isblkxder(1,is),isblkxder(3,is))
          itry2=iptoh(isblkxder(2,is),isblkxder(1,is),isblkxder(4,is))
          write(6,*)('try1,try2 '),itry1,itry2
          call dws_sync
          call dws_finalize
          stop
         end if
        end if
       end do
      end if                                                            8d17s24
      call dws_gsumf(bc(i4od(1)),n4o3o)
      if(idwsdeb.gt.10.and.ionexd(1).gt.0)then                          8d17s24
       if(nsblkxder1.gt.0)then                                           10d3s16
        write(6,*)('for full (ooox)'' ')
       else
        write(6,*)('for full (ooox)" ')
       end if
       do is=1,nsblkxder
        nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
        ncol=nocc(isblkxder(3,is))*nvirtc(isblkxder(4,is))
        write(6,*)('consider is = '),is,('type '),
     $       (isblkxder(j,is),j=1,4),nrow,ncol
        if(min(nrow,ncol).gt.0)then
         write(6,12)(isblkxder(j,is),j=1,4)
         write(6,*)('args: '),ionexd(is),nrow,ncol,bc(ilookat)
          call prntm2(bc(ionexd(is)),nrow,ncol,nrow)
        end if
       end do
      end if                                                            3d21s16
      if(idwsdeb.gt.10)write(6,*)('now form (oooo)" '),bc(igoal)
      im0=ibcoff
      do is=1,nsblkder
       nrow=nocc(isblkder(1,is))*nocc(isblkder(2,is))
       ncol=nocc(isblkder(3,is))*nocc(isblkder(4,is))
       if(nrow*ncol.gt.0)then
        imt=ibcoff
        ibcoff=imt+nrow*ncol
        call enough('parajkfromhd1.  5',bc,ibc)
        do i=0,nrow*ncol-1
         bc(imt+i)=0d0
        end do
        ihitit(1)=0                                                     6d29s22
        ihitit(2)=0                                                     6d29s22
        ihitit(3)=0                                                     7d1s22
        ihitit(4)=0                                                     6d29s22
        do is2=1,nsblkxder
         if(isblkder(1,is).eq.isblkxder(4,is2).and.
     $        isblkder(2,is).eq.isblkxder(3,is2))then
          if(isblkder(3,is).eq.isblkxder(1,is2))then
           ihitit(1)=1
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i4+nocc(isblkxder(4,is2))*(i3
     $              +nocc(isblkxder(3,is2))*(i1
     $              +nocc(isblkxder(1,is2))*i2))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(3,is).eq.isblkxder(2,is2))then               6d29s22
           ihitit(1)=1
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i4+nocc(isblkxder(4,is2))*(i3
     $              +nocc(isblkxder(3,is2))*(i2
     $              +nocc(isblkxder(2,is2))*i1))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
         if(isblkder(2,is).eq.isblkxder(4,is2).and.
     $        isblkder(1,is).eq.isblkxder(3,is2))then
          if(isblkder(3,is).eq.isblkxder(1,is2))then
           iad1=i4od(is2)
           ihitit(2)=1
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i3+nocc(isblkxder(3,is2))*(i4
     $              +nocc(isblkxder(4,is2))*(i1
     $              +nocc(isblkxder(1,is2))*i2))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(3,is).eq.isblkxder(2,is2))then               6d29s22
           iad1=i4od(is2)
           ihitit(2)=1
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i3+nocc(isblkxder(3,is2))*(i4
     $              +nocc(isblkxder(4,is2))*(i2                         6d29s22
     $              +nocc(isblkxder(2,is2))*i1))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
         if(isblkder(3,is).eq.isblkxder(4,is2).and.
     $        isblkder(4,is).eq.isblkxder(3,is2))then
          if(isblkder(1,is).eq.isblkxder(1,is2))then
           iad1=i4od(is2)
           ihitit(3)=1
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i1+nocc(isblkxder(1,is2))*(i2
     $              +nocc(isblkxder(2,is2))*(i4
     $              +nocc(isblkxder(4,is2))*i3))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(1,is).eq.isblkxder(2,is2))then
           iad1=i4od(is2)
           ihitit(3)=1
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i2+nocc(isblkxder(2,is2))*(i1
     $              +nocc(isblkxder(1,is2))*(i4
     $              +nocc(isblkxder(4,is2))*i3))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
         if(isblkder(4,is).eq.isblkxder(4,is2).and.
     $        isblkder(3,is).eq.isblkxder(3,is2))then
          if(isblkder(1,is).eq.isblkxder(1,is2))then
           iad1=i4od(is2)
           ihitit(4)=1
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i1+nocc(isblkxder(1,is2))*(i2
     $              +nocc(isblkxder(2,is2))*(i3
     $              +nocc(isblkxder(3,is2))*i4))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(1,is).eq.isblkxder(2,is2))then
           iad1=i4od(is2)
           ihitit(4)=1
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i2+nocc(isblkxder(2,is2))*(i1
     $              +nocc(isblkxder(1,is2))*(i3
     $              +nocc(isblkxder(3,is2))*i4))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
        end do
        nhitit=ihitit(1)+ihitit(2)+ihitit(3)+ihitit(4)                  6d29s22
        if(nhitit.ne.4)then                                             6d29s22
         write(6,*)('we missed an index!!! '),ihitit                    6d29s22
         write(6,*)('for oooo in parajkfromhd1 ')                       6d29s22
         call dws_synca                                                 6d29s22
         call dws_finalize                                              6d29s22
         stop                                                           6d29s22
        end if                                                          6d29s22
        if(nsblkxder1.lt.-1)then                                        5d9s22
         if(idwsdeb.gt.10)then
          write(6,12)(isblkder(j,is),j=1,4)
          call prntm2(bc(imt),nrow,ncol,nrow)
          write(6,*)('adding in whats under i4od2b ')                    10d3s16
         end if                                                         12d2s16
         call dws_gsumf(bc(i4od2b(is)),nrow*ncol)                       10d3s16
         if(idwsdeb.gt.10)then
          call prntm2(bc(i4od2b(is)),nrow,ncol,nrow)
          end if                                                        12d2s16
         do j=0,nrow*ncol-1                                             10d3s16
          bc(imt+j)=bc(imt+j)+bc(i4od2b(is)+j)
         end do                                                         10d3s16
        else
        end if                                                          10d3s16
        if(idwsdeb.gt.10)then
         if(nrow*ncol.gt.0)then                                         4d19s16
          write(6,12)(isblkder(j,is),j=1,4)
          call prntm2(bc(imt),nrow,ncol,nrow)
          if(nsymb.eq.1)call printa(bc(imt),nocc,0,nocc,0,nocc,0,
     $         nocc,0,bc(ibcoff))
         end if                                                         4d19s16
        end if
       end if
      end do
      imt=im0
c
c     i want to copy imt to i4od, however i4od is pointed to by
c     nsblkxder and isblkxder while imt is point to by nsblkder and
c     isblkder. now the storage initially allocated for i4od will
c     not be smaller than that allocated for imt, but the symmetries
c     (and sizes) pointed to by isblkxder could be different from
c     those used in imt. however since memory for both were allocated
c     contiguously, we just need to reset the i4od(i),i gt 1.
c
      i4od0=i4od(1)                                                     12d8s16
      if(idwsdeb.gt.10)write(6,*)('next part '),bc(igoal)
      do is=1,nsblkder                                                  12d8s16
       nrow=nocc(isblkder(1,is))*nocc(isblkder(2,is))                   12d8s16
       ncol=nocc(isblkder(3,is))*nocc(isblkder(4,is))                   12d8s16
       nwds=nrow*ncol                                                   12d8s16
       if(nwds.gt.0)then                                                12d8s16
        i4od(is)=i4od0                                                  12d8s16
        do i=0,nwds-1                                                   12d8s16
         bc(i4od(is)+i)=bc(imt+i)                                       12d8s16
        end do                                                          12d8s16
        imt=imt+nwds                                                    12d8s16
        i4od0=i4od0+nwds                                                12d8s16
       end if                                                           12d8s16
      end do                                                            12d8s16
      ibcoff=im0                                                        12d8s16
      der2e=0d0
      e2j=0d0
      e2k=0d0
      do is=1,nsblkder
       imt=i4od(is)                                                     12d8s16
       nrow=nocc(isblkder(1,is))*nocc(isblkder(2,is))
       ncol=nocc(isblkder(3,is))*nocc(isblkder(4,is))
       if(min(nrow,ncol).gt.0)then                                      4d28s22
        if(idwsdeb.gt.10)then                                           4d28s22
         if(ndercode.eq.1)then                                           4d28s22
          write(6,*)('final (oooo)'' for sym '),(isblkder(j,is),j=1,4)
         else if(ndercode.eq.2)then                                      4d28s22
          write(6,*)('final (oooo)" for sym '),(isblkder(j,is),j=1,4)
         end if                                                          4d28s22
         call prntm2(bc(imt),nrow,ncol,nrow)
         if(nsymb.eq.1)call printa(bc(imt),nocc,0,nocc,0,nocc,0,nocc,0,
     $        bc(ibcoff))
        end if                                                          4d28s22
        if(isblkder(1,is).eq.isblkder(2,is).and.                         4d4s16
     $      isblkder(3,is).eq.isblkder(4,is))then                       4d4s16
         fmul=2d0
         do i1=0,nocc(isblkder(1,is))-1
          do i2=0,nocc(isblkder(3,is))-1
           iad=imt+i1+nocc(isblkder(1,is))*(i1+nocc(isblkder(1,is))
     $         *(i2+nocc(isblkder(3,is))*i2))
           der2e=der2e+fmul*bc(iad)
          end do
         end do
        end if
        if(isblkder(1,is).eq.isblkder(3,is).and.                         4d4s16
     $      isblkder(2,is).eq.isblkder(4,is))then                       4d4s16
         if(isblkder(1,is).ne.isblkder(2,is))then                        4d4s16
          fmul=2d0                                                       4d4s16
         else                                                            4d4s16
          fmul=1d0                                                       4d4s16
         end if                                                          4d4s16
         do i1=0,nocc(isblkder(1,is))-1
          do i2=0,nocc(isblkder(2,is))-1
           iad=imt+i1+nocc(isblkder(1,is))*(i2+nocc(isblkder(2,is))
     $         *(i1+nocc(isblkder(1,is))*i2))
           der2e=der2e-bc(iad)*fmul                                      4d4s16
          end do
         end do
        end if
       end if                                                           4d28s22
       imt=imt+nrow*ncol
      end do
   12 format('integral type ',4i2,5x,i5)                                3d4s13\
c                                                                       9d16s16
c     now consider ders of j,k and mixed 2nd ders of 4o and onex.       9d16s16
c
c
c     for j need
c     (o'o|vv), (oo'|vv), (oo|v'v) and (oo|vv')
c     recall we have (v'v|ov)
c                     a b dc  stored under iptoh(d,c,b),
c     thus presently we can only get the last two
c     recall k_nm^ab=(nb|ma)
c     for k need
c     (o'v|ov), (ov'|ov), (ov|o'v) and (ov|ov').
c     we can get them all.
c     for mixed 4o we need
c     (o'a|oo), (o'o|ao), (o'o|oa) and (a'o|oo)
c     (ao'|oo), (oo'|ao), (oo'|oa) and (oa'|oo)
c     (ao|o'o), (oa|o'o), (oo|o'a) and (oo|a'o)
c     (ao|oo'), (oa|oo'), (oo|ao') and (oo|oa')
c     we can get them all.
c     for mixed onex we need
c     (o'a|ov), (o'o|av), (o'o|oa) and (a'o|ov)
c     (ao'|ov), (oo'|av), (oo'|oa) and (oa'|ov)
c     (ao|o'v), (oa|o'v), (oo|o'a) and (oo|a'v)
c     (ao|ov'), (oa|ov'), (oo|av') and (oo|oa')
c     we can get all except ('oo|av) and (oo'|av).
c                                                                       9d16s16
      do isb=1,nsblkkder                                                9d16s16
       nrow=nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb))               8d22s16
       ncol=nvirtc(isblkkder(3,isb))*nvirtc(isblkkder(4,isb))            11d30s16
       nall=nrow*ncol                                                   8d4s16
       call ilimts(nvirtc(isblkkder(3,isb)),nvirtc(isblkkder(4,isb)),   8d22s16
     $      mynprocg,mynowprog,il,ih,i1s,i1e,is,i2e)                    8d22s16
       jmatda=ibcoff
       ibcoff=jmatda+nall
       kmatda=ibcoff                                                    9d16s16
       ibcoff=kmatda+nall                                               9d16s16
       ntot=nall*2                                                      9d16s16
       call enough('parajkfromhd1.  6',bc,ibc)
       do i=0,ntot-1                                                    8d4s16
        bc(jmatda+i)=0d0                                                8d4s16
       end do                                                           8d4s16
c
c     (oo|v'v) part ok
c
       itry=iptoh(isblkkder(1,isb),isblkkder(2,isb),isblkkder(4,isb))   9d16s16
       itry2=iptoh(isblkkder(2,isb),isblkkder(1,isb),isblkkder(4,isb))  9d16s16
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(1,isb)),nbasdwsc(isblkkder(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        do i2=i2s,i2e
         i2m=i2-1
         if(i2.eq.i2e)i1n=i1e
         do i1=i10,i1n
          i1m=i1-1
          if(i2.le.nocc(isblkkder(2,isb)))then
           do i4=0,nvirtc(isblkkder(4,isb))-1
            i4p=i4+nocc(isblkkder(4,isb))
            do i3=0,nvirtc(isblkkder(3,isb))-1
             i3p=i3+nocc(isblkkder(3,isb))
             iadh=ii0+i3p+nbasdwsc(isblkkder(3,isb))*i4p
             iadj=jmatda+i1m+nocc(isblkkder(1,isb))*(i2m
     $         +nocc(isblkkder(2,isb))*(i3+nvirtc(isblkkder(3,isb))*i4))
             bc(iadj)=bc(iadj)+bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d16s16
        call ilimts(nocc(isblkkder(2,isb)),nbasdwsc(isblkkder(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        do i2=i2s,i2e
         i2m=i2-1
         if(i2.eq.i2e)i1n=i1e
         do i1=i10,i1n
          i1m=i1-1
          if(i2.le.nocc(isblkkder(1,isb)))then
           do i4=0,nvirtc(isblkkder(4,isb))-1
            i4p=i4+nocc(isblkkder(4,isb))
            do i3=0,nvirtc(isblkkder(3,isb))-1
             i3p=i3+nocc(isblkkder(3,isb))
             iadh=ii0+i3p+nbasdwsc(isblkkder(3,isb))*i4p
             iadj=jmatda+i2m+nocc(isblkkder(1,isb))*(i1m
     $         +nocc(isblkkder(2,isb))*(i3+nvirtc(isblkkder(3,isb))*i4))
             bc(iadj)=bc(iadj)+bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find h for (oo|v''v) '),isblkkder(1,isb),
     $       isblkkder(2,isb),isblkkder(3,isb),isblkkder(4,isb)
        itry3=iptoh(isblkkder(1,isb),isblkkder(2,isb),isblkkder(3,isb))   9d16s16
        itry4=iptoh(isblkkder(2,isb),isblkkder(1,isb),isblkkder(3,isb))  9d16s16
        do isx1=1,nsymb
         do isx2=1,nsymb
          do isx3=1,nsymb
           write(6,*)isx1,isx2,isx3,iptoh(isx1,isx2,isx3)
          end do
         end do
        end do
        write(6,*)('itry3,itry4 '),itry3,itry4
        call dws_sync
        call dws_finalize
        stop
       end if                                                           9d16s16
c
c     (oo|vv') part ok
c
       itry=iptoh(isblkkder(1,isb),isblkkder(2,isb),isblkkder(3,isb))   9d16s16
       itry2=iptoh(isblkkder(2,isb),isblkkder(1,isb),isblkkder(3,isb))  9d16s16
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(1,isb)),nbasdwsc(isblkkder(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        do i2=i2s,i2e
         i2m=i2-1
         if(i2.eq.i2e)i1n=i1e
         do i1=i10,i1n
          i1m=i1-1
          if(i2.le.nocc(isblkkder(2,isb)))then
           do i4=0,nvirtc(isblkkder(4,isb))-1
            i4p=i4+nocc(isblkkder(4,isb))
            do i3=0,nvirtc(isblkkder(3,isb))-1
             i3p=i3+nocc(isblkkder(3,isb))
             iadh=ii0+i4p+nbasdwsc(isblkkder(4,isb))*i3p
             iadj=jmatda+i1m+nocc(isblkkder(1,isb))*(i2m
     $         +nocc(isblkkder(2,isb))*(i3+nvirtc(isblkkder(3,isb))*i4))
             bc(iadj)=bc(iadj)+bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d16s16
        call ilimts(nocc(isblkkder(2,isb)),nbasdwsc(isblkkder(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))
        do i2=i2s,i2e
         i2m=i2-1
         if(i2.eq.i2e)i1n=i1e
         do i1=i10,i1n
          i1m=i1-1
          if(i2.le.nocc(isblkkder(1,isb)))then
           do i4=0,nvirtc(isblkkder(4,isb))-1
            i4p=i4+nocc(isblkkder(4,isb))
            do i3=0,nvirtc(isblkkder(3,isb))-1
             i3p=i3+nocc(isblkkder(3,isb))
             iadh=ii0+i4p+nbasdwsc(isblkkder(4,isb))*i3p
             iadj=jmatda+i2m+nocc(isblkkder(1,isb))*(i1m
     $         +nocc(isblkkder(2,isb))*(i3+nvirtc(isblkkder(3,isb))*i4))
             bc(iadj)=bc(iadj)+bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find h for (oo|vv'') ')
        write(6,12)(isblkkder(j,isb),j=1,4)
        write(6,*)(nocc(isblkkder(j,isb)),j=1,2),
     $       (nvirtc(isblkkder(j,isb)),j=3,4)
        itest=multh(isblkkder(4,isb),iprop)
        write(6,12)(isblkkder(j,isb),j=1,3),itest
     $
        write(6,*)(nocc(isblkkder(j,isb)),j=1,2),
     $       nvirtc(isblkkder(3,isb)),nvirtc(itest)
        call dws_sync
        call dws_finalize
        stop
       end if                                                           9d16s16
c
c     (o'v|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       itry=iptoh(isblkkder(2,isb),isblkkder(3,isb),isblkkder(4,isb))   9d16s16
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(2,isb)),nbasdwsc(isblkkder(3,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(2,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkkder(1,isb))*nbasdwsc(isblkkder(4,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkkder(3,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(icm.ge.0)then
           do ib=0,nvirtc(isblkkder(4,isb))-1
            ibp=ib+nocc(isblkkder(4,isb))
            do ia=0,nocc(isblkkder(1,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkkder(1,isb))*ibp
             iadk=kmatda+ia+nocc(isblkkder(1,isb))*(idm
     $        +nocc(isblkkder(2,isb))*(icm+nvirtc(isblkkder(3,isb))*ib))
             bc(iadk)=bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find h for (o''v|ov) ')
        call dws_sync
        call dws_finalize
        stop
       end if                                                           9d16s16
c
c     (ov'|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       itry=iptoh(isblkkder(2,isb),isblkkder(3,isb),isblkkder(1,isb))   9d16s16
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(2,isb)),nbasdwsc(isblkkder(3,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(2,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkkder(1,isb))*nbasdwsc(isblkkder(4,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkkder(3,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(icm.ge.0)then
           do ib=0,nocc(isblkkder(1,isb))-1
            do ia=0,nvirtc(isblkkder(4,isb))-1
             iap=ia+nocc(isblkkder(4,isb))
             iadh=ii0+iap+nbasdwsc(isblkkder(4,isb))*ib
             iadk=kmatda+ib+nocc(isblkkder(1,isb))*(idm
     $        +nocc(isblkkder(2,isb))*(icm+nvirtc(isblkkder(3,isb))*ia))
             bc(iadk)=bc(iadk)+bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find h for (ov''|ov) ')
        call dws_sync
        call dws_finalize
        stop
       end if                                                           9d16s16
c
c     (ov|o'v) part ok
c     recall k_nm^ab=(nb|ma)
c
       itry=iptoh(isblkkder(1,isb),isblkkder(4,isb),isblkkder(3,isb))   9d16s16
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(1,isb)),nbasdwsc(isblkkder(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkkder(2,isb))*nbasdwsc(isblkkder(3,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkkder(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(icm.ge.0)then
           do ib=0,nvirtc(isblkkder(3,isb))-1
            ibp=ib+nocc(isblkkder(3,isb))
            do ia=0,nocc(isblkkder(2,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkkder(2,isb))*ibp
             iadk=kmatda+idm+nocc(isblkkder(1,isb))*(ia
     $        +nocc(isblkkder(2,isb))*(ib+nvirtc(isblkkder(3,isb))*icm))
             bc(iadk)=bc(iadk)+bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find h for (ov|o''v) ')
        call dws_sync
        call dws_finalize
        stop
       end if                                                           9d16s16
c
c     (ov|ov') part ok
c     recall k_nm^ab=(nb|ma)
c
       itry=iptoh(isblkkder(1,isb),isblkkder(4,isb),isblkkder(2,isb))   9d16s16
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkkder(1,isb)),nbasdwsc(isblkkder(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkkder(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkkder(2,isb))*nbasdwsc(isblkkder(3,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkkder(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(icm.ge.0)then
           do ib=0,nocc(isblkkder(2,isb))-1
            do ia=0,nvirtc(isblkkder(3,isb))-1
             iap=ia+nocc(isblkkder(3,isb))
             iadh=ii0+iap+nbasdwsc(isblkkder(3,isb))*ib
             iadk=kmatda+idm+nocc(isblkkder(1,isb))*(ib
     $        +nocc(isblkkder(2,isb))*(ia+nvirtc(isblkkder(3,isb))*icm))
             bc(iadk)=bc(iadk)+bc(iadh)
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find h for (ov|o''v) ')
        call dws_sync
        call dws_finalize
        stop
       end if                                                           9d16s16
       call dws_gsumf(bc(jmatda),ntot)
       call ilimts(nvirtc(isblkkder(3,isb)),nvirtc(isblkkder(4,isb)),   8d22s16
     $      mynprocg,mynowprog,il,ih,i1s,i1e,is,i2e)                    8d22s16
       if(idwsdeb.gt.10)then
       write(6,*)('what we are retrieving for k ')
       call prntm2(bc(kmatd(isb)),nrow,ih+1-il,nrow)
       end if
       do i12=il,ih                                                     8d22s16
        i12m=i12-1                                                      8d22s16
        i12k=i12-il                                                     8d22s16
        iad1=jmatd(isb)+i12k*nrow                                       8d22s16
        iad2=jmatda+i12m*nrow                                           8d22s16
        do i34=0,nrow-1                                                 8d22s16
         bc(iad1+i34)=bc(iad1+i34)+bc(iad2+i34)                         9d16s16
        end do                                                          8d22s16
        iad1=kmatd(isb)+i12k*nrow                                       8d22s16
        iad2=kmatda+i12m*nrow                                           8d22s16
        do i34=0,nrow-1                                                 8d22s16
         bc(iad1+i34)=bc(iad1+i34)+bc(iad2+i34)                         9d16s16
        end do                                                          8d22s16
       end do
       ncolj=nvirtc(isblkkder(3,isb))*nvirtc(isblkkder(4,isb))
       nrowj=nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb))
       if(idwsdeb.gt.10)then
        write(6,*)('final ao der contribution to k'),bc(igoal)
        call prntm2(bc(kmatda),nrowj,ncolj,nrowj)
        write(6,*)('and to j ')
        call prntm2(bc(jmatda),nrowj,ncolj,nrowj)
        write(6,*)('final kmat der '),ntot,il,ih,jmatda,kmatda,nall
        do i=0,ntot-1                                                   10d31s16
         bc(jmatda+i)=0d0                                               10d31s16
        end do                                                          10d31s16
        nhere=ih+1-il                                                   10d31s16
        do i=0,nhere-1                                                  10d31s16
         ip=il+i-1                                                      10d31s16
         do j=0,nrow-1                                                  10d31s16
          iad1=jmatda+j+nrow*ip                                         10d31s16
          iad2=jmatd(isb)+j+nrow*i                                      10d31s16
          bc(iad1)=bc(iad2)                                             10d31s16
          iad1=kmatda+j+nrow*ip                                         10d31s16
          iad2=kmatd(isb)+j+nrow*i                                      10d31s16
          bc(iad1)=bc(iad2)                                             10d31s16
         end do                                                         10d31s16
        end do                                                          10d31s16
        write(6,*)('kmatda before global summing ')
        call prntm2(bc(kmatda),nrow,ncol,nrow)
        call dws_gsumf(bc(jmatda),ntot)                                 10d31s16
        write(6,12)(isblkkder(j,isb),j=1,4)
        ncolk=ih+1-il
        write(6,*)('distributed '),kmatd(isb),kmatd(isb)+201*16,
     $       igoal,bc(igoal),bc(kmatd(isb)+201*16)
        call prntm2(bc(kmatd(isb)),nrow,nhere,nrow)
        write(6,*)('full ')
        call prntm2(bc(kmatda),nrow,ncol,nrow)
       end if
       if(idwsdeb.gt.10)then
        write(6,*)('jmat ders so far ')
        write(6,12)(isblkkder(j,isb),j=1,4)
        write(6,*)('distributed ')
        call prntm2(bc(jmatd(isb)),nrow,nhere,nrow)
        write(6,*)('full ')
        call prntm2(bc(jmatda),nrow,ncol,nrow)
       end if
       ibcoff=jmatda                                                    8d22s16
      end do                                                            9d16s16
      do isb=1,nsblkder1                                                9d2s16
       nrow=nocc(isblkder1 (1,isb))*nocc(isblkder1 (2,isb))
       ncol=nocc(isblkder1 (3,isb))*nocc(isblkder1 (4,isb))
       ntot=nrow*ncol
       if(idwsdeb.gt.10)then
       write(6,*)('gsummed starting mixed 4od2b')
       write(6,12)(isblkder1(j,isb),j=1,4)
       itmp=ibcoff
       ibcoff=itmp+ntot
       call enough('parajkfromhd1.  7',bc,ibc)
       do i=0,ntot-1
        bc(itmp+i)=bc(i4od2b(isb)+i)
       end do
       call dws_gsumf(bc(itmp),ntot)
       call prntm2(bc(itmp),nrow,ncol,nrow)
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkder1 (1,isb)),0,
     $        nocc(isblkder1 (2,isb)),0,nocc(isblkder1 (3,isb)),0,
     $        nocc(isblkder1 (4,isb)),0,
     $        bc(ibcoff))
        end if
       ibcoff=itmp
       end if
       isw1=multh(isblkder1(1,isb),iprop)
       isw2=multh(isblkder1(2,isb),iprop)
       isw3=multh(isblkder1(3,isb),iprop)
       isw4=multh(isblkder1(4,isb),iprop)
c
c     (o'`o|oo): (a'o|oo) part ok
c     recall we have (v'v|ov)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkder1(3,isb),isblkder1(4,isb),isblkder1(2,isb))
       itry2=iptoh(isblkder1(4,isb),isblkder1(3,isb),isblkder1(2,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isblkder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw1)*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(4,isb)))then
           do ib=0,nocc(isblkder1(2,isb))-1
            do ia=0,nbasdwsc(isw1)-1
             iadh=ii0+ia+nbasdwsc(isw1)*ib
             ff=bc(iadh)*2d0
             do m1=0,nocc(isblkder1(1,isb))-1
              iad1=itrans(isw1)+ia+nbasdwsc(isw1)*m1
              iadk=i4od2b(isb)+m1+nocc(isblkder1(1,isb))*(ib
     $        +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*icm))
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
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isblkder1(3,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw1)*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(3,isb)))then
           do ib=0,nocc(isblkder1(2,isb))-1
            do ia=0,nbasdwsc(isw1)-1
             iadh=ii0+ia+nbasdwsc(isw1)*ib
             ff=bc(iadh)
             do m1=0,nocc(isblkder1(1,isb))-1
              iad1=itrans(isw1)+ia+nbasdwsc(isw1)*m1
              iadk=i4od2b(isb)+m1+nocc(isblkder1(1,isb))*(ib
     $        +nocc(isblkder1(2,isb))*(icm+nocc(isblkder1(3,isb))*idm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (a''o|oo) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o'o`|oo): (o'a|oo) part ok
c
       itry=iptoh(isblkder1(3,isb),isblkder1(4,isb),isw2)
       itry2=iptoh(isblkder1(4,isb),isblkder1(3,isb),isw2)
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isblkder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw2)*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(4,isb)))then
           do ib=0,nbasdwsc(isw2)-1
            do ia=0,nocc(isblkder1(1,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(1,isb))*ib
             ff=bc(iadh)*2d0
             do m2=0,nocc(isblkder1(2,isb))-1
              iad1=itrans(isw2)+ib+nbasdwsc(isw2)*m2
              iadk=i4od2b(isb)+ia+nocc(isblkder1(1,isb))*(m2
     $        +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*icm))
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
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isblkder1(3,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw2)*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(3,isb)))then
           do ib=0,nbasdwsc(isw2)-1
            do ia=0,nocc(isblkder1(1,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(1,isb))*ib
             ff=bc(iadh)*2d0
             do m2=0,nocc(isblkder1(2,isb))-1
              iad1=itrans(isw2)+ib+nbasdwsc(isw2)*m2
              iadk=i4od2b(isb)+ia+nocc(isblkder1(1,isb))*(m2
     $        +nocc(isblkder1(2,isb))*(icm+nocc(isblkder1(3,isb))*idm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''a|oo) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o'o|o`o): (o'o|ao) part ok
c
       itry=iptoh(isblkder1(4,isb),isw3,isblkder1(2,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isw3),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(2,isb))*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(2,isb))-1
           do ia=0,nocc(isblkder1(1,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(1,isb))*ib
            ff=bc(iadh)*2d0
            do m3=0,nocc(isblkder1(3,isb))-1
             iad1=itrans(isw3)+icm+nbasdwsc(isw3)*m3
             iadk=i4od2b(isb)+ia+nocc(isblkder1(1,isb))*(ib
     $        +nocc(isblkder1(2,isb))*(m3+nocc(isblkder1(3,isb))*idm))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|ao) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o'o|oo`): (o'o|oa) part ok
c
       itry=iptoh(isblkder1(3,isb),isw4,isblkder1(2,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isw4),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(2,isb))*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(2,isb))-1
           do ia=0,nocc(isblkder1(1,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(1,isb))*ib
            ff=bc(iadh)*2d0
            do m4=0,nocc(isblkder1(4,isb))-1
             iad1=itrans(isw4)+icm+nbasdwsc(isw4)*m4
             iadk=i4od2b(isb)+ia+nocc(isblkder1(1,isb))*(ib
     $        +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*m4))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|oa) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o`o'|oo): (ao'|oo) part ok
c
       itry=iptoh(isblkder1(3,isb),isblkder1(4,isb),isw1)
       itry2=iptoh(isblkder1(4,isb),isblkder1(3,isb),isw1)
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isblkder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw1)*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(4,isb)))then
           do ib=0,nbasdwsc(isw1)-1
            do ia=0,nocc(isblkder1(2,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(2,isb))*ib
             ff=bc(iadh)*2d0
             do m1=0,nocc(isblkder1(1,isb))-1
              iad1=itrans(isw1)+ib+nbasdwsc(isw1)*m1
              iadk=i4od2b(isb)+m1+nocc(isblkder1(1,isb))*(ia
     $        +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*icm))
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
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isblkder1(3,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw1)*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(3,isb)))then
           do ib=0,nbasdwsc(isw1)-1
            do ia=0,nocc(isblkder1(2,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(2,isb))*ib
             ff=bc(iadh)*2d0
             do m1=0,nocc(isblkder1(1,isb))-1
              iad1=itrans(isw1)+ib+nbasdwsc(isw1)*m1
              iadk=i4od2b(isb)+m1+nocc(isblkder1(1,isb))*(ia
     $        +nocc(isblkder1(2,isb))*(icm+nocc(isblkder1(3,isb))*idm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (ao''|oo) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'`|oo): (oa'|oo) part ok
c
       itry=iptoh(isblkder1(3,isb),isblkder1(4,isb),isblkder1(1,isb))
       itry2=iptoh(isblkder1(4,isb),isblkder1(3,isb),isblkder1(1,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isblkder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw2)*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(4,isb)))then
           do ib=0,nocc(isblkder1(1,isb))-1
            do ia=0,nbasdwsc(isw2)-1
             iadh=ii0+ia+nbasdwsc(isw2)*ib
             ff=bc(iadh)*2d0
             do m2=0,nocc(isblkder1(2,isb))-1
              iad1=itrans(isw2)+ia+nbasdwsc(isw2)*m2
              iadk=i4od2b(isb)+ib+nocc(isblkder1(1,isb))*(m2
     $        +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*icm))
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
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isblkder1(3,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw2)*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(3,isb)))then
           do ib=0,nocc(isblkder1(1,isb))-1
            do ia=0,nocc(isw2)-1
             iadh=ii0+ia+nbasdwsc(isw2)*ib
             ff=bc(iadh)*2d0
             do m2=0,nocc(isblkder1(2,isb))-1
              iad1=itrans(isw2)+ia+nbasdwsc(isw2)*m2
              iadk=i4od2b(isb)+ib+nocc(isblkder1(1,isb))*(m2
     $        +nocc(isblkder1(2,isb))*(icm+nocc(isblkder1(3,isb))*idm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oa''|oo) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|o`o): (oo'|ao) part ok
c
       itry=iptoh(isblkder1(4,isb),isw3,isblkder1(1,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isw3),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(2,isb))*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(1,isb))-1
           do ia=0,nocc(isblkder1(2,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(2,isb))*ib
            ff=bc(iadh)*2d0
            do m3=0,nocc(isblkder1(3,isb))-1
             iad1=itrans(isw3)+icm+nbasdwsc(isw3)*m3
             iadk=i4od2b(isb)+ib+nocc(isblkder1(1,isb))*(ia
     $        +nocc(isblkder1(2,isb))*(m3+nocc(isblkder1(3,isb))*idm))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|ao) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|oo`): (oo'|oa) part ok
c
       itry=iptoh(isblkder1(3,isb),isw4,isblkder1(1,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isw4),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(2,isb))*nbasdwsc(isblkder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(1,isb))-1
           do ia=0,nocc(isblkder1(2,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(2,isb))*ib
            ff=bc(iadh)*2d0
            do m4=0,nocc(isblkder1(4,isb))-1
             iad1=itrans(isw4)+icm+nbasdwsc(isw4)*m4
             iadk=i4od2b(isb)+ib+nocc(isblkder1(1,isb))*(ia
     $        +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*m4))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|oa) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o`o|o'o): (ao|o'o) part ok
c
       itry=iptoh(isblkder1(2,isb),isw1,isblkder1(4,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isw1),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(2,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isblkder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(4,isb))-1
           do ia=0,nocc(isblkder1(3,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(3,isb))*ib
            ff=bc(iadh)*2d0
            do m1=0,nocc(isblkder1(1,isb))-1
             iad1=itrans(isw1)+icm+nbasdwsc(isw1)*m1
             iadk=i4od2b(isb)+m1+nocc(isblkder1(1,isb))*(idm
     $        +nocc(isblkder1(2,isb))*(ia+nocc(isblkder1(3,isb))*ib))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (ao|o''o) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo`|o'o): (oa|o'o) part ok
c
       itry=iptoh(isblkder1(1,isb),isw2,isblkder1(4,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isw2),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isblkder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(4,isb))-1
           do ia=0,nocc(isblkder1(3,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(3,isb))*ib
            ff=bc(iadh)*2d0
            do m2=0,nocc(isblkder1(2,isb))-1
             iad1=itrans(isw2)+icm+nbasdwsc(isw2)*m2
             iadk=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(m2
     $        +nocc(isblkder1(2,isb))*(ia+nocc(isblkder1(3,isb))*ib))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oa|o''o) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|o'`o): (oo|a'o) part ok
c
       itry=iptoh(isblkder1(2,isb),isblkder1(1,isb),isblkder1(4,isb))
       itry2=iptoh(isblkder1(1,isb),isblkder1(2,isb),isblkder1(4,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isblkder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(2,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw3)*nbasdwsc(isblkder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkder1(1,isb)))then
           idm=id-1
           do ib=0,nocc(isblkder1(4,isb))-1
            do ia=0,nbasdwsc(isw3)-1
             iadh=ii0+ia+nbasdwsc(isw3)*ib
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(icm
     $         +nocc(isblkder1(2,isb))*(m3+nocc(isblkder1(3,isb))*ib))
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
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isblkder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(1,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw3)*nbasdwsc(isblkder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkder1(2,isb)))then
           idm=id-1
           do ib=0,nocc(isblkder1(4,isb))-1
            do ia=0,nbasdwsc(isw3)-1
             iadh=ii0+ia+nbasdwsc(isw3)*ib
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=i4od2b(isb)+icm+nocc(isblkder1(1,isb))*(idm
     $         +nocc(isblkder1(2,isb))*(m3+nocc(isblkder1(3,isb))*ib))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo|a''o) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|o'o`): (oo|o'a) part ok
c
       itry=iptoh(isblkder1(1,isb),isblkder1(2,isb),isw4)
       itry2=iptoh(isblkder1(2,isb),isblkder1(1,isb),isw4)
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isblkder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(2,isb)))then
           do ib=0,nbasdwsc(isw4)-1
            do ia=0,nocc(isblkder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(3,isb))*ib
             ff=bc(iadh)*2d0
             do m4=0,nocc(isblkder1(4,isb))-1
              iad1=itrans(isw4)+ib+nbasdwsc(isw4)*m4
              iadk=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(icm
     $          +nocc(isblkder1(2,isb))*(ia+nocc(isblkder1(3,isb))*m4))
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
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isblkder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(1,isb)))then
           do ib=0,nbasdwsc(isw4)-1
            do ia=0,nocc(isblkder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(3,isb))*ib
             ff=bc(iadh)*2d0
             do m4=0,nocc(isblkder1(4,isb))-1
              iad1=itrans(isw4)+ib+nbasdwsc(isw4)*m4
              iadk=i4od2b(isb)+icm+nocc(isblkder1(1,isb))*(idm
     $          +nocc(isblkder1(2,isb))*(ia+nocc(isblkder1(3,isb))*m4))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (ao|o''o) for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o`o|oo'): (ao|oo') part ok
c
       itry=iptoh(isblkder1(2,isb),isw1,isblkder1(3,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isw1),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(2,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isblkder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(3,isb))-1
           do ia=0,nocc(isblkder1(4,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(4,isb))*ib
            ff=bc(iadh)*2d0
            do m1=0,nocc(isblkder1(1,isb))-1
             iad1=itrans(isw1)+icm+nbasdwsc(isw1)*m1
             iadk=i4od2b(isb)+m1+nocc(isblkder1(1,isb))*(idm
     $        +nocc(isblkder1(2,isb))*(ib+nocc(isblkder1(3,isb))*ia))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (ao|oo'') for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo`|oo'): (oa|oo') part ok
c
       itry=iptoh(isblkder1(1,isb),isw2,isblkder1(3,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isw2),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isblkder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkder1(3,isb))-1
           do ia=0,nocc(isblkder1(4,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkder1(4,isb))*ib
            ff=bc(iadh)*2d0
            do m2=0,nocc(isblkder1(2,isb))-1
             iad1=itrans(isw2)+icm+nbasdwsc(isw2)*m2
             iadk=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(m2
     $        +nocc(isblkder1(2,isb))*(ib+nocc(isblkder1(3,isb))*ia))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oa|oo'') for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|o`o'): (oo|ao') part ok
c
       itry=iptoh(isblkder1(1,isb),isblkder1(2,isb),isw3)
       itry2=iptoh(isblkder1(2,isb),isblkder1(1,isb),isw3)
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isblkder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(4,isb))*nbasdwsc(isw3)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(2,isb)))then
           do ib=0,nbasdwsc(isw3)-1
            do ia=0,nocc(isblkder1(4,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(4,isb))*ib
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkder1(3,isb))-1
              iad1=itrans(isw3)+ib+nbasdwsc(isw3)*m3
              iadk=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(icm
     $          +nocc(isblkder1(2,isb))*(m3+nocc(isblkder1(3,isb))*ia))
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
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isblkder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkder1(4,isb))*nbasdwsc(isw3)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(1,isb)))then
           do ib=0,nbasdwsc(isw3)-1
            do ia=0,nocc(isblkder1(4,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(4,isb))*ib
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkder1(3,isb))-1
              iad1=itrans(isw3)+ib+nbasdwsc(isw3)*m3
              iadk=i4od2b(isb)+icm+nocc(isblkder1(1,isb))*(idm
     $          +nocc(isblkder1(2,isb))*(m3+nocc(isblkder1(3,isb))*ia))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo|ao'') for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|o'`): (oo|oa') part ok
c
       itry=iptoh(isblkder1(1,isb),isblkder1(2,isb),isblkder1(3,isb))
       itry2=iptoh(isblkder1(2,isb),isblkder1(1,isb),isblkder1(3,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isblkder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(2,isb)))then
           do ib=0,nocc(isblkder1(3,isb))-1
            do ia=0,nbasdwsc(isw4)-1
             iadh=ii0+ia+nbasdwsc(isw4)*ib
             ff=bc(iadh)*2d0
             do m4=0,nocc(isblkder1(4,isb))-1
              iad1=itrans(isw4)+ia+nbasdwsc(isw4)*m4
              iadk=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(icm
     $          +nocc(isblkder1(2,isb))*(ib+nocc(isblkder1(3,isb))*m4))
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
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isblkder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(1,isb)))then
           do ib=0,nocc(isblkder1(3,isb))-1
            do ia=0,nbasdwsc(isw4)-1
             iadh=ii0+ia+nbasdwsc(isw4)*ib
             ff=bc(iadh)*2d0
             do m4=0,nocc(isblkder1(4,isb))-1
              iad1=itrans(isw4)+ia+nbasdwsc(isw4)*m4
              iadk=i4od2b(isb)+icm+nocc(isblkder1(1,isb))*(idm
     $          +nocc(isblkder1(2,isb))*(ib+nocc(isblkder1(3,isb))*m4))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo|ao'') for i4od2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
       if(idwsdeb.gt.10)then
        write(6,12)(isblkder1(j,isb),j=1,4)
        write(6,*)('un gsummed 4od2b ')
        call prntm2(bc(i4od2b(isb)),nrow,ncol,nrow)
        itmp=ibcoff
        ibcoff=itmp+nrow*ncol
        call enough('parajkfromhd1.  8',bc,ibc)
        do i=0,nrow*ncol-1
         bc(itmp+i)=bc(i4od2b(isb)+i)
        end do
        write(6,*)('gsummed ')
        call dws_gsumf(bc(itmp),nrow*ncol)
        call prntm2(bc(itmp),nrow,ncol,nrow)
        if(nsymb.eq.1)call printa(bc(itmp),nocc,0,nocc,0,nocc,0,nocc,0,
     $       bc(ibcoff))
       ibcoff=itmp
       end if
      end do
      do isb=1,nsblkxder1                                                9d2s16
       nrow=nocc(isblkxder1(1,isb))*nocc(isblkxder1(2,isb))
       ncol=nocc(isblkxder1(3,isb))*nvirtc(isblkxder1(4,isb))
       ntot=nrow*ncol
       if(idwsdeb.gt.10)then
       write(6,*)('starting onexd2 '),bc(ionexd2(isb)),ionexd2(isb),
     $      loc(bc(ionexd2(isb)))
       write(6,12)(isblkxder1(j,isb),j=1,4)
       write(6,*)('gsummed ')
       itmp=ibcoff
       ibcoff=itmp+ntot
       call enough('parajkfromhd1.  9',bc,ibc)
       do i=0,ntot-1
        bc(itmp+i)=bc(ionexd2(isb)+i)
       end do
       call dws_gsumf(bc(itmp),ntot)
       call prntm2(bc(itmp),nrow,ncol,nrow)
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkxder1 (1,isb)),0,
     $        nocc(isblkxder1 (2,isb)),0,nocc(isblkxder1 (3,isb)),0,
     $        nvirt(isblkxder1(4,isb)),nocc(isblkxder1 (4,isb)),
     $        bc(ibcoff))
        end if
       ibcoff=itmp
       end if
       isw1=multh(isblkxder1(1,isb),iprop)
       isw2=multh(isblkxder1(2,isb),iprop)
       isw3=multh(isblkxder1(3,isb),iprop)
       isw4=multh(isblkxder1(4,isb),iprop)
c
c     (o'`o|ov): (a'o|ov) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(3,isb),isblkxder1(4,isb),isblkxder1(2,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isblkxder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw1)*nbasdwsc(isblkxder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkxder1(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.gt.nocc(isblkxder1(4,isb)))then
           do ib=0,nocc(isblkxder1(2,isb))-1
            do ia=0,nbasdwsc(isw1)-1
             iadh=ii0+ia+nbasdwsc(isw1)*ib
             ff=bc(iadh)*2d0
             do m1=0,nocc(isblkxder1(1,isb))-1
              iad1=itrans(isw1)+ia+nbasdwsc(isw1)*m1
              iadk=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(ib
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*icm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (a''o|ov) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o'o`|ov): (o'a|ov) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(3,isb),isblkxder1(4,isb),isw2)
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isblkxder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw2)*nbasdwsc(isblkxder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkxder1(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.gt.nocc(isblkxder1(4,isb)))then
           do ib=0,nbasdwsc(isw2)-1
            do ia=0,nocc(isblkxder1(1,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(1,isb))*ib
             ff=bc(iadh)*2d0
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+ib+nbasdwsc(isw2)*m2
              iadk=ionexd2(isb)+ia+nocc(isblkxder1(1,isb))*(m2
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*icm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''a|ov) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o'o|o`v): requires (o'o|av) which we dont currently have
c
c     (o'o|ov`): (o'o|oa) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(3,isb),isw4,isblkxder1(2,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isw4),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(2,isb))*nbasdwsc(isblkxder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkxder1(2,isb))-1
           do ia=0,nocc(isblkxder1(1,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkxder1(1,isb))*ib
            ff=bc(iadh)*2d0
            do m4=0,nvirtc(isblkxder1(4,isb))-1
             m4p=m4+nocc(isblkxder1(4,isb))
             iad1=itrans(isw4)+icm+nbasdwsc(isw4)*m4p
             iadk=ionexd2(isb)+ia+nocc(isblkxder1(1,isb))*(ib
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*m4))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|oa) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o`o'|ov): (ao'|ov) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(3,isb),isblkxder1(4,isb),isw1)
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isblkxder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw1)*nbasdwsc(isblkxder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkxder1(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.gt.nocc(isblkxder1(4,isb)))then
           do ib=0,nbasdwsc(isw1)-1
            do ia=0,nocc(isblkxder1(2,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(2,isb))*ib
             ff=bc(iadh)*2d0
             do m1=0,nocc(isblkxder1(1,isb))-1
              iad1=itrans(isw1)+ib+nbasdwsc(isw1)*m1
              iadk=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(ia
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*icm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (ao''|ov) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'`|ov): (oa'|ov) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(3,isb),isblkxder1(4,isb),isblkxder1(1,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isblkxder1(4,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw2)*nbasdwsc(isblkxder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkxder1(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.gt.nocc(isblkxder1(4,isb)))then
           idm=id-1
           do ib=0,nocc(isblkxder1(1,isb))-1
            do ia=0,nbasdwsc(isw2)-1
             iadh=ii0+ia+nbasdwsc(isw2)*ib
             ff=bc(iadh)*2d0
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+ia+nbasdwsc(isw2)*m2
              iadk=ionexd2(isb)+ib+nocc(isblkxder1(1,isb))*(m2
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*icm))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oa''|ov) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|o`v): (oo'|av) part but we dont have this yet
c
c     (oo'|ov`): (oo'|oa) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(3,isb),isw4,isblkxder1(1,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isw4),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(2,isb))*nbasdwsc(isblkxder1(1,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkxder1(1,isb))-1
           do ia=0,nocc(isblkxder1(2,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkxder1(2,isb))*ib
            ff=bc(iadh)*2d0
            do m4=0,nvirtc(isblkxder1(4,isb))-1
             m4p=m4+nocc(isblkxder1(4,isb))
             iad1=itrans(isw4)+icm+nbasdwsc(isw4)*m4p
             iadk=ionexd2(isb)+ib+nocc(isblkxder1(1,isb))*(ia
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*m4))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|oa) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (`oo|o'v): (ao|o'v) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(2,isb),isw1,isblkxder1(4,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isw1),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nvirtc(isblkxder1(4,isb))-1
           ibp=ib+nocc(isblkxder1(4,isb))
           do ia=0,nocc(isblkxder1(3,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkxder1(3,isb))*ibp
            ff=bc(iadh)*2d0
            do m1=0,nocc(isblkxder1(1,isb))-1
             iad1=itrans(isw1)+icm+nbasdwsc(isw1)*m1
             iadk=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(idm
     $       +nocc(isblkxder1(2,isb))*(ia+nocc(isblkxder1(3,isb))*ib))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (ao|o''v) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo`|o'v): (oa|o'v) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(1,isb),isw2,isblkxder1(4,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isw2),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nvirtc(isblkxder1(4,isb))-1
           ibp=ib+nocc(isblkxder1(4,isb))
           do ia=0,nocc(isblkxder1(3,isb))-1
            iadh=ii0+ia+nbasdwsc(isblkxder1(3,isb))*ibp
            ff=bc(iadh)*2d0
            do m2=0,nocc(isblkxder1(2,isb))-1
             iad1=itrans(isw2)+icm+nbasdwsc(isw2)*m2
             iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(m2
     $       +nocc(isblkxder1(2,isb))*(ia+nocc(isblkxder1(3,isb))*ib))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oa|o''v) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|o'`v): (oo|a'v) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isblkxder1(4,isb))
       itry2=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),
     $      isblkxder1(4,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isblkxder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isw3)*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(2,isb)))then
           idm=id-1
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nbasdwsc(isw3)-1
             iadh=ii0+ia+nbasdwsc(isw3)*ibp
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $         +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ib))
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
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isblkxder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isw3)*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(1,isb)))then
           idm=id-1
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nbasdwsc(isw3)-1
             iadh=ii0+ia+nbasdwsc(isw3)*ibp
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ia+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $         +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ib))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo|a''v) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|o'v`): (oo|o'a) part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isw4)
       itry2=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),isw4)
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isblkxder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(2,isb)))then
           idm=id-1
           do ib=0,nbasdwsc(isw4)-1                                     12d2s16
            do ia=0,nocc(isblkxder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(3,isb))*ib                 12d2s16
             ff=bc(iadh)*2d0
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+ib+nbasdwsc(isw4)*m4p                   12d2s16
              iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $         +nocc(isblkxder1(2,isb))*(ia+nocc(isblkxder1(3,isb))*m4))
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
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isblkxder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(1,isb)))then
           idm=id-1
           do ib=0,nbasdwsc(isw4)-1                                     12d2s16
            do ia=0,nocc(isblkxder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(3,isb))*ib                 12d2s16
             ff=bc(iadh)*2d0
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+ib+nbasdwsc(isw4)*m4p                   12d2s16
              iadk=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $         +nocc(isblkxder1(2,isb))*(ia+nocc(isblkxder1(3,isb))*m4))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo|o''a) for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (`oo|ov'): (ao|ov') part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(2,isb),isw1,isblkxder1(3,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isw1),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkxder1(3,isb))-1
           do ia=0,nvirtc(isblkxder1(4,isb))-1
            iap=ia+nocc(isblkxder1(4,isb))
            iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
            ff=bc(iadh)*2d0
            do m1=0,nocc(isblkxder1(1,isb))-1
             iad1=itrans(isw1)+icm+nbasdwsc(isw1)*m1
             iadk=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(idm
     $       +nocc(isblkxder1(2,isb))*(ib+nocc(isblkxder1(3,isb))*ia))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (ao|ov'') for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o`o|ov'): (oa|ov') part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(1,isb),isw2,isblkxder1(3,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isw2),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isblkxder1(4,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          do ib=0,nocc(isblkxder1(3,isb))-1
           do ia=0,nvirtc(isblkxder1(4,isb))-1
            iap=ia+nocc(isblkxder1(4,isb))
            iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
            ff=bc(iadh)*2d0
            do m2=0,nocc(isblkxder1(2,isb))-1
             iad1=itrans(isw2)+icm+nbasdwsc(isw2)*m2
             iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(m2
     $       +nocc(isblkxder1(2,isb))*(ib+nocc(isblkxder1(3,isb))*ia))
             bc(iadk)=bc(iadk)+ff*bc(iad1)
            end do
           end do
          end do
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oa|ov'') for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|o`v'): (oo|av') part ok
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isw3)
       itry2=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),isw3)
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isblkxder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(4,isb))*nbasdwsc(isw3)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(2,isb)))then
           idm=id-1
           do ib=0,nbasdwsc(isw3)-1
            do ia=0,nvirtc(isblkxder1(4,isb))-1
             iap=ia+nocc(isblkxder1(4,isb))
             iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ib+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $         +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ia))
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
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isblkxder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkxder1(4,isb))*nbasdwsc(isw3)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(1,isb)))then
           idm=id-1
           do ib=0,nbasdwsc(isw3)-1
            do ia=0,nvirtc(isblkxder1(4,isb))-1
             iap=ia+nocc(isblkxder1(4,isb))
             iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
             ff=bc(iadh)*2d0
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+ib+nbasdwsc(isw3)*m3
              iadk=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $         +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*ia))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo|av'') for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo|ov'`): (oo|oa') part
c     recall we have (a'a|oa)
c                     a b dc  stored under iptoh(d,c,b),
c
       itry=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isblkxder1(3,isb))
       itry2=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),
     $      isblkxder1(3,isb))
       if(itry.gt.0)then
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isblkxder1(2,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(1,isb))
        ii0=ihcol(itry)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(2,isb)))then
           idm=id-1
           do ib=0,nocc(isblkxder1(3,isb))-1
            do ia=0,nbasdwsc(isw4)-1
             iadh=ii0+ia+nbasdwsc(isw4)*ib
             ff=bc(iadh)*2d0
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+ia+nbasdwsc(isw4)*m4p
              iadk=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $         +nocc(isblkxder1(2,isb))*(ib+nocc(isblkxder1(3,isb))*m4))
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
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isblkxder1(1,isb)),
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(2,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkxder1(3,isb))*nbasdwsc(isw4)
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          if(ic.le.nocc(isblkxder1(1,isb)))then
           idm=id-1
           do ib=0,nocc(isblkxder1(3,isb))-1
            do ia=0,nbasdwsc(isw4)-1
             iadh=ii0+ia+nbasdwsc(isw4)*ib
             ff=bc(iadh)*2d0
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+ia+nbasdwsc(isw4)*m4p
              iadk=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $         +nocc(isblkxder1(2,isb))*(ib+nocc(isblkxder1(3,isb))*m4))
              bc(iadk)=bc(iadk)+ff*bc(iad1)
             end do
            end do
           end do
          end if
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo|oa'') for ionexd2 ')
        call dws_sync
        call dws_finalize
        stop
       end if
       if(idwsdeb.gt.10.and.min(nrow,ncol).gt.0)then
        write(6,*)('finishing onexd2 ooox '),bc(ionexd2(isb)),
     $       ionexd2(isb),loc(ionexd2(isb))
       write(6,12)(isblkxder1(j,isb),j=1,4)
       write(6,*)('unglobal summed ')
       call prntm2(bc(ionexd2(isb)),nrow,ncol,nrow)
       if(nsymb.eq.1)call printa(bc(ionexd2(1)),nocc,0,nocc,0,
     $      nocc,0,nvirt,nocc,bc(ibcoff))
       itmp=ibcoff
       ibcoff=itmp+nrow*ncol
       call enough('parajkfromhd1. 10',bc,ibc)
       do i=0,nrow*ncol-1
        bc(itmp+i)=bc(ionexd2(isb)+i)
       end do
       call dws_gsumf(bc(itmp),nrow*ncol)
       ibcoff=itmp
        write(6,*)('finishing onex '),bc(ionexd2(isb)),
     $       ionexd2(isb),loc(ionexd2(isb))
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkxder1 (1,isb)),0,
     $        nocc(isblkxder1 (2,isb)),0,nocc(isblkxder1 (3,isb)),0,
     $        nvirt(isblkxder1(4,isb)),nocc(isblkxder1 (4,isb)),
     $        bc(ibcoff))
        end if
       end if
      end do
      return
      end
