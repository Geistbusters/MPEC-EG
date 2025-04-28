c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parajkfromhd2c(ipair,ngaus,ibdat,nhcolt,ihmat,nocc,icol9d20s16
     $     ,iorbx,iapair,isstor,ibstor,multh,iptoh,imsg,                9d28s16
     $     iter,idwsdeb,ncomp,nvirtc,i4od2b,ionexd2,                    9d28s16
     $     der2e,nbasdwsc,                                              9d28s16
     $     isblkder,isblkxder,nsblkder,nsblkxder,isblkxder1,nsblkxder1, 9d16s16
     $     isblkder1,nsblkder1,isblkkder,nsblkkder,i4od,nbasisp,xgoal,  11d10s22
     $     bc,ibc,idorel)                                               3d7s23
      implicit real*8 (a-h,o-z)
c
c     like jkfromhd1, working with (ab'|d'c) instead of (a'b|dc)
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
      dimension multh(8,8),iorbn(8),itrans(1),iorbx(8),                 9d16s16
     $     ioooo(1),isblkder(4,idbk),isblkxder(4,idbk),                 8d24s16
     $     ionex(idbk),myh(idbk),nvirtc(8),isblkkder(4,idbk),           9d28s16
     $     i4od2b(idbk),nbasisp(*),                                     4d7s22
     $     ionexd2(idbk),isblkxder1(4,idbk),isblkder1(4,idbk)           9d16s16
      call second(time1)                                                11d27s12
      if(iter.eq.1.and.mynowprog.eq.0)write(6,*)('in parajkfromhd2c ')   3d3s17
c
c     in this version, do transpose of hmat, then transform full hmat
c     to mo basis, then store appropriate parts of the result.
c
c
      ntmpmat=2                                                         5d31s22
      if(idorel.eq.0.or.idorel.eq.1)then                                5d31s22
       ntmpmat=1                                                        5d31s22
      end if                                                            5d31s22
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
      nn=ngaus*ngaus
      jbdat=ibdat-1
      if(mynprocg.gt.idp)then                                           11d30s12
       write(6,*)('in parajkfromhd2c, mynprocg = '),mynprocg
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
          isa=multh(isd,isbc)                                           9d28s16
          ii=iptoh(isd,isc,isb)                                         11d30s12
          ncd=nbasdwsc(isc)*nocc(isd)                                   9d28s16
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
        bx=0d0
        do i=0,msend-1
         rms=rms+bc(ibufs+i)**2
         if(abs(bc(ibufs+i)).gt.bx)then
          bx=abs(bc(ibufs+i))
          ibx=i
         end if
        end do
        rms=sqrt(rms/dfloat(msend))                                     3d3s23
        if(rms.gt.1d10)then
         write(6,*)('rms size of integrals in parajkfromhd2c is '),
     $        rms
         write(6,*)('while really doesn''t look right!!')
         call dws_sync
         call dws_finalize
         stop 'parajkfromhd2c'
        end if
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
        do isb=1,nsymb                                                  2d11s13
         nbasdwsb=nbasisp(isb)                                          3d7s23
         do isc=1,nsymb                                                 9d27s16
          isbc=multh(isc,isb)                                           9d27s16
          do isd=1,nsymb                                                2d11s13
           isa=multh(isd,isbc)                                          9d28s16
           nbasdwsa=nbasisp(isa)
           ntotrun=0
           ncd=nbasdwsc(isc)*nocc(isd)                                  9d28s16
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
 3353          format(3i8,5x,3i8,5x,i8)
 5653          format(4i8)
                if(issb.eq.isb.and.issa.eq.isa)then                      9d20s16
                 if(ipass.eq.1)then                                      2d11s13
                 irecv(ip+1)=irecv(ip+1)+ntmpmat*nhere                  3d7s23
                  myh(iib)=myh(iib)+ntmpmat*nhere                       3d7s23
                  nsumx=nsumx+1
                 else                                                    2d11s13
                  ix=iaaa-1+nbasdwsa*(ibbb-1)                             3d1s16
                  nxn=nbasdwsb*nbasdwsa                                      3d1s16
                  do in=0,ntmpmat-1                                     3d7s23
                   do j=0,nhere-1                                         3d2s16
                    iad=myh(iib)+ix+nxn*(j+nhere*in)                    3d7s23
                    bc(iad)=bc(ipr(ip+1))                                 2d11s13
                    ntotrun=ntotrun+1
                    ipr(ip+1)=ipr(ip+1)+1                                 2d11s13
                   end do                                               3d7s23
                  end do                                                 2d11s13
                 end if
                end if                                                  2d11s13
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
        if(iter.eq.1.and.mynowprog.eq.0)then                             3d3s17
         write(6,*)('no. of words to receive:  '),nrecv                   11d30s12
         write(6,*)('no. of words under ihmat: '),ngot,ibcoff                   11d30s12
        end if
        if(nrecv.gt.ngot)then                                           11d30s12
         more=nrecv-ngot                                                11d30s12
         if(iter.eq.1.and.mynowprog.eq.0)then                            3d3s17
         write(6,*)('need to make some space for the receive ')         11d30s12
         write(6,*)('need '),more,(' additional words ')                11d30s12
         end if                                                          3d3s17
         ibcoff=ihmat+nrecv                                             12d3s12
         call enough('parajkfromhd2c.  1',bc,ibc)
        end if                                                          11d30s12
        ibufs=ibcoff                                                    11d28s12
        ibcoff=ibufs+max(nsend,nrecv)                                   11d30s12
        call enough('parajkfromhd2c.  2',bc,ibc)
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
       if(iter.eq.1.and.ipass.eq.1.and.mynowprog.eq.0)then
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
       nbasdwsb=nbasisp(isb)
       do isc=1,nsymb                                                   11d30s12
        isbc=multh(isc,isb)                                             3d4s16
        do isd=1,nsymb                                                  11d30s12
         isa=multh(isd,isbc)                                            9d28s16
         nbasdwsa=nbasisp(isa)
         nxn=nbasdwsa*nbasdwsb                                          3d1s16
         nxnc=nbasdws(isa)*nbasdws(isb)                                 4d12s22
         iis=iptoh(isd,isc,isb)                                         12d21s12
         ncd=nbasdwsc(isc)*nocc(isd)                                    9d28s16
         if(ncd.gt.0.and.iis.gt.0)then                                  12d21s12
          nhere=ncd/mynprocg                                            11d30s12
          nleft=ncd-nhere*mynprocg                                      11d30s12
          if(mynowprog.lt.nleft)nhere=nhere+1                           11d30s12
          if(nhere.gt.0)then                                            11d30s12
            itmp=ibcoff                                                  11d30s12
            itmp2=itmp+nxn                                               3d2s16
            ibcoff=itmp2+nxn                                             3d2s16
            call enough('parajkfromhd2c.  3',bc,ibc)
            ihcol(iis)=ihlast                                            3d5s13
            ihlast=ihlast+nxnc*nhere                                    4d12s22
            if(ihlast.gt.ibufs)then                                      3d5s13
             write(6,*)('ihlast now exceeds ibufs!!! '),ihlast,ibufs     3d5s13
             call dws_sync                                               3d5s13
             call dws_finalize                                           3d5s13
             stop                                                        3d5s13
            end if                                                       3d5s13
            do icol=1,nhere                                              11d30s12
             iads=ihcol(iis)+nxnc*(icol-1)                              4d12s22
             factr=0d0                                                   10d30s20
             jorbb=iorbx(isb)                                            10d30s20
             jorba=iorbx(isa)                                            10d30s20
             do in=0,ntmpmat-1                                          3d7s23
              iad=myh(iis)+nxn*(icol-1+nhere*in)                        3d7s23
              rms=0d0
              do i=0,nbasdwsa*nbasdwsb-1
               rms=rms+bc(iad+i)**2
              end do
              if(rms.ne.rms)then
               write(6,*)('input is nan for col '),icol
               call dws_sync
               call dws_finalize
               stop
              end if
             call dgemm('n','n',nbasdwsa,nbasdwsc(isb),nbasdwsb,        4d7s22
     $            1d0,bc(iad),nbasdwsa,bc(jorbb),nbasdwsb*ncomp,        3d7s23
     $            0d0,bc(itmp),nbasdwsa,                                8d31s15
     d' parajkfromhd2c.  1')
             rms=0d0
             do i=1,nbasdwsa                                             8d31s15
              do j=1,nbasdwsc(isb)                                      4d7s22
               iad1=iad+j-1+nbasdwsc(isb)*(i-1)                         4d7s22
               iad2=itmp+i-1+nbasdwsa*(j-1)                              8d31s15
               bc(iad1)=bc(iad2)                                         11d30s12
               rms=rms+bc(iad2)**2
              end do                                                     11d30s12
             end do                                                      11d30s12
             if(rms.ne.rms)then
              write(6,*)('half transformed is nan for col '),icol
              call dws_sync
              call dws_finalize
              stop
             end if
             call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isa),nbasdwsa,   4d7s22
     $            2d0,bc(iad),nbasdwsc(isb),bc(jorba),nbasdwsa*ncomp,   3d7s23
     $            factr,bc(itmp2),nbasdwsc(isb),                        3d7s23
     d' parajkfromhd2c.  2')
              factr=1d0                                                  10d30s20
              jorbb=jorbb+nbasdwsb                                       10d30s20
              jorba=jorba+nbasdwsa                                       10d30s20
             end do                                                       10d30s20
             rms=0d0
             do i=1,nbasdwsc(isb)                                       4d7s22
              do j=1,nbasdwsc(isa)                                      4d7s22
               iad1=iads+j-1+nbasdwsc(isa)*(i-1)                        4d7s22
               iad2=itmp2+i-1+nbasdwsc(isb)*(j-1)                       3d7s23
               bc(iad1)=bc(iad2)                                         11d30s12
               rms=rms+bc(iad2)**2
              end do                                                     11d30s12
             end do                                                      11d30s12
             if(rms.ne.rms)then
              write(6,*)('fully transformed is nan for col '),icol
              call dws_sync
              call dws_finalize
              stop
             end if
            end do                                                       11d30s12
            ibcoff=itmp                                                  11d30s12
          end if                                                        11d30s12
         end if                                                         11d30s12
        end do                                                          11d30s12
       end do                                                           11d30s12
      end do                                                            11d30s12
c
c     for (oo|oo) go for (o'o|o'o), (oo'|o'o), (o'o|oo'), and (oo'|oo')
c     for (oo|ov) go for (o'o|o'v), (oo'|o'v), (o'o|ov'), and (oo'|ov')
c     recall we have 2(vv'|o'v)
c                      a b  dc  stored under iptoh(d,c,b),
c                                                                       9d16s16
      do isd=1,nsymb
       do isc=1,nsymb
        iscd=multh(isc,isd)
        do isb=1,nsymb
         isa=multh(isb,iscd)
         if(iptoh(isd,isc,isb).gt.0)then
          nrow=nbasdwsc(isb)*nbasdwsc(isa)
          ncol=nocc(isd)*nbasdwsc(isc)
         end if
        end do
       end do
      end do
      igoal=i4od2b(1)+1
      do isb=1,nsblkder1                                                9d28s16
       nrow=nocc(isblkder1(1,isb))*nocc(isblkder1(2,isb))               9d28s16
       ncol=nocc(isblkder1(3,isb))*nocc(isblkder1(4,isb))               9d28s16
       nall=nrow*ncol                                                   8d4s16
c
c     (o'o|o'o) part ok
c
       itry=iptoh(isblkder1(1,isb),isblkder1(2,isb),isblkder1(3,isb))
       itry2=iptoh(isblkder1(3,isb),isblkder1(4,isb),isblkder1(1,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isblkder1(2,isb)),  9d28s16
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
          if(ic.le.nocc(isblkder1(2,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(3,isb))-1
            do ia=0,nocc(isblkder1(4,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(4,isb))*ib
             iadj=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(icm
     $       +nocc(isblkder1(2,isb))*(ib+nocc(isblkder1(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d28s16
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isblkder1(4,isb)),  9d28s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkder1(1,isb))*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(4,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(1,isb))-1
            do ia=0,nocc(isblkder1(2,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(2,isb))*ib
             iadj=i4od2b(isb)+ib+nocc(isblkder1(1,isb))*(ia
     $       +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*icm))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|o''o) ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|o'o) part ok
c
       itry=iptoh(isblkder1(2,isb),isblkder1(1,isb),isblkder1(3,isb))
       itry2=iptoh(isblkder1(3,isb),isblkder1(4,isb),isblkder1(2,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isblkder1(1,isb)),  9d28s16
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
          if(ic.le.nocc(isblkder1(1,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(3,isb))-1
            do ia=0,nocc(isblkder1(4,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(4,isb))*ib
             iadj=i4od2b(isb)+icm+nocc(isblkder1(1,isb))*(idm
     $       +nocc(isblkder1(2,isb))*(ib+nocc(isblkder1(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d28s16
        call ilimts(nocc(isblkder1(3,isb)),nbasdwsc(isblkder1(4,isb)),  9d28s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(3,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkder1(1,isb))*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(4,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(2,isb))-1
            do ia=0,nocc(isblkder1(1,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(1,isb))*ib
             iadj=i4od2b(isb)+ia+nocc(isblkder1(1,isb))*(ib
     $       +nocc(isblkder1(2,isb))*(idm+nocc(isblkder1(3,isb))*icm))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|o''o) ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o'o|oo') part ok
c     recall we have (vv'|o'v)
c                     a b  dc  stored under iptoh(d,c,b)
c
       itry=iptoh(isblkder1(1,isb),isblkder1(2,isb),isblkder1(4,isb))
       itry2=iptoh(isblkder1(4,isb),isblkder1(3,isb),isblkder1(1,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkder1(1,isb)),nbasdwsc(isblkder1(2,isb)),  9d28s16
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
          if(ic.le.nocc(isblkder1(2,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(4,isb))-1
            do ia=0,nocc(isblkder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(3,isb))*ib
             iadj=i4od2b(isb)+idm+nocc(isblkder1(1,isb))*(icm
     $       +nocc(isblkder1(2,isb))*(ia+nocc(isblkder1(3,isb))*ib))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d28s16
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isblkder1(3,isb)),  9d28s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkder1(1,isb))*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(3,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(1,isb))-1
            do ia=0,nocc(isblkder1(2,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(2,isb))*ib
             iadj=i4od2b(isb)+ib+nocc(isblkder1(1,isb))*(ia
     $       +nocc(isblkder1(2,isb))*(icm+nocc(isblkder1(3,isb))*idm))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|oo'') ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|oo') part ok
c     recall we have (vv'|o'v)
c                     a b  dc  stored under iptoh(d,c,b)
c
       itry=iptoh(isblkder1(2,isb),isblkder1(1,isb),isblkder1(4,isb))
       itry2=iptoh(isblkder1(4,isb),isblkder1(3,isb),isblkder1(2,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkder1(2,isb)),nbasdwsc(isblkder1(1,isb)),  9d28s16
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
          if(ic.le.nocc(isblkder1(1,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(4,isb))-1
            do ia=0,nocc(isblkder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(3,isb))*ib
             iadj=i4od2b(isb)+icm+nocc(isblkder1(1,isb))*(idm
     $       +nocc(isblkder1(2,isb))*(ia+nocc(isblkder1(3,isb))*ib))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d28s16
        call ilimts(nocc(isblkder1(4,isb)),nbasdwsc(isblkder1(3,isb)),  9d28s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkder1(4,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkder1(1,isb))*nbasdwsc(isblkder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.le.nocc(isblkder1(3,isb)))then                          9d28s16
           do ib=0,nocc(isblkder1(2,isb))-1
            do ia=0,nocc(isblkder1(1,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkder1(1,isb))*ib
             iadj=i4od2b(isb)+ia+nocc(isblkder1(1,isb))*(ib
     $       +nocc(isblkder1(2,isb))*(icm+nocc(isblkder1(3,isb))*idm))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|oo'') ')
        call dws_sync
        call dws_finalize
        stop
       end if
       if(idwsdeb.gt.10)then
        write(6,*)('for 4od2 at the end of parajkfromhd2c')
        write(6,12)(isblkder1(j,isb),j=1,4)
        call prntm2(bc(i4od2b(isb)),nrow,ncol,nrow)
        itmp=ibcoff                                                     10d31s16
        ibcoff=itmp+nrow*ncol                                           10d31s16
        do i=0,nrow*ncol-1                                              10d31s16
         bc(itmp+i)=bc(i4od2b(isb)+i)                                   10d31s16
        end do
        call dws_gsumf(bc(itmp),nrow*ncol)
        write(6,*)('global summed ')
        call prntm2(bc(itmp),nrow,ncol,nrow)
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkder1 (1,isb)),0,
     $        nocc(isblkder1 (2,isb)),0,nocc(isblkder1 (3,isb)),0,
     $        nocc(isblkder1 (4,isb)),0,
     $        bc(ibcoff))
        end if
        ibcoff=itmp                                                     10d31s16
       end if
      end do                                                            9d16s16
   12 format('integral type ',4i2,5x,i5)                                3d4s13\
      ilookat=6110931
      do isb=1,nsblkxder1
       nrow=nocc(isblkxder1(1,isb))*nocc(isblkxder1(2,isb))               9d28s16
       ncol=nocc(isblkxder1(3,isb))*nvirtc(isblkxder1(4,isb))           9d28s16
       nall=nrow*ncol                                                   8d4s16
c
c     (o'o|o'v) part ok
c     recall we have 2(vv'|o'v)
c                      a b  dc  stored under iptoh(d,c,b)
c
       itry=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isblkxder1(3,isb))
       itry2=iptoh(isblkxder1(3,isb),isblkxder1(4,isb),
     $      isblkxder1(1,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isblkxder1(2,isb)),  9d28s16
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
          if(ic.le.nocc(isblkxder1(2,isb)))then                          9d28s16
           do ib=0,nocc(isblkxder1(3,isb))-1
            do ia=0,nvirtc(isblkxder1(4,isb))-1
             iap=ia+nocc(isblkxder1(4,isb))
             iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
             iadj=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $       +nocc(isblkxder1(2,isb))*(ib+nocc(isblkxder1(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d28s16
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isblkxder1(4,isb)),  9d28s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkxder1(1,isb))*nbasdwsc(isblkxder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkxder1(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.gt.nocc(isblkxder1(4,isb)))then                          9d28s16
           do ib=0,nocc(isblkxder1(1,isb))-1
            do ia=0,nocc(isblkxder1(2,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(2,isb))*ib
             iadj=ionexd2(isb)+ib+nocc(isblkxder1(1,isb))*(ia
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*icm))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|o''v) ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|o'v) part ok
c     recall we have 2(vv'|o'v)
c                      a b  dc  stored under iptoh(d,c,b)
c
       itry=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),isblkxder1(3,isb))
       itry2=iptoh(isblkxder1(3,isb),isblkxder1(4,isb),
     $      isblkxder1(2,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isblkxder1(1,isb)),  9d28s16
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
          if(ic.le.nocc(isblkxder1(1,isb)))then                          9d28s16
           do ib=0,nocc(isblkxder1(3,isb))-1
            do ia=0,nvirtc(isblkxder1(4,isb))-1
             iap=ia+nocc(isblkxder1(4,isb))
             iadh=ii0+iap+nbasdwsc(isblkxder1(4,isb))*ib
             iadj=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $       +nocc(isblkxder1(2,isb))*(ib+nocc(isblkxder1(3,isb))*ia))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else if(itry2.gt.0)then                                          9d28s16
        call ilimts(nocc(isblkxder1(3,isb)),nbasdwsc(isblkxder1(4,isb)),  9d28s16
     $      mynprocg,mynowprog,ilh,ihh,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nocc(isblkxder1(3,isb))
        ii0=ihcol(itry2)
        nrowh=nbasdwsc(isblkxder1(1,isb))*nbasdwsc(isblkxder1(2,isb))
        do ic=i2s,i2e
         icm=ic-1-nocc(isblkxder1(4,isb))
         if(ic.eq.i2e)i1n=i1e
         do id=i10,i1n
          idm=id-1
          if(ic.gt.nocc(isblkxder1(4,isb)))then                          9d28s16
           do ib=0,nocc(isblkxder1(2,isb))-1
            do ia=0,nocc(isblkxder1(1,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(1,isb))*ib
             iadj=ionexd2(isb)+ia+nocc(isblkxder1(1,isb))*(ib
     $       +nocc(isblkxder1(2,isb))*(idm+nocc(isblkxder1(3,isb))*icm))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (oo''|o''v) ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (o'o|ov') part ok
c     recall we have 2(vv'|o'v)
c                      a b  dc  stored under iptoh(d,c,b)
c
       itry=iptoh(isblkxder1(1,isb),isblkxder1(2,isb),isblkxder1(4,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkxder1(1,isb)),nbasdwsc(isblkxder1(2,isb)),  9d28s16
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
          if(ic.le.nocc(isblkxder1(2,isb)))then                          9d28s16
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nocc(isblkxder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(3,isb))*ibp
             iadj=ionexd2(isb)+idm+nocc(isblkxder1(1,isb))*(icm
     $       +nocc(isblkxder1(2,isb))*(ia+nocc(isblkxder1(3,isb))*ib))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|ov'') ')
        call dws_sync
        call dws_finalize
        stop
       end if
c
c     (oo'|ov') part ok
c     recall we have (vv'|o'v)
c                     a b  dc  stored under iptoh(d,c,b)
c
       itry=iptoh(isblkxder1(2,isb),isblkxder1(1,isb),isblkxder1(4,isb))
       if(itry.gt.0)then                                                9d16s16
        call ilimts(nocc(isblkxder1(2,isb)),nbasdwsc(isblkxder1(1,isb)),  9d28s16
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
          if(ic.le.nocc(isblkxder1(1,isb)))then                          9d28s16
           do ib=0,nvirtc(isblkxder1(4,isb))-1
            ibp=ib+nocc(isblkxder1(4,isb))
            do ia=0,nocc(isblkxder1(3,isb))-1
             iadh=ii0+ia+nbasdwsc(isblkxder1(3,isb))*ibp
             iadj=ionexd2(isb)+icm+nocc(isblkxder1(1,isb))*(idm
     $       +nocc(isblkxder1(2,isb))*(ia+nocc(isblkxder1(3,isb))*ib))
             bc(iadj)=bc(iadj)+bc(iadh)                                 10d31s16
            end do
           end do
          end if                                                        9d28s16
          ii0=ii0+nrowh
         end do
         i10=1
        end do
       else
        write(6,*)('could not find (o''o|ov'') ')
        call dws_sync
        call dws_finalize
        stop
       end if
       if(idwsdeb.gt.10.and.min(nrow,ncol).gt.0)then
        write(6,*)('for onexd2 ooox at the bottom of parajkfromhd2c')
        write(6,12)(isblkxder1(j,isb),j=1,4)
        call prntm2(bc(ionexd2(isb)),nrow,ncol,nrow)
        itmp=ibcoff
        ibcoff=itmp+nrow*ncol
        do i=0,nrow*ncol-1
         bc(itmp+i)=bc(ionexd2(isb)+i)
        end do
        call dws_gsumf(bc(itmp),nrow*ncol)
        write(6,*)('global summed ')
        call prntm2(bc(itmp),nrow,ncol,nrow)
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkxder1 (1,isb)),0,
     $        nocc(isblkxder1 (2,isb)),0,nocc(isblkxder1 (3,isb)),0,
     $        nvirt(isblkxder1(4,isb)),nocc(isblkxder1 (4,isb)),
     $        bc(ibcoff))
        end if
       ibcoff=itmp
       end if
      end do
      return
      end
